#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <mpi.h>
#include <math.h>
#include "pgm.h"
#include "cshift.h"

#define PI 3.14159265358979323846

typedef double complex cplx;

// Function to convert a matrix in form of a vector
// mat = matrix
// width = number of columns
// height = number of rows
// return = vector
cplx* mat2vet(cplx** mat, int width, int height){
	cplx *v = (cplx*)malloc(height * width * sizeof(cplx));
	for(int i = 0; i < height; i++){
        for (int j = 0; j < width; j++){
            v[i * width + j] = mat[i][j];
        }
    }
	return v;
}

// Function to convert a vector in form of a matrix
// v = vector
// width = number of columns
// height = number of rows
// return = matrix
cplx** vet2mat(cplx* v, int width, int height){
	cplx **mat = (cplx**)malloc(height * sizeof(cplx*));
	for(int i = 0; i < height; i++){
		mat[i] = (cplx*)malloc(width * sizeof(cplx));
		for (int j = 0; j < width; j++){
			mat[i][j] = v[i * width + j];
		}
	}
	return mat;
}

// Function to calculate the next power of 2
// num = number to calculate the next power of 2
// return = next power of 2
int nextPowerOf2(int num) {
    int power = 1;
    while (power < num) {
        power *= 2;
    }
    return power;
}

// Function to check if a number is a power of 2
// x = number to check
// return = 1 if x is a power of 2, 0 otherwise
int is_power_of_two(int x)
{
    return (x != 0) && ((x & (x - 1)) == 0);
}

// Function to perform zero padding on an image
// image = image to pad
// return = padded image
pgm_t zeroPadding(const pgm_t image) {

    pgm_t paddedImage;
    int newWidth = nextPowerOf2(image.width);
    int newHeight = nextPowerOf2(image.height);

    // If the image is not square, pad the image to make it square
    if (newWidth != newHeight) {
        newWidth = newHeight = (newWidth > newHeight) ? newWidth : newHeight;
    }

    // Allocate memory for the padded image
    paddedImage.data = (cplx**)malloc(newHeight * sizeof(cplx*));
    for (int i = 0; i < newHeight; i++) {
        paddedImage.data[i] = (cplx*)calloc(newWidth, sizeof(cplx));
    }

    // Copy the image into the padded image
    for (int i = 0; i < image.height; i++) {
        for (int j = 0; j < image.width; j++) {
            paddedImage.data[i][j] = image.data[i][j];
        }
    }

    // Copy the info of the image into the padded image
    paddedImage.width = newWidth;
    paddedImage.height = newHeight;
    paddedImage.max = image.max;
    strcpy(paddedImage.type, image.type);

    return paddedImage;
}

// Function to perform Cooley-Tukey FFT
// x = input vector
// N = size of the input vector
// Forwards if inverse = 0, backwards if inverse = 1
void cooley_tukey_fft(cplx x[], int N, int inverse) {
    // Bit-reversal permutation
    int i, j, k;
    for (i = 1, j = N / 2; i < N - 1; i++) {
        if (i < j) {
            cplx temp = x[i];
            x[i] = x[j];
            x[j] = temp;
        }
        k = N / 2;
        while (k <= j) {
            j -= k;
            k /= 2;
        }
        j += k;
    }

    // Iterative FFT or IFFT
    double sign = (inverse) ? 1.0 : -1.0; // Sign for IFFT
    for (int s = 1; s <= log2(N); s++) {
        int m = 1 << s; // Subproblem size
        cplx omega_m = cexp(sign * I * 2.0 * PI / m);

        for (int k = 0; k < N; k += m) {
            cplx omega = 1.0;

            for (int j = 0; j < m / 2; j++) {
                cplx t = omega * x[k + j + m / 2];
                cplx u = x[k + j];
                x[k + j] = u + t;
                x[k + j + m / 2] = u - t;
                omega *= omega_m;
            }
        }
    }
}

// Function to transpose a matrix in form of a vector
// v = vector
// width = number of columns
// height = number of rows
// return = transposed vector
cplx* transpose(cplx* v, int width, int height){
	cplx *tmp = (cplx*)malloc(height * width * sizeof(cplx));

	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			tmp[j * height + i] = v[i * width + j];
		}
	}

	return tmp;
}


int main(int argc, char** argv) {

    double start_time = MPI_Wtime();

    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    pgm_t img;
    cplx *v_send;
    cplx *v_revc;
    int len_info[5];

    if(rank == 0){

        img = pgm_read(argv[1]);

        len_info[3] = img.width;
        len_info[4] = img.height;

        // Check if padding is needed
        if(!is_power_of_two(img.width*img.height) || img.width != img.height){
            img = zeroPadding(img);
        }
       
        len_info[0] = img.width;
        len_info[1] = img.height;
        len_info[2] = img.max;

        // Allocate memory for the send vector
        v_send = (cplx*)malloc(img.width * img.height * sizeof(cplx));
        // Convert image to vector
        v_send = mat2vet(img.data, img.width, img.height);

    }

    // Broadcast the usefull length information
    MPI_Bcast(len_info, 5, MPI_INT, 0, MPI_COMM_WORLD);

    // Find the number of rows per processor
    int rows_per_processor = len_info[1] / size;
    int remainder = len_info[1] % size;
    
    // Find the number of processors that will receive an extra row
    int processors_with_extra_rows = (rank < remainder) ? 1 : 0;
    
    // Find the number of rows to send to this processor
    int my_num_rows = rows_per_processor + processors_with_extra_rows;
    
    // Calculate the displacements for MPI_Scatterv
    int* displacements = (int*)malloc(size * sizeof(int));
    int* recvcounts = (int*)malloc(size * sizeof(int));

    // Calculate the displacements and the recvcounts
    int displacement = 0;
    for (int i = 0; i < size; i++) {
        recvcounts[i] = rows_per_processor * len_info[0];
        if (i < remainder) {
            recvcounts[i] += len_info[0]; // Distribute remaining rows
        }
        displacements[i] = displacement;
        displacement += recvcounts[i];
    }

    // Allocate memory for the received vector 
    v_revc = (cplx*)malloc(my_num_rows * len_info[0] * sizeof(cplx));

    // Scatter the data
    MPI_Scatterv(v_send, recvcounts, displacements, MPI_C_DOUBLE_COMPLEX, v_revc, my_num_rows * len_info[0], MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    
    //#################### Start 2D FFT ####################
    // Perform 1D FFT
    for(int i = 0; i < my_num_rows; i++)
	{
        cooley_tukey_fft(v_revc + i * len_info[0], len_info[0], 0);
    }

    // Gather the data
    MPI_Gatherv(v_revc, my_num_rows * len_info[0], MPI_C_DOUBLE_COMPLEX, v_send, recvcounts, displacements, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    if(rank == 0){
        // Transpose
        v_send = transpose(v_send, len_info[0], len_info[1]);

    }

    // Scatter the data
    MPI_Scatterv(v_send, recvcounts, displacements, MPI_C_DOUBLE_COMPLEX, v_revc, my_num_rows * len_info[0], MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    // Perform 1D FFT
    for(int i = 0; i < my_num_rows; i++)
	{
        cooley_tukey_fft(v_revc + i * len_info[0], len_info[0], 0);
    }

    //#################### End 2D FFT ####################

    // Gather the data
    MPI_Gatherv(v_revc, my_num_rows * len_info[0], MPI_C_DOUBLE_COMPLEX, v_send, recvcounts, displacements, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    if(rank == 0){
        // Print the FFT image
        img.data = vet2mat(fftshift(v_send, len_info[0], len_info[1]), len_info[0], len_info[1]);
        
        // Write FFT image
        pgm_write_fft(img, "fft.pgm", "");

        v_send = ifftshift(mat2vet(img.data, len_info[0], len_info[1]), len_info[0], len_info[1]);

    }

    // Scatter the data
    MPI_Scatterv(v_send, recvcounts, displacements, MPI_C_DOUBLE_COMPLEX, v_revc, my_num_rows * len_info[0], MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    //#################### Start 2D FFT ####################
    // Perform 1D FFT
    for(int i = 0; i < my_num_rows; i++)
	{
        cooley_tukey_fft(v_revc + i * len_info[0], len_info[0], 1);
    }

    // Gather the data
    MPI_Gatherv(v_revc, my_num_rows * len_info[0], MPI_C_DOUBLE_COMPLEX, v_send, recvcounts, displacements, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    if(rank == 0){
        // Transpose
        v_send = transpose(v_send, len_info[0], len_info[1]);
       
    }

    // Scatter the data
    MPI_Scatterv(v_send, recvcounts, displacements, MPI_C_DOUBLE_COMPLEX, v_revc, my_num_rows * len_info[0], MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
   
    // Perform 1D FFT
    for(int i = 0; i < my_num_rows; i++)
	{
        cooley_tukey_fft(v_revc + i * len_info[0], len_info[0], 1);
    }

    //#################### End 2D FFT ####################
    //(missing the division by the number of elements, we will do it after the gather)

    // Gather the data
    MPI_Gatherv(v_revc, my_num_rows * len_info[0], MPI_C_DOUBLE_COMPLEX, v_send, recvcounts, displacements, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);


    if(rank == 0){

        // Divide by the number of elements
        for(int i=0; i<len_info[0]*len_info[1]; i++){
            v_send[i] /= (double)(len_info[0]*len_info[1]);
        }

        img.data = vet2mat(v_send, len_info[0], len_info[1]);

        free(v_send);

        img.data = realloc(img.data, len_info[4] * sizeof(cplx*));
        for (int i = 0; i < len_info[4]; i++) {
            img.data[i] = realloc(img.data[i], len_info[3] * sizeof(cplx));
        }

        img.width = len_info[3];
        img.height = len_info[4];

        // Write the ifft
        pgm_write(img, "ifft.pgm", "");

        free(img.data);
    }

    free(v_revc);
    free(displacements);
    free(recvcounts);

    double end_time = MPI_Wtime();

    if(rank == 0){
        printf("Time: %.10lf\n", end_time - start_time);
    }

    MPI_Finalize();
    return 0;
}



