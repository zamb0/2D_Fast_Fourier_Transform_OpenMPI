#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include "pgm.h"
#include "cshift.h"
#include <time.h>

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


int main(int argc, char** argv){

    int start = clock();

    pgm_t img;
    cplx* v_data;
    int o_width, o_height;

    // Read image
    img = pgm_read(argv[1]);

    o_height = img.height;
    o_width = img.width;

    // Check if padding is needed
    if(!is_power_of_two(img.width*img.height) || img.width != img.height){
        img = zeroPadding(img);
    }

    // Convert image to vector
    v_data = mat2vet(img.data, img.width, img.height);

    //################# START 2D FFT #################
    // Perform 1D FFT
    for(int i=0; i < img.height; i++){
        cooley_tukey_fft(v_data + i*img.width, img.width, 0);
    }
    //

    // Transpose vector
    v_data = transpose(v_data, img.width, img.height);

    // Perform 1D FFT
    for(int i=0; i < img.height; i++){
        cooley_tukey_fft(v_data + i*img.width, img.width, 0);
    }
    //################# END 2D FFT #################

    // Print the FFT image
    img.data = vet2mat(fftshift(v_data, img.width, img.height), img.width, img.height);

    // Write FFT image
    pgm_write_fft(img, "fft.pgm", "");

    v_data = ifftshift(mat2vet(img.data, img.width, img.height), img.width, img.height);


    //################# START 2D iFFT #################
    // Perform 1D iFFT
    for(int i=0; i < img.height; i++){
        cooley_tukey_fft(v_data + i*img.width, img.width, 1);
    }

    // Transpose vector
    v_data = transpose(v_data, img.width, img.height);

    // Perform 1D iFFT
    for(int i=0; i < img.height; i++){
        cooley_tukey_fft(v_data + i*img.width, img.width, 1);
    }

    // Divide by the number of pixels   
    for(int i=0; i < img.width*img.height; i++){
        v_data[i] /= (img.height*img.width);
    }
    //################# END 2D iFFT #################

    // Convert vector to matrix
    img.data = vet2mat(v_data, img.width, img.height);

    free(v_data);

    // Realloc
    img.data = (cplx**)realloc(img.data, o_height * sizeof(cplx*));
    for (int i = 0; i < o_height; i++) {
        img.data[i] = (cplx*)realloc(img.data[i], o_width * sizeof(cplx));
    }

    img.width = o_width;
    img.height = o_height;

    // Write inverse FFT image
    pgm_write(img, "ifft.pgm", "");
    
    free(img.data);

    int end = clock();

    printf("Time: %.10lf\n", (double)(end - start) / CLOCKS_PER_SEC);
    return 0;

}