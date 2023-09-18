![C-Language](https://img.shields.io/badge/C%20Language-red)
![OpenMPI](https://img.shields.io/badge/OpenMPI-blue)

# 2D Fast Fourier Transform (2D FFT) OpenMPI Implementation

![fft_cover](https://github.com/zamb0/2D_Fast_Fourier_Transform_OpenMPI/assets/69969487/5c5e032d-43c7-44fc-8456-6a172bfda06f)

This repository contains the code for the Advanced Computer Architecture project (2022/2023).

## Aim of the Project

The aim of the project is to implement a 2D Fast Fourier Transform (FFT) algorithm using OpenMPI. The algorithm is implemented in C and it is based on the Cooley-Tukey iterative algorithm. The algorithm is implemented in two versions: a serial version and a parallel version. The parallel version is implemented using OpenMPI.
The input of the algorithm is a 2D PGM image and the output is the 2D FFT of the image and the 2D iFFT image.

## Execution Instructions

### 1 - Setup the Git Repository

If not already done, download and install ```git``` from [git-scm.com](https://git-scm.com). Then clone the repository with the following command:

```bash
git clone https://github.com/zamb0/ACA.git
```

### 2 - Install OpenMPI

Download Open MPI 4.1.1 from the official [Open MPI website](https://www.open-mpi.org/software/ompi/v4.1/) and install it.

### 3 - Compile

To compile and run the serial application, run the following commands:

```bash
gcc -Wall -o executableFile -lm -std=c99 *.c
```

you need to compile using all the .c files in the folder but ```fft_parallel.c``` and ```img_generator.c```.

If you want to compile and run the parallel application instead,
then compile in the following way:

```bash
mpicc -Wall -o executableFile -lm -std=c99 *.c
```

also in this case you need to compile using all the .c files in the folder but ```fft_serial.c``` and ```img_generator.c```.

### 3b - Generate the Images (if needed)

If you do not have images to run the 2D FFT algorithm on, you can compile the ```img_generator.c``` to generate sample images of your choice. To compile the image generator, run the following command:

```bash
gcc -Wall -o img_gen img_generator.c pgm.c
```

Then, to generate an image, run the following command:

```bash
./img_gen [outputFile] [width] [height] [pgmType] [maxValue]
```

where ```outputFile``` is the name of the file where the image will be saved, ```width``` is the width of the image, ```height``` is the height of the image, ```pgmType``` is the type of pgm image (i.e. P2 is the one allowed to be used in our project), ```maxValue``` is the maximum value of the pixels in the image.

### 4 - Run the Program

To compile and run the serial application, run the following commands:

```bash
./executableFile [inputImage]
```

If you want to run the parallel application instead:

```bash
mpirun -np N executableFile [inputImage]
```

where ```N``` is the number of processes you want to use.
