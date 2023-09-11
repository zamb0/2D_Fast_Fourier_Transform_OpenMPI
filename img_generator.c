#include "pgm.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h> 

int main(int argc, char *argv[]){

    pgm_t image;

    image.width = atoi(argv[2]);
    image.height = atoi(argv[3]);
    strcpy(image.type, argv[4]);
    image.max = atoi(argv[5]);

    image.data = (cplx**)malloc(image.height * sizeof(cplx*));
    for (int i = 0; i < image.height; i++){
        image.data[i] = (cplx*)malloc(image.width * sizeof(cplx));
    }

    for(int i = 0; i < image.height; i++){
        for(int j = 0; j < image.width; j++){
            image.data[i][j] = rand() % image.max;
        }
    }

    pgm_write(image, argv[1], "");

    return 0;
}