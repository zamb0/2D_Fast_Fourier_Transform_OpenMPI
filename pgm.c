#include "pgm.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>

typedef double complex cplx;

pgm_t log_scale(pgm_t img){

    double c = 255/log(1 + img.max);

    for(int i = 0; i < img.height; i++){
        for(int j = 0; j < img.width; j++){
            img.data[i][j] = c * log(1 + cabs(img.data[i][j]));
        }
    }

    return img;
}

void pgm_write(pgm_t img, char *fabs, char *farg){

    if(strcmp(farg, "") != 0){

        FILE *fpabs;
        FILE *fparg;

        fpabs = fopen(fabs, "wb");
        fparg = fopen(farg, "wb");

        if(fpabs == NULL || fparg == NULL){
            printf("Error opening file\n");
            exit(1);
        }

        //printf("writing header\n");

        fprintf(fpabs, "%s\n", img.type);
        fprintf(fparg, "%s\n", img.type);

        fprintf(fpabs, "%d %d\n", img.width, img.height);
        fprintf(fparg, "%d %d\n", img.width, img.height);

        fprintf(fpabs, "%d\n", img.max);
        fprintf(fparg, "%d\n", img.max);

        //printf("writing data\n");
        for(int i = 0; i < img.height; i++){
            for(int j = 0; j < img.width; j++){

                if(isinf(cabs(img.data[i][j]))){
                    printf("inf\n");
                    exit(1);
                }

                fprintf(fpabs, "%0.lf ", cabs(img.data[i][j]));
                fprintf(fparg, "%0.lf ", carg(img.data[i][j]));
            }
            fprintf(fpabs, "\n");
            fprintf(fparg, "\n");
        }

        fclose(fpabs);
        fclose(fparg);
    } else {

        FILE *fpabs;

        fpabs = fopen(fabs, "wb");

        if(fpabs == NULL){
            printf("Error opening file\n");
            exit(1);
        }

        //printf("writing header\n");

        fprintf(fpabs, "%s\n", img.type);

        fprintf(fpabs, "%d %d\n", img.width, img.height);

        fprintf(fpabs, "%d\n", img.max);

        //printf("writing data\n");
        for(int i = 0; i < img.height; i++){
            for(int j = 0; j < img.width; j++){

                if(isinf(cabs(img.data[i][j]))){
                    printf("inf\n");
                    exit(1);
                }

                fprintf(fpabs, "%0.lf ", cabs(img.data[i][j]));
            }
            fprintf(fpabs, "\n");
        }

        fclose(fpabs);
    }

} 

void pgm_write_fft(pgm_t img, char *fabs, char *farg){

    if (strcmp(farg, "") != 0)
    {
        FILE *fpabs;
        FILE *fparg;

        fpabs = fopen(fabs, "wb");
        fparg = fopen(farg, "wb");

        if (fpabs == NULL || fparg == NULL)
        {
            printf("Error opening file\n");
            exit(1);
        }

        fprintf(fpabs, "%s\n", img.type);
        fprintf(fparg, "%s\n", img.type);

        fprintf(fpabs, "%d %d\n", img.width, img.height);
        fprintf(fparg, "%d %d\n", img.width, img.height);

        fprintf(fpabs, "%d\n", img.max);
        fprintf(fparg, "%d\n", img.max);

        double img_max = 0;
        for(int i = 0; i < img.height; i++){
            for(int j = 0; j < img.width; j++){
                if(cabs(img.data[i][j]) > img_max){
                    img_max = cabs(img.data[i][j]);
                }
            }
        }

        double c = 255/log(1 + img_max);

        for(int i = 0; i < img.height; i++){
            for(int j = 0; j < img.width; j++){

                double img_tmp = c * log(1 + cabs(img.data[i][j]));

                fprintf(fpabs, "%0.lf ", img_tmp);
                fprintf(fparg, "%0.lf ", carg(img.data[i][j]));
            }
            fprintf(fpabs, "\n");
            fprintf(fparg, "\n");
        }

        fclose(fpabs);
        fclose(fparg);
    } else {

        FILE *fpabs;

        fpabs = fopen(fabs, "wb");

        if(fpabs == NULL){
            printf("Error opening file\n");
            exit(1);
        }

        fprintf(fpabs, "%s\n", img.type);

        fprintf(fpabs, "%d %d\n", img.width, img.height);

        fprintf(fpabs, "%d\n", img.max);

        double img_max = 0;
        for(int i = 0; i < img.height; i++){
            for(int j = 0; j < img.width; j++){
                if(cabs(img.data[i][j]) > img_max){
                    img_max = cabs(img.data[i][j]);
                }
            }
        }

        double c = 255/log(1 + img_max);

        for(int i = 0; i < img.height; i++){
            for(int j = 0; j < img.width; j++){

                double img_tmp = c * log(1 + cabs(img.data[i][j]));

                fprintf(fpabs, "%0.lf ", img_tmp);
            }
            fprintf(fpabs, "\n");
        }

        fclose(fpabs);
    }
    
    

} 

pgm_t pgm_read(char *filename){
    FILE *fp;
    pgm_t img;

    fp = fopen(filename, "rb");

    if(fp == NULL){
        printf("Error opening file\n");
        exit(1);
    }

    char buff[2048];
    int tmp;

    fgets(buff, sizeof(buff), fp);
    sscanf(buff, "%s", img.type);

    fgets(buff, sizeof(buff), fp);
    sscanf(buff, "%d %d", &img.width, &img.height);

    fgets(buff, sizeof(buff), fp);
    sscanf(buff, "%d", &img.max);

    img.data = (cplx**)malloc(img.height * sizeof(cplx*));

    for(int i = 0; i < img.height; i++){
        img.data[i] = (cplx*)malloc(img.width * sizeof(cplx));

        for (int j = 0; j < img.width; j++){
            fscanf(fp, "%d", &tmp);
            img.data[i][j] = (cplx)tmp;
        }
        
    }

    fclose(fp);
    return img;
}