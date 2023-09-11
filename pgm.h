#ifndef PGM_H_
#define PGM_H_

#include <complex.h>

typedef double complex cplx;

typedef struct pgm{
    char type[3];
    int width;
    int height;
    int max;
    cplx **data;
}pgm_t;

void pgm_write(pgm_t, char *, char *);

void pgm_write_fft(pgm_t, char *, char *);

pgm_t pgm_read(char *);

#endif