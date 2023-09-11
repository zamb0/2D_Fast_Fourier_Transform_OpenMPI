#include <complex.h>
#include "pgm.h"
#include <stdlib.h>

//https://stackoverflow.com/questions/5915125/fftshift-ifftshift-c-c-source-code

cplx* _circshift(const cplx* in, int xdim, int ydim, int xshift, int yshift)
{
    cplx* out = (cplx*)malloc(xdim * ydim * sizeof(cplx));
    
    for (int i = 0; i < xdim; i++) {
        int ii = (i + xshift) % xdim;
        for (int j = 0; j < ydim; j++) {
        int jj = (j + yshift) % ydim;
        out[ii * ydim + jj] = in[i * ydim + j];
        }
    }

    return out;
}

cplx* fftshift(cplx* in, int x, int y){
    return _circshift(in, x, y, (x/2), (y/2));
} 

cplx* ifftshift(cplx* in, int x, int y){
    return _circshift(in, x, y, ((x+1)/2), ((y+1)/2));
}