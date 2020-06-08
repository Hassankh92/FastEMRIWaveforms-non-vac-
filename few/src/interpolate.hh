#ifndef __INTERP_H__
#define __INTERP_H__


#include "kernel.hh"
#include "global.h"
#include "cusparse.h"

void interpolate_arrays(double *t_arr, double *y_all, double *c1, double *c2, double *c3, int ninterps, int length, double *B, double *upper_diag, double *diag, double *lower_diag);

void get_waveform(cmplx *d_waveform, double *y_vals, double *c1, double *c2, double *c3,
              int *d_m, int *d_n, int init_len, int out_len, int num_teuk_modes, cmplx *d_Ylms,
              double delta_t, double *h_t);

/*
CuSparse error checking
*/
#define ERR_NE(X,Y) do { if ((X) != (Y)) { \
                            fprintf(stderr,"Error in %s at %s:%d\n",__func__,__FILE__,__LINE__); \
                            exit(-1);}} while(0)

#define CUDA_CALL(X) ERR_NE((X),cudaSuccess)
#define CUSPARSE_CALL(X) ERR_NE((X),CUSPARSE_STATUS_SUCCESS)


#endif // __INTERP_H__
