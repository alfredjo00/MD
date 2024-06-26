/*
fft.h

Header file for fft.c
 
Created by Martin Gren on 2014-10-22.
*/

#ifndef _fft_h
#define _fft_h

extern void powerspectrum(double *, double *, int);

extern void fft_freq(double *, double, int);

extern void powerspectrum_shift(double *, int);

extern void fft_freq_shift(double *, double, int);

void fast_correlation(double *fft_dt, double *ifft_dt, double vel[][3], int n, int n_atoms);

#endif
