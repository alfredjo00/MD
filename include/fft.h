#ifndef _fft_h
#define _fft_h

extern void powerspectrum(double *, double *, int);

extern void fft_freq(double *, double, int);

extern void powerspectrum_shift(double *, int);

extern void fft_freq_shift(double *, double, int);

void fft_and_ifft(double *ifft_dt, double *fft_dt, double *vel, int n, int n_atoms, int n_zeros);

#endif
