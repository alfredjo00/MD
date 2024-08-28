#ifndef _fft_h
#define _fft_h

void fft_fast_vel_correlation(double vel[][3], int n_time, int n_atoms, double dt);

void temp_pressure_equilibrium_task(double T_eq, double deviation, char *fname);

#endif
