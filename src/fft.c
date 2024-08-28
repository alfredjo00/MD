#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <complex.h>
#include <omp.h>


#define SQ(X) ((X) * (X))


/* Makes fft of data and returns the powerspectrum in powspec_data */
void powerspectrum(
	double *data, 
	double *powspec_data, 
	int n) /* input data, output powspec_data, number of timesteps */
{
	/* Declaration of variables */
	int i;
	double complex_coefficient[2*n]; // array for the complex fft data
	double data_cp[n]; 

	/*make copy of data to avoid messing with data in the transform*/
	for (i = 0; i < n; i++)	{
		data_cp[i] = data[i];
	}

	/*Declare wavetable and workspace for fft*/
	gsl_fft_real_wavetable *real;
	gsl_fft_real_workspace *work;

	/*Allocate space for wavetable and workspace for fft*/
	work = gsl_fft_real_workspace_alloc(n);
	real = gsl_fft_real_wavetable_alloc(n);

	/*Do the fft*/
	gsl_fft_real_transform(data_cp, 1, n, real, work);	

	/*Unpack the output into array with alternating real and imaginary part*/	
	gsl_fft_halfcomplex_unpack(data_cp, complex_coefficient,1,n);
	
	/*fill the output powspec_data with the powerspectrum */
	for (i = 0; i < n; i++)	{
		powspec_data[i] = (complex_coefficient[2*i]*complex_coefficient[2*i]+complex_coefficient[2*i+1]*complex_coefficient[2*i+1])/n; 
	}
	
	/*Free memory of wavetable and workspace*/
	gsl_fft_real_wavetable_free(real);
	gsl_fft_real_workspace_free(work);
}

/* Shifts the input powspec_data to center the 0 frequency */
void powerspectrum_shift(
	double *powspec_data, 
	int n) /* input data, timestep, number of timesteps */
{
	/* Declaration of variables */
	int i;
	
	/* make copy of fft_data as reference for the shift */ 
	double powspec_cp[n];
	for (i = 0; i < n; i++)	{
		powspec_cp[i] = powspec_data[i];
	}

	/* make shift */
	for (i = 0; i < n; i++)	{
		if (n % 2) /*if n odd*/	{ 
			if (i<=(n-2)/2)	{
				powspec_data[i] = powspec_cp[(i+(n+1)/2)];
			}
			else {
				powspec_data[i] = powspec_cp[(i+(n+1)/2)%(n)];
			}			
		}
		else {
			if (i<n/2) {
				powspec_data[i] = powspec_cp[i+n/2];
			}
			else {
				powspec_data[i] = powspec_cp[(i+n/2)%(n)];
			}			
		}
	}
}

/* Makes a frequency array fft_freq with frequency interval 1/(dt*n) */
void fft_freq(
	double *fft_freq, 
	double dt, 
	int n) /* output frequency array, timestep, number of timesteps */
{
	/* Declaration of variables */
    	int i;
	/* Fill the output aaray with frequencies */
	for (i = 0; i < n; i++)	{
		fft_freq[i] = i/dt/n;
	}
}

/* Makes a frequency array fft_freq with frequency interval 1/(dt*n) with a centered O frequency */
void fft_freq_shift(
	double *fft_freq,
	double dt,
	int n) /* output frequency array, timestep, number of timesteps */
{
	/* Declaration of variables */
    	int i;
	/* Fill the output aaray with shifted frequencies */
	for (i = 0; i < n; i++)	{
		if (n % 2) /*if n odd*/	{ 
			if (i<=(n-2)/2)	{
				fft_freq[i] = (-(n-1)/2+i)/dt/n;
			}
			else {
				fft_freq[i] = (i-(n-1)/2)/dt/n;
			}			
		}
		else {
			if (i<n/2) {
				fft_freq[i] = (-n/2+i)/dt/n;
			}
			else {
				fft_freq[i] = (i-n/2)/dt/n;
			}			
		}
	}
}

/*
 * Does FFT of velocities, squares the FFT, does the IFFT,
 * returns the IFFT data and FFT data.
 *
 * ifft_dt  => store IFFT data
 * fft_dt   => store FFT data
 * vel[][3] => velocity matrix to FFT [n_atoms * timesteps][x y z]
 * n        => number of time steps
 * n_atoms  => number of atoms
 *
 */
void fft_and_ifft(
	double *ifft_dt, 
	double *fft_dt, 
	double *vel, 
	int n, 
	int n_atoms, 
	int n_zeros)
{
	

	int N_total = n + 2 * n_zeros;

	/* Declaration of variables */
	double cc[2*N_total]; // array for the complex fft data
	double cc2[2*N_total]; // array for the sq complex fft data
	double cc2_s[2*N_total]; // array for the sq complex fft data
	double icc2[2*N_total]; // array for the inv sq complex fft data


	/*Declare wavetable and workspace for fft*/

	gsl_fft_halfcomplex_wavetable * hc;
	gsl_fft_real_wavetable *real;
	gsl_fft_real_workspace *work;

	/*Allocate space for wavetable and workspace for fft*/
	work = gsl_fft_real_workspace_alloc(N_total);
	real = gsl_fft_real_wavetable_alloc(N_total);
	hc   = gsl_fft_halfcomplex_wavetable_alloc(N_total);

	/*Do the fft*/
	double *v_sum = malloc(sizeof(double[N_total * n_atoms * 3]));

	// Summing up all atoms to FFT
	/*
	*	v_sum = sum_x0 sum_y0 sum_z0 sum_x1 sum_y1 sum_z1 ....
	*/
	for (int j = 0; j < N_total; ++j){
		fft_dt[j] 		= 0.;
		ifft_dt[j] 		= 0.;
		cc2_s[2*j] 		= 0.;
		cc2_s[2*j + 1] 	= 0.;
	}

	for (int k = 0; k < 3; ++k){
		for (int i = 0; i < n_atoms; ++i){
			for (int j = 0; j < n; ++j){
				v_sum[j + n_zeros + i * N_total + k * N_total * n_atoms] = vel[j * n_atoms + i + k * n_atoms * n];	
			}
		}
	}
	for (int i = 0; i < 3; ++i){
		for (int k = 0; k < n_atoms; ++k){
			// FFT on v_sum
			gsl_fft_real_transform(&v_sum[k * N_total + i * N_total * n_atoms], 1, N_total, real, work);
			gsl_fft_halfcomplex_unpack(&v_sum[k * N_total + i * N_total * n_atoms], cc, 1, N_total);

			// Computes the square of FFT of v_sum
			// cc[2*j] = REAL, cc[2*j + 1] = IMAGINARY
			#pragma omp parallel for
			for (int j = 0; j < N_total; ++j){
				cc2[2*j] = SQ( cc[2*j]) - SQ( cc[2*j + 1]); // Real = a^2 - b^2
				cc2[2*j + 1] = 2 *  cc[2*j] * cc[2*j + 1]; // Complex = 2 a b i

				cc2_s[2*j] += cc2[2*j];
				cc2_s[2*j + 1] += cc2[2*j + 1];
				fft_dt[j] += SQ(cc2[2*j]) + SQ(cc2[2*j + 1]);
			}	
		}
	}	


	for (int j = 0; j < N_total; ++j){
		fft_dt[j] /= (n_atoms * N_total * 3);
		cc2_s[2*j] /= n_atoms * 3;
		cc2_s[2*j + 1] /= n_atoms * 3;
	}

	gsl_fft_real_wavetable_free(real);

	double abs_cc2_s[N_total];
	
	for (int i = 0; i < N_total; ++i){
		abs_cc2_s[i] = sqrt(SQ(cc2_s[2*i]) +  SQ(cc2_s[2*i + 1]));
	} 

	gsl_fft_halfcomplex_transform(abs_cc2_s, 1, N_total, hc, work);
	gsl_fft_halfcomplex_unpack(abs_cc2_s, icc2, 1, N_total);

	for (int i = 0; i < N_total; ++i){
		ifft_dt[i] = sqrt(SQ(icc2[2*i]) + SQ(icc2[2*i + 1]))/ N_total;
	}

	/*Free memory of wavetable and workspace*/
	gsl_fft_halfcomplex_wavetable_free(hc);
	gsl_fft_real_workspace_free(work);
	free(v_sum); v_sum = NULL;
}