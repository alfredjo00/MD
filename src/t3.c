#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"
#include "verlet.h"
#include "fft.h"

#define SQ(X) ((X) * (X))

/*
 * Computes the velocity correlation function by
 * doing FFT on velocity matrix and then IFFT of the
 * squared FFT.
 * Write the IFFT data and squared FFT data to a txt file
 */
void fft_fast_vel_correlation(double vel[][3], int n_time, int n_atoms, double dt)
{
  int n_zeros = 200;
  int n_z = n_time + n_zeros * 2;

	double fft_dt[n_z],	ifft_dt[n_z], frequencies[n_z], time_range[n_z];  
  double *vel_ccat = malloc(sizeof(double[n_time * n_atoms * 3]));

  for (int k = 0; k < 3; ++k)
  {
    for (int i = 0; i < n_time; ++i)
    {
      for (int j = 0; j < n_atoms; ++j)
      {
        vel_ccat[i * n_atoms + j + k * n_atoms * n_time ] = vel[i * n_atoms + j][k];
      }
    }
  }

	arange(time_range, 0, n_z, dt);
	/*
	 * Does the FFT and IFFT and saves it to the first 2 arrays
	 */
	int f_data_cols = 4; // Number of columns to write
	double **f_data = (double**) malloc( f_data_cols * sizeof(double*));
  double *f_fft   = calloc(n_z, sizeof(double));
  double *f_ifft  = calloc(n_z, sizeof(double));
	
  fft_and_ifft(ifft_dt, fft_dt, vel_ccat, n_time, n_atoms, n_zeros);
  powerspectrum_shift(fft_dt, n_z);
		
  for (int j = 0; j < n_z; ++j){
    f_fft[j]  = (fft_dt[j]);
    f_ifft[j] = (ifft_dt[j]);
  }
  
  for (int j = 0; j < n_z; ++j){
    fft_dt[j]   = 0.;
    ifft_dt[j]  = 0.;
  }

	fft_freq_shift(frequencies, dt, n_z);

	f_data[0] = time_range;
	f_data[1] = f_ifft;
	f_data[2] = frequencies;
	f_data[3] = f_fft;

	write_to_file("data/fft_data.txt", f_data, n_z, f_data_cols);
	free(f_data);	    f_data   = NULL;
	free(f_fft);	    f_fft    = NULL;
	free(f_ifft);	    f_ifft   = NULL;
    free(vel_ccat);   vel_ccat = NULL;
}

void temp_pressure_equilibrium_task(double T_eq, double deviation, char *fname)
{
  char file_name_f[100];
  char file_name_pos[100];
  char file_name_vel[100];

  sprintf(file_name_f, "data/task4_%s.txt", fname);
  sprintf(file_name_pos, "data/task4_%s_pos.txt", fname);
  sprintf(file_name_vel, "data/task4_%s_vel.txt", fname);
  /*
   * Temperature and pressure to equilibrate to.
   * T_eq = (Celcius) + 273.15 >> in K
   * P_eq = (bar / 10 000) / const. >> in eV / Å^3
   * double deviation // 0.35 for liquid, 0.065 for solid.
   */
  double P_eq = (1 * 1e-4) / 160.21766208; // eV / Å^3

	/*
	 * Converts time units later when plotting
	 */
  double t_conversion = 3.335641e-7;
  double dt = 5; // fs
  dt *= 1e-3/ t_conversion;
  double a0 = 4.0479; // Ångström
  int N = 4; // Unit cells in each direction
  int n_atoms = 4*N*N*N;
  int n_timesteps = 25000;
  int N_eq = 15000;
  int time_its = n_timesteps - N_eq;

  double X[n_atoms][3];	
	
  // Mass of aluminium
  long double m = 2.51512636522e10; // eV   
  init_fcc(X, N, a0);
  //double E_pot = get_energy_AL(X , N * a0 , n_atoms);
  rand_deviations(X, n_atoms, a0, deviation);
	
  /*
   * Arrays to store position, velocity
   */
  double (*pos)[3]    	= malloc(sizeof(double[n_atoms * n_timesteps][3]));
  double (*vel)[3]    	= malloc(sizeof(double[n_atoms * n_timesteps][3]));

  /*
   * Temp arrays to store one time step of positon, velocity
   */
  double (*pos_sub)[3] 	= malloc(sizeof(double[n_atoms][3]));
  double (*vel_sub)[3] 	= malloc(sizeof(double[n_atoms][3]));
	
  /*
   * Arrrays to Pressure, Temperature, lattice const,
   * mean squared disp.  velocity correlation fn.
   */
  double *P_array 		= malloc(n_timesteps * sizeof(double));
  double *T_array 		= malloc(n_timesteps * sizeof(double));
  double *a0_array 		= malloc(n_timesteps * sizeof(double));
  double *msd_array 	= calloc(n_timesteps, sizeof(double));
  double *vfc_array 	= calloc(n_timesteps, sizeof(double));
	
  for (int i = 0; i < n_atoms; ++i){
    pos[i][0] = X[i][0];
    pos[i][1] = X[i][1];
    pos[i][2] = X[i][2];
  }

  /*
   * Time evolve system, more info in verlet.c
   */
  verlet_equilibration_scaling(n_timesteps, n_atoms,
			       N_eq, 	 N, 		pos, 	vel,
				   dt,		 m, 		a0, 	T_eq,
				   P_eq, 	 P_array,	T_array,
			       a0_array);

  
  double (*pos_eq)[3]    	= malloc(sizeof(double[n_atoms * time_its][3]));
  double (*vel_eq)[3]    	= malloc(sizeof(double[n_atoms * time_its][3]));

  for (int i = 0; i < time_its; ++i){
    for (int j = 0; j < n_atoms; ++j){
      pos_eq[j + i * n_atoms][0] = pos[j + (i + N_eq) * n_atoms][0];
      pos_eq[j + i * n_atoms][1] = pos[j + (i + N_eq) * n_atoms][1];
      pos_eq[j + i * n_atoms][2] = pos[j + (i + N_eq) * n_atoms][2];

      vel_eq[j + i * n_atoms][0] = vel[j + (i + N_eq) * n_atoms][0];
      vel_eq[j + i * n_atoms][1] = vel[j + (i + N_eq) * n_atoms][1];
      vel_eq[j + i * n_atoms][2] = vel[j + (i + N_eq) * n_atoms][2];
    }
  }
  
  mean_squared_displacement(msd_array, pos_eq, time_its, n_atoms);

  /*
   * Computes velocity correlation function
   */
  velocity_correlation(vfc_array, vel_eq, time_its, n_atoms);

  double time_range[n_timesteps];
  arange(time_range, 0, n_timesteps, dt );

  /*
   * Create array to store all values to later write to .txt
   */
  double **f_array; int fcols = 8;
  create_2D_array(&f_array, fcols, n_timesteps);
		

  #pragma omp parallel for
  for (int i = 0; i < n_timesteps; ++i){
    for (int j = 0; j < n_atoms; ++j){
      vel_sub[j][0] = vel[j + i * n_atoms][0];
      vel_sub[j][1] = vel[j + i * n_atoms][1];
      vel_sub[j][2] = vel[j + i * n_atoms][2];
			
      pos_sub[j][0] = pos[j + i * n_atoms][0];
      pos_sub[j][1] = pos[j + i * n_atoms][1];
      pos_sub[j][2] = pos[j + i * n_atoms][2];
    }
		
    f_array[0][i] = time_range[i];	
    f_array[1][i] = get_energy_AL(pos_sub, N * a0, n_atoms);
    f_array[2][i] = kinetic_energy(vel_sub, n_atoms, m);
    f_array[3][i] = P_array[i];
    f_array[4][i] = T_array[i];
    f_array[5][i] = a0_array[i];
    if ( i <= time_its){
      f_array[6][i] = msd_array[i];
      f_array[7][i] = vfc_array[i];
    }
    else {
      f_array[6][i] = 0;
      f_array[7][i] = 0;
    }
  }


  fft_fast_vel_correlation(vel_eq, time_its, n_atoms, dt * t_conversion);

  write_to_file(file_name_f, f_array, time_its, fcols);	

  write_positions(file_name_pos, pos, time_its, n_atoms, time_range);

  write_positions(file_name_vel, vel_eq, time_its, n_atoms, time_range);
    
  destroy_2D_array(f_array, fcols); f_array = NULL;
  free(pos); 		pos       = NULL;
  free(vel); 		vel       = NULL;
  free(pos_sub); 	pos_sub   = NULL;
  free(vel_sub); 	vel_sub   = NULL;
  free(pos_eq); 	pos_eq    = NULL;
  free(vel_eq); 	vel_eq    = NULL;
  free(P_array); 	P_array   = NULL;
  free(T_array); 	T_array   = NULL;
  free(a0_array);   a0_array  = NULL;
  free(msd_array);  msd_array = NULL;
  free(vfc_array);  vfc_array = NULL;

}