#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"
#include "verlet.h"
#include "fft.h"

#define SQ(X) ((X) * (X))

// Data for plot of Energy/Unit cell - Volume of Al: Task 1
void write_Al_E_pots()
{
  int n_points = 15;
  double a0[n_points];
  linspace(a0, 4.0, 4.082, n_points);
	
  int N = 4;
  int n_atoms = 4 * N * N * N;
  double X[n_atoms][3];
	
  double E_potentials[n_points]; // Energy pot / unit cell
  double volume_Al[n_points];
  for (int i = 0; i < n_points; ++i){
    init_fcc(X, N, a0[i]);		
    E_potentials[i] = get_energy_AL(X, N * a0[i], n_atoms)
      / (N * N * N);
    volume_Al[i] = (a0[i] * a0[i] * a0[i]);
  }	
	
  double **write_array = malloc(2. * sizeof(double*));
	
  write_array[0] = volume_Al; 
  write_array[1] = E_potentials;
	
  write_to_file("data/volume_Al-energy_pot.txt",
		/**/		  write_array, n_points, 2);
	
  free(write_array); write_array = NULL;
}

/*
 * Computes the instantaneous temperature
 */
double temperature(int n_atoms, double kinetic_energy)
{
  double k_b = 8.617333262e-5;	
  return kinetic_energy * 2. / (3. * n_atoms * k_b);
}

/*
 * Used for task 2
 * not relevant for other tasks
 */
void task2()
{
	double t_conversion = 3.335641e-7;
	double dt = 5; // fs
	dt *= 1e-3/ t_conversion;
	double a0 = 4.05; // Ångström
	double deviation = 0.065;
	int N = 4; // Unit cells in each direction
	int n_atoms = 4*N*N*N;
	int n_timesteps = 20000;
	int skip = 1; // Must skip % n_timesteps == 0
	int time_its = n_timesteps/skip;

	double X[n_atoms][3];


	// Mass of aluminium
	long double m = 2.51512636522e10; // eV
	init_fcc(X, N, a0);
	//double E_pot = get_energy_AL(X , N * a0 , n_atoms);
	rand_deviations(X, n_atoms, a0, deviation);
	
	double (*pos)[3] = malloc(sizeof(double[n_atoms * time_its][3]));
	double (*vel)[3] = malloc(sizeof(double[n_atoms * time_its][3]));
	double (*pos_sub)[3] = malloc(sizeof(double[n_atoms][3]));
	double (*vel_sub)[3] = malloc(sizeof(double[n_atoms][3]));
	
	for (int i = 0; i < n_atoms; ++i){
		pos[i][0] = X[i][0];
		pos[i][1] = X[i][1];
		pos[i][2] = X[i][2];
	}
	
	verlet_vel_pos(n_timesteps, n_atoms, skip,
	 pos, vel, dt, m, N * a0);
	
	double time_range[time_its];
	arange(time_range, 0, time_its, dt * skip * t_conversion);
	
	double **f_array; int fcols = 4;
	create_2D_array(&f_array, fcols, time_its);
	
	for (int i = 0; i < time_its; ++i){
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
		f_array[2][i] = get_virial_AL(pos_sub, N * a0, n_atoms);
		f_array[3][i] = kinetic_energy(vel_sub, n_atoms, m);
	}
	
	write_to_file("data/task2_0.txt", f_array,
	time_its, fcols);
	
	double average_ke = 0.;

	for (int i = 500; i < n_timesteps; ++i){
		average_ke += f_array[3][i] / (n_timesteps - 500);
	}

	double T = temperature(n_atoms, average_ke);
	printf("Kinetic energy %f eV\nTemperature %f C\n", average_ke, T - 273.15);
	
	destroy_2D_array(f_array, fcols); f_array = NULL;
	free(pos); 		pos = NULL;
	free(vel); 		vel = NULL;
	free(pos_sub); 	pos_sub = NULL;
	free(vel_sub); 	vel_sub = NULL;
}


/*
 * Write trajectories of three particles to txt file
 */
void write_positions(char *fname, double pos[][3],
		     int n_time,
		     int n_atoms,
		     double *timerange)
{
  double **f_array;
  create_2D_array(&f_array, 10, n_time);

  for (int i = 0; i < n_time; ++i){
    f_array[0][i] = timerange[i];
    f_array[1][i] = pos[0 +   i * n_atoms][0];
    f_array[2][i] = pos[0 +   i * n_atoms][1];
    f_array[3][i] = pos[0 +   i * n_atoms][2];
    f_array[4][i] = pos[49 +  i * n_atoms][1];
    f_array[5][i] = pos[49 +  i * n_atoms][2];
    f_array[6][i] = pos[49 +  i * n_atoms][3];
    f_array[7][i] = pos[127 + i * n_atoms][1];
    f_array[8][i] = pos[127 + i * n_atoms][2];
    f_array[9][i] = pos[127 + i * n_atoms][3];
  }

  write_to_file(fname, f_array, n_time, 10);

  destroy_2D_array(f_array, 10);
}


/*
 * Computes the velocity correlation function by
 * doing FFT on velocity matrix and then IFFT of the
 * squared FFT.
 * Write the IFFT data and squared FFT data to a txt file
 */
void fft_task6(double vel[][3], int n_time, int n_atoms, double dt, int n_eq)
{

	int add_zeros = 0; // Not utilized at the moment
	int N = add_zeros + n_time - n_eq;

	double fft_dt[N],		ifft_dt[N],
		   frequencies[N], time_range[N];

	double (*vel_add)[3]    	= malloc(sizeof(double[n_atoms * N][3]));

	for (int i = 0; i < N; ++i ){
		for (int j = 0; j < n_atoms; ++j){
			if (i < n_time){
				vel_add[j + i * n_atoms][0] = vel[j + (n_eq + i) * n_atoms][0];
				vel_add[j + i * n_atoms][1] = vel[j + (n_eq + i) * n_atoms][1];
				vel_add[j + i * n_atoms][2] = vel[j + (n_eq + i) * n_atoms][2];
			}
			else {
				vel_add[j + i * n_atoms][0] = 0.0;
				vel_add[j + i * n_atoms][1] = 0.0;
				vel_add[j + i * n_atoms][2] = 0.0;
			}
		}
	}

	arange(time_range, 0, N, dt);
	/*
	 * Does the FFT and IFFT and saves it to the first 2 arrays
	 */
	fast_correlation(ifft_dt, fft_dt, vel_add, N, n_atoms);

	int f_data_cols = 4; // Number of columns to write
	double **f_data = (double**) malloc( f_data_cols * sizeof(double*));

	powerspectrum_shift(fft_dt, N);
	fft_freq_shift(frequencies, dt, N);

	f_data[0] = time_range;
	f_data[1] = ifft_dt;
	f_data[2] = frequencies;
	f_data[3] = fft_dt;


	write_to_file("data/fft_data.txt", f_data, N, f_data_cols);
	free(f_data);	f_data = NULL;
	free(vel_add);  vel_add = NULL;
}



void run_later_tasks()
{

	/*
	 * Converts time units later when plotting
	 */
  double t_conversion = 3.335641e-7;
  double dt = 5; // fs
  dt *= 1e-3/ t_conversion;
  double a0 = 4.0479; // Ångström
  double deviation = 0.35; // 0.35 for liquid, 0.065 for solid.
  int N = 4; // Unit cells in each direction
  int n_atoms = 4*N*N*N;
  int n_timesteps = 20000;
  int skip = 1; // Must skip % n_timesteps == 0z
  int time_its = n_timesteps/skip;
  int N_eq = 10000;

  double X[n_atoms][3];		
	
  // Mass of aluminium
  long double m = 2.51512636522e10; // eV   
  init_fcc(X, N, a0);
  //double E_pot = get_energy_AL(X , N * a0 , n_atoms);
  rand_deviations(X, n_atoms, a0, deviation);
	
  /*
   * Arrays to store position, velocity
   */
  double (*pos)[3]    	= malloc(sizeof(double[n_atoms * time_its][3]));
  double (*vel)[3]    	= malloc(sizeof(double[n_atoms * time_its][3]));

  /*
   * Temp arrays to store one time step of positon, velocity
   */
  double (*pos_sub)[3] 	= malloc(sizeof(double[n_atoms][3]));
  double (*vel_sub)[3] 	= malloc(sizeof(double[n_atoms][3]));
	
  /*
   * Arrrays to Pressure, Temperature, lattice const,
   * mean squared disp.  velocity correlation fn.
   */
  double *P_array 		= malloc(time_its * sizeof(double));
  double *T_array 		= malloc(time_its * sizeof(double));
  double *a0_array 		= malloc(time_its * sizeof(double));
  double *msd_array 	= calloc(time_its, sizeof(double));
  double *vfc_array 	= calloc(time_its, sizeof(double));
	
  for (int i = 0; i < n_atoms; ++i){
    pos[i][0] = X[i][0];
    pos[i][1] = X[i][1];
    pos[i][2] = X[i][2];
  }

  /*
   * Temperature and pressure to equilibrate to.
   * T_eq = (Celcius) + const >> in K
   * P_eq = (bar / 10 000) / const. >> in eV / Å^3
   */
  double T_eq = 700 + 273.15; // K
  double P_eq = 0.0001 / 160.21766208; // eV / Å^3
  printf("P_eq %f\n",P_eq);
  /*
   * Time evolve system, more info in verlet.c
   */
  verlet_equilibration_scaling(n_timesteps, n_atoms,
			       skip, 	 N, 		pos, 	vel,
				   dt,		 m, 		a0, 	T_eq,
				   P_eq, 	 P_array,	T_array,
			       a0_array);

  /*
   * Computes Mean squared displacement
   */
  mean_squared_displacement(msd_array, pos, time_its, n_atoms, N_eq);

  /*
   * Computes velocity correlation function
   */
  velocity_correlation(vfc_array, vel, time_its, n_atoms, N_eq);

  double time_range[time_its];
  arange(time_range, 0, time_its, dt * skip);

  /*
   * Create array to store all values to later write to .txt
   */
  double **f_array; int fcols = 8;
  create_2D_array(&f_array, fcols, time_its);
		
  for (int i = 0; i < time_its; ++i){
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
    f_array[6][i] = msd_array[i];
    f_array[7][i] = vfc_array[i];
  }
  write_to_file("data/task4_liq.txt", f_array,
		time_its, fcols);	



  write_positions("data/task4_liq_pos.txt", pos, time_its,
		  n_atoms, time_range);
  

  //double T = temperature(n_atoms, 25.0);
  //printf("Temperature %f K\n", T);


  fft_task6(vel, time_its, n_atoms, dt * skip * t_conversion, N_eq);


  destroy_2D_array(f_array, fcols); f_array = NULL;
  free(pos); 		pos       = NULL;
  free(vel); 		vel       = NULL;
  free(pos_sub); 	pos_sub   = NULL;
  free(vel_sub); 	vel_sub   = NULL;
  free(P_array); 	P_array   = NULL;
  free(T_array); 	T_array   = NULL;
  free(a0_array);   a0_array  = NULL;
  free(msd_array);  msd_array = NULL;
  free(vfc_array);  vfc_array = NULL;

}

int run(int argc, char *argv[])
{
	run_later_tasks();

	return 0.;
}
























