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
 * Computes the instantaneous temperature
 */
double get_T(int n_atoms, double kinetic_energy)
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
	
	write_to_file("data/task2_0.txt", f_array, time_its, fcols);
	
	double average_ke = 0.;

	for (int i = 500; i < n_timesteps; ++i){
		average_ke += f_array[3][i] / (n_timesteps - 500);
	}

	double T = get_T(n_atoms, average_ke);
	printf("Kinetic energy %f eV\nTemperature %f C\n", average_ke, T - 273.15);
	
	destroy_2D_array(f_array, fcols); f_array = NULL;
	free(pos); 		pos = NULL;
	free(vel); 		vel = NULL;
	free(pos_sub); 	pos_sub = NULL;
	free(vel_sub); 	vel_sub = NULL;
}