#include <stdlib.h>
#include <stdio.h>
#include "potential.h"
#include "tools.h"
#include <math.h>

#define k_b 8.617333262e-5
#define SQ(X) ((X) * (X))

void verlet_vel_pos(
		int n_timesteps, 
		int n_atoms,
		int skip,
		double positions[][3],
		double velocities[][3],
		double dt,
		double m,
		double cell_length)
{			
	
	double (*a)[3] = malloc(sizeof(double[n_atoms][3]));
	double (*u)[3] = malloc(sizeof(double[n_atoms][3]));
	double (*w)[3] = malloc(sizeof(double[n_atoms][3]));
	
	for (int i = 0; i < n_atoms; ++i){
			u[i][0] = positions[i][0];	
			u[i][1] = positions[i][1];	
			u[i][2] = positions[i][2];	
			
			w[i][0] = 0.;
			w[i][1] = 0.;
			w[i][2] = 0.;
	}
	
	// calc acc
	get_forces_AL(a, u, cell_length, n_atoms);	
	int h = 0;
    for (int i = 1; i < n_timesteps + 1; i++) {
        /* v(t+dt/2) */
        for (int j = 0; j < n_atoms; j++) {
			w[j][0] += dt * 0.5 * a[j][0] / m;
			w[j][1] += dt * 0.5 * a[j][1] / m;
			w[j][2] += dt * 0.5 * a[j][2] / m;
        }
        
        /* q(t+dt) */
        for (int j = 0; j < n_atoms; j++) {
			u[j][0] += dt * w[j][0];
			u[j][1] += dt * w[j][1];
			u[j][2] += dt * w[j][2];
        }
        
        /* a(t+dt) */
		get_forces_AL(a, u, cell_length, n_atoms);	
        
        /* v(t+dt) */
        for (int j = 0; j < n_atoms; j++) {
			w[j][0] += dt * 0.5 * a[j][0] / m;
			w[j][1] += dt * 0.5 * a[j][1] / m;
			w[j][2] += dt * 0.5 * a[j][2] / m;
        }
		if (i % skip == 0){
			for (int j = 0; j < n_atoms; j++) {
				positions[j + h * n_atoms][0] = u[j][0];
				positions[j + h * n_atoms][1] = u[j][1];
				positions[j + h * n_atoms][2] = u[j][2];
							
				velocities[j + h * n_atoms][0] = w[j][0];
				velocities[j + h * n_atoms][1] = w[j][1];
				velocities[j + h * n_atoms][2] = w[j][2];
			}
			h++;
		}
	}	
	free(a); a = NULL;
	free(u); u = NULL;
	free(w); w = NULL;	
}


double instant_temperature(double v[][3], double m, int n_atoms)
{
	double kinetic_energy = 0.;	
	
	for (int i = 0; i < n_atoms; ++i){
		kinetic_energy += (m / 2.) * ( SQ(v[i][0]) +
	    /* */						   SQ(v[i][1]) +
	    /* */						   SQ(v[i][2]) );
	}

	return kinetic_energy * 2. / (3. * n_atoms * k_b);
}

double instant_pressure(
		double v[][3], 
		double virial, 
		double cell_V,
		double m,
		int n_atoms)
{
	double T = instant_temperature(v, m, n_atoms);
	return (n_atoms * k_b * T + virial)/cell_V;
}
/*
 		int n_timesteps => Total number of time steps
		int n_atoms     => Number of atoms
		int skip        => Only save every (int) timestep
		int n_cells     => Number of cells
		double positions[][3],  => Position matrix to store positions
		double velocities[][3], => Store velocities
		double dt       => time step
		double m        => mass
		double a0       => Lattice constant
		double T_eq     => Equilibration constant Temperature
		double P_eq     => Equilibration constant Pressure
		double *Pressure     => Store presssure time series
		double *Temperature  => Store temperature time series
		double *a0_array     => Store lattice constant time series

 */
void verlet_equilibration_scaling(
		int n_timesteps, 
		int n_atoms,
		int skip,
		int n_cells,
		double positions[][3],
		double velocities[][3],
		double dt,
		double m,
		double a0,
		double T_eq,
		double P_eq,
		double *Pressure,
		double *Temperature,
		double *a0_array)
{			
	double T_inst, P_inst, alpha_T, alpha_P;
	double virial;
	double (*F)[3] 		= malloc(sizeof(double[n_atoms][3]));
	double (*u)[3] 		= malloc(sizeof(double[n_atoms][3]));
	double (*w)[3] 		= malloc(sizeof(double[n_atoms][3]));
	int N_eq = 10000;
	int a0_N = 1500;
	double V;
	
	double tau_T = 200 * dt;
	double tau_P = 200 * dt;
	/*
	kappa_t = 0.01385 GPa-1
	1 GPa-1 = 160.21766208 Å^3 eV^-1
	*/
	double kappa_t = 0.0161290 * 160.21766208;
	
	for (int i = 0; i < n_atoms; ++i){
		u[i][0] 	= positions[i][0];
		u[i][1] 	= positions[i][1];
		u[i][2] 	= positions[i][2];

		w[i][0] 	= 0.;
		w[i][1] 	= 0.;
		w[i][2] 	= 0.;
	}
	

	// calc acc
	get_forces_AL(F, u, n_cells * a0, n_atoms);	
	int h = 0;
	for (int i = 1; i < n_timesteps + 1; i++) {
        /* v(t+dt/2) */
        for (int j = 0; j < n_atoms; j++) {
			w[j][0] += dt * 0.5 * F[j][0] / m;
			w[j][1] += dt * 0.5 * F[j][1] / m;
			w[j][2] += dt * 0.5 * F[j][2] / m;
        }
        
        /* q(t+dt) */
        for (int j = 0; j < n_atoms; j++) {
			u[j][0] += dt * w[j][0];
			u[j][1] += dt * w[j][1];
			u[j][2] += dt * w[j][2];
        }
        
        /* a(t+dt) */
		get_forces_AL(F, u, n_cells * a0, n_atoms);	
        
        /* v(t+dt) */
        for (int j = 0; j < n_atoms; j++) {
			w[j][0] += dt * 0.5 * F[j][0] / m;
			w[j][1] += dt * 0.5 * F[j][1] / m;
			w[j][2] += dt * 0.5 * F[j][2] / m;
        }

		V = (n_cells * a0) * (n_cells * a0) * (n_cells * a0);
		virial = get_virial_AL(u, n_cells * a0, n_atoms);
		T_inst = instant_temperature(w, m, n_atoms);
		P_inst = instant_pressure(w, virial, V, m, n_atoms);
		
		if (i <= N_eq){

			alpha_T = 1 + (2. * dt / (tau_T)) * ((T_eq - T_inst) / T_inst);
			alpha_P = 1 - kappa_t * (dt / tau_P) * (P_eq - P_inst);

			for (int j = 0; j < n_atoms; j++) {
				w[j][0] *= sqrt(alpha_T);
				w[j][1] *= sqrt(alpha_T);
				w[j][2] *= sqrt(alpha_T);

				u[j][0] *= cbrt(alpha_P);
				u[j][1] *= cbrt(alpha_P);
				u[j][2] *= cbrt(alpha_P);
			}
			a0 *= cbrt(alpha_P);

			if (i == N_eq){
				double a0_mean = 0.;
				for (int u = a0_N; u < N_eq; ++u){
					a0_mean += a0_array[u] / (N_eq - a0_N);
				}
				a0 = a0_mean;
			}
		}


		if (i % skip == 0){
			// Printing to see progression.
			printf("T %.3f \t P %.3f \t x0 %.4f \t %d / %d\n",
					T_inst - 273.15,
					P_inst * 160.21766208 / 0.0001,
					u[0][0],
					i/skip,
					n_timesteps/skip);
			// Saving data to arrays
			for (int j = 0; j < n_atoms; j++) {
				positions[j + h * n_atoms][0] 	= u[j][0];
				positions[j + h * n_atoms][1] 	= u[j][1];
				positions[j + h * n_atoms][2] 	= u[j][2];
							
				velocities[j + h * n_atoms][0] 	= w[j][0];
				velocities[j + h * n_atoms][1] 	= w[j][1];
				velocities[j + h * n_atoms][2] 	= w[j][2];
				
				Temperature[h] 	= T_inst;
				Pressure[h] 	= P_inst;
				a0_array[h]		= a0;
			}
			h++;
		}
	}	
	free(F); 		F = NULL;
	free(u); 		u = NULL;
	free(w); 		w = NULL;
}
