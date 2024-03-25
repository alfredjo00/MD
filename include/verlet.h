#pragma once

void verlet_vel_pos(
		int n_timesteps, 
		int n_atoms,
		int skip,
		double positions[][3],
		double velocities[][3],
		double dt,
		double m,
		double a0);		

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
		double *a0_array);
