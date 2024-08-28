#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"
#include "verlet.h"

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
    E_potentials[i] = get_energy_AL(X, N * a0[i], n_atoms) / (N * N * N);
    volume_Al[i] = a0[i] * a0[i] * a0[i];
  }	
	
  double **write_array = malloc(2. * sizeof(double*));
	
  write_array[0] = volume_Al; 
  write_array[1] = E_potentials;
	
  write_to_file("data/volume_Al-energy_pot.txt",
		/**/		  write_array, n_points, 2);
	
  free(write_array); write_array = NULL;
}