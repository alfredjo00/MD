/*
H1lattice.c
Program that arranges atoms on a fcc lattice. 
Created by Anders Lindman on 2013-03-15.
*/

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* Function takes a matrix of size [4*N*N*N][3] as input and stores a fcc lattice in it. N is the number of unit cells in each dimension and lattice_param is the lattice parameter. */
void init_fcc(double positions[][3], int N, double lattice_param)
{
    int i, j, k;
    int xor_value;
    
    for (i = 0; i < 2 * N; i++){
        for (j = 0; j < 2 * N; j++){
            for (k = 0; k < N; k++){
                if (j % 2 == i % 2 ){
                    xor_value = 0;
                }
                else {
                    xor_value = 1;
                }
                positions[i * N * 2 * N + j * N + k][0] = lattice_param * (0.5 * xor_value + k);
                positions[i * N * 2 * N + j * N + k][1] = lattice_param * (j * 0.5);
                positions[i * N * 2 * N + j * N + k][2] = lattice_param * (i * 0.5);
            }
        }
    }
}

void rand_deviations(double positions[][3], int n_atoms, double lattice_param, double dev)
{
	double lim_0 = lattice_param * (1.0 - dev);
	double lim_1 = lattice_param * (1.0 + dev);
		
	double C = (lim_1 - lim_0);
	double D = -(lim_1 - lim_0) / 2.0; 
	
	const gsl_rng_type *T;
	int seed = 42;

	gsl_rng *r;
	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	gsl_rng_set(r, seed); 	
	
	for (int i = 0; i < n_atoms; ++i){
		positions[i][0] += gsl_rng_uniform(r) * C + D;
		positions[i][1] += gsl_rng_uniform(r) * C + D;
		positions[i][2] += gsl_rng_uniform(r) * C + D;
	}
	
	gsl_rng_free(r);  	
}
