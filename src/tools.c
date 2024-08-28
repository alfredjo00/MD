#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools.h"
#include <omp.h>


#define SQ(X) ((X) * (X))
void elementwise_addition(
		     double *res,
		     double *v1,
		     double *v2,
		     unsigned int len)
{
  for (int i = 0; i < len; ++i){
    res[i] = v1[i] + v2[i];    
  }
}

void elementwise_multiplication(
        double *res,
        double *v1,
        double *v2,
        unsigned int len)
{
  for (int i = 0; i < len; ++i){
    res[i] = v1[i] * v2[i];  
  }
}

double dot_product(
        double *v1,
        double *v2,
        unsigned int len)
{
  double result = 0.;

  for (int i = 0; i < len; ++i){
    result += v1[i] * v2[i];  
  }
  return result;
}

void create_2D_array(
        double ***array,
        unsigned int column_size,
        unsigned int row_size)
{
  *array = (double**) calloc(column_size, sizeof(double*));
  
  for (int i = 0; i < column_size; ++i){
    *(*array + i) = (double*) calloc(row_size, sizeof(double));
  }  
}

void destroy_2D_array(
        double **array,
        unsigned int column_size)
{
  for (int i = 0; i < column_size; ++i){
    free(array[i]);
  }
  free(array);
}

void matrix_multiplication(
        double **result,
        double **v1,
        double **v2,
        unsigned int m,
        unsigned int n)
{
  for (int i = 0; i < n; ++i){
    for (int j = 0; j < n; ++j){
      for (int k = 0; k < m; ++k){
        result[i][j] += v1[i][k] * v2[k][j];
      }    
    }  
  }
}

double vector_norm(
        double *v1,
        unsigned int len)
{
  double result = 0.;
  
  for (int i = 0; i < len; ++i){
    result += v1[i] * v1[i];  
  }
  return sqrt(result);
}


void normalize_vector(
        double *v1,
        unsigned int len)
{
  double norm = vector_norm(v1, len);
  
  for (int i = 0; i < len; ++i){
    v1[i] = v1[i] / norm;  
  }
}

double average(
        double *v1,
        unsigned int len)
{
  double sum = 0.;
  for (int i = 0; i < len; ++i){
    sum += v1[i];
  }
  return (sum / len);
}


double standard_deviation(
        double *v1,
        unsigned int len)
{
  double ave = average(v1, len);
  
  double *v_dev = malloc(len * sizeof(double));
  
  for (int i = 0; i < len; ++i){
    v_dev[i] = pow(v1[i] - ave , 2);
  }
  double variance =  average(v_dev, len);
  free(v_dev);
  return sqrt(variance);
}

double distance_between_vectors(
			 double *v1,
			 double *v2,
			 unsigned int len)
{
  double *v_del = malloc(len * sizeof(double));

  for (int i = 0; i < len; ++i){
    v_del[i] = v1[i] - v2[i];
  }
  double dist = vector_norm(v_del, len);
  free(v_del);
  return dist;
}

void arange(double *array, double start, int len_t, double dt){
    for(int i = 0; i < len_t; i++){
		array[i] = start + i*dt;
    }
} 

void write_to_file(char *fname, double **arr, int n_points, int n_cols)
{
	FILE *fp = fopen(fname, "w");
	for(int i = 0; i < n_points; ++i){
		for (int j = 0; j < n_cols; ++j){
      if (isnan(arr[j][i]))
        printf("%f\n", arr[j][i]);
			if (fprintf(fp, "%.5e\t\t", arr[j][i]) < 0 ) {
        perror("Error writing file.");
        fclose(fp);
        exit(EXIT_FAILURE);
      }
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void read_data(char *fname, double **array, int cols, int rows)
{
    FILE *f = fopen(fname, "r");
	int i, j;
	
    for(j=0; j<rows; j++)
        for(i=0; i<cols; i++)
            if(fscanf(f, "%lf", &array[i][j]) != 1)
                exit(1);
    fclose(f);
}

void array_ones(double *arr, int n)
{
	for (int i = 0; i < n; ++i){
		arr[i] = 1.;		
	}	
}

void print_array(double *arr, int n){
	for (int i = 0; i < n; ++i){
		printf("%f\n", arr[i]);		
	}	
}

double trapezoidal_int(double *f, double dt, int n)
{
	double sum = 0.;
	for (int i = 1; i < n - 1; ++i){
		sum += f[i];
	}
	return dt * (sum + (f[0] + f[n-1])/2.);
}

void linspace(double *arr, double a, double b, int n)
{
	double delta = (b - a)/(n - 1.);
	for (int i = 0; i < n; ++i){
		arr[i] = a + delta * i;		
	}	
}

double kinetic_energy(double vel[][3], int n_atoms, double m)
{
	double ke = 0.;
	
	for (int i = 0; i < n_atoms; ++i){		
		ke += (m / 2.) * ((vel[i][0] * vel[i][0]) + 
	/* */                 (vel[i][1] * vel[i][1]) + 
	/* */                 (vel[i][2] * vel[i][2]));
	}	
	return ke;
}

void velocity_correlation(double *vfc, double vel[][3], int n_times, int n_atoms) {
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < n_times; ++i) {
        for (int h = 0; h < n_times - i; ++h) {
            for (int j = 0; j < n_atoms; ++j) {
                vfc[i] += (vel[j + (h + i) * n_atoms][0] * vel[j + h * n_atoms][0] +
                           vel[j + (h + i) * n_atoms][1] * vel[j + h * n_atoms][1] +
                           vel[j + (h + i) * n_atoms][2] * vel[j + h * n_atoms][2]) / (n_times - h);
            }
        }
    }
  #pragma omp parallel for
  for (int i = 0; i < n_times; ++i) {
    vfc[i] /= n_atoms;
  }
}

void mean_squared_displacement(double *msd, double pos[][3], int n_times, int n_atoms) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < n_times; ++i) {
        for (int j = 0; j < n_atoms; ++j) {
            msd[i] += (SQ(pos[j + i * n_atoms][0] - pos[j][0]) +
                        SQ(pos[j + i * n_atoms][1] - pos[j][1]) +
                        SQ(pos[j + i * n_atoms][2] - pos[j][2])) / n_atoms;
        }
    }
}

void write_positions(char *fname, double pos[][3], int n_time, int n_atoms, double *timerange)
{
  double **f_array;
  create_2D_array(&f_array, 10, n_time);

  for (int i = 0; i < n_time; ++i){
    f_array[0][i] = timerange[i];
    f_array[1][i] = pos[0 +   i * n_atoms][0];
    f_array[2][i] = pos[0 +   i * n_atoms][1];
    f_array[3][i] = pos[0 +   i * n_atoms][2];
    f_array[4][i] = pos[49 +  i * n_atoms][0];
    f_array[5][i] = pos[49 +  i * n_atoms][1];
    f_array[6][i] = pos[49 +  i * n_atoms][2];
    f_array[7][i] = pos[127 + i * n_atoms][0];
    f_array[8][i] = pos[127 + i * n_atoms][1];
    f_array[9][i] = pos[127 + i * n_atoms][2];
  }

  write_to_file(fname, f_array, n_time, 10);
  destroy_2D_array(f_array, 10);
}
