/* **********************************************
*
* Add v1 and v2 elementwise
* results is stored in res
*
* res should be properly initialised to zero
* for this function to work correctly
*
* **********************************************/
void elementwise_addition(
	double *res,
	double *v1,
	double *v2,
	unsigned int len);

/* **********************************************
*
* Multiply v1 and v2 elementwise
* results is stored in res
*
* res should be properly initialised to zero
* for this function to work correctly
*
* **********************************************/
void elementwise_multiplication(
	double *res,
	double *v1,
	double *v2,
	unsigned int len);

/* **********************************************
*
* Calculate the dot product between
* v1 and v2
*
* the result is returned as a double
*
* **********************************************/
double dot_product(
	double *v1,
	double *v2,
	unsigned int len);

/* **********************************************
*
* Allocate the memory to a 2D array
*
* **********************************************/
void create_2D_array(
	double ***array,
	unsigned int column_size,
	unsigned int row_size);

void destroy_2D_array(
	double **array,
	unsigned int column_size);

void matrix_multiplication(
	double **result,
	double **v1,
	double **v2,
	unsigned int m,
	unsigned int n);

double vector_norm(
	double *v1,
	unsigned int len);

void normalize_vector(
	double *v1,
	unsigned int len);

double average(
	double *v1,
	unsigned int len);

double standard_deviation(
	double *v1,
	unsigned int len);

double distance_between_vectors(
	double *v1,
	double *v2,
	unsigned int len);

void arange(
	double *array,
	double start, 
	int len_t, 
	double dt);

void write_to_file(
	char *fname, 
	double **arr,
	int n_points,
	int n_cols);

void read_data(
	char *fname,
	double **array,
	int cols,
	int rows);

void array_ones(
	double *arr, 
	int n);

void print_array(
	double *arr,
	int n);


double trapezoidal_int(
	double *f,
	double dt,
	int n);
	

void linspace(
	double *arr, 
	double a, 
	double b, 
	int n);
	
	
double kinetic_energy(
	double vel[][3], 
	int n_atoms,
	double m);


void velocity_correlation(
	double *vfc,
	double vel[][3],
	int n_times,
	int n_atoms);


void mean_squared_displacement(
	double *msd,
	double pos[][3],
	int n_times,
	int n_atoms);

void write_positions(
	char *fname, 
	double pos[][3],
	int n_time, 
	int n_atoms,
	double *timerange);
