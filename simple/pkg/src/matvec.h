/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating host populations.
*/

#ifndef NAN
#define NAN log(0)
#endif

#ifndef __MATVEC_H
#define __MATVEC_H

/*
   ==================
   === STRUCTURES ===
   ==================
*/


typedef struct{
	int *values, length;
} vec_int;


typedef struct{
	vec_int ** rows;
	int n;
} mat_int;




typedef struct{
	int length;
	double *values;
} vec_double;


typedef struct{
	vec_double ** rows;
	int n;
} mat_double;




/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/



vec_int * create_vec_int(int n);

mat_int * create_mat_int(int n);

vec_double * create_vec_double(int n);

mat_double * create_mat_double(int n);


/*
   ===================
   === DESTRUCTORS ===
   ===================
*/


void free_mat_int(mat_int *in);

void free_vec_int(vec_int *in);

void free_mat_double(mat_double *in);

void free_vec_double(vec_double *in);



/*
   ==========================
   === EXTERNAL FUNCTIONS ===
   ==========================
*/
int vecint_i(vec_int *in, int i);

int matint_ij(mat_int *in, int i, int j);

double vecdouble_i(vec_double *in, int i);

double matdouble_ij(mat_double *in, int i, int j);

void print_vec_int(vec_int *in);

void print_mat_int(mat_int *in);

void print_vec_double(vec_double *in);

void print_mat_double(mat_double *in);

void print_gsl_vector(gsl_vector *in, char format[256]);

int max_vec_int(vec_int *vec);

int min_vec_int(vec_int *vec);

void permut_vec_int(vec_int *in, gsl_rng * rng);

void sample_vec_int(vec_int *in, vec_int *out, bool replace, gsl_rng * rng);

void sort_vec_int(vec_int *in, vec_int *out, vec_int *idx);


#endif