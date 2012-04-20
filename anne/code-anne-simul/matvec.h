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





/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/



vec_int * create_vec_int(int n);

mat_int * create_mat_int(int n);



/*
   ===================
   === DESTRUCTORS ===
   ===================
*/


void free_mat_int(mat_int *in);

void free_vec_int(vec_int *in);




/*
   ==========================
   === EXTERNAL FUNCTIONS ===
   ==========================
*/
int vecint_i(vec_int *in, int i);

int matint_ij(mat_int *in, int i, int j);

void print_vec_int(vec_int *in);

void print_mat_int(mat_int *in);

void print_gsl_vector(gsl_vector *in, char format[256]);

int max_vec_int(vec_int *vec);

int min_vec_int(vec_int *vec);

void permut_vec_int(vec_int *in, gsl_rng * rng);

void sample_vec_int(vec_int *in, vec_int *out, bool replace, gsl_rng * rng);




#endif
