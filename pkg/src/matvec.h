/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating host populations.
*/

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
	struct vec_int ** rows;
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



#endif
