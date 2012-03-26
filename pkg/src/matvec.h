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


struct vec_int{
	int *values, length;
};


struct mat_int{
	struct vec_int ** rows;
	int n;
};





/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/



struct vec_int * create_vec_int(int n);

struct mat_int * create_mat_int(int n);

struct vec_int * create_vec_int_zero(int n);



/*
   ===================
   === DESTRUCTORS ===
   ===================
*/


void free_mat_int(struct mat_int *in);

void free_vec_int(struct vec_int *in);




/*
   ==========================
   === EXTERNAL FUNCTIONS ===
   ==========================
*/

void print_vec_int(struct vec_int *in);

void print_mat_int(struct mat_int *in);



#endif
