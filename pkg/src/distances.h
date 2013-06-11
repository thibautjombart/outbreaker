/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating host populations.
*/
#ifndef __MATVEC_H
#include "matvec.h"
#endif

#ifndef __GENCLASSES_H
#include "genclasses.h"
#endif


#ifndef __DISTANCES_H
#define  __DISTANCES_H

/*
   ==================
   === STRUCTURES ===
   ==================
*/


typedef struct{
	mat_int *mutation1, *mutation2, *nbcommon;
	int n;
} dna_dist;



/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/

dna_dist * create_dna_dist(int n);



/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

void free_dna_dist(dna_dist * in);





/*
   =================
   === AUXILIARY ===
   =================
*/

bool is_atgc(char in);

int get_mutation1(dna_dist * in, int i, int j);

int get_mutation2(dna_dist * in, int i, int j);

int get_nbcommon(dna_dist * in, int i, int j);





/*
   ==========================
   === EXTERNAL FUNCTIONS ===
   ==========================
*/

void print_dna_dist(dna_dist *in);

dna_dist * compute_dna_distances(list_dnaseq *in);


#endif
