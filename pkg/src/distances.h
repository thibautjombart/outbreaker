/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating host populations.
*/

#ifndef __DISTANCES_H
#define  __DISTANCES_H

/*
   ==================
   === STRUCTURES ===
   ==================
*/


struct dna_dist{
	struct mat_int *transi, *transv, *nbcommon;
	int n;
};



/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/

struct dna_dist * create_dna_dist(int n);



/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

void free_dna_dist(struct dna_dist * in);





/*
   =================
   === AUXILIARY ===
   =================
*/

bool is_atgc(char in);

int get_transi(struct dna_dist * in, int i, int j);

int get_transv(struct dna_dist * in, int i, int j);

int get_nbcommon(struct dna_dist * in, int i, int j);





/*
   ==========================
   === EXTERNAL FUNCTIONS ===
   ==========================
*/

void print_dna_dist(struct dna_dist *in);

struct dna_dist * compute_dna_distances(struct list_dnaseq *in);


#endif
