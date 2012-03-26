/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating host populations.
*/

#ifndef __GENLIKE_H
#define __GENLIKE_H


/*
   =================
   === AUXILIARY ===
   =================
*/





/*
   ==========================
   === EXTERNAL FUNCTIONS ===
   ==========================
*/
double compute_genlike(int i, int j, double ti, double tj, double nu1, double nu2, double alpha, double tau, struct dna_dist *dnainfo, struct param *par);


#endif
