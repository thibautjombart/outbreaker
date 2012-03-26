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
double genlike_ij(int *s_i, int *s_j, double *t_i, double *t_j, int m_i, int m_j, double nu1, double nu2, double alpha, double tau, struct dna_dist *dnainfo, struct param *par);


#endif
