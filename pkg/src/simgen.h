/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

  These functions are basic routines for simulating host populations.
*/

#ifndef __SIMGEN_H
#define __SIMGEN_H



/*
  ===================
  AUXILIARY FUNCTIONS
  ===================
*/

char transi(char in);

char transv(char in, gsl_rng *rng);




/*
   ==========================
   === EXTERNAL FUNCTIONS ===
   ==========================
*/

dnaseq *create_haplo(int length, gsl_rng *rng);

void evolve_haplo(dnaseq *in, double nu1, double nu2, double t, gsl_rng *rng);

dnaseq *replicate_haplo(dnaseq *in, double nu1, double nu2, double t, gsl_rng *rng);



#endif
