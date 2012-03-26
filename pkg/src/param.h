/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), January 2012.
  Licence: GPL >=2.

*/


#ifndef __PARAM_H
#define __PARAM_H



/*
   =======================
   === DATA STRUCTURES ===
   =======================
*/

/* LEGEND OF THE PARAMETERS */
/*
  - nu1: rate of transitions
  - nu2: rate of transversions; nu2 = kappa * nu1 with kappa a constant
  - tau: offset for the coalescent time
  - alpha: proba that two sampled isolates belong to the same lineage, i.e. the oldest is the MRCA
  - weightNaGen: weight to be used when genetic likelihood can't be computed
*/

struct param{
	double nu1, nu2, kappa, tau, alpha, weightNaGen;
};




/*
   ====================
   === CONSTRUCTORS ===
   ====================
*/

struct param * create_param();





/*
   ===================
   === DESTRUCTORS ===
   ===================
*/

void free_param(struct param * in);


#endif
