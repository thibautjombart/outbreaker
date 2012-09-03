
#ifndef __MCMC_H
#define __MCMC_H



/*
  ===================
  AUXILIARY FUNCTIONS
  ===================
*/

void fprint_chains(FILE *file, data *dat, dna_dist *dnainfo, gentime *gen, param *par, int step, gsl_rng *rng, bool quiet);

void fprint_mcmc_param(FILE *file, mcmc_param *mcmcPar, int step);




/*
   ================
   TUNING FUNCTIONS
   ================
*/

void tune_mu1(mcmc_param * in, gsl_rng *rng);

void tune_gamma(mcmc_param * in, gsl_rng *rng);

/* void tune_Tinf(mcmc_param * in, gsl_rng *rng); */




/*
   ===============================================
   METROPOLIS-HASTING ALGORITHM FOR ALL PARAMETERS
   ===============================================
*/

void mcmc(int nIter, int outEvery, char outputFile[256], char mcmcOutputFile[256], int tuneEvery, bool quiet, param *par, data *dat, dna_dist *dnainfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng);

#endif
