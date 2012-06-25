
#ifndef __MCMC_H
#define __MCMC_H

void fprint_param(FILE *file, param *par, int step, bool quiet);

void mcmc(int nIter, int outEvery, char outputFile[256], bool quiet, param *par, data *dat, dna_dist *dnainfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng);

#endif
