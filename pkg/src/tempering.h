
#ifndef __TEMPERING_H
#define __TEMPERING_H


/* COMPUTE TEMPERED LOG-POSTERIOR */
double temper(double *logPost, int temperature, mcmc_param *mcmcPar);

/* COMPUTE TEMPERATURE PRIOR */
double logprior_temperature(int temperature, mcmc_param *mcmcPar);

/* MOVE TEMPERATURE */
void move_temperature(data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, param *par, mcmc_param *mcmcPar, gsl_rng *rng);



#endif
