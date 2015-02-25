
#ifndef __TEMPERING_H
#define __TEMPERING_H


/* COMPUTE TEMPERED LOG-POSTERIOR */
double temper(double *logPost, int temperature, mcmc_param *mcmcPar);

/* COMPUTE TEMPERATURE PRIOR */
double logprior_temperature(int temperature, mcmc_param *mcmcPar);


#endif
