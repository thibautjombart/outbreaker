
#ifndef __TEMPERING_H
#define __TEMPERING_H


/* compute tempered log-posterior */
double temper(double *logPost, int temperature, mcmc_param *mcmcPar);




#endif
