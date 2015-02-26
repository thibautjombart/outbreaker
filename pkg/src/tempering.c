#include "common.h"
#include "structures.h"
#include "matvec.h"
#include "genclasses.h"
#include "distances.h"
#include "init.h"
#include "prior.h"
#include "likelihood.h"
#include "moves.h"
#include "tempering.h"


/* note: the function 'temper' has to be in moves.c to avoid circular referencing */




/* COMPUTE TEMPERATURE PRIOR */
double logprior_temperature(int temperature, mcmc_param *mcmcPar){
  double out=0.0;
  out = log(gsl_ran_poisson_pdf((unsigned int) temperature, (double) mcmcPar->prior_temperature));

  /* filter and return */
  filter_logprob(&out);
  return out;
}



/* MOVE TEMPERATURE */
void move_temperature(data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, param *par, mcmc_param *mcmcPar, gsl_rng *rng){
  double logPost=0.0, logRatio=0.0, correcRatio=0.0;
  int newTemperature =  mcmcPar->current_temperature;

  Rprintf("\n\n-- temperature movement --");
  /* PROPOSE NEW TEMPERATURE */
  if(mcmcPar->current_temperature==1) {
    newTemperature = 2;
    correcRatio += log(0.5);
    Rprintf("\nAutomatic proposal 1->2");
  } else if(mcmcPar->current_temperature==mcmcPar->max_temperature){
    newTemperature=mcmcPar->max_temperature-1;
    correcRatio += log(0.5);
    Rprintf("\nAutomatic proposal %d->%d", mcmcPar->max_temperature, mcmcPar->max_temperature-1);
  } else {
    newTemperature += (gsl_rng_uniform(rng) >= 0.5 ? 1 : -1);
    Rprintf("\nRandom proposal %d->%d",mcmcPar->current_temperature, newTemperature);
  }

  /* compute current logposterior value */
  logPost = logposterior_all(dat, dnaInfo, spaInfo, gen, par, rng);
  Rprintf("\nlogPost=%.5f    proba=%.5f", logPost, exp(logPost));

  /* compute acceptance ratio */
  Rprintf("\nj->i: hj(x)=%.5f  hi(x)=%.5f  p(j)=%.5f  p(i)=%.5f", temper(&logPost, newTemperature), temper(&logPost, mcmcPar->current_temperature), logprior_temperature(newTemperature, mcmcPar), logprior_temperature(mcmcPar->current_temperature, mcmcPar) );

  logRatio = temper(&logPost, newTemperature) - temper(&logPost, mcmcPar->current_temperature) +
    logprior_temperature(newTemperature, mcmcPar) - logprior_temperature(mcmcPar->current_temperature, mcmcPar);

  Rprintf("\nacceptance log-ratio: %.5f   ratio: %.5f", logRatio, exp(logRatio));
  Rprintf("\ncorrection log-ratio: %.5f   ratio: %.5f", correcRatio, exp(correcRatio));

  /* correction factor */
  logRatio += correcRatio;
  Rprintf("\nfinal acceptance log-ratio: %.5f   ratio: %.5f", logRatio, exp(logRatio));

  /* METROPOLIS-HASTING ACCEPTANCE */
  /* if p(new/old) > 1, accept new */
  if(logRatio>=0.0) {
    Rprintf("\nACCEPTED %d->%d", mcmcPar->current_temperature, newTemperature);
    mcmcPar->current_temperature = newTemperature;
    mcmcPar->n_accept_temperature += 1;
  } else { /* else accept new with proba (new/old) */
    if(log(gsl_rng_uniform(rng)) <= logRatio){ /* accept */
      Rprintf("\nACCEPTED %d->%d", mcmcPar->current_temperature, newTemperature);
      mcmcPar->current_temperature = newTemperature;
      mcmcPar->n_accept_temperature += 1;
    } else { /* reject */
      Rprintf("\nREJECTED %d->%d", mcmcPar->current_temperature, newTemperature);
      mcmcPar->n_reject_temperature += 1;
    }
  }
} /* end move_temperature */
