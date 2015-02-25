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





/* COMPUTE TEMPERED LOG-POSTERIOR */
/* this one is only used when moving temperature */
/* temperature is otherwise embedded in parameter movement functions  */
double temper(double *logPost, int temperature, mcmc_param *mcmcPar){
  double out=0.0;
  /* checks */
  if(temperature<1) {
    warning("requested temperature (%d) <1 - setting it to one", temperature);
    temperature = 1;
  } else if(temperature>mcmcPar->max_temperature){
    warning("requested temperature (%d) > max_temperature (%d); setting it to %d", temperature, mcmcPar->max_temperature, mcmcPar->max_temperature);
    temperature = mcmcPar->max_temperature;
  }

  /* filter and return */
  out =  (*logPost)/temperature;
  filter_logprob(&out);
  return out;
}




/* COMPUTE TEMPERATURE PRIOR */
double logprior_temperature(int temperature, mcmc_param *mcmcPar){
  double out=0.0;
  out = gsl_ran_poisson_pdf((unsigned int) temperature, mcmcPar->prior_temperature);

  /* filter and return */
  filter_logprob(&out);
  return out;
}



/* MOVE TEMPERATURE */
void move_temperature(data *dat, dna_dist *dnaInfo, spatial_dist *spaInfo, gentime *gen, param *par, mcmc_param *mcmcPar, gsl_rng *rng){
  double logPost, logRatio=0.0, correcRatio=0.0;
  int newTemperature =  mcmcPar->current_temperature;

  /* PROPOSE NEW TEMPERATURE */
  if(mcmcPar->current_temperature==1) {
    newTemperature = 2;
    correcRatio += log(0.5);
  } else if(mcmcPar->current_temperature==mcmcPar->max_temperature){
    newTemperature=mcmcPar->max_temperature-1;
    correcRatio += log(0.5);
  } else {
    newTemperature += (gsl_rng_uniform(rng) >= 0.5 ? 1 : -1);
  }

  /* compute current logposterior value */
  logPost = logposterior_all(dat, dnaInfo, spaInfo, gen, par, rng);

  /* compute acceptance ratio */
  logRatio += logPost/newTemperature - logPost/mcmcPar->current_temperature +
    logprior_temperature(newTemperature, mcmcPar) - logprior_temperature(mcmcPar->current_temperature, mcmcPar);

  /* correction factor */
  logRatio += correcRatio;

  /* METROPOLIS-HASTING ACCEPTANCE */
  /* if p(new/old) > 1, accept new */
  if(logRatio>=0.0) {
    mcmcPar->current_temperature = newTemperature;
    mcmcPar->n_accept_temperature += 1;
  } else { /* else accept new with proba (new/old) */
    if(log(gsl_rng_uniform(rng)) <= logRatio){ /* accept */
      mcmcPar->current_temperature = newTemperature;
      mcmcPar->n_accept_temperature += 1;
    } else { /* reject */
      mcmcPar->n_reject_temperature += 1;
    }
  }
} /* end move_temperature */
