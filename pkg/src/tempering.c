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



/* /\* MOVE TEMPERATURE *\/ */
/* void move_temperature(param *par, data *dat, mcmc_param *mcmcPar, gsl_rng *rng){ */
/*   double oldLogPost, newLogPost; */

/* /\* compute current logposterior value *\/ */
/*   double logPost = logposterior_all(dat, dnaInfo, spatialInfo, gen, par, rng); */
  
/*   /\* *\/ */
/* } */
