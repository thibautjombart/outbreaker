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





/* compute tempered log-posterior */
double temper(double *logPost, int temperature, mcmc_param *mcmcPar){
  if(temperature<1) {
    warning("requested temperature (%d) <1 - setting it to one", temperature);
    temperature = 1;
  } else if(temperature>mcmcPar->max_temperature){
    warning("requested temperature (%d) > max_temperature (%d); setting it to %d", temperature, mcmcPar->max_temperature, mcmcPar->max_temperature);
    temperature = mcmcPar->max_temperature;
  }
  
  return (*logPost)/temperature;
}
