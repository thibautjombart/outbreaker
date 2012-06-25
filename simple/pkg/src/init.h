#ifndef __COMMON_H
#include "common.h"
#endif

#ifndef __INIT_H
#define __INIT_H

gsl_rng * create_gsl_rng(time_t t);

void init_gentime(gentime *in, int type, double param1, double param2, double param3);

void init_param(param *par, data *dat,  gentime *gen, int *ances);

void init_mcmc_param(mcmc_param *in, data *dat);


/* void CalculIsInHosp(nb_data *, raw_data *); */

/* void InitAugData(parameters *, nb_data * , raw_data *, aug_data *); */

/* void InitMCMCSettings(mcmcInternals *); */

/* void InitParam(parameters *param); */

#endif

