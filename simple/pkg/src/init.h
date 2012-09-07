#ifndef __COMMON_H
#include "common.h"
#endif

#ifndef __INIT_H
#define __INIT_H

gsl_rng * create_gsl_rng(time_t t);

void init_gentime(gentime *in, double *values);

int find_maxLike_kappa_i(int T, gentime *gen);

void init_param(param *par, data *dat,  gentime *gen, int *ances, double pi_param1, double pi_param2, double phi_param1, double phi_param2, double init_mu1, double init_gamma, gsl_rng *rng);

void init_mcmc_param(mcmc_param *in, data *dat, bool move_mut, int *move_alpha, bool move_kappa, bool move_Tinf, bool move_pi, bool move_phi, bool find_import, int burnin, int find_import_at);


/* void CalculIsInHosp(nb_data *, raw_data *); */

/* void InitAugData(parameters *, nb_data * , raw_data *, aug_data *); */

/* void InitMCMCSettings(mcmcInternals *); */

/* void InitParam(parameters *param); */

#endif

