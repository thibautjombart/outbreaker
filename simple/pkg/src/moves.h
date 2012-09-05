
#ifndef __MOVES_H
#define __MOVES_H



/*
  ===================
  AUXILIARY FUNCTIONS
  ===================
*/

int choose_kappa_i(int T, gentime *gen, gsl_rng *rng);

int choose_alpha_i(int i, data *dat, dna_dist *dnainfo, param *currentPar, mcmc_param *mcmcPar, gsl_rng *rng);



/*
  =====
  MOVES
  =====
*/

void move_mu1(param *currentPar, param *tempPar, data *dat, dna_dist *dnainfo, mcmc_param *mcmcPar, gsl_rng *rng);

void move_gamma(param *currentPar, param *tempPar, data *dat, dna_dist *dnainfo, mcmc_param *mcmcPar, gsl_rng *rng);

void move_pi(param *currentPar, param *tempPar, data *dat, mcmc_param *mcmcPar, gsl_rng *rng);

void move_phi(param *currentPar, param *tempPar, data *dat, mcmc_param *mcmcPar, gsl_rng *rng);

void move_Tinf(param *currentPar, param *tempPar, data *dat, dna_dist *dnainfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng);

void move_alpha(param *currentPar, param *tempPar, data *dat, dna_dist *dnainfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng);

void move_kappa(param *currentPar, param *tempPar, data *dat, dna_dist *dnainfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng);

/* void move_alpha_kappa(param *currentPar, param *tempPar, data *dat, dna_dist *dnainfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng); */


#endif
