#ifndef __LIKELIHOOD_H
#define __LIKELIHOOD_H



/*
  ===================
  AUXILIARY FUNCTIONS
  ===================
*/
int find_sequenced_ancestor(int i, data *dat, dna_dist *dnainfo, param *par);

int mutation1_ij(int i, int j, data *dat, dna_dist *dnainfo);

int mutation2_ij(int i, int j, data *dat, dna_dist *dnainfo);

int com_nucl_ij(int i, int j, data *dat, dna_dist *dnainfo);

double gsl_ran_poisson_pdf_fixed(unsigned int k, double mu);


/*
  ====================
  LIKELIHOOD FUNCTIONS
  ====================
*/

double loglikelihood_i(int i, data *dat, dna_dist *dnainfo, gentime *gen, param *par, gsl_rng *rng);

double loglikelihood_gen_i(int i, data *dat, dna_dist *dnainfo, param *par, gsl_rng *rng);

double loglikelihood_all(data *dat, dna_dist *dnainfo, gentime *gen, param *par, gsl_rng *rng);

double loglikelihood_gen_all(data *dat, dna_dist *dnainfo, param *par, gsl_rng *rng);

double loglike_kappa_all(param *par);

/* double loglike_alpha_all(param *par); */

double logposterior_all(data *dat, dna_dist *dnainfo, gentime *gen, param *par, gsl_rng *rng);

double sim_loglike_gen(data *dat, param *par, gsl_rng *rng);

bool check_loglikelihood_all(data *dat, dna_dist *dnainfo, gentime *gen, param *par, gsl_rng *rng);

#endif
