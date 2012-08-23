#ifndef __LIKELIHOOD_H
#define __LIKELIHOOD_H


double loglikelihood_i(int i, data *dat, dna_dist *dnainfo, gentime *gen, param *par);

double loglikelihood_gen_i(int i, dna_dist *dnainfo, param *par);

double loglikelihood_all(data *dat, dna_dist *dnainfo, gentime *gen, param *par);

double loglikelihood_gen_all(data *dat, dna_dist *dnainfo, param *par);

double loglike_kappa_all(param *par);

double loglike_alpha_all(param *par);

double logposterior_all(data *dat, dna_dist *dnainfo, gentime *gen, param *par);

bool check_loglikelihood_all(data *dat, dna_dist *dnainfo, gentime *gen, param *par);

#endif
