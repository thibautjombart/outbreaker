#ifndef __LIKELIHOOD_H
#define __LIKELIHOOD_H


double loglikelihood_i(int i, data *dat, dna_dist *dnainfo, gentime *gen, param *par);

double loglikelihood_gen_i(int i, dna_dist *dnainfo, param *par);

double loglikelihood_all(data *dat, dna_dist *dnainfo, gentime *gen, param *par);

double loglikelihood_gen_all(data *dat, dna_dist *dnainfo, param *par);

double logposterior_all(data *dat, dna_dist *dnainfo, gentime *gen, param *par);

void check_loglikelihood_all(data *dat, dna_dist *dnainfo, gentime *gen, param *par);

#endif
