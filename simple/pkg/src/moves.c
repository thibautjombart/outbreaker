
#include "common.h"
#include "common.h"
#include "structures.h"
#include "matvec.h"
#include "genclasses.h"
#include "init.h"
#include "prior.h"
#include "likelihood.h"



/* MOVE VALUES OF MU1 */
void move_mu1(param *currentPar, param *tempPar, data *dat, dna_dist *dnainfo, gentime *gen, param *par, mcmc_param *mcmcPar, gsl_rng *rng){
    double logRatio=0.0;

    /* GENERATE CANDIDATE VALUE FOR MU1 */
    while(tempPar->mu1 < 0){ /* avoid negative values */
	tempPar->mu1 += gsl_ran_gaussian(rng, mcmcPar->sigma_mu1);
    }


    /* ACCEPT / REJECT */
    /* only likelihood as priors are flat for mu1 */
    /* compute only genetic part as the epi part is unchanged */
    logRatio += loglikelihood_gen_all(dat, dnainfo, tempPar);
    logRatio -= loglikelihood_gen_all(dat, dnainfo, currentpar);

    /* if p(new/old) > 1, accept new */
    if(logRatio>=0) {
	currentPar->mu1 = tempPar->mu1;
	mcmcPar->n_accept_mu1 += 1;
    } else { /* else accept new with proba (new/old) */
	if(log(gsl_rng_uniform(rng)) <= logRatio){ /* accept */
	    currentPar->mu1 = tempPar->mu1;
	    mcmcPar->n_accept_mu1 += 1;
	} else { /* reject */
	    tempPar->mu1 = currentPar->mu1;
	    mcmcPar->n_reject_mu1 += 1;
	}
    }

} /* move_mu1 */







/* MOVE VALUES OF GAMMA */
void move_gamma(param *currentPar, param *tempPar, data *dat, dna_dist *dnainfo, gentime *gen, param *par, mcmc_param *mcmcPar, gsl_rng *rng){
    double logRatio=0.0;

    /* GENERATE CANDIDATE VALUE FOR GAMMA */
    while(tempPar->gamma < 0){ /* avoid negative values */
	tempPar->gamma += gsl_ran_gaussian(rng, mcmcPar->sigma_gamma);
    }

    /* ACCEPT / REJECT */
    /* compute only genetic part as the epi part is unchanged */
    logRatio += loglikelihood_gen_all(dat, dnainfo, tempPar);
    logRatio -= loglikelihood_gen_all(dat, dnainfo, currentpar);

    /* compute the priors */
    logRatio += logprior_gamma(tempPar);
    logRatio -= logprior_gamma(currentpar);


    /* if p(new/old) > 1, accept new */
    if(logRatio>=0) {
	currentPar->gamma = tempPar->gamma;
	mcmcPar->n_accept_gamma += 1;
    } else { /* else accept new with proba (new/old) */
	if(log(gsl_rng_uniform(rng)) <= logRatio){ /* accept */
	    currentPar->gamma = tempPar->gamma;
	    mcmcPar->n_accept_gamma += 1;
	} else { /* reject */
	    tempPar->gamma = currentPar->gamma;
	    mcmcPar->n_reject_gamma += 1;
	}
    }

} /* end move_gamma*/






/* move values of kappa */
void move_kappa(param *old_par, param *new_par, data *dat, dna_dist *dnainfo, gentime *gen, param *par, gsl_rng *rng){
    /* DETERMINE THE NUMBER OF kappa_i VALUES TO MOVE */

    /* DETERMINE WHICH kappa_i TO MOVE */

    /* MOVE EACH kappa_i */

}




/*

MODEL CODE FROM ANNE 

*/
/* void movePi(parameters * curParam, parameters * newParam, raw_data * data, aug_data *augData, acceptance *accept){ */
/*     /\* updating Pi with a Gibbs sampler *\/ */
/*     double ColonisedAtFirstAdmission=0; */
/*     double NonColonisedAtFirstAdmission=0; */
/*     int i; */

/*     for(i=0;i<data->NbPatients;i++){ */
/* 	if(augData->C[i]<gsl_vector_get(data->A[i],0)){ */
/* 	    ColonisedAtFirstAdmission++; */
/* 	} else { */
/* 	    NonColonisedAtFirstAdmission++; */
/* 	} */
/*     } */

/*     curParam->Pi = gsl_ran_beta (data->rng, ColonisedAtFirstAdmission+1, NonColonisedAtFirstAdmission+1); */
/* } */






/* void moveNu(mcmcInternals * MCMCSettings, parameters * curParam, parameters *newParam, raw_data * data, nb_data *nb, aug_data *augData, dna_dist *dnainfo, acceptance *accept, NbProposals *NbProp){ */
/*     /\* updating nu1 with Metropolis algorithm *\/ */
/*     /\* proposal distribution = lognormal *\/ */
/*     double *newVal, *curVal; */
/*     double * nbAccept; */
/*     double * nbPropos; */
/*     double sigmaProp; */
/*     double r,z; */
/*     double pAccept = 0; */
/*     /\* double (* logPrior) (parameters *); *\/ */

/*     curVal = &curParam->nu1; */
/*     newVal = &newParam->nu1; */
/*     sigmaProp = MCMCSettings->Sigma_nu1; */
/*     nbAccept = &accept->PourcAcc_nu1; */
/*     nbPropos = &NbProp->NbProp_nu1; */

/*     *(newVal) = *(curVal)*gsl_ran_lognormal(data->rng,0,sigmaProp); */

/*     pAccept += Colon(data, nb, augData, dnainfo, newParam); */
/*     pAccept -= Colon(data, nb, augData, dnainfo, curParam); */

/*     pAccept +=  logpriorNu1(newParam) - logpriorNu1(curParam); */

/*     pAccept +=  log(*(newVal)) - log(*(curVal)); /\* correction for lognormal *\/ */

/*     if (pAccept>0) r=0; else r=pAccept; /\* r=log(min(1,ratioOfMHalgo))*\/ */
/*     z=gsl_rng_uniform(data->rng); */
/*     if (log(z)<=r) { */
/* 	*curVal = *newVal; */
/* 	*nbAccept +=1; */
/*     }else{ */
/* 	*newVal = *curVal; */
/*     } */
/*     *nbPropos+=1; */

/* } */








