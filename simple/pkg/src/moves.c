
#include "common.h"
#include "structures.h"
#include "matvec.h"
#include "genclasses.h"
#include "distances.h"
#include "init.h"
#include "prior.h"
#include "likelihood.h"
#include "moves.h"


/* MOVE VALUES OF MU1 */
void move_mu1(param *currentPar, param *tempPar, data *dat, dna_dist *dnainfo, mcmc_param *mcmcPar, gsl_rng *rng){
    double logRatio=0.0;

    /* GENERATE CANDIDATE VALUE FOR MU1 */
    do{
	tempPar->mu1 += gsl_ran_gaussian(rng, mcmcPar->sigma_mu1);
    } while(tempPar->mu1 < 0); /* avoid negative values */


    /* ACCEPT / REJECT */
    /* only likelihood as priors are flat for mu1 */
    /* compute only genetic part as the epi part is unchanged */
    logRatio += loglikelihood_gen_all(dat, dnainfo, tempPar);
    logRatio -= loglikelihood_gen_all(dat, dnainfo, currentPar);

    /* printf("\nCurrent mu1: %.10f   New mu1: %.10f\n", currentPar->mu1, tempPar->mu1); */
    /* printf("\nLL old: %.5f  LL new:%.5f  logRatio: %.5f\n", loglikelihood_gen_all(dat, dnainfo, currentPar), loglikelihood_gen_all(dat, dnainfo, tempPar), logRatio); */

    /* if p(new/old) > 1, accept new */
    if(logRatio>=0) {
	currentPar->mu1 = tempPar->mu1;
	mcmcPar->n_accept_mu1 += 1;
	/* printf("\nAccepting new value\n"); */
    } else { /* else accept new with proba (new/old) */
	if(log(gsl_rng_uniform(rng)) <= logRatio){ /* accept */
	    currentPar->mu1 = tempPar->mu1;
	    mcmcPar->n_accept_mu1 += 1;
	    /* printf("\nAccepting new value\n"); */
	} else { /* reject */
	    tempPar->mu1 = currentPar->mu1;
	    mcmcPar->n_reject_mu1 += 1;
	    /* printf("\nRejecting new value\n"); */
	}
    }

} /* end move_mu1 */







/* MOVE VALUES OF GAMMA */
void move_gamma(param *currentPar, param *tempPar, data *dat, dna_dist *dnainfo, mcmc_param *mcmcPar, gsl_rng *rng){
    double logRatio=0.0;

    /* GENERATE CANDIDATE VALUE FOR GAMMA */
    do{
	tempPar->gamma += gsl_ran_gaussian(rng, mcmcPar->sigma_gamma);
    } while(tempPar->gamma < 0); /* avoid negative values */
    

    /* ACCEPT / REJECT */
    /* compute only genetic part as the epi part is unchanged */
    logRatio += loglikelihood_gen_all(dat, dnainfo, tempPar);
    logRatio -= loglikelihood_gen_all(dat, dnainfo, currentPar);

    /* compute the priors */
    logRatio += logprior_gamma(tempPar);
    logRatio -= logprior_gamma(currentPar);


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





void move_Tinf(param *currentPar, param *tempPar, data *dat, dna_dist *dnainfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng){
    double logRatio=0.0;
    int i, toMove = 0;

    /* DETERMINE WHICH Tinf_i TO MOVE */
    sample_vec_int(mcmcPar->all_idx, mcmcPar->idx_move_Tinf, FALSE, rng);

    /* MOVE EACH Tinf_i IN TURN */
    for(i=0;i<mcmcPar->idx_move_Tinf->length;i++){
	toMove = vec_int_i(mcmcPar->idx_move_Tinf,i);

	/* move i-th Tinf */
	do{
	    tempPar->Tinf->values[toMove] += gsl_rng_uniform(rng) >= 0.5 ? 1 : -1; /* move : +/- 1 unit time */
	} while(vec_int_i(tempPar->Tinf,toMove) > vec_int_i(dat->dates,toMove)); /* ensure Tinf_i <= T_i*/


	/* ACCEPT/REJECT STEP */
	/* compute only genetic part as the epi part is unchanged */
	logRatio = loglikelihood_all(dat, dnainfo, gen, tempPar) - loglikelihood_all(dat, dnainfo, gen, currentPar);

	/* compute the priors */
	logRatio += logprior_gamma(tempPar) - logprior_gamma(currentPar);

	/* if p(new/old) > 1, accept new */
	if(logRatio>=0) {
	    currentPar->Tinf->values[toMove] = vec_int_i(tempPar->Tinf,toMove);
	    mcmcPar->n_accept_Tinf += 1;
	} else { /* else accept new with proba (new/old) */
	    if(log(gsl_rng_uniform(rng)) <= logRatio){ /* accept */
		currentPar->Tinf->values[toMove] = vec_int_i(tempPar->Tinf,toMove);
		mcmcPar->n_accept_Tinf += 1;
	    } else { /* reject */
		tempPar->Tinf->values[toMove] = vec_int_i(currentPar->Tinf,toMove);
		mcmcPar->n_reject_Tinf += 1;
	    }
	}

    }
} /* end move_Tinf*/





/* MOVE VALUES OF ALPHA AND KAPPA */
void move_alpha(param *old_par, param *new_par, data *dat, dna_dist *dnainfo, gentime *gen, gsl_rng *rng){
    /* DETERMINE WHICH alpha_i TO MOVE */


    /* find relevant candidates (we need T^inf_{alpha_i} < T^inf_i) */

    
    /* MOVE EACH alpha_i */

}



/* MOVE VALUES OF KAPPA */
void move_kappa(param *old_par, param *new_par, data *dat, dna_dist *dnainfo, gentime *gen, gsl_rng *rng){
    /* DETERMINE WHICH kappa_i TO MOVE */
    /* find relevant candidates (*/
    
    /* MOVE EACH kappa_i */

}




/*

CODE MODEL FROM ANNE 

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










/*
>>>> TESTING <<<<
*/

int main(){
  /* DECLARATIONS */
    int TIMESPAN;
    data *dat;
    gentime *gen;
    param *par, *tempPar;
    dna_dist * dnainfo;
    mcmc_param * mcmcPar;

    double logPrior, logLike, logPost;

    /* INITIALIZE RNG */
    gsl_rng *rng = create_gsl_rng(time(NULL));


    /* CONVERT DATA */
    dat = alloc_data(3,10);
    dat->dates->values[0] = 0;
    dat->dates->values[1] = 2;
    dat->dates->values[1] = 3;
    dat->dna->list[0]->seq[0] = 'a';
    dat->dna->list[1]->seq[0] = 'a';
    dat->dna->list[2]->seq[0] = 't';
    printf("\n>>> Data <<<\n");
    print_data(dat);


    /* GET TIME SPAN */
    TIMESPAN = max_vec_int(dat->dates) - min_vec_int(dat->dates);
    printf("\nTimespan is %d\n",TIMESPAN);


    /* CREATE AND INIT GENERATION TIME */
    gen = alloc_gentime(TIMESPAN, 5);
    init_gentime(gen, 1, 1.0, 0.0, 0.0);
    printf("\n>>> gentime info <<<\n");
    print_gentime(gen);


     /* CREATE AND INIT PARAMETERS */
    par = alloc_param(3);
    par->alpha->values[0] = -1;
    par->alpha->values[1] = 0;
    par->alpha->values[2] = 0;
    par->kappa->values[0] = 1;
    par->kappa->values[1] = 1;
    par->kappa->values[2] = 1;
    par->mu1 = 0.0001;
    par->gamma = 1.0;
    par->pi = 0.5;
    printf("\nParameters (par)\n");
    print_param(par);

    tempPar = alloc_param(3);
    copy_param(par,tempPar);
    printf("\nParameters (tempPar)\n");
    print_param(par);


    /* ALLOCATE MCMCPAR */
    mcmcPar = alloc_mcmc_param(3);
    init_mcmc_param(mcmcPar, dat);
    printf("\nMCMC parameters (mcmcPar)\n");
    print_mcmc_param(mcmcPar);


    /* COMPUTE GENETIC DISTANCES */
    dnainfo = compute_dna_distances(dat->dna);
    printf("\n>>> DNA info <<<\n");
    print_dna_dist(dnainfo);


    /* COMPUTE PRIORS */
    logPrior = logprior_all(par);
    printf("\nPrior value (log): %.10f\n", logPrior);

   /* COMPUTE LIKELIHOOD */
    logLike = loglikelihood_all(dat, dnainfo, gen, par);
    printf("\nLog-likelihood value: %.10f\n", logLike);

    /* COMPUTE POSTERIOR */
    logPost = logposterior_all(dat, dnainfo, gen, par);
    printf("\nLog-posterior value: %.10f\n", logPost);


    /* MOVE MU1 */
    int i;
    mcmcPar->sigma_mu1 = 0.001;
    for(i=0;i<500;i++){
	move_mu1(par, tempPar, dat, dnainfo, mcmcPar, rng);
	printf("\nmu1: %.10f (reject: %d  accept: %d  ratio: %.3f)", par->mu1, mcmcPar->n_reject_mu1, mcmcPar->n_accept_mu1, (double) mcmcPar->n_reject_mu1 / mcmcPar->n_accept_mu1);
    }
    printf("\n");
    fflush(stdout);

    /* MOVE GAMMA */
    for(i=0;i<500;i++){
	move_gamma(par, tempPar, dat, dnainfo, mcmcPar, rng);
	printf("\ngamma: %.10f (reject: %d  accept: %d  ratio: %.3f)", par->gamma, mcmcPar->n_reject_gamma, mcmcPar->n_accept_gamma, (double) mcmcPar->n_reject_gamma / mcmcPar->n_accept_gamma);
    }
    printf("\n");
    fflush(stdout);

    /* MOVE TINF */
    for(i=0;i<500;i++){
	move_Tinf(par, tempPar, dat, dnainfo, gen, mcmcPar, rng);
	printf("\nTinf:");
	print_vec_int(par->Tinf);
	printf(" (reject: %d  accept: %d  ratio: %.3f)", mcmcPar->n_reject_Tinf, mcmcPar->n_accept_Tinf, (double) mcmcPar->n_reject_Tinf / mcmcPar->n_accept_Tinf);
    }
    printf("\n");
    fflush(stdout);

    /* /\* RUNTIME TEST *\/ */
    /* int ITER=10e6, i; */
    /* time_t t1, t2; */
    /* time(&t1); */
    /* printf("\nRuntime (%d computations of posterior): \n", ITER); */
    /* for(i=0;i<ITER;i++){ */
    /* 	logPost = logposterior_all(dat, dnainfo, gen, par); */
    /* } */
    /* time(&t2); */
    /* printf("\nellapsed time: %d seconds\n", (int) t2 - (int) t1); */


    /* FREE / RETURN */
    gsl_rng_free(rng);
    free_data(dat);
    free_gentime(gen);
    free_dna_dist(dnainfo);
    free_param(par);
    free_param(tempPar);
    free_mcmc_param(mcmcPar);

    return 0;
}





/*
  gcc instructions

  gcc -o moves matvec.c genclasses.c structures.c init.c distances.c prior.c likelihood.c moves.c -lgsl -lgslcblas -Wall -g


  gcc -o moves matvec.c genclasses.c structures.c init.c distances.c prior.c likelihood.c moves.c -lgsl -lgslcblas -Wall -O3

 ./moves

  valgrind --leak-check=full -v moves

*/

