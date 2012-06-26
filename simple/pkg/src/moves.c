
#include "common.h"
#include "structures.h"
#include "matvec.h"
#include "genclasses.h"
#include "distances.h"
#include "init.h"
#include "prior.h"
#include "likelihood.h"
#include "moves.h"




/*
  ===================
  AUXILIARY FUNCTIONS
  ===================
*/

/* FIND MOST LIKELY KAPPA_I (given time to infection 'T') */
int find_maxLike_kappa_i(int T, gentime *gen){
    int i, out=1;
    double temp=0.0, currentMax=0.0;

    for(i=1;i<gen->maxK;i++){
	temp = gentime_dens(gen, T, i);
	if(currentMax < temp) {
	    currentMax = temp;
	    out = i+1;
	}
    }

    return out;
} /* end find_maxLike_kappa_i */







/*
  =====
  MOVES
  =====
*/

/* MOVE VALUES OF MU1 */
void move_mu1(param *currentPar, param *tempPar, data *dat, dna_dist *dnainfo, mcmc_param *mcmcPar, gsl_rng *rng){
    double logRatio=0.0;

    /* GENERATE CANDIDATE VALUE FOR MU1 */
    /* do{ */
    /* 	tempPar->mu1 += gsl_ran_gaussian(rng, mcmcPar->sigma_mu1); */
    /* } while(tempPar->mu1 < 0.0); /\* avoid negative values *\/ */
    tempPar->mu1 += gsl_ran_gaussian(rng, mcmcPar->sigma_mu1);
    if(tempPar->mu1 < 0.0) tempPar->mu1 = 0.0;


    /* ACCEPT / REJECT */
    /* only likelihood as priors are flat for mu1 */
    /* compute only genetic part as the epi part is unchanged */
    logRatio += loglikelihood_gen_all(dat, dnainfo, tempPar);
    logRatio -= loglikelihood_gen_all(dat, dnainfo, currentPar);

    /* if p(new/old) > 1, accept new */
    if(logRatio>=0.0) {
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







/* /\* MOVE VALUES OF PI *\/ */
/* void move_pi(param *currentPar, param *tempPar, data *dat, dna_dist *dnainfo, mcmc_param *mcmcPar, gsl_rng *rng){ */
/*     double logRatio=0.0; */

/*     /\* GENERATE CANDIDATE VALUE FOR PI *\/ */
/*     /\* do{ *\/ */
/*     /\* 	tempPar->pi += gsl_ran_gaussian(rng, mcmcPar->sigma_pi); *\/ */
/*     /\* } while(tempPar->pi < 0.0); /\\* avoid negative values *\\/ *\/ */
/*     tempPar->pi += gsl_ran_gaussian(rng, mcmcPar->sigma_pi); */
/*     if(tempPar->pi < 0.0) { */
/* 	tempPar->pi = 0.0; */
/*     } else if(tempPar->pi > 1.0) tempPar->pi = 1.0; */


/*     /\* ACCEPT / REJECT *\/ */
/*     logRatio += loglikelihood_all(dat, dnainfo, tempPar); */
/*     logRatio -= loglikelihood_all(dat, dnainfo, currentPar); */
/*     logRatio += logprior_pi(tempPar); */
/*     logRatio -= logprior_pi(currentPar); */


/*     /\* if p(new/old) > 1, accept new *\/ */
/*     if(logRatio>=0.0) { */
/* 	currentPar->pi = tempPar->pi; */
/* 	mcmcPar->n_accept_pi += 1; */
/* 	/\* printf("\nAccepting new value\n"); *\/ */
/*     } else { /\* else accept new with proba (new/old) *\/ */
/* 	if(log(gsl_rng_uniform(rng)) <= logRatio){ /\* accept *\/ */
/* 	    currentPar->pi = tempPar->pi; */
/* 	    mcmcPar->n_accept_pi += 1; */
/* 	    /\* printf("\nAccepting new value\n"); *\/ */
/* 	} else { /\* reject *\/ */
/* 	    tempPar->pi = currentPar->pi; */
/* 	    mcmcPar->n_reject_pi += 1; */
/* 	    /\* printf("\nRejecting new value\n"); *\/ */
/* 	} */
/*     } */

/* } /\* end move_pi *\/ */







/* MOVE VALUES OF GAMMA */
void move_gamma(param *currentPar, param *tempPar, data *dat, dna_dist *dnainfo, mcmc_param *mcmcPar, gsl_rng *rng){
    double logRatio=0.0;

    /* GENERATE CANDIDATE VALUE FOR GAMMA */
    /* do{ */
    /* 	tempPar->gamma += gsl_ran_gaussian(rng, mcmcPar->sigma_gamma); */
    /* } while(tempPar->gamma < 0); /\* avoid negative values *\/ */
    tempPar->gamma += gsl_ran_gaussian(rng, mcmcPar->sigma_gamma);
    if(tempPar->gamma < 0.0) tempPar->gamma = 0.0;


    /* ACCEPT / REJECT */
    /* compute only genetic part as the epi part is unchanged */
    logRatio += loglikelihood_gen_all(dat, dnainfo, tempPar);
    logRatio -= loglikelihood_gen_all(dat, dnainfo, currentPar);

    /* compute the priors */
    logRatio += logprior_gamma(tempPar);
    logRatio -= logprior_gamma(currentPar);


    /* if p(new/old) > 1, accept new */
    if(logRatio>=0.0) {
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
	/* do{ */
	/*     tempPar->Tinf->values[toMove] += gsl_rng_uniform(rng) >= 0.5 ? 1 : -1; /\* move : +/- 1 unit time *\/ */
	/* } while(vec_int_i(tempPar->Tinf,toMove) > vec_int_i(dat->dates,toMove)); /\* ensure Tinf_i <= T_i*\/ */
	tempPar->Tinf->values[toMove] += gsl_rng_uniform(rng) >= 0.5 ? 1 : -1;
	if(vec_int_i(tempPar->Tinf,toMove) > vec_int_i(dat->dates,toMove)) tempPar->Tinf->values[toMove] = vec_int_i(dat->dates,toMove);

	/* ACCEPT/REJECT STEP */
	/* compute the likelihood */
	logRatio = loglikelihood_all(dat, dnainfo, gen, tempPar) - loglikelihood_all(dat, dnainfo, gen, currentPar);

	/* compute the priors */
	logRatio += logprior_gamma(tempPar) - logprior_gamma(currentPar);

	/* if p(new/old) > 1, accept new */
	if(logRatio>=0.0) {
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
void move_alpha_kappa(param *currentPar, param *tempPar, data *dat, dna_dist *dnainfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng){
    int i, j, oldestDate = 0, nCandidates=0, toMove=0, T, ances;
    double logRatio = 0.0;


    /* DETERMINE WHICH alpha_i/kappa_i TO MOVE */
    /*  - we change same kappa_i as alpha_i - */
    sample_vec_int(mcmcPar->all_idx, mcmcPar->idx_move_alpha, FALSE, rng);

    /* FIND OLDEST INFECTION (it has no ancestor) */
    oldestDate = min_vec_int(currentPar->Tinf);

    /* MOVE EACH alpha_i TO MOVE */
    for(i=0;i<mcmcPar->idx_move_alpha->length;i++){
	toMove = vec_int_i(mcmcPar->idx_move_alpha,i);

	/* move only isolates with an ancestor in the sample */
	if(vec_int_i(currentPar->Tinf, toMove) > oldestDate){
	    /* find candidate ancestors ('alpha_i' so that T^inf_{alpha_i} < T^inf_i) */
	    nCandidates=0;
	    for(j=0;j<dat->n;j++){
		if(vec_int_i(currentPar->Tinf,j) < vec_int_i(currentPar->Tinf,toMove))
		    mcmcPar->candid_ances->values[nCandidates++] =  j;
	    }

	    /* check that there is at least a candidate or issue error */
	    if(nCandidates==0){
		fprintf(stderr, "\n[in: moves.c->move_alpha_kappa]\nNo candidate ancestor for 'i' but still trying to move 'alpha_i' (i: %d).\n", toMove);
		exit(1);
	    }

	    /* GET PROPOSED ALPHA_I */
	    tempPar->alpha->values[toMove] = vec_int_i(mcmcPar->candid_ances, gsl_rng_uniform_int(rng, nCandidates));

	    /* /\* GET PROPOSED (MOST PROBABLE GIVEN W) KAPPA_I *\/ */
	    /* ances: proposed ancestor for 'toMove' */
	    ances = vec_int_i(tempPar->alpha,toMove);

	    /* time between ancestor and descendent */
	    T = vec_int_i(tempPar->Tinf, toMove) - vec_int_i(tempPar->Tinf, ances);

	    /* most likely value of kappa */
	    /* tempPar->kappa->values[toMove] = find_maxLike_kappa_i(T, gen); */


	    /* ACCEPT/REJECT STEP */
	    /* compute the likelihood */
	    logRatio = loglikelihood_all(dat, dnainfo, gen, tempPar) - loglikelihood_all(dat, dnainfo, gen, currentPar);

	    /* compute the priors */
	    logRatio += logprior_gamma(tempPar) - logprior_gamma(currentPar);

	    /* if p(new/old) > 1, accept new */
	    if(logRatio>=0.0) {
		currentPar->alpha->values[toMove] = vec_int_i(tempPar->alpha,toMove);
		currentPar->kappa->values[toMove] = vec_int_i(tempPar->kappa,toMove);
		mcmcPar->n_accept_alpha += 1;
		mcmcPar->n_accept_kappa += 1;
	    } else { /* else accept new with proba (new/old) */
		if(log(gsl_rng_uniform(rng)) <= logRatio){ /* accept */
		    currentPar->alpha->values[toMove] = vec_int_i(tempPar->alpha,toMove);
		    currentPar->kappa->values[toMove] = vec_int_i(tempPar->kappa,toMove);
		    mcmcPar->n_accept_alpha += 1;
		    mcmcPar->n_accept_kappa += 1;
		} else { /* reject */
		    tempPar->alpha->values[toMove] = vec_int_i(currentPar->alpha,toMove);
		    tempPar->kappa->values[toMove] = vec_int_i(currentPar->kappa,toMove);
		    mcmcPar->n_reject_alpha += 1;
		    mcmcPar->n_reject_kappa += 1;
		}
	    } /* end  ACCEPT/REJECT STEP */

	} /* end if isolate is not the oldest one */
    } /* end for loop (for all 'i' to move) */
} /* end move_alpha */







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

/* int main(){ */
/*   /\* DECLARATIONS *\/ */
/*     int TIMESPAN, i; */
/*     data *dat; */
/*     gentime *gen; */
/*     param *par, *tempPar; */
/*     dna_dist * dnainfo; */
/*     mcmc_param * mcmcPar; */

/*     double logPrior, logLike, logPost; */

/*     /\* INITIALIZE RNG *\/ */
/*     gsl_rng *rng = create_gsl_rng(time(NULL)); */


/*     /\* CONVERT DATA *\/ */
/*     dat = alloc_data(3,10); */
/*     dat->dates->values[0] = 0; */
/*     dat->dates->values[1] = 2; */
/*     dat->dates->values[1] = 3; */
/*     dat->dna->list[0]->seq[0] = 'a'; */
/*     dat->dna->list[1]->seq[0] = 'a'; */
/*     dat->dna->list[2]->seq[0] = 't'; */
/*     printf("\n>>> Data <<<\n"); */
/*     print_data(dat); */


/*     /\* GET TIME SPAN *\/ */
/*     TIMESPAN = max_vec_int(dat->dates) - min_vec_int(dat->dates); */
/*     printf("\nTimespan is %d\n",TIMESPAN); */


/*     /\* CREATE AND INIT GENERATION TIME *\/ */
/*     gen = alloc_gentime(TIMESPAN, 5); */
/*     init_gentime(gen, 1, 1.0, 0.0, 0.0); */
/*     printf("\n>>> gentime info <<<\n"); */
/*     print_gentime(gen); */
/*     printf("sizes of rows in gen: "); */
/*     for(i=0;i<gen->dens->n;i++) printf("%d ", gen->dens->rows[i]->length); */


/*      /\* CREATE AND INIT PARAMETERS *\/ */
/*     par = alloc_param(3); */
/*     par->alpha->values[0] = -1; */
/*     par->alpha->values[1] = 0; */
/*     par->alpha->values[2] = 0; */
/*     par->kappa->values[0] = 1; */
/*     par->kappa->values[1] = 1; */
/*     par->kappa->values[2] = 1; */
/*     par->mu1 = 0.0001; */
/*     par->gamma = 1.0; */
/*     par->pi = 0.5; */
/*     printf("\nParameters (par)\n"); */
/*     print_param(par); */

/*     tempPar = alloc_param(3); */
/*     copy_param(par,tempPar); */
/*     printf("\nParameters (tempPar)\n"); */
/*     print_param(par); */


/*     /\* ALLOCATE MCMCPAR *\/ */
/*     mcmcPar = alloc_mcmc_param(3); */
/*     init_mcmc_param(mcmcPar, dat); */
/*     printf("\nMCMC parameters (mcmcPar)\n"); */
/*     print_mcmc_param(mcmcPar); */


/*     /\* COMPUTE GENETIC DISTANCES *\/ */
/*     dnainfo = compute_dna_distances(dat->dna); */
/*     printf("\n>>> DNA info <<<\n"); */
/*     print_dna_dist(dnainfo); */


/*     /\* COMPUTE PRIORS *\/ */
/*     logPrior = logprior_all(par); */
/*     printf("\nPrior value (log): %.10f\n", logPrior); */

/*    /\* COMPUTE LIKELIHOOD *\/ */
/*     logLike = loglikelihood_all(dat, dnainfo, gen, par); */
/*     printf("\nLog-likelihood value: %.10f\n", logLike); */

/*     /\* COMPUTE POSTERIOR *\/ */
/*     logPost = logposterior_all(dat, dnainfo, gen, par); */
/*     printf("\nLog-posterior value: %.10f\n", logPost); */


/*     /\* MOVE MU1 *\/ */
/*     mcmcPar->sigma_mu1 = 0.001; */
/*     for(i=0;i<500;i++){ */
/* 	move_mu1(par, tempPar, dat, dnainfo, mcmcPar, rng); */
/* 	printf("\nmu1: %.10f (reject: %d  accept: %d  ratio: %.3f)", par->mu1, mcmcPar->n_reject_mu1, mcmcPar->n_accept_mu1, (double) mcmcPar->n_reject_mu1 / mcmcPar->n_accept_mu1); */
/*     } */
/*     printf("\n"); */
/*     fflush(stdout); */

/*     /\* MOVE GAMMA *\/ */
/*     for(i=0;i<500;i++){ */
/* 	move_gamma(par, tempPar, dat, dnainfo, mcmcPar, rng); */
/* 	printf("\ngamma: %.10f (reject: %d  accept: %d  ratio: %.3f)", par->gamma, mcmcPar->n_reject_gamma, mcmcPar->n_accept_gamma, (double) mcmcPar->n_reject_gamma / mcmcPar->n_accept_gamma); */
/*     } */
/*     printf("\n"); */
/*     fflush(stdout); */

/*     /\* MOVE TINF *\/ */
/*     for(i=0;i<500;i++){ */
/* 	move_Tinf(par, tempPar, dat, dnainfo, gen, mcmcPar, rng); */
/* 	printf("\nTinf:"); */
/* 	print_vec_int(par->Tinf); */
/* 	printf(" (reject: %d  accept: %d  ratio: %.3f)", mcmcPar->n_reject_Tinf, mcmcPar->n_accept_Tinf, (double) mcmcPar->n_reject_Tinf / mcmcPar->n_accept_Tinf); */
/*     } */
/*     printf("\n"); */
/*     fflush(stdout); */



/*     /\* MOVE ALPHA AND KAPPA *\/ */
/*     for(i=0;i<500;i++){ */
/* 	move_alpha_kappa(par, tempPar, dat, dnainfo, gen, mcmcPar, rng); */
/* 	printf("\nAlpha:"); */
/* 	print_vec_int(par->alpha); */
/* 	printf("\nKappa:"); */
/* 	print_vec_int(par->kappa); */
/* 	printf(" (reject: %d  accept: %d  ratio: %.3f)", mcmcPar->n_reject_alpha, mcmcPar->n_accept_alpha, (double) mcmcPar->n_reject_alpha / mcmcPar->n_accept_alpha); */
/*     } */
/*     printf("\n"); */
/*     fflush(stdout); */


/*     /\* /\\* RUNTIME TEST *\\/ *\/ */
/*     /\* int ITER=10e6, i; *\/ */
/*     /\* time_t t1, t2; *\/ */
/*     /\* time(&t1); *\/ */
/*     /\* printf("\nRuntime (%d computations of posterior): \n", ITER); *\/ */
/*     /\* for(i=0;i<ITER;i++){ *\/ */
/*     /\* 	logPost = logposterior_all(dat, dnainfo, gen, par); *\/ */
/*     /\* } *\/ */
/*     /\* time(&t2); *\/ */
/*     /\* printf("\nellapsed time: %d seconds\n", (int) t2 - (int) t1); *\/ */


/*     /\* FREE / RETURN *\/ */
/*     gsl_rng_free(rng); */
/*     free_data(dat); */
/*     free_gentime(gen); */
/*     free_dna_dist(dnainfo); */
/*     free_param(par); */
/*     free_param(tempPar); */
/*     free_mcmc_param(mcmcPar); */

/*     return 0; */
/* } */





/*
  gcc instructions

  gcc -o moves matvec.c genclasses.c structures.c init.c distances.c prior.c likelihood.c moves.c -lgsl -lgslcblas -Wall -g


  gcc -o moves matvec.c genclasses.c structures.c init.c distances.c prior.c likelihood.c moves.c -lgsl -lgslcblas -Wall -O3

 ./moves

  valgrind --leak-check=full -v moves

*/

