
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

/* /\* FIND MOST LIKELY KAPPA_I (given time to infection 'T') *\/ */
/* int find_maxLike_kappa_i(int T, gentime *gen){ */
/*     int i, out=1; */
/*     double temp=0.0, currentMax=0.0; */

/*     for(i=1;i<gen->maxK;i++){ */
/* 	temp = gentime_dens(gen, T, i); */
/* 	if(currentMax < temp) { */
/* 	    currentMax = temp; */
/* 	    out = i+1; */
/* 	} */
/*     } */

/*     return out; */
/* } /\* end find_maxLike_kappa_i *\/ */





/* SAMPLE KAPPA_I USING PROB */
/* sample a value of kappa_i according to respective proba */
/* of kappa_1, kappa_2, ..., kappa_maxK */
int choose_kappa_i(int T, gentime *gen, gsl_rng *rng){
    int i;
    double probVec[gen->maxK], sumKappa=0.0, cumSum=0.0, rnd;

    /* STANDARDIZE P(KAPPA_I) TO GET PROBA */
    for(i=1;i<gen->maxK;i++){
	probVec[i-1] = gentime_dens(gen, T, i);
	sumKappa += probVec[i-1];
    }

    for(i=0;i<gen->maxK;i++){
	probVec[i] = probVec[i]/sumKappa;
    }


    /* CHOOSE KAPPA_I WITH PROBA PROBVEC[I] */
    rnd = gsl_rng_uniform(rng);
    i=0;
    do{
	cumSum += probVec[i++];
    } while(rnd>=cumSum);

    if(i>gen->maxK){
	fprintf(stderr, "\n[in: moves.c->choose_kappa_i]\nInvalid value of kappa_i returned (%d, but maxK=%d).\n", i, gen->maxK);
	exit(1);
    }

    return i;
} /* end choose_kappa_i */







/*
  =====
  MOVES
  =====
*/

/* MOVE VALUES OF MU1 */
void move_mu1(param *currentPar, param *tempPar, data *dat, dna_dist *dnainfo, mcmc_param *mcmcPar, gsl_rng *rng){
    double logRatio=0.0, QCur, QTemp;

    /* GENERATE CANDIDATE VALUE FOR MU1 ACCORDING TO LOGNORMAL DISTRIBUTION */
    do{
	tempPar->mu1 = gsl_ran_lognormal(rng,log(currentPar->mu1),mcmcPar->sigma_mu1);
    } while(tempPar->mu1>1.0);
    /* which should be the same as:
    tempPar->mu1 = currentPar->mu1 * gsl_ran_lognormal(rng,0,mcmcPar->sigma_mu1); */

    /* other possibility for proposal */
    /* tempPar->mu1 += gsl_ran_gaussian(rng, mcmcPar->sigma_mu1); */
    /* if(tempPar->mu1 < 0.0) tempPar->mu1 = 0.0;*/

    /* ACCEPT / REJECT */
    /* only likelihood as priors are flat for mu1 */
    /* compute only genetic part as the epi part is unchanged */
    logRatio += loglikelihood_gen_all(dat, dnainfo, tempPar);
    logRatio -= loglikelihood_gen_all(dat, dnainfo, currentPar);

    /* ADD CORRECTION FOR MH truncated lognormal */
    QCur = gsl_cdf_gaussian_P(-log(currentPar->mu1),mcmcPar->sigma_mu1);
    QTemp = gsl_cdf_gaussian_P(-log(tempPar->mu1),mcmcPar->sigma_mu1);
    logRatio +=  log(tempPar->mu1) - log(currentPar->mu1); /* correction for lognormal */
    logRatio +=   log(QCur) - log(QTemp); /* correction for truncation (no values >1) */


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





/* MOVE VALUES OF GAMMA */
void move_gamma(param *currentPar, param *tempPar, data *dat, dna_dist *dnainfo, mcmc_param *mcmcPar, gsl_rng *rng){
    double logRatio=0.0;

    /* GENERATE CANDIDATE VALUE FOR GAMMA ACCORDING TO LOGNORMAL DISTRIBUTION */
    tempPar->gamma = gsl_ran_lognormal(rng,log(currentPar->gamma),mcmcPar->sigma_gamma);
    /* which should be the same as:
    tempPar->gamma = currentPar->gamma * gsl_ran_lognormal(rng,0,mcmcPar->sigma_gamma); */

      /* other possibility for proposal */
      /* tempPar->gamma += gsl_ran_gaussian(rng, mcmcPar->sigma_gamma);
    if(tempPar->gamma < 0.0) tempPar->gamma = 0.0;*/

    /* ACCEPT / REJECT */
    /* compute only genetic part as the epi part is unchanged */
    logRatio += loglikelihood_gen_all(dat, dnainfo, tempPar);
    logRatio -= loglikelihood_gen_all(dat, dnainfo, currentPar);

   /* add correction (MH) for lognormal proposal */
    logRatio += log(tempPar->gamma) - log(currentPar->gamma);

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
    double logRatio=0.0, ll1, ll2;
    int i, toMove = 0;

    /* DETERMINE WHICH Tinf_i TO MOVE */
    sample_vec_int(mcmcPar->all_idx, mcmcPar->idx_move_Tinf, FALSE, rng);

    /* MOVE EACH Tinf_i IN TURN */
    for(i=0;i<mcmcPar->idx_move_Tinf->length;i++){
	toMove = vec_int_i(mcmcPar->idx_move_Tinf,i);

	/* move i-th Tinf */
	tempPar->Tinf->values[toMove] += (gsl_rng_uniform(rng) >= 0.5 ? 1 : -1);

	/* MAY NEED TO CHANGE THIS AND ADD CORRECTION */
	/* constraint: Tinf_i <= t_i */
	if(vec_int_i(tempPar->Tinf,toMove) > vec_int_i(dat->dates,toMove)) tempPar->Tinf->values[toMove] = vec_int_i(dat->dates,toMove);
	/* constraint: Tinf_i >= -trunc */
	if(vec_int_i(tempPar->Tinf,toMove) < -gen->trunc) tempPar->Tinf->values[toMove] = -gen->trunc;

	/* PROCEED TO ACCEPT/REJECT ONLY IF TINF HAS CHANGED */
	if(vec_int_i(tempPar->Tinf,toMove) != vec_int_i(currentPar->Tinf,toMove)){
	    /* ACCEPT/REJECT STEP */
	    /* compute the likelihood (no priors for Tinf) */
	    logRatio = loglikelihood_all(dat, dnainfo, gen, tempPar) - loglikelihood_all(dat, dnainfo, gen, currentPar);

	    ll1 = loglikelihood_all(dat, dnainfo, gen, currentPar);
	    ll2 = loglikelihood_all(dat, dnainfo, gen, tempPar);

	    /* if p(new/old) > 1, accept new */
	    if(logRatio>=0.0) {
		/* printf("\nTinf_%d: accepting automatically move from %d to %d (respective loglike:%f and %f)\n",toMove+1, vec_int_i(currentPar->Tinf,toMove), vec_int_i(tempPar->Tinf,toMove), ll1, ll2); */
		/* fflush(stdout); */

		currentPar->Tinf->values[toMove] = vec_int_i(tempPar->Tinf,toMove);
		mcmcPar->n_accept_Tinf += 1;
	    } else { /* else accept new with proba (new/old) */
		if(log(gsl_rng_uniform(rng)) <= logRatio){ /* accept */
		/*     printf("\nTinf_%d: accepting move from %d to %d (respective loglike:%f and %f)\n",toMove+1, vec_int_i(currentPar->Tinf,toMove), vec_int_i(tempPar->Tinf,toMove), ll1, ll2); */
		/* fflush(stdout); */

		    currentPar->Tinf->values[toMove] = vec_int_i(tempPar->Tinf,toMove);
		    mcmcPar->n_accept_Tinf += 1;
		} else { /* reject */
		    tempPar->Tinf->values[toMove] = vec_int_i(currentPar->Tinf,toMove);
		    mcmcPar->n_reject_Tinf += 1;
		}
	    }
	} /* end if Tinf has changed */
    } /* end for each indiv to move */
} /* end move_Tinf*/






/* MOVE VALUES OF ALPHA */
void move_alpha(param *currentPar, param *tempPar, data *dat, dna_dist *dnainfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng){
    int i, j, oldestDate = 0, nCandidates=0, toMove=0, T;
    double logRatio = 0.0, ll1, ll2;


    /* DETERMINE WHICH alpha_i TO MOVE */
    sample_vec_int(mcmcPar->all_idx, mcmcPar->idx_move_alpha, FALSE, rng);

    /* FIND OLDEST INFECTION (it has no ancestor) */
    oldestDate = min_vec_int(currentPar->Tinf);

    /* MOVE EACH alpha_i TO MOVE */
    /* !!! need to propose a new kappa at the same time */
    for(i=0;i<mcmcPar->idx_move_alpha->length;i++){
	toMove = vec_int_i(mcmcPar->idx_move_alpha,i);

	/* move only isolates with an ancestor in the sample */
	if(vec_int_i(currentPar->Tinf, toMove) > oldestDate){
	    /* find candidate ancestors ('alpha_i' so that T^inf_{alpha_i} < T^inf_i) */
	    /* nCandidates=0; */
	    nCandidates=1;
	    mcmcPar->candid_ances->values[0] = -1;
	    for(j=0;j<dat->n;j++){
		if(vec_int_i(currentPar->Tinf,j) < vec_int_i(currentPar->Tinf,toMove))
		    mcmcPar->candid_ances->values[nCandidates++] =  j;
	    }

	    /* check that there is at least a candidate or issue error */
	    if(nCandidates==0){
		fprintf(stderr, "\n[in: moves.c->move_alpha]\nNo candidate ancestor for 'i' but still trying to move 'alpha_i' (i: %d).\n", toMove);
		exit(1);
	    }


	    /* GET PROPOSED ALPHA_I */
	    tempPar->alpha->values[toMove] = vec_int_i(mcmcPar->candid_ances, gsl_rng_uniform_int(rng, nCandidates));

	    /* PROCEED ONLY IF ALPHA HAS CHANGED */
	    if(vec_int_i(tempPar->alpha,toMove) != vec_int_i(currentPar->alpha,toMove)){

		if(vec_int_i(tempPar->alpha,toMove)>=0){
		    /* MOVE KAPPA_I - MULTINOMIAL */
		    /* T: Tinf_i - Tinf_ances */
		    T = vec_int_i(tempPar->Tinf,toMove) - vec_int_i(tempPar->Tinf, tempPar->alpha->values[toMove]);
		    tempPar->kappa->values[toMove] = choose_kappa_i(T, gen, rng);
		} else {
		    tempPar->kappa->values[toMove] = 1;
		}

		/* ACCEPT/REJECT STEP */
		/* compute the likelihood */
		logRatio = loglikelihood_all(dat, dnainfo, gen, tempPar) - loglikelihood_all(dat, dnainfo, gen, currentPar);
		/* debugging */
		ll1=loglikelihood_all(dat, dnainfo, gen, currentPar);
		ll2=loglikelihood_all(dat, dnainfo, gen, tempPar);

		/* /\* compute the priors *\/ */
		/* logRatio += logprior_kappa_i(toMove,tempPar) - logprior_kappa_i(toMove,currentPar); */

		/* if p(new/old) > 1, accept new */
		if(logRatio>=0.0) {
		    /* /\* debugging *\/ */
		    /* printf("\naccepting automatically move from %d->%d to %d->%d (respective loglike:%f and %f)\n",vec_int_i(currentPar->alpha,toMove), toMove+1, vec_int_i(tempPar->alpha,toMove), toMove+1, ll1, ll2); */
		    /* fflush(stdout); */

		    currentPar->alpha->values[toMove] = vec_int_i(tempPar->alpha,toMove);
		    currentPar->kappa->values[toMove] = vec_int_i(tempPar->kappa,toMove);
		    mcmcPar->n_accept_alpha += 1;
		} else { /* else accept new with proba (new/old) */
		    if(log(gsl_rng_uniform(rng)) <= logRatio){ /* accept */
			/* /\* debugging *\/ */
			/* printf("\naccepting move from %d->%d to %d->%d (respective loglike:%f and %f)\n",vec_int_i(currentPar->alpha,toMove), toMove+1, vec_int_i(tempPar->alpha,toMove), toMove+1, ll1, ll2); */
			/* fflush(stdout); */

			currentPar->alpha->values[toMove] = vec_int_i(tempPar->alpha,toMove);
			currentPar->kappa->values[toMove] = vec_int_i(tempPar->kappa,toMove);
			mcmcPar->n_accept_alpha += 1;
		    } else { /* reject */
			tempPar->alpha->values[toMove] = vec_int_i(currentPar->alpha,toMove);
			tempPar->kappa->values[toMove] = vec_int_i(currentPar->kappa,toMove);
			mcmcPar->n_reject_alpha += 1;
		    }
		} /* end  ACCEPT/REJECT STEP */
	    } /* end if ancestor has changed */
	} /* end if isolate is not the oldest one */
    } /* end for loop (for all 'i' to move) */
} /* end move_alpha */






/* MOVE VALUES OF KAPPA */
void move_kappa(param *currentPar, param *tempPar, data *dat, dna_dist *dnainfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng){
    int i, j, oldestDate = 0, toMove=0, temp;
    double logRatio = 0.0;


    /* DETERMINE WHICH kappa_i TO MOVE */
    sample_vec_int(mcmcPar->all_idx, mcmcPar->idx_move_kappa, FALSE, rng);

    /* FIND OLDEST INFECTION (it has no ancestor) */
    oldestDate = min_vec_int(currentPar->Tinf);

    /* MOVE EACH kappa_i TO MOVE */
    for(i=0;i<mcmcPar->idx_move_kappa->length;i++){
	toMove = vec_int_i(mcmcPar->idx_move_kappa,i);

	/* move only isolates with an ancestor in the sample */
	/* if(vec_int_i(currentPar->Tinf, toMove) > oldestDate){ */
	if(vec_int_i(currentPar->alpha, toMove) >= 0){
	    /* GET PROPOSED KAPPA_I */
	    /* needs to be on [1;maxK]*/
	    temp = tempPar->kappa->values[toMove] + (gsl_rng_uniform(rng) >= 0.5 ? 1 : -1);
	    if(temp < 1) {
		temp = 1;
	    } else if(temp>gen->maxK){
		temp = gen->maxK;
	    }
	    tempPar->kappa->values[toMove] = temp;

	    /* PROCEED TO ACCEPT/REJECT STEP ONLY IF ALPHA HAS CHANGED */
	    if(vec_int_i(tempPar->kappa,toMove) != vec_int_i(currentPar->kappa,toMove)){
		/* ACCEPT/REJECT STEP */
		/* compute the likelihood */
		logRatio = loglikelihood_all(dat, dnainfo, gen, tempPar) - loglikelihood_all(dat, dnainfo, gen, currentPar);

		/* /\* compute the priors *\/ */
		/* logRatio += logprior_kappa_i(toMove,tempPar) - logprior_kappa_i(toMove,currentPar); */

		/* if p(new/old) > 1, accept new */
		if(logRatio>=0.0) {
		    currentPar->kappa->values[toMove] = vec_int_i(tempPar->kappa,toMove);
		    mcmcPar->n_accept_kappa += 1;
		} else { /* else accept new with proba (new/old) */
		    if(log(gsl_rng_uniform(rng)) <= logRatio){ /* accept */
			currentPar->kappa->values[toMove] = vec_int_i(tempPar->kappa,toMove);
			mcmcPar->n_accept_kappa += 1;
		    } else { /* reject */
			tempPar->kappa->values[toMove] = vec_int_i(currentPar->kappa,toMove);
			mcmcPar->n_reject_kappa += 1;
		    }
		} /* end  ACCEPT/REJECT STEP */
	    }
	} /* end if isolate is not the oldest one */
    } /* end for loop (for all 'i' to move) */
} /* end move_kappa */










/* MOVE VALUES OF PI */
void move_pi(param *currentPar, param *tempPar, data *dat, mcmc_param *mcmcPar, gsl_rng *rng){
    double logRatio=0.0;
    double QCur, QTemp;

    /* GENERATE CANDIDATE VALUE FOR PI */
    /* HERE REPLACE WITH TRUNCATED LOGNORMAL (no values >1) )*/
    do
    {
    	tempPar->pi = gsl_ran_lognormal(rng,log(currentPar->pi),mcmcPar->sigma_pi);
    	/* which should be the same as: */
	/* tempPar->pi = currentPar->pi*gsl_ran_lognormal(rng,0,mcmcPar->sigma_pi); */
    } while(tempPar->pi>1.0);


    /* ACCEPT / REJECT */
    /* likelihood */
    logRatio += loglike_kappa_all(tempPar) - loglike_kappa_all(currentPar);

    /* prior */
    logRatio += logprior_pi(tempPar) - logprior_pi(currentPar);

    /* ADD CORRECTION FOR MH truncated lognormal */
    QCur = gsl_cdf_gaussian_P(-log(currentPar->pi),mcmcPar->sigma_pi);
    QTemp = gsl_cdf_gaussian_P(-log(tempPar->pi),mcmcPar->sigma_pi);
    logRatio +=  log(tempPar->pi) - log(currentPar->pi); /* correction for lognormal */
    logRatio +=   log(QCur) - log(QTemp); /* correction for truncation (no values >1) */

    /* if p(new/old) > 1, accept new */
    if(logRatio>=0.0) {
	currentPar->pi = tempPar->pi;
	mcmcPar->n_accept_pi += 1;
	/* printf("\nAccepting new value\n"); */
    } else { /* else accept new with proba (new/old) */
	if(log(gsl_rng_uniform(rng)) <= logRatio){ /* accept */
	    currentPar->pi = tempPar->pi;
	    mcmcPar->n_accept_pi += 1;
	    /* printf("\nAccepting new value\n"); */
	} else { /* reject */
	    tempPar->pi = currentPar->pi;
	    mcmcPar->n_reject_pi += 1;
	    /* printf("\nRejecting new value\n"); */
	}
    }

} /* end move_pi */








/* MOVE VALUES OF PHI */
void move_phi(param *currentPar, param *tempPar, data *dat, mcmc_param *mcmcPar, gsl_rng *rng){
    int i;
    double logRatio=0.0;
    double QCur, QTemp;

    /* GENERATE CANDIDATE VALUE FOR PHI */
    /* HERE REPLACE WITH TRUNCATED LOGNORMAL (no values >1) )*/
    do
    {
    	tempPar->phi = gsl_ran_lognormal(rng,log(currentPar->phi),mcmcPar->sigma_phi);
    	/* which should be the same as: */
	/* tempPar->phi = currentPar->phi*gsl_ran_lognormal(rng,0,mcmcPar->sigma_phi); */
    } while(tempPar->phi>1.0);


    /* ACCEPT / REJECT */
   /* likelihood */
    logRatio += loglike_alpha_all(tempPar) - loglike_alpha_all(currentPar);

    /* prior */
    logRatio += logprior_phi(tempPar) - logprior_phi(currentPar);

    /* ADD CORRECTION FOR MH truncated lognormal */
    QCur = gsl_cdf_gaussian_P(-log(currentPar->phi),mcmcPar->sigma_phi);
    QTemp = gsl_cdf_gaussian_P(-log(tempPar->phi),mcmcPar->sigma_phi);
    logRatio +=  log(tempPar->phi) - log(currentPar->phi); /* correction for lognormal */
    logRatio +=   log(QCur) - log(QTemp); /* correction for truncation (no values >1) */

    /* if p(new/old) > 1, accept new */
    if(logRatio>=0.0) {
	currentPar->phi = tempPar->phi;
	mcmcPar->n_accept_phi += 1;
	/* printf("\nAccepting new value\n"); */
    } else { /* else accept new with proba (new/old) */
	if(log(gsl_rng_uniform(rng)) <= logRatio){ /* accept */
	    currentPar->phi = tempPar->phi;
	    mcmcPar->n_accept_phi += 1;
	    /* printf("\nAccepting new value\n"); */
	} else { /* reject */
	    tempPar->phi = currentPar->phi;
	    mcmcPar->n_reject_phi += 1;
	    /* printf("\nRejecting new value\n"); */
	}
    }

} /* end move_phi */














/* /\* MOVE VALUES OF ALPHA AND KAPPA *\/ */
/* void move_alpha_kappa(param *currentPar, param *tempPar, data *dat, dna_dist *dnainfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng){ */
/*     int i, j, oldestDate = 0, nCandidates=0, toMove=0, T, ances; */
/*     double logRatio = 0.0; */


/*     /\* DETERMINE WHICH alpha_i/kappa_i TO MOVE *\/ */
/*     /\*  - we change same kappa_i as alpha_i - *\/ */
/*     sample_vec_int(mcmcPar->all_idx, mcmcPar->idx_move_alpha, FALSE, rng); */

/*     /\* FIND OLDEST INFECTION (it has no ancestor) *\/ */
/*     oldestDate = min_vec_int(currentPar->Tinf); */

/*     /\* MOVE EACH alpha_i TO MOVE *\/ */
/*     for(i=0;i<mcmcPar->idx_move_alpha->length;i++){ */
/* 	toMove = vec_int_i(mcmcPar->idx_move_alpha,i); */

/* 	/\* move only isolates with an ancestor in the sample *\/ */
/* 	if(vec_int_i(currentPar->Tinf, toMove) > oldestDate){ */
/* 	    /\* find candidate ancestors ('alpha_i' so that T^inf_{alpha_i} < T^inf_i) *\/ */
/* 	    nCandidates=0; */
/* 	    for(j=0;j<dat->n;j++){ */
/* 		if(vec_int_i(currentPar->Tinf,j) < vec_int_i(currentPar->Tinf,toMove)) */
/* 		    mcmcPar->candid_ances->values[nCandidates++] =  j; */
/* 	    } */

/* 	    /\* check that there is at least a candidate or issue error *\/ */
/* 	    if(nCandidates==0){ */
/* 		fprintf(stderr, "\n[in: moves.c->move_alpha_kappa]\nNo candidate ancestor for 'i' but still trying to move 'alpha_i' (i: %d).\n", toMove); */
/* 		exit(1); */
/* 	    } */

/* 	    /\* GET PROPOSED ALPHA_I *\/ */
/* 	    tempPar->alpha->values[toMove] = vec_int_i(mcmcPar->candid_ances, gsl_rng_uniform_int(rng, nCandidates)); */

/* 	    /\* GET PROPOSED KAPPA_I *\/ */
/* 	    /\*  (GIBBS AS DISTRIBUTION IS NEG BIN) *\/ */
/* 	    /\* ances: proposed ancestor for 'toMove' *\/ */
/* 	    ances = vec_int_i(tempPar->alpha,toMove); */

/* 	    /\* time between ancestor and descendent *\/ */
/* 	    T = vec_int_i(tempPar->Tinf, toMove) - vec_int_i(tempPar->Tinf, ances); */

/* 	    /\* most likely value of kappa *\/ */
/* 	    /\* tempPar->kappa->values[toMove] = find_maxLike_kappa_i(T, gen); *\/ */


/* 	    /\* ACCEPT/REJECT STEP *\/ */
/* 	    /\* compute the likelihood *\/ */
/* 	    logRatio = loglikelihood_all(dat, dnainfo, gen, tempPar) - loglikelihood_all(dat, dnainfo, gen, currentPar); */

/* 	    /\* compute the priors *\/ */
/* 	    logRatio += logprior_alpha_i(toMove,tempPar) - logprior_alpha_i(toMove,currentPar); */
/* 	    /\* logRatio += logprior_kappa_i(toMove,tempPar) - logprior_kappa_i(toMove,currentPar); *\/ */

/* 	    /\* if p(new/old) > 1, accept new *\/ */
/* 	    if(logRatio>=0.0) { */
/* 		currentPar->alpha->values[toMove] = vec_int_i(tempPar->alpha,toMove); */
/* 		currentPar->kappa->values[toMove] = vec_int_i(tempPar->kappa,toMove); */
/* 		mcmcPar->n_accept_alpha += 1; */
/* 		mcmcPar->n_accept_kappa += 1; */
/* 	    } else { /\* else accept new with proba (new/old) *\/ */
/* 		if(log(gsl_rng_uniform(rng)) <= logRatio){ /\* accept *\/ */
/* 		    currentPar->alpha->values[toMove] = vec_int_i(tempPar->alpha,toMove); */
/* 		    currentPar->kappa->values[toMove] = vec_int_i(tempPar->kappa,toMove); */
/* 		    mcmcPar->n_accept_alpha += 1; */
/* 		    mcmcPar->n_accept_kappa += 1; */
/* 		} else { /\* reject *\/ */
/* 		    tempPar->alpha->values[toMove] = vec_int_i(currentPar->alpha,toMove); */
/* 		    tempPar->kappa->values[toMove] = vec_int_i(currentPar->kappa,toMove); */
/* 		    mcmcPar->n_reject_alpha += 1; */
/* 		    mcmcPar->n_reject_kappa += 1; */
/* 		} */
/* 	    } /\* end  ACCEPT/REJECT STEP *\/ */

/* 	} /\* end if isolate is not the oldest one *\/ */
/*     } /\* end for loop (for all 'i' to move) *\/ */
/* } /\* end move_alpha *\/ */






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

