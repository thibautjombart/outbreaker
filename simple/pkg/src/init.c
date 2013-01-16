
#include "common.h"
#include "structures.h"

/* #include "init.h" */
/* #include "InputOutput.h" */
/* #include "logL.h" */
/* #include "mcmc.h" */
/* #include "moves.h" */
/* #include "tuneVariances.h" */



/* create a random number generator */
/* time_t: time in seconds, used to change the seed of the random generator */
gsl_rng * create_gsl_rng(time_t t){
    /* time_t t = time(NULL); /\* time in seconds, used to change the seed of the random generator *\/ */
    gsl_rng_env_setup();
    gsl_rng *out = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(out,t); /* changes the seed of the random generator */
    return out;
}



/* initialize and pre-compute generation time */
/* W: generation time distribution */
/* F: time to collection distribution */
void init_gentime(gentime *in, double *Wvalues, double *Fvalues){
    double sumDens=0.0;
    int i,j;


    /* PRE-COMPUTE DENSITIES */
    /* for kappa=1 */
    /* copy densities provided by user (from R) */
    for(j=0;j<in->truncW;j++){
	in->dens->rows[0]->values[j] = Wvalues[j];
    }

    /* normalize this density */
    sumDens = sum_vec_double(in->dens->rows[0]);
    for(j=0;j<in->truncW;j++){
	in->dens->rows[0]->values[j] = in->dens->rows[0]->values[j]/sumDens;
    }

    /* compute convolutions for kappa>1 */
    for(i=1;i<in->dens->n;i++){
	convol_vec_double(in->dens->rows[0], in->dens->rows[i-1], in->dens->rows[i]);
    }

    /* fill in values for f */
    /* copy densities provided by user (from R) */
    for(j=0;j<in->truncF;j++){
	in->collTime->values[j] = Fvalues[j];
    }

} /* end init_gentime */





/* FIND MOST LIKELY KAPPA_I (given time to infection 'T') */
int find_maxLike_kappa_i(int T, gentime *gen){
    int i, out=1;
    double temp=0.0, currentMax=0.0;

    /* printf("\nKappa_i values (T=%d):",T);fflush(stdout); */
    for(i=1;i<gen->maxK;i++){
	temp = gentime_dens(gen, T, i);
	/* printf("\nkappa_i=%d: %f",i,temp);fflush(stdout); */
	if(currentMax < temp) {
	    currentMax = temp;
	    out = i;
	}
    }

    /* printf("\nReturned value of kappa: %d\n",out);fflush(stdout); */

    return out;
} /* end find_maxLike_kappa_i */






/* INITIALIZE PARAMETERS */
void init_param(param *par, data *dat,  gentime *gen, int *ances, int *init_kappa, double pi_param1, double pi_param2, double init_mu1, double init_gamma, double outlier_threshold, gsl_rng *rng){
    int i, ancesId, T, TmaxLike;

    /* Tinf */
    TmaxLike = which_max_vec_double(gen->dens->rows[0]);
    for(i=0;i<dat->n;i++){
	par->Tinf->values[i] = vec_int_i(dat->dates,i) - TmaxLike;
    }

    /* alpha */
    for(i=0;i<dat->n;i++){
	par->alpha->values[i] = ances[i];
    }

    /* kappa */
    /* printf("\nGeneration time data\n");fflush(stdout); */
    /* print_gentime(gen); */
    for(i=0;i<dat->n;i++){
	/* par->kappa->values[i] = 1; */
	ancesId = vec_int_i(par->alpha,i);
	/* printf("\nInitial kappa_%d: %d\n",i,init_kappa[i]);fflush(stdout); */
	if(ancesId>-1){
	    if(init_kappa[i]<1){ /* value < 1 => find ML kappa */
		T = vec_int_i(par->Tinf, i) - vec_int_i(par->Tinf, ancesId);
		par->kappa->values[i] = find_maxLike_kappa_i(T, gen);
	    } else {
		par->kappa->values[i] = init_kappa[i]; /* otherwise, use specified kappa */
	    }
	} else { /* kappa = 1 by convention for imported cases */
	    par->kappa->values[i] = 1;
	}
	/* printf("\nInitialized kappa_%d: %d\n",i,par->kappa->values[i]);fflush(stdout); */
    }


    /* doubles*/
    par->mu1 = init_mu1;
    par->mu1_prior = init_mu1;
    par->gamma = init_gamma;
    par->pi = gsl_ran_beta (rng,pi_param1,pi_param2);
    par->pi_param1 = pi_param1;
    par->pi_param2 = pi_param2;
    par->outlier_threshold = (outlier_threshold>10.0) ? outlier_threshold : 10.0;
    /* par->phi = gsl_ran_beta (rng,phi_param1,phi_param2); */
    /* par->phi_param1 = phi_param1; */
    /* par->phi_param2 = phi_param2; */
}





void init_mcmc_param(mcmc_param *in, data *dat, bool move_mut, int *move_alpha, int *move_kappa, bool move_Tinf, bool move_pi, bool find_import, int burnin, int find_import_at){
    int i, N = dat->n;

    /* INITIALIZE COUNTERS */
    /* the first set of parameters is accepted by definition */
    /* param accepted: mu1 (1), gamma (1), some kappa (n_move_kappa), some alpha (n_move_alpha) */
    in->n_accept_mu1 = 1;
    in->n_reject_mu1 = 0;
    in->n_accept_gamma = 1;
    in->n_reject_gamma = 0;
    in->n_accept_Tinf = in->n_move_Tinf;
    in->n_reject_Tinf = 0;
    in->n_accept_alpha = in->n_move_alpha;
    in->n_reject_alpha = 0;
    in->n_accept_kappa = in->n_move_kappa;
    in->n_reject_kappa = 0;


    /* INITIALIZE MCMC PARAMETERS */
    in->sigma_mu1 = 0.0001;
    /* in->sigma_mu1 = 0.001; */
    in->sigma_gamma = 1;
    in->sigma_pi = 0.01;
    /* in->sigma_phi = 0.01; */
    in->lambda_Tinf = 1;
    in->step_notune = 0;
    /* in->Pmove_alpha_old = 1.0; */
    /* in->Pmove_alpha_new = 1.0; */
    in->burnin = burnin;
    in->find_import_at = find_import_at;

    /* FILL IN VECTORS */
    for(i=0;i<N;i++) {
	/* vector of all indices */
	in->all_idx->values[i] = i;
	/* vector of moved alpha_i*/
	in->move_alpha->values[i] = move_alpha[i] > 0.0 ? 1.0 : 0.0;
  	/* vector of moved kappa_i*/
	in->move_kappa->values[i] = move_kappa[i] > 0.0 ? 1.0 : 0.0;
    }

    /* FILL IN BOOLEANS */
    in->move_mut = move_mut;
    in->move_Tinf = move_Tinf;
    in->move_pi = move_pi;
    /* in->move_phi = move_phi; */
    in->find_import = find_import;

    /* ENSURE MOVE-TUNING CONSISTENCY */
    if(!move_mut){
	in->tune_mu1 = FALSE;
	in->tune_gamma = FALSE;
    }
    if(!move_pi) in->tune_pi = FALSE;
    /* if(!move_phi) in->tune_phi = FALSE; */
    in->tune_all = in->tune_mu1 || in->tune_gamma || in->tune_pi;

} /* end init_mcmc_param */



/* void InitMCMCSettings(mcmcInternals *MCMCSettings){ */
/*     MCMCSettings->NbSimul = 110000; */
/*     MCMCSettings->SubSample = 10; */
/*     MCMCSettings->BurnIn = 10000; */

/*     gsl_matrix_set(MCMCSettings->Sigma_beta,0,0,0.1); */
/*     gsl_matrix_set(MCMCSettings->Sigma_beta,0,1,0.1); */
/*     gsl_matrix_set(MCMCSettings->Sigma_beta,1,0,0.1); */
/*     gsl_matrix_set(MCMCSettings->Sigma_beta,1,1,0.1); */

/*     MCMCSettings->Sigma_betaWardOut=0.1; */
/*     MCMCSettings->Sigma_betaOutOut=0.1; */
/*     MCMCSettings->Sigma_mu=5; */
/*     MCMCSettings->Sigma_sigma=1; */
/*     MCMCSettings->Sigma_nu1=0.005; */
/*     MCMCSettings->Sigma_kappa=0.005; */
/*     MCMCSettings->Sigma_tau=10; */
/*     MCMCSettings->Sigma_alpha=0.1; */
/* } */




