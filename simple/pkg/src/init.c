

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
void init_gentime(gentime *in, double *values){
    double sumDens=0.0;
    int i,j;


    /* PRE-COMPUTE DENSITIES */
    /* for kappa=1 */
    /* copy densities provided by user (from R) */
    for(j=0;j<in->trunc;j++){
	in->dens->rows[0]->values[j] = values[j];
    }

    /* normalize this density */
    sumDens = sum_vec_double(in->dens->rows[0]);
    for(j=0;j<in->trunc;j++){
	in->dens->rows[0]->values[j] = in->dens->rows[0]->values[j]/sumDens;
    }

    /* compute convolutions for kappa>1 */
    for(i=1;i<in->dens->n;i++){
	convol_vec_double(in->dens->rows[0], in->dens->rows[i-1], in->dens->rows[i]);
    }
} /* end init_gentime */




void init_param(param *par, data *dat,  gentime *gen, int *ances){
    int i, TmaxLike;

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
    for(i=0;i<dat->n;i++){
	par->kappa->values[i] = 1;
    }

    /* doubles*/
    par->mu1 = 1e-5;
    par->gamma = 1.0;
    par->pi = 0.5;
}





void init_mcmc_param(mcmc_param *in, data *dat){
    int i, N = dat->n;

    /* INITIALIZE COUNTERS */
    /* the first set of parameters is accepted by definition */
    /* param accepted: mu1 (1), gamma (1), some kappa (n_move_kappa), some alpha (n_move_alpha) */
    in->n_reject = 0;
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
    in->n_accept = in->n_accept_mu1 + in->n_accept_gamma + in->n_accept_Tinf + in->n_accept_alpha + in->n_accept_kappa;


    /* INITIALIZE MCMC PARAMETERS */
    /* in->sigma_mu1 = 0.000001; */
    in->sigma_mu1 = 0.0001;
    in->sigma_gamma = 1;
    in->lambda_Tinf = 1;

    /* FILL IN VECTOR OF ALL INDICES */
    for(i=0;i<N;i++) in->all_idx->values[i] = i;
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




