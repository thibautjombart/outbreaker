

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
void init_gentime(gentime *in, int type, double param1, double param2, double param3){
    double sumDens=0.0;

    /* UPDATE OBJECT'S CONTENT */
    in->type = type;
    in->param1 = param1;
    in->param2 = param2;
    in->param3 = param3;

    /* PRE-COMPUTE DENSITIES */
    /* ! need to incorporate cases with Ki>1 */
    int i;
    switch(in->type){
    case 1: /* Poisson */
	for(i=0;i<in->dens->length;i++){
	    in->dens->values[i] =  gsl_ran_poisson_pdf((unsigned int) i, in->param1);
	}
	break;
    default:
	fprintf(stderr, "\n[in: init.c->init_gentime]\nMethod %d is unknown. Exiting.\n", in->type);
	exit(1);
    }

    /* NORMALIZE DENSITIES */
    sumDens = sum_vec_double(in->dens);
    for(i=0;i<in->dens->length;i++){
	in->dens->values[i] =in->dens->values[i]/sumDens;
    }
}




void init_param(param *par, data *dat,  gentime *gen, int *ances){
    int i, TmaxLike;

    /* Tinf */
    TmaxLike = which_max_vec_double(gen->dens);
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




