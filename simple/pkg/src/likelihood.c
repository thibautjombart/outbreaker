
#include "common.h"
#include "matvec.h"
#include "genclasses.h"
#include "distances.h"
#include "structures.h"
#include "prior.h"
#include "likelihood.h"






/* LOG-LIKELIHOOD FOR INDIVIDUAL 'i' */
double loglikelihood_i(int i, data *dat, dna_dist *dnainfo, gentime *gen, param *par, gsl_rng *rng){
    int ances=vec_int_i(par->alpha,i);
    double out=0.0;


    /* = EXTERNAL CASES = */
    if(ances < 0){
	/* PROBA OF SAMPLING TIME */
	out = log(gentime_dens(gen, vec_int_i(dat->dates,i) - vec_int_i(par->Tinf,i), 1));

	/* PROBA OF EXTERNAL CASE */
	out += log(par->phi);

	/* PROBA OF INFECTION TIME (UNIFORM OVER TIMESPAN) */
	out -= log((double) dat->timespan);

	/* /\* SIMULATED GENETIC PROBA *\/ */
	/* out += sim_loglike_gen(dat, par, rng); */

	/* FILTER AND RETURN */
	filter_logprob(&out);
	/* printf("\nlikelihood (imported case): %f\n", out);fflush(stdout); */
	return out;
    }


    /* = INTERNAL CASES = */
    /* GENETIC LIKELIHOOD */
    /* dpois(0,0) returns -NaN, not 1! */
    if(mat_int_ij(dnainfo->nbcommon, i, ances)>0){
	/* transitions */
	out += log(gsl_ran_poisson_pdf((unsigned int) mat_int_ij(dnainfo->transi, i, ances), (double) mat_int_ij(dnainfo->nbcommon, i, ances) * (double) vec_int_i(par->kappa,i) * par->mu1));

	/* transversions */
	out += log(gsl_ran_poisson_pdf((unsigned int) mat_int_ij(dnainfo->transv, i, ances), (double) mat_int_ij(dnainfo->nbcommon, i, ances) * (double) vec_int_i(par->kappa,i) * par->gamma *par->mu1));
    }


    /* EPIDEMIOLOGICAL LIKELIHOOD */
    /* LIKELIHOOD OF COLLECTION DATE */
    out += log(gentime_dens(gen, vec_int_i(dat->dates,i) - vec_int_i(par->Tinf,i), 1));

    /* LIKELIHOOD OF INFECTION TIME */
    /* printf("\ninfection date: %.10f\n", log(gentime_dens(gen, vec_int_i(par->Tinf,i) - vec_int_i(par->Tinf,ances), vec_int_i(par->kappa,i)))); */
    out += log(gentime_dens(gen, vec_int_i(par->Tinf,i) - vec_int_i(par->Tinf,ances), vec_int_i(par->kappa,i)));

    /* PROBA OF NON-EXTERNAL INFECTION */
    out += log(1.0 - par->phi);

    /* PROBA OF (KAPPA_I-1) UNOBSERVED CASES */
    out += log(gsl_ran_negative_binomial_pdf((unsigned int) vec_int_i(par->kappa,i)-1, par->pi, 1.0));

    /* FILTER AND RETURN */
    filter_logprob(&out);

    /* printf("\nlikelihood (internal case): %f\n", out);fflush(stdout); */

    return out;
} /* end loglikelihood_i */





/* GENETIC LOG-LIKELIHOOD FOR INDIVIDUAL 'i' */
double loglikelihood_gen_i(int i, data *dat, dna_dist *dnainfo, param *par, gsl_rng *rng){
    int ances=vec_int_i(par->alpha,i);
    double out=0.0;

    /* IMPORTED CASES */
    if(ances < 0) {
	/* return sim_loglike_gen(dat, par, rng); */
	return 0.0;
    }


    /* NON-IMPORTED CASES */
    /* dpois(0,0) returns -NaN, not 1! */
    if(mat_int_ij(dnainfo->nbcommon, i, ances)>0){
	/* TRANSITIONS */
	/* printf("\ntransitions: %.10f\n", log(gsl_ran_poisson_pdf((unsigned int) mat_int_ij(dnainfo->transi, i, ances), (double) mat_int_ij(dnainfo->nbcommon, i, ances) * (double) vec_int_i(par->kappa,i) * par->mu1))); */

	out += log(gsl_ran_poisson_pdf((unsigned int) mat_int_ij(dnainfo->transi, i, ances), (double) mat_int_ij(dnainfo->nbcommon, i, ances) * (double) vec_int_i(par->kappa,i) * par->mu1));

	/* TRANSVERSIONS */
	/* printf("\ntransversions: %.10f\n",log(gsl_ran_poisson_pdf((unsigned int) mat_int_ij(dnainfo->transv, i, ances), (double) mat_int_ij(dnainfo->nbcommon, i, ances) * (double) vec_int_i(par->kappa,i) * par->gamma *par->mu1))); */

	out += log(gsl_ran_poisson_pdf((unsigned int) mat_int_ij(dnainfo->transv, i, ances), (double) mat_int_ij(dnainfo->nbcommon, i, ances) * (double) vec_int_i(par->kappa,i) * par->gamma *par->mu1));
    }


    /* RETURN */
    filter_logprob(&out);

    return out;
} /* end loglikelihood_i */





/* LOG-LIKELIHOOD FOR ALL INDIVIDUALS */
double loglikelihood_all(data *dat, dna_dist *dnainfo, gentime *gen, param *par, gsl_rng *rng){
    int i;
    double out=0.0;

    for(i=0;i<dat->n;i++){
	out += loglikelihood_i(i, dat, dnainfo, gen, par, rng);
    }

    filter_logprob(&out);

    return out;
}





/* GENETIC LOG-LIKELIHOOD FOR ALL INDIVIDUALS */
double loglikelihood_gen_all(data *dat, dna_dist *dnainfo, param *par, gsl_rng *rng){
    int i;
    double out=0.0;

    for(i=0;i<dat->n;i++){
	out += loglikelihood_gen_i(i, dat, dnainfo, par, rng);
    }

    filter_logprob(&out);

    return out;
}





/* LOG-LIKELIHOOD KAPPA FOR ALL INDIVIDUALS */
double loglike_kappa_all(param *par){
    double out = 0.0;
    int i;
    for(i=0;i<par->n;i++){
	out += log(gsl_ran_negative_binomial_pdf((unsigned int) vec_int_i(par->kappa,i)-1, par->pi, 1.0));
    }
    filter_logprob(&out);
    return out;
}






/* LOG-LIKELIHOOD ALPHA FOR ALL INDIVIDUALS */
double loglike_alpha_all(param *par){
    double out = 0.0;
    int i;
    for(i=0;i<par->n;i++){
	out += (vec_int_i(par->alpha,i)<0) ? log(par->phi) : log(1.0-par->phi);
    }
    filter_logprob(&out);
    return out;
}





/* LOG-POSTERIOR FOR ALL INDIVIDUALS */
double logposterior_all(data *dat, dna_dist *dnainfo, gentime *gen, param *par, gsl_rng *rng){
    double out = logprior_all(par) + loglikelihood_all(dat, dnainfo, gen, par, rng);

    filter_logprob(&out);

    return out;
}





/* SIMULATE GENETIC LOG-LIKELIHOOD */
double sim_loglike_gen(data *dat, param *par, gsl_rng *rng){
    double out = 0.0, lambda1, lambda2;

    /* compute param of the Poisson distributions */
    lambda1 = (double) dat->length * par->mu1;
    lambda2 = (double) dat->length * par->mu1 * par->gamma;

    /* compute log-likelihood */
    out += log(gsl_ran_poisson_pdf(gsl_ran_poisson(rng, lambda1) , lambda1));
    out += log(gsl_ran_poisson_pdf(gsl_ran_poisson(rng, lambda2) , lambda2));

    /* filter and return */
    filter_logprob(&out);

    return out;
} /* end sim_loglike_gen */





/* CHECK LOG-LIKELIHOOD FOR ALL INDIVIDUALS */
/* returns TRUE if all is fine, FALSE if likelihood is zero */
bool check_loglikelihood_all(data *dat, dna_dist *dnainfo, gentime *gen, param *par, gsl_rng *rng){
    int i, ances;
    double temp;
    bool out=TRUE;

    for(i=0;i<dat->n;i++){
	temp = loglikelihood_i(i, dat, dnainfo, gen, par, rng);
	filter_logprob(&temp);

	if(temp <= NEARMINUSINF){
	    out = FALSE;
	    printf("\nlikelihood for ancestry of %d is zero", i+1);
	    fflush(stdout);

	    /* display genetic likelihood */
	    temp = loglikelihood_gen_i(i, dat, dnainfo, par, rng);
	    filter_logprob(&temp);
	    printf("\ni=%d: genetic like is: %f", i+1, temp);
	    fflush(stdout);
	    if(temp <= NEARMINUSINF) printf(" (i.e., zero)");
	    fflush(stdout);

	    /* display epi likelihood */
	    ances=vec_int_i(par->alpha,i);

	    /* likelihood of collection date */
	    temp = log(gentime_dens(gen, vec_int_i(dat->dates,i) - vec_int_i(par->Tinf,i), 1));
	    filter_logprob(&temp);
	    printf("\ni=%d: collection date (t_%d=%d,Tinf_%d=%d) like is: %f", i+1, i+1, vec_int_i(dat->dates,i), i+1, vec_int_i(par->Tinf,i), temp);
	    fflush(stdout);
	    if(temp <= NEARMINUSINF) printf(" (i.e., zero)");
	    fflush(stdout);

	    /* likelihood of infection time */
	    temp = log(gentime_dens(gen, vec_int_i(par->Tinf,i) - vec_int_i(par->Tinf,ances), vec_int_i(par->kappa,i)));
	    filter_logprob(&temp);
	    printf("\ni=%d: infection time (Tinf=%d,Tances=%d) like is: %f", i+1, vec_int_i(par->Tinf,i),vec_int_i(par->Tinf,ances), temp);
	    fflush(stdout);
	    if(temp <= NEARMINUSINF) printf(" (i.e., zero)");
	    fflush(stdout);

	}
    }

    return out;
} /* end check_loglikelihood_all */



/*
>>>> TESTING <<<<
*/

/* #include "init.h" */
/* int main(){ */
/*   /\* DECLARATIONS *\/ */
/*     int TIMESPAN; */
/*     data *dat; */
/*     gentime *gen; */
/*     param *par; */
/*     dna_dist * dnainfo; */

/*     double logPrior, logLike, logPost; */



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
/*     printf("\nParameters\n"); */
/*     print_param(par); */


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
/*     free_data(dat); */
/*     free_gentime(gen); */
/*     free_dna_dist(dnainfo); */
/*     free_param(par); */

/*     return 0; */
/* } */



/*
  gcc instructions

  gcc -o likelihood matvec.c genclasses.c structures.c init.c distances.c prior.c likelihood.c -lgsl -lgslcblas -Wall -g


  gcc -o likelihood matvec.c genclasses.c structures.c init.c distances.c prior.c likelihood.c -lgsl -lgslcblas -Wall -O3

 ./likelihood

  valgrind --leak-check=full -v likelihood

*/

