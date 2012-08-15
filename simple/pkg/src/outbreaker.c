#include "common.h"
#include "structures.h"
#include "genclasses.h"
#include "distances.h"
#include "init.h"
#include "prior.h"
#include "likelihood.h"
#include "moves.h"
#include "mcmc.h"




/*
  ======================
  MAIN EXTERNAL FUNCTION
  ======================
*/

void R_outbreaker(unsigned char *DNAbinInput, int *Tcollec, int *n, int *length, 
		  double *gentimeDens, int *wTrunc, 
		  int *ances, int *nIter, int *outputEvery, int *tuneEvery, 
		  double *pi_param1, double *pi_param2, int *quiet, int *vecDist, int *stepStopTune){
    /* DECLARATIONS */
    int N = *n, TIMESPAN;
    gsl_rng *rng;
    data *dat;
    gentime *gen;
    param *par;
    dna_dist * dnainfo;
    mcmc_param * mcmcPar;
    int i,j, counter;

    double logPrior, logLike, logPost;

    bool checkLike;


    /* INITIALIZE RNG */
    rng = create_gsl_rng(time(NULL));


    /* CONVERT DATA */
    dat = Rinput2data(DNAbinInput, Tcollec, n, length);
    /* printf("\n>>> Data <<<\n"); */
    /* print_data(dat); */


    /* GET TIME SPAN */
    TIMESPAN = max_vec_int(dat->dates) - min_vec_int(dat->dates);
    printf("\nTimespan is %d\n",TIMESPAN);


    /* CREATE AND INIT GENERATION TIME */
    gen = alloc_gentime(TIMESPAN, *wTrunc);
    init_gentime(gen, gentimeDens);
    /* printf("\n>>> gentime info <<<\n"); */
    /* print_gentime(gen); */


    /* CREATE AND INIT PARAMETERS */
    par = alloc_param(N);
    init_param(par, dat,  gen, ances, *pi_param1, *pi_param2, rng);
    print_param(par);


    /* COMPUTE GENETIC DISTANCES */
    dnainfo = compute_dna_distances(dat->dna);
    /* printf("\n>>> DNA info <<<\n"); */
    /* print_dna_dist(dnainfo); */


   /*  /\* COMPUTE PRIORS *\/ */
   /*  logPrior = logprior_all(par); */
   /*  printf("\nPrior value (log): %.10f\n", logPrior); */

   /* /\* COMPUTE LIKELIHOOD *\/ */
   /*  logLike = loglikelihood_all(dat, dnainfo, gen, par); */
   /*  printf("\n\n = Initial Log-likelihood value: %f\n", logLike); */

   /*  /\* COMPUTE POSTERIOR *\/ */
   /*  logPost = logposterior_all(dat, dnainfo, gen, par); */
   /*  printf("\nLog-posterior value: %.10f\n", logPost); */

    /* ALLOCATE AND INITIALIZE MCMC PARAMETERS */
    mcmcPar = alloc_mcmc_param(dat->n);
    init_mcmc_param(mcmcPar, dat);

    /* SET MCMC_PARAM */
    /* mcmcPar->sigma_mu1 = *sigma_mu1; */

    /* OPTIONAL - fix some parameters */

    /* /\* COMPUTE LIKELIHOOD *\/ */
    /* logLike = loglikelihood_all(dat, dnainfo, gen, par); */
    /* printf("\n\n = Initial Log-likelihood value (before mcmc call): %f\n", logLike); */
    /* fflush(stdout); */

    /* CHECK THAT INITIAL STATE HAS A NON-NULL LIKELIHOOD */
    checkLike = check_loglikelihood_all(dat, dnainfo, gen, par);
    if(!checkLike){
	fprintf(stderr, "\n\n!WARNING! Initial state of the chain has a likelihood of zero. The chain may never converge. Please consider using a different initial tree.\n");
	fflush(stdout);
    }

    /* RUN MCMC */
    mcmc(*nIter, *outputEvery, "output.txt", "mcmcOutput.txt", *tuneEvery, (bool) *quiet, par, dat, dnainfo, gen, mcmcPar, rng);


    /* FILL IN GENETIC DISTANCE VECTOR */
    counter = 0;
    for(i=0;i<(dat->n-1);i++){
	for(j=i+1;j<dat->n;j++){
	    vecDist[counter++] = mat_int_ij(dnainfo->transi,i,j) + mat_int_ij(dnainfo->transv,i,j);
	}
    }

    /* STORE STEP AT WHICH TUNING STOPPED */
    *stepStopTune = mcmcPar->step_notune;

    /* FREE MEMORY */
    free_data(dat);
    free_gentime(gen);
    free_param(par);
    free_dna_dist(dnainfo);
    free_mcmc_param (mcmcPar);
    gsl_rng_free(rng);
} /* end R_outbreaker */





/* 

   Compilation instructions: 

   gcc -o outbreaker -Wall -g alloc.c matvec.c genclasses.c distances.c genlike.c logL.c prior.c moves.c mcmc.c init.c InputOutput.c tuneVariances.c outbreaker.c -lgsl -lgslcblas

   valgrind --leak-check=full outbreaker 

   valgrind -v --leak-check=full --track-origins=yes --show-reachable=yes outbreaker 


*/


