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
		  int *wType, double *wParam1, double *wParam2, double *wParam3, int *wTrunc, 
		  int *ances, int *nIter, int *outputEvery, int *quiet){
    /* DECLARATIONS */
    int N = *n, TIMESPAN;
    gsl_rng *rng;
    data *dat;
    gentime *gen;
    param *par;
    dna_dist * dnainfo;
    mcmc_param * mcmcPar;

    double logPrior, logLike, logPost;

    /* INITIALIZE RNG */
    rng = create_gsl_rng(time(NULL));


    /* CONVERT DATA */
    dat = Rinput2data(DNAbinInput, Tcollec, n, length);
    printf("\n>>> Data <<<\n");
    print_data(dat);


    /* GET TIME SPAN */
    TIMESPAN = max_vec_int(dat->dates) - min_vec_int(dat->dates);
    printf("\nTimespan is %d\n",TIMESPAN);


    /* CREATE AND INIT GENERATION TIME */
    gen = alloc_gentime(TIMESPAN, *wTrunc);
    init_gentime(gen, *wType, *wParam1, *wParam2, *wParam3);
    printf("\n>>> gentime info <<<\n");
    print_gentime(gen);


    /* CREATE AND INIT PARAMETERS */
    par = alloc_param(N);
    init_param(par, dat,  gen, ances);
    print_param(par);


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

    /* ALLOCATE AND INITIALIZE MCMC PARAMETERS */
    mcmcPar = alloc_mcmc_param(dat->n);
    init_mcmc_param(mcmcPar, dat);

    /* RUN MCMC */
    mcmc(*nIter, *outputEvery, "output.txt", (bool) *quiet, par, dat, dnainfo, gen, mcmcPar, rng);


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


