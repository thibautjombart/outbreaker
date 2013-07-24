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

void R_outbreaker(unsigned char *DNAbinInput, int *Tcollec, int *n, int *nSeq, int *length, 
		  int *idxCasesInDna, int *mutModel, double *gentimeDens, int *wTrunc, 
		  double *colltimeDens, int *fTrunc,
		  double *distMat, int *spaModel,
		  int *ances, int *init_kappa, int *nIter, int *outputEvery, int *tuneEvery, 
		  double *pi_param1, double *pi_param2, 
		  double *init_mu1, double *init_gamma, 
		  double *init_spa1, double *init_spa2, 
		  double *spa1_prior, double *spa2_prior,
		  int *move_mut, int *move_alpha, int *move_kappa, int *move_Tinf, 
		  int *move_pi, int *move_spa,
		  int *find_import, int *burnin, int *find_import_at, 
		  double *outlier_threshold,
		  int *quiet, int *vecDist, int *stepStopTune,
		  char **res_file_name, char **tune_file_name, int *seed){
    /* DECLARATIONS */
    int N = *n;
    gsl_rng *rng;
    data *dat;
    gentime *gen;
    param *par;
    dna_dist * dnainfo;
    spatial_dist * spatialinfo;
    mcmc_param * mcmcPar;
    int i,j, counter;

    bool checkLike;


    /* INITIALIZE RNG */
    /* rng = create_gsl_rng((time_t) time(NULL)); */
    rng = create_gsl_rng((time_t) *seed);


    /* CONVERT DATA */
    dat = Rinput2data(DNAbinInput, Tcollec, n, nSeq, length, idxCasesInDna);
    /* Rprintf("\n>>> Data <<<\n"); */
    /* print_data(dat); */


    /* GET TIME SPAN */
    /* TIMESPAN = max_vec_int(dat->dates) - min_vec_int(dat->dates); */
    /* Rprintf("\nTimespan is %d\n",TIMESPAN); */


    /* CREATE AND INIT GENERATION TIME */
    gen = alloc_gentime(dat->timespan, *wTrunc, *fTrunc);
    init_gentime(gen, gentimeDens, colltimeDens);
    /* Rprintf("\n>>> gentime info <<<\n"); */
    /* print_gentime(gen); */


    /* CREATE AND INIT PARAMETERS */
    par = alloc_param(N);
    init_param(par, dat,  gen, ances, init_kappa, *pi_param1, *pi_param2, *init_mu1, *init_gamma, *init_spa1, *init_spa2, *spa1_prior, *spa2_prior, *outlier_threshold, *mutModel, *spaModel, rng);
    print_param(par);


    /* COMPUTE GENETIC DISTANCES */
    dnainfo = compute_dna_distances(dat->dna, *mutModel);
    /* Rprintf("\n>>> DNA info <<<\n"); */
    /* print_dna_dist(dnainfo); */


    /* CONVERT AND STORE SPATIAL DISTANCES */
    spatialinfo = doublevec2spatial_dist(distMat, n);
    Rprintf("\n>>> SPATIAL info <<<\n");
    print_spatial_dist(spatialinfo);


   /*  /\* COMPUTE PRIORS *\/ */
   /*  double logPrior = logprior_all(par); */
   /*  Rprintf("\nPrior value (log): %.10f\n", logPrior);fflush(stdout); */

   /* /\* COMPUTE LIKELIHOOD *\/ */
   /*  double logLike = loglikelihood_all(dat, dnainfo, gen, par, rng); */
   /*  Rprintf("\n\n = Initial Log-likelihood value: %f\n", logLike); */

   /*  /\* COMPUTE POSTERIOR *\/ */
   /*  double logPost = logposterior_all(dat, dnainfo, gen, par, rng); */
   /*  Rprintf("\nLog-posterior value: %.10f\n", logPost);fflush(stdout); */

   /*  /\* ALLOCATE AND INITIALIZE MCMC PARAMETERS *\/ */
   /*  Rprintf("\nBefore check init LL\n");fflush(stdout);fflush(stdout); */

    mcmcPar = alloc_mcmc_param(dat->n);
    init_mcmc_param(mcmcPar, dat, (bool) *move_mut, move_alpha, move_kappa, (bool) *move_Tinf, 
		    (bool) *move_pi, (bool) *move_spa, (bool) *find_import, *burnin, *find_import_at);
    /* Rprintf("\nMCMC parameters\n");fflush(stdout); */
    /* print_mcmc_param(mcmcPar); */

    /* CHECK THAT INITIAL STATE HAS A NON-NULL LIKELIHOOD */
    checkLike = check_loglikelihood_all(dat, dnainfo, spatialinfo, gen, par, rng);
    if(!checkLike){
      warning("\n\n!WARNING! Initial state of the chain has a likelihood of zero. The chain may never converge. Please consider using a different initial tree.\n");
      /* fprintf(stderr, "\n\n!WARNING! Initial state of the chain has a likelihood of zero. The chain may never converge. Please consider using a different initial tree.\n"); */
      /* fflush(stdout); */
    }

    /* Rprintf("\nAfter check init LL\n");fflush(stdout); */
    /* Rprintf("\nBefore MCMC\n");fflush(stdout); */

    /* RUN MCMC */
    mcmc(*nIter, *outputEvery, *res_file_name, *tune_file_name, *tuneEvery,
	 (bool) *quiet, par, dat, dnainfo, spatialinfo, gen, mcmcPar, rng);

    /* Rprintf("\nAfter MCMC\n");fflush(stdout); */

    /* FILL IN GENETIC DISTANCE VECTOR */
    counter = 0;
    for(i=0;i<(N-1);i++){
	for(j=i+1;j<N;j++){
	    vecDist[counter++] = mutation1_ij(i,j,dat,dnainfo) + mutation2_ij(i,j,dat,dnainfo);
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


