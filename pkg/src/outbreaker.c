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
		  double *distMat, int *locations, int *spaModel,
		  int *ances, int *init_kappa, int *nIter, int *outputEvery, int *tuneEvery, 
		  double *piParam1, double *piParam2, 
		  double *phiParam1, double *phiParam2, 
		  double *initMu1, double *initGamma, 
		  double *initSpa1, double *initSpa2, 
		  double *spa1Prior, double *spa2Prior,
		  int *moveMut, int *moveAlpha, int *moveKappa, int *moveTinf, 
		  int *movePi, int *movePhi, int *moveSpa,
		  int *importMethod, int *findImportAt, int *burnin, 
		  double *outlierThreshold, int *maxK,
		  int *quiet, int *vecDist, int *stepStopTune,
		  char **resFileName, char **tuneFileName, int *seed, int *l, int *group_vec){
    /* DECLARATIONS */
    int N = *n;
    gsl_rng *rng;
    data *dat;
    gentime *gen;
    param *par;
    dna_dist * dnaInfo;
    spatial_dist * spatialInfo;
    mcmc_param * mcmcPar;
    int i,j, counter;
    mat_double *trans_mat;
    int num_of_groups = *l;

    bool checkLike;
    bool findImport = (bool) *importMethod>0;
    

    /* INITIALIZE RNG */
    /* rng = create_gsl_rng((time_t) time(NULL)); */
    rng = create_gsl_rng((time_t) *seed);


    /* CONVERT DATA */
    dat = Rinput2data(DNAbinInput, Tcollec, n, nSeq, length, idxCasesInDna, locations, group_vec, l);
    /* Rprintf("\n>>> Data <<<\n"); */
    /* print_data(dat); */


    /* /\* GET TIME SPAN *\/ */
    /* int TIMESPAN = max_vec_int(dat->dates) - min_vec_int(dat->dates); */
    /* Rprintf("\nTimespan is %d\n",TIMESPAN); */


    /* CREATE AND INIT GENERATION TIME */
    gen = alloc_gentime(*maxK, *wTrunc, *fTrunc);
    init_gentime(gen, gentimeDens, colltimeDens);
    /* Rprintf("\n>>> gentime info <<<\n"); */
    /* print_gentime(gen); */


    /* CREATE AND INIT PARAMETERS */
    par = alloc_param(N,num_of_groups);
    init_param(par, dat,  gen, ances, init_kappa, *piParam1, *piParam2, *phiParam1, *phiParam2, *initMu1, *initGamma, *initSpa1, *initSpa2, *spa1Prior, *spa2Prior, *outlierThreshold, *mutModel, *spaModel, *importMethod, rng, num_of_groups);
    /* Rprintf("\n>>> param <<<\n"); */
    /* print_param(par); */


    /* COMPUTE GENETIC DISTANCES */
    dnaInfo = compute_dna_distances(dat->dna, *mutModel);
    /* Rprintf("\n>>> DNA info <<<\n"); */
    /* print_dna_dist(dnaInfo); */


    /* CONVERT AND STORE SPATIAL DISTANCES */
    spatialInfo = doublevec2spatial_dist(distMat, n);
    /* Rprintf("\n>>> SPATIAL info <<<\n"); */
    /* print_spatial_dist(spatialInfo); */


    /* COMPUTE PRIORS */
    double logPrior = logprior_all(par);
    Rprintf("\nPrior value (log): %.10f\n", logPrior);/* fflush(stdout); */

   /* COMPUTE LIKELIHOOD */
    double logLike = loglikelihood_all(dat, dnaInfo, spatialInfo, gen, par, rng);
    Rprintf("\n\n = Initial Log-likelihood value: %f\n", logLike);

    /* COMPUTE POSTERIOR */
    double logPost = logposterior_all(dat, dnaInfo, spatialInfo, gen, par, rng);
    Rprintf("\nLog-posterior value: %.10f\n", logPost);/* fflush(stdout); */

    /* ALLOCATE AND INITIALIZE MCMC PARAMETERS */
    Rprintf("\nBefore check init LL\n");/* fflush(stdout); */

    mcmcPar = alloc_mcmc_param(N);
    init_mcmc_param(mcmcPar, par, dat, (bool) *moveMut, moveAlpha, moveKappa, (bool) *moveTinf, 
		    (bool) *movePi, (bool) *movePhi, (bool) *moveSpa, findImport, *burnin, *findImportAt);
    /* Rprintf("\nMCMC parameters\n");fflush(stdout); */
    /* print_mcmc_param(mcmcPar); */

    /* CHECK THAT INITIAL STATE HAS A NON-NULL LIKELIHOOD */
    checkLike = check_loglikelihood_all(dat, dnaInfo, spatialInfo, gen, par, rng);
    if(!checkLike){
      warning("\n\n!WARNING! Initial state of the chain has a likelihood of zero. The chain may never converge. Please consider using a different initial tree.\n");
    }

    Rprintf("\nAfter check init LL\n");/* fflush(stdout); */
    Rprintf("\nBefore MCMC\n");/* fflush(stdout); */

    /* RUN MCMC */
    mcmc(*nIter, *outputEvery, *resFileName, *tuneFileName, *tuneEvery,
	 (bool) *quiet, par, dat, dnaInfo, spatialInfo, gen, mcmcPar, rng);

    /* Rprintf("\nAfter MCMC\n");fflush(stdout); */

    /* FILL IN GENETIC DISTANCE VECTOR */
    counter = 0;
    for(i=0;i<(N-1);i++){
	for(j=i+1;j<N;j++){
	    vecDist[counter++] = mutation1_ij(i,j,dat,dnaInfo) + mutation2_ij(i,j,dat,dnaInfo);
	}
    }

    /* STORE STEP AT WHICH TUNING STOPPED */
    *stepStopTune = mcmcPar->step_notune;

    /* FREE MEMORY */
    free_data(dat);
    free_gentime(gen);
    free_param(par);
    free_dna_dist(dnaInfo);
    free_mcmc_param (mcmcPar);
    gsl_rng_free(rng);
} /* end R_outbreaker */


/* TESTING FUNCTION - the purpose of this function is so that I can pass data in from R, check that all structs are initialised and created correctly and then do simple likelihood calculations without running a whole mcmc*/
/* to run this: install the package, open R and run source("tester.R") which will define and then run an interface function with parameters defined at the bottom of tester.R */
void test_R(unsigned char *DNAbinInput, int *Tcollec, int *n, int *nSeq, int *length, 
		  int *idxCasesInDna, int *mutModel, double *gentimeDens, int *wTrunc, 
		  double *colltimeDens, int *fTrunc,
		  double *distMat, int *locations, int *spaModel,
		  int *ances, int *init_kappa, int *nIter, int *outputEvery, int *tuneEvery, 
		  double *piParam1, double *piParam2, 
		  double *phiParam1, double *phiParam2, 
		  double *initMu1, double *initGamma, 
		  double *initSpa1, double *initSpa2, 
		  double *spa1Prior, double *spa2Prior,
		  int *moveMut, int *moveAlpha, int *moveKappa, int *moveTinf, 
		  int *movePi, int *movePhi, int *moveSpa,
		  int *importMethod, int *findImportAt, int *burnin, 
		  double *outlierThreshold, int *maxK,
		  int *quiet, int *vecDist, int *stepStopTune,
		  char **resFileName, char **tuneFileName, int *seed, int *l, int *group_vec){
    /* DECLARATIONS */
    int N = *n;
    int num_of_groups = *l;
    gsl_rng *rng;
    data *dat;
    gentime *gen;
    param *par;
    dna_dist * dnaInfo;
    spatial_dist * spatialInfo;
    mcmc_param * mcmcPar;
    int i,j, counter;
    mat_double *trans_mat;

    bool checkLike;
    bool findImport = (bool) *importMethod>0;
    

    /* INITIALIZE RNG */
    rng = create_gsl_rng((time_t) *seed);


    /* CONVERT DATA */
    dat = Rinput2data(DNAbinInput, Tcollec, n, nSeq, length, idxCasesInDna, locations, group_vec, l);
    /* Rprintf("\n>>> Data <<<\n"); 
    print_data(dat); */

/* CREATE AND INIT GENERATION TIME */
    gen = alloc_gentime(*maxK, *wTrunc, *fTrunc);
    init_gentime(gen, gentimeDens, colltimeDens);
    /* Rprintf("\n>>> gentime info <<<\n");
    	print_gentime(gen); */


    /* CREATE AND INIT PARAMETERS */
    par = alloc_param(N,num_of_groups);
   init_param(par, dat,  gen, ances, init_kappa, *piParam1, *piParam2, *phiParam1, *phiParam2, *initMu1, *initGamma, *initSpa1, *initSpa2, *spa1Prior, *spa2Prior, *outlierThreshold, *mutModel, *spaModel, *importMethod, rng, num_of_groups);
    /* Rprintf("\n>>> param <<<\n");
    print_param(par);*/

	/* COMPUTE GENETIC DISTANCES */
    dnaInfo = compute_dna_distances(dat->dna, *mutModel);
    /* Rprintf("\n>>> DNA info <<<\n");
    print_dna_dist(dnaInfo); */


    /* CONVERT AND STORE SPATIAL DISTANCES */
    spatialInfo = doublevec2spatial_dist(distMat, n);
    /* Rprintf("\n>>> SPATIAL info <<<\n"); 
     print_spatial_dist(spatialInfo); */

	mcmcPar = alloc_mcmc_param(N);
    init_mcmc_param(mcmcPar, par, dat, (bool) *moveMut, moveAlpha, moveKappa, (bool) *moveTinf, 
		    (bool) *movePi, (bool) *movePhi, (bool) *moveSpa, findImport, *burnin, *findImportAt);
    /* Rprintf("\nMCMC parameters\n");
    print_mcmc_param(mcmcPar); */


   param *temp_par;
   temp_par = alloc_param(N,num_of_groups);
   copy_param(par,temp_par);
   /* mcmcPar->sigma_trans_mat = 1; */
   /* jiggle_trans_mat(par,temp_par,dat,spatialInfo,mcmcPar,rng,num_of_groups); */


   /* Rprintf("\n loglikelihood_i:%f \n", loglikelihood_i(2,dat, dnaInfo, spatialInfo, gen, temp_par, rng));
   Rprintf("\n loglikelihood_grp_all: %f \n", loglikelihood_grp_all(dat, temp_par, rng));
   Rprintf("\n loglikelihood_all: %f \n", loglikelihood_all(dat, dnaInfo, spatialInfo, gen, temp_par, rng));
   Rprintf("\n loglikelihood_local_i: %f \n", loglikelihood_local_i(2,dat, dnaInfo, spatialInfo, gen, temp_par, rng)); */

   /* int caseno = 4;
    int ances2 = find_sequenced_ancestor(caseno, dat, dnaInfo, par);
    Rprintf("\n ancestor: %d\n",ances2);
    Rprintf("\n kappa_temp: %d\n",par->kappa_temp);
Rprintf("\n mu1:%f\n",par->mu1);
Rprintf("\n mutation_ij:%d\n",mutation1_ij(caseno, ances2, dat, dnaInfo));
Rprintf("\n matrix value: %d\n",mat_int_ij(dnaInfo->mutation1, 0, 0));
Rprintf("\n com_nucl_ij:%d\n",com_nucl_ij(caseno, ances2, dat, dnaInfo));
Rprintf("\n proba_mut: %f\n", proba_mut(1, 9999, 1, 0.5));
print_dna_dist(dnaInfo);

    double genlike = loglikelihood_gen_i(caseno,dat,dnaInfo,par,rng);
    Rprintf("\n Genetic likelihood for i: %f\n",genlike); */

    /* COMPUTE PRIORS */
    double logPrior = logprior_all(par);
    Rprintf("\nPrior value (log): %.10f\n", logPrior);

   /* COMPUTE LIKELIHOOD */
    double logLike = loglikelihood_all(dat, dnaInfo, spatialInfo, gen, par, rng);
    Rprintf("\n\n = Initial Log-likelihood value: %f\n", logLike);

    /* COMPUTE POSTERIOR */
    double logPost = logposterior_all(dat, dnaInfo, spatialInfo, gen, par, rng);
    Rprintf("\nLog-posterior value: %.10f\n", logPost);

    /* ALLOCATE AND INITIALIZE MCMC PARAMETERS */
    Rprintf("\nBefore check init LL\n");

    mcmcPar = alloc_mcmc_param(N);
    init_mcmc_param(mcmcPar, par, dat, (bool) *moveMut, moveAlpha, moveKappa, (bool) *moveTinf, 
		    (bool) *movePi, (bool) *movePhi, (bool) *moveSpa, findImport, *burnin, *findImportAt);
    /* Rprintf("\nMCMC parameters\n");fflush(stdout); */
    /* print_mcmc_param(mcmcPar); */

    /* CHECK THAT INITIAL STATE HAS A NON-NULL LIKELIHOOD */
    checkLike = check_loglikelihood_all(dat, dnaInfo, spatialInfo, gen, par, rng);
    if(!checkLike){
      warning("\n\n!WARNING! Initial state of the chain has a likelihood of zero. The chain may never converge. Please consider using a different initial tree.\n");
    }

    /* RUN MCMC */
    /* mcmc(*nIter, *outputEvery, *resFileName, *tuneFileName, *tuneEvery,
	 (bool) *quiet, par, dat, dnaInfo, spatialInfo, gen, mcmcPar, rng); */

    
    Rprintf("%f",loglikelihood_grp_all(dat, par, rng));
}






/* 

   Compilation instructions: 

   gcc -o outbreaker -Wall -g alloc.c matvec.c genclasses.c distances.c genlike.c logL.c prior.c moves.c mcmc.c init.c InputOutput.c tuneVariances.c outbreaker.c -lgsl -lgslcblas

   valgrind --leak-check=full outbreaker 

   valgrind -v --leak-check=full --track-origins=yes --show-reachable=yes outbreaker 


*/


