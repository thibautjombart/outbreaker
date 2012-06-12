#include "common.h"
#include "structures.h"
#include "genclasses.h"
#include "distances.h"
#include "init.h"
/* #include "InputOutput.h" */
/* #include "logL.h" */
/* #include "mcmc.h" */
/* #include "moves.h" */
/* #include "prior.h" */
/* #include "tuneVariances.h" */




/*
  ======================
  MAIN EXTERNAL FUNCTION
  ======================
*/

void R_outbreaker(unsigned char *DNAbinInput, int *Tcollec, int *n, int *length, 
		  int *wType, double *wParam1, double *wParam2, double *wParam3, int *wTrunc, 
		  int *ances){
    /* DECLARATIONS */
    int N = *n;
    gsl_rng *rng;
    data *dat;
    gentime *gen;
    param *par;
    dna_dist * dnainfo;


    /* INITIALIZE RNG */
    rng = create_gsl_rng(time(NULL));


    /* CONVERT DATA */
    dat = Rinput2data(DNAbinInput, Tcollec, n, length);
    printf("\n>>> Data <<<\n");
    print_data(dat);


    /* /\* GET TIME SPAN *\/ */
    /* TIMESPAN = max_vec_int(dat->dates) - min_vec_int(dat->dates); */
    /* printf("\nTimespan is %d\n",TIMESPAN); */


    /* CREATE AND INIT GENERATION TIME */
    gen = alloc_gentime(*wTrunc);
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


    /* FREE MEMORY */
    gsl_rng_free(rng);
    free_data(dat);
    free_gentime(gen);
    free_param(par);
    free_dna_dist(dnainfo);
} /* end R_outbreaker */





/* 

   Compilation instructions: 

   gcc -o outbreaker -Wall -g alloc.c matvec.c genclasses.c distances.c genlike.c logL.c prior.c moves.c mcmc.c init.c InputOutput.c tuneVariances.c outbreaker.c -lgsl -lgslcblas

   valgrind --leak-check=full outbreaker 

   valgrind -v --leak-check=full --track-origins=yes --show-reachable=yes outbreaker 


*/


