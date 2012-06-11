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
		  int *wType, double *wParam1, double *wParam2, double *wParam3, int *wTrunc){
    /* DECLARATIONS */
    int N, TIMESPAN;
    data *dat;
    gsl_rng *rng;
    gentime *gen;
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


    /* COMPUTE GENETIC DISTANCES */
    dnainfo = compute_dna_distances(dat->dna);
    printf("\n>>> DNA info <<<\n");
    print_dna_dist(dnainfo);


    /* FREE MEMORY */
    free_data(dat);
    gsl_rng_free(rng);
    free_dna_dist(dnainfo);
} /* end R_outbreaker */





/* 

   Compilation instructions: 

   gcc -o outbreaker -Wall -g alloc.c matvec.c genclasses.c distances.c genlike.c logL.c prior.c moves.c mcmc.c init.c InputOutput.c tuneVariances.c outbreaker.c -lgsl -lgslcblas

   valgrind --leak-check=full outbreaker 

   valgrind -v --leak-check=full --track-origins=yes --show-reachable=yes outbreaker 


*/


