#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
/* #include <conio.h> */

#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sf_gamma.h> /* for combinations calculation */
#include <gsl/gsl_eigen.h>

#include <omp.h>

#include <time.h>

#include "common.h"
#include "init.h"
#include "InputOutput.h"
#include "logL.h"
#include "mcmc.h"
#include "moves.h"
#include "prior.h"
#include "alloc.h"
#include "tuneVariances.h"

/******************************************************************************/
/* DECLARATION OF GLOBAL VARIABLES                                            */
/******************************************************************************/

gsl_rng * rng;

/******************************************************************************/
/* MAIN                                                                       */
/******************************************************************************/

#define OPTIONS "awrh"

int main(int argc, char *argv[]) 
{
	clock_t start = clock() ;
	int where = 1;
	char workspace[500];

	/*****************************************************/
	/***      CHANGE THE WORKSPACE IF NEEDED           ***/
	/*****************************************************/

	if(where==1)
	{
		strcpy(workspace, "");
	}else
	{
		strcpy(workspace, "\\\\fi--san01\\homes\\acori\\workspace\\EpiGenet\\");
	}

	/*****************************************************/

	/*****************************************************/
	/***     DEFINITION OF A RANDOM GENERATOR          ***/
	/*****************************************************/

	time_t t;
	t = 1;//time(NULL); /* time in seconds, used to change the seed of the random generator */
	const gsl_rng_type *typ;
	gsl_rng_env_setup();
	typ=gsl_rng_default;
	rng=gsl_rng_alloc(typ);
	gsl_rng_set(rng,t); /* changes the seed of the random generator */

	/*****************************************************/

	/*****************************************************/
	/***           DECLARE ALL STRUCTURES              ***/
	/*****************************************************/

	/*MCMC */
	mcmcInternals *MCMCSettings = createMcmcInternals();
	
	/*RAW DATA*/
	nb_data *nb = createNbData();
	/* reading Nbdata */
	/***************************************
	 *************** TO WRITE ***************
	 ****************************************/
	readFakeNbData(nb);
	raw_data *data = createRawData(nb);
	
	/*AUG DATA*/
	aug_data *augData = createAugData();

	/*parameters */
	parameters *param = createParam();

	/* Output files */
	output_files *Files = createFILES(workspace);

	/* Acceptance */
	acceptance *accept = createAcceptance(); /* accept is initialized to zero */
	isAcceptOK *acceptOK = createIsAcceptOK();
	NbProposals *nbProp = createNbProposals();

	/*****************************************************/

	/*****************************************************/
	/***         READ DATA AND MCMC SETTINGS           ***/
	/***        + INIT AUGDATA AND PARAMETERS          ***/
	/*****************************************************/

	/* reading data */
	/***************************************
	 *************** TO WRITE ***************
	 ****************************************/
	readFakeData(nb, data);

	/* filling in data->IsInHosp */
	CalculIsInHosp(nb, data);

	/* fill in MCMC settings */
	InitMCMCSettings(MCMCSettings);

	/* initialize parameters */
	/***************************************
	 *************** TO WRITE ***************
	 ****************************************/
	InitParam(param);

	/* initialize augmented data */
	InitAugData(param, nb, data, augData);

	/*****************************************************/



	/*****************************************************/
	/***   Preparing output files (writing headers)    ***/
	/*****************************************************/

	prepAllFiles(Files);

	/*****************************************************/

	/*****************************************************/
	/***                 Launch MCMC                   ***/
	/*****************************************************/

	printf("Starting estimation\n");
   	fflush(stdout);
   	metro (MCMCSettings, param, data, nb, augData, accept, acceptOK, nbProp, Files);


   	/*****************************************************/
   	/***        Close files and Free memory            ***/
   	/*****************************************************/

   	freeMcmcInternals(MCMCSettings);
   	freeNbData(nb);
   	freeRawData(data);
   	freeAugData(augData);
   	freeParam(param);
   	freeFILES(Files);
   	freeAcceptance(accept);
   	freeIsAcceptOK(acceptOK);
   	freeNbProposals(nbProp);
   	gsl_rng_free(rng);

   	clock_t end = clock() ;
	double TimeElapsed = (end-start)/(double)CLOCKS_PER_SEC ;

	printf("Computing time: %lf seconds\n",TimeElapsed);

	/* getchar(); */
	return 0;

}
