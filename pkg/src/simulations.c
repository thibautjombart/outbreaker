#include "common.h"
#include "alloc.h"
#include "genclasses.h"
#include "InputOutput.h"
#include "simepid.h"
#include "simgen.h"
#include "init.h"
#include "logL.h"
#include "mcmc.h"
#include "moves.h"
#include "prior.h"
#include "tuneVariances.h"



/* TEST OF SIMULATIONS */
int main(){
    /* time_t t = time(NULL); /\* time in seconds, used to change the seed of the random generator *\/ */
    time_t t = 1; /* time in seconds, used to change the seed of the random generator */
    const gsl_rng_type *typ;
    gsl_rng_env_setup();
    typ=gsl_rng_default;
    gsl_rng * rng=gsl_rng_alloc(typ);
    gsl_rng_set(rng,t); /* changes the seed of the random generator */

    int where = 1;
    int i;

    char workspace[500];
    if(where==1)
	{
	    strcpy(workspace, "");
	}else
	{
	    strcpy(workspace, "\\\\fi--san01\\homes\\acori\\workspace\\SimulateStaphEpid\\");
	}

    /* for(i=0;i<115;i++)
       {
       printf("%c", workspace[i]);
       fflush(stdout);
       }
       printf("\n");
       fflush(stdout); */

    /* === SIMULATE EPIDEMIC (WITHOUT GENETIC) FIRST === */
    /* PARAMETERS */
    parameters *param = createParam();
    hospDurationParam *paramHosp = createHospDurationParam();
    readParameters(workspace,param,paramHosp);

    /* DATA */
    int Tmax=100;
    /* int NbPatientsMax = Tmax*(SizeWard0+SizeWard1); */
    int NbPatientsMax = 150;
    int NbColonisedPatients = 0;
    int NbSequences = 0;
    nb_data *nbData = createNbData(NbPatientsMax, Tmax, NbSequences, NbColonisedPatients);
    for(i=0 ; i<NbPatientsMax ; i++)
	{
	    nbData->NbAdmissions[i]=1;
	    nbData->NbPosSwabs[i]=Tmax;
	    nbData->NbNegSwabs[i]=Tmax;
	}
    raw_data *data = createRawData(nbData);
    for(i=0 ; i<NbPatientsMax ; i++)
	{
	    nbData->NbPosSwabs[i]=0;
	    nbData->NbNegSwabs[i]=0;
	}
    aug_data *augData = createAugData(NbPatientsMax, Tmax);
    int *indexInfector = (int *) calloc(nbData->NbPatients, sizeof(int)); /* -1 if infected from outside ; -100 never infected */

    /* SIMULATION ALGORITHM */
    int NbCases = SimulEpid(param, paramHosp, nbData, data, augData,indexInfector);
    writeData(NbCases, workspace, nbData, data, augData,indexInfector);


    /* === SIMULATE GENETIC DATA === */

    /* printf("\n>>> Initial nb data: \n"); */
    /* print_nbData(nbData); */
    /* fflush(stdout); */

    /* printf("\n\n>>> Initial raw data: \n"); */
    /* print_rawData(data); */
    /* fflush(stdout); */

    /* printf("\n\n>>> Initial augmented data: \n"); */
    /* print_augData(augData); */
    /* fflush(stdout); */

    int max_nb_lineages=5, haplo_length=100;
    double mu_dist=3.0, sigma_dist=0.1, lambda_nlin=1;
    double nu1=5e-5, nu2=1e-4, lambdaNseq=1.0;

    epid_dna *alldna = create_epid_dna(NbCases, max_nb_lineages, haplo_length);

    evolve_epid_dna(alldna, indexInfector, mu_dist, sigma_dist, lambda_nlin, nu1, nu2, augData->C, rng);
    /* print_epid_dna(alldna); */

    /* SAMPLING PROCEDURE */
    list_dnaseq *dnasample;
    dnasample = sample_epid_dna(alldna, nbData, data, lambdaNseq, nu1, nu2, augData->C, rng);

    /* printf("\n>>> Sampled DNA:\n"); */
    /* print_list_dnaseq(dnasample); */

    /* GET GENETIC DISTANCES */
    dna_dist *dnainfo = compute_dna_distances(dnasample); /* genetic distances - all DNA info needed */

    /* ESTIMATION */
    mcmcInternals *MCMCSettings = createMcmcInternals();
    param->weightNaGen = 0.001;

    /* OUTPUT FILES */
    output_files *Files = createFILES(workspace);

    /* ACCEPTANCE */
    acceptance *accept = createAcceptance(); /* accept is initialized to zero */
    isAcceptOK *acceptOK = createIsAcceptOK();
    NbProposals *nbProp = createNbProposals();

    /* INITIALIZE IN MCMC SETTINGS */
    InitMCMCSettings(MCMCSettings);
    /* printf("\nMCMC initialized\n"); */

    /* INITIALIZE PARAMETERS */
    InitParam(param);
    /* printf("\nParam initialized\n"); */

    /* INITIALIZE AUGMENTED DATA */
    InitAugData(param, nbData, data, augData);
    /* printf("\nAugmented data initialized\n"); */


    /* OUTPUT FILE PREPARATION  */
    prepAllFiles(Files, data->NbPatients);
    /* printf("\nOutput files prepared\n"); */


    /*****************************************************/
    /***                 Launch MCMC                   ***/
    /*****************************************************/
    MCMCSettings->BurnIn = 100;
    MCMCSettings->NbSimul = 100;
    MCMCSettings->SubSample = 10;

    printf("\nStarting estimation\n");
    fflush(stdout);
    metro(MCMCSettings, param, data, nbData, augData, dnainfo, accept, acceptOK, nbProp, Files);


    /* Closing files and freeing memory */
    freeMcmcInternals(MCMCSettings);
    freeParam(param);
    freeAugData(augData);
    freeNbData(nbData);
    freeRawData(data);
    freeHospDurationParam(paramHosp);
    free_epid_dna(alldna);
    free_list_dnaseq(dnasample);
    free(indexInfector);
    gsl_rng_free(rng);

    //getchar();
    return 0;
}



/*
  gcc instructions:

  gcc -o simulation matvec.c alloc.c genclasses.c distances.c InputOutput.c init.c prior.c genlike.c logL.c moves.c mcmc.c tuneVariances.c simepid.c simgen.c simulations.c -Wall -g -lgsl -lgslcblas

  ./simulation

  valgrind --leak-check=full -v simulation

*/
