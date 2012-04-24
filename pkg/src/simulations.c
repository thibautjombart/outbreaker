#include "common.h"
#include "alloc.h"
#include "genclasses.h"
#include "InputOutput.h"
#include "simepid.h"
#include "simgen.h"




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

    printf("\nInitial nb data: \n");
    print_nbData(nbData);
    fflush(stdout);

    printf("\n\nInitial raw data: \n");
    print_rawData(data);
    fflush(stdout);

    printf("\n\nInitial augmented data: \n");
    print_augData(augData);
    fflush(stdout);

    int max_nb_lineages=5, haplo_length=10000;
    double mu_dist=3.0, sigma_dist=0.1, lambda_nlin=2;
    double nu1=5e-5, nu2=1e-4, lambdaNseq=1.0;

    epid_dna *alldna = create_epid_dna(NbCases, max_nb_lineages, haplo_length);

    evolve_epid_dna(alldna, indexInfector, mu_dist, sigma_dist, lambda_nlin, nu1, nu2, augData->C, rng);
    print_epid_dna(alldna);

    return 0;

    /* SAMPLING PROCEDURE */
    list_dnaseq *dnasample;
    dnasample = sample_epid_dna(alldna, nbData, data, lambdaNseq, nu1, nu2, augData->C, rng);



    /* Closing files and freeing memory */
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

  gcc -o simulation matvec.c alloc.c genclasses.c InputOutput.c prior.c genlike.c logL.c simepid.c simgen.c simulations.c -Wall -g -lgsl -lgslcblas

  ./simulation

  valgrind --leak-check=full -v simulation

*/
