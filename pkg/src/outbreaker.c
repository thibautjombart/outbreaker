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
/* MAIN                                                                       */
/******************************************************************************/

void R_outbreaker(int *nbPatients, int *duration, int *nbAdmVec, int *nbPosSwab, int *nbNegSwab, int *nbColPatients, int *idxColPatients, int *nbSeqPat, int *wardVec, int *tAdmVec, int *tDisVec, int *tPosSwab, int *tNegSwab, int *hospPres, int *idxSeqVec, double *tCollecVec, unsigned char *DNAbinInput, int *totNbSeq, int *seqLength, double *weightNaGen){

    int NbPatients=*nbPatients, T = *duration, NbSeq = *totNbSeq;
    char workspace[500]= "";

    /* === CREATE AND FILL OBJECTS === */
    /* MCMC */
    mcmcInternals *MCMCSettings = createMcmcInternals();

    /* DATA */
    /* allocate / fill in */
    nb_data *nb = createNbData(NbPatients, T, NbSeq); /* nb_data*/
    importNbData(nbAdmVec, nbPosSwab, nbNegSwab, nbColPatients, nbPatients, duration, idxColPatients, nbSeqPat, totNbSeq, nb);
    printf("\n== Nb data ==\n");
    print_nbData(nb);

    raw_data *data = createRawData(nb); /* epi data */
    printf("\ntPosSwab 10 first values:\n");
    int i;
    for(i=0;i<10;i++) printf("%d ", tPosSwab[i]);
    importRawData(wardVec, tAdmVec, tDisVec, tPosSwab, tNegSwab, hospPres, idxSeqVec, totNbSeq, tCollecVec, nb, data);
    printf("\n\n== Raw data ==\n");
    print_rawData(data);

    list_dnaseq *dna = DNAbin2list_dnaseq(DNAbinInput, totNbSeq, seqLength); /* DNA sequences */
    dna_dist *dnainfo = compute_dna_distances(dna); /* genetic distances - all DNA info needed */
    printf("\n\n== Dna dist info ==\n");
    print_dna_dist(dnainfo);

    aug_data *augData = createAugData(NbPatients, T); /* augmented data */
    printf("\n\n== Augmented data ==\n");
    print_augData(augData);

    /* PARAMETERS */
    parameters *param = createParam();
    printf("\n\n== Parameters ==\n");
    print_param(param);
    param->weightNaGen = *weightNaGen;


    /* OUTPUT FILES */
    output_files *Files = createFILES(workspace);


    /* ACCEPTANCE */
    acceptance *accept = createAcceptance(); /* accept is initialized to zero */
    isAcceptOK *acceptOK = createIsAcceptOK();
    NbProposals *nbProp = createNbProposals();

    printf("\nAcceptance done.\n");

    /* AUGMENTED DATA AND PARAMETER INITIALIZATION */
    /* Initialize in MCMC settings */
    InitMCMCSettings(MCMCSettings);
    printf("\nMCMC initialized\n");

    /* !!! TO WRITE !!! */
    InitParam(param);
    printf("\nParam initialized\n");

    /* initialize augmented data */
    InitAugData(param, nb, data, augData);
    printf("\nAugmented data initialized\n");


    /* OUTPUT FILE PREPARATION  */
    prepAllFiles(Files, NbPatients);
    printf("\nOutput files prepared\n");


    /*****************************************************/
    /***                 Launch MCMC                   ***/
    /*****************************************************/
    MCMCSettings->BurnIn = 100;
    MCMCSettings->NbSimul = 100;
    MCMCSettings->SubSample = 10;

    printf("\nStarting estimation\n");
    fflush(stdout);
    metro(MCMCSettings, param, data, nb, augData, dnainfo, accept, acceptOK, nbProp, Files);


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

} /* end R_outbreaker */






/******************************************************************************/
/* MAIN                                                                       */
/******************************************************************************/

int main(int argc, char *argv[]){
    time_t start, end;
    int TimeElapsed, NbPatients=142, T=100, NbSeq=10;
    char workspace[500];

    /* SET TIMER */
    time(&start);

    /**********************************/
    /***      SET THE WORKSPACE     ***/
    /**********************************/
    strcpy(workspace, "");


    /*****************************************************/
    /***           DECLARE ALL STRUCTURES              ***/
    /*****************************************************/

    /*MCMC */
    mcmcInternals *MCMCSettings = createMcmcInternals();

    /*RAW DATA*/
    nb_data *nb = createNbData(NbPatients, T, NbSeq);
    /* printf("\ncreated nb data\n"); */

    /* reading Nbdata */
    readFakeNbData(nb);
    /* printf("\nread fake nb data\n"); */

    raw_data *data = createRawData(nb);
    /* printf("\ncreated raw data\n"); */

    /*AUG DATA*/
    aug_data *augData = createAugData(NbPatients, T);
    /* printf("\ncreated aug data\n"); */

    /*parameters */
    parameters *param = createParam();
    param->weightNaGen = 0.01;
    /* printf("\ncreated param\n"); */

    /* Output files */
    output_files *Files = createFILES(workspace);

    /* Acceptance */
    acceptance *accept = createAcceptance(); /* accept is initialized to zero */
    isAcceptOK *acceptOK = createIsAcceptOK();
    NbProposals *nbProp = createNbProposals();


    /* DNA DISTANCES */
    dna_dist *dnainfo = NULL;

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

    prepAllFiles(Files, NbPatients);

    /*****************************************************/

    /*****************************************************/
    /***                 Launch MCMC                   ***/
    /*****************************************************/
    MCMCSettings->BurnIn = 100;
    MCMCSettings->NbSimul = 100;
    MCMCSettings->SubSample = 10;

    printf("Starting estimation\n");
    fflush(stdout);
    metro(MCMCSettings, param, data, nb, augData, dnainfo, accept, acceptOK, nbProp, Files);


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

    time(&end);
    TimeElapsed = (int) end - (int) start;

    printf("Computing time: %i seconds\n",TimeElapsed);

    /* getchar(); */
    return 0;

}




/* 

   Compilation instructions: 

   gcc -o outbreaker -Wall -g alloc.c matvec.c genclasses.c distances.c genlike.c logL.c prior.c moves.c mcmc.c init.c InputOutput.c tuneVariances.c outbreaker.c -lgsl -lgslcblas

   valgrind --leak-check=full outbreaker 

   valgrind -v --leak-check=full --track-origins=yes --show-reachable=yes outbreaker 


*/
