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

int main(int argc, char *argv[]){
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

    clock_t end = clock() ;
    double TimeElapsed = (end-start)/(double)CLOCKS_PER_SEC ;

    printf("Computing time: %lf seconds\n",TimeElapsed);

    /* getchar(); */
    return 0;

}




/* 

   Compilation instructions: 

   gcc -o epilike -Wall -g alloc.c logL.c prior.c moves.c mcmc.c init.c InputOutput.c tuneVariances.c epilike.c -lgsl -lgslcblas && ./epilike

   valgrind epilike -v --leak-check=full --track-origins=yes

*/
