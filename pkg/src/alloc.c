#include "common.h"
#include "init.h"
#include "InputOutput.h"
#include "logL.h"
#include "mcmc.h"
#include "moves.h"
#include "prior.h"
#include "alloc.h"
#include "genclasses.h"
#include "matvec.h"
#include "tuneVariances.h"





/**************** nb_data ****************/
nb_data * createNbData(int NbPatients, int T){
    nb_data *nb = (nb_data *) malloc(sizeof(nb_data));
    if(nb == NULL){
	fprintf(stderr, "\n[in: alloc.c->createNbData]\nNo memory left for creating nbData. Exiting.\n");
	exit(1);
    }

    nb->NbAdmissions = (int *) calloc(NbPatients, sizeof(int));
    if(nb->NbAdmissions == NULL){
	fprintf(stderr, "\n[in: alloc.c->createNbData]\nNo memory left for creating nbData. Exiting.\n");
	exit(1);
    }

    nb->NbPosSwabs = (int *) calloc(NbPatients, sizeof(int));
    if(nb->NbPosSwabs == NULL){
	fprintf(stderr, "\n[in: alloc.c->createNbData]\nNo memory left for creating nbData. Exiting.\n");
	exit(1);
    }

    nb->NbNegSwabs = (int *) calloc(NbPatients, sizeof(int));
    if(nb->NbNegSwabs == NULL){
	fprintf(stderr, "\n[in: alloc.c->createNbData]\nNo memory left for creating nbData. Exiting.\n");
	exit(1);
    }

    nb->indexColonisedPatients = (int *) calloc(NbPatients, sizeof(int));
    if(nb->indexColonisedPatients == NULL){
	fprintf(stderr, "\n[in: alloc.c->createNbData]\nNo memory left for creating nbData. Exiting.\n");
	exit(1);
    }

    nb->M = (int *) calloc(NbPatients, sizeof(int));
    if(nb->M == NULL){
	fprintf(stderr, "\n[in: alloc.c->createNbData]\nNo memory left for creating nbData. Exiting.\n");
	exit(1);
    }

    nb->NbColonisedPatients = 0;
    nb->NbPatients = NbPatients;
    nb->T = T;

    return nb;
}





void freeNbData(nb_data *nb){
    free(nb->NbAdmissions);
    free(nb->NbPosSwabs);
    free(nb->NbNegSwabs);
    free(nb->indexColonisedPatients);
    free(nb->M);
    free(nb);
}






/**************** raw_data ****************/
raw_data *createRawData(nb_data *nb){
    int i;
    raw_data *data = (raw_data *) malloc(sizeof(raw_data));

    /* EPI DATA */
    if(data == NULL){
	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	exit(1);
    }

    data->ward = (int *) calloc(nb->NbPatients, sizeof(int));
    if(data->ward == NULL){
	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	exit(1);
    }

    data->PatientIndex = (int *) calloc(nb->NbPatients, sizeof(int));
    if(data->PatientIndex == NULL){
	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	exit(1);
    }

    /* data->timeSeq = (int *) calloc(nb->NbPatients, sizeof(int)); */
    /* if(data->timeSeq == NULL){ */
    /* 	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n"); */
    /* 	exit(1); */
    /* } */

    data->A = (gsl_vector **) calloc(nb->NbPatients, sizeof(gsl_vector *));
    if(data->A == NULL){
	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	exit(1);
    }

    data->D = (gsl_vector **) calloc(nb->NbPatients, sizeof(gsl_vector *));
    if(data->D == NULL){
	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	exit(1);
    }

    data->P = (gsl_vector **) calloc(nb->NbPatients, sizeof(gsl_vector *));
    if(data->P == NULL){
	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	exit(1);
    }

    data->N = (gsl_vector **) calloc(nb->NbPatients, sizeof(gsl_vector *));
    if(data->N == NULL){
	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	exit(1);
    }

    data->IsInHosp = (gsl_vector **) calloc(nb->NbPatients, sizeof(gsl_vector *));
    if(data->IsInHosp == NULL){
	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	exit(1);
    }

    for(i=0;i<nb->NbPatients;i++){
	data->A[i] = gsl_vector_calloc(nb->NbAdmissions[i]);
	data->D[i] = gsl_vector_calloc(nb->NbAdmissions[i]);
	if(nb->NbPosSwabs[i]>0)
	    {
		data->P[i] = gsl_vector_calloc(nb->NbPosSwabs[i]);
		nb->indexColonisedPatients[nb->NbColonisedPatients] = i;
		nb->NbColonisedPatients++;
	    }
	if(nb->NbNegSwabs[i]>0){data->N[i] = gsl_vector_calloc(nb->NbNegSwabs[i]);}
	data->IsInHosp[i] = gsl_vector_calloc(nb->T + 1); /* +1 because time goes from 0 to T */
    }

    /* simple integers */
    data->NbPatients = nb->NbPatients;
    data->T = nb->T;


    /* GENETIC DATA */
    /* S: indices of sequences collected for each patient */
    data->S = (int **) malloc(nb->NbPatients*sizeof(int *));
    if(data->S == NULL){
	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	exit(1);
    }
    for(i=0;i<nb->NbPatients;i++){
	data->S[i] = (int *) calloc(nb->M[i], sizeof(int));
	if(data->S[i] == NULL){
	    fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	    exit(1);
	}
    }

    /* Tcollec: collection times for each sequence */
    data->Tcollec = (double *) calloc(nb->NbPatients, sizeof(double));
    if(data->Tcollec == NULL){
	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	exit(1);
    }

    /* M: number of sequences collected for each patient */
    data->M = (int *) calloc(nb->NbPatients, sizeof(int));
    if(data->M == NULL){
	fprintf(stderr, "\n[in: alloc.c->createRawData]\nNo memory left for creating rawData. Exiting.\n");
	exit(1);
    }
    for(i=0;i<nb->NbPatients;i++){
	data->M[i] = nb->M[i];
    }


    /* RANDOM NUMBER GENERATOR */
    time_t t = time(NULL); // time in seconds, used to change the seed of the random generator
    const gsl_rng_type *typ;
    gsl_rng_env_setup();
    typ = gsl_rng_default;
    data->rng = gsl_rng_alloc(typ);
    gsl_rng_set(data->rng,t); // changes the seed of the random generator

    return data;
}





void freeRawData(raw_data *data){
    int i;

    for(i=0 ; i<data->NbPatients ; i++){
	gsl_vector_free(data->A[i]);
	gsl_vector_free(data->D[i]);
	gsl_vector_free(data->P[i]);
	gsl_vector_free(data->N[i]);
	gsl_vector_free(data->IsInHosp[i]);
	free(data->S[i]);
    }

    free(data->ward);
    free(data->PatientIndex);
    /* free(data->timeSeq); */
    free(data->A);
    free(data->D);
    free(data->P);
    free(data->N);
    free(data->IsInHosp);
    free(data->S);
    free(data->Tcollec);
    free(data->M);
    gsl_rng_free(data->rng);
    free(data);
}






/**************** aug_data ****************/
aug_data *createAugData(int NbPatients, int T){
    aug_data *augData = (aug_data *) malloc(sizeof(aug_data));
    if(augData == NULL){
	fprintf(stderr, "\n[in: alloc.c->createAugData]\nNo memory left for creating augData. Exiting.\n");
	exit(1);
    }


    augData->C = (int *) calloc(NbPatients, sizeof(int));
    if(augData->C == NULL){
	fprintf(stderr, "\n[in: alloc.c->createAugData]\nNo memory left for creating augData. Exiting.\n");
	exit(1);
    }

    augData->E = (int *) calloc(NbPatients, sizeof(int));
    if(augData->E == NULL){
	fprintf(stderr, "\n[in: alloc.c->createAugData]\nNo memory left for creating augData. Exiting.\n");
	exit(1);
    }

    augData->I0 = (int *) calloc(T, sizeof(int));
    if(augData->I0 == NULL){
	fprintf(stderr, "\n[in: alloc.c->createAugData]\nNo memory left for creating augData. Exiting.\n");
	exit(1);
    }

    augData->I1 = (int *) calloc(T, sizeof(int));
    if(augData->I1 == NULL){
	fprintf(stderr, "\n[in: alloc.c->createAugData]\nNo memory left for creating augData. Exiting.\n");
	exit(1);
    }

    augData->NbPatients = NbPatients;
    augData->T = T;

    return augData;
}





void freeAugData(aug_data *augData){
    free(augData->C);
    free(augData->E);
    free(augData->I0);
    free(augData->I1);
    free(augData);
}





void copyAugData(aug_data *augDataDest, aug_data *augDataSource){
    augDataDest->NbPatients = augDataSource->NbPatients;
    augDataDest->T = augDataSource->T;
    memcpy(augDataDest->C, augDataSource->C, augDataDest->NbPatients*sizeof(int));
    memcpy(augDataDest->E, augDataSource->E, augDataDest->NbPatients*sizeof(int));
    memcpy(augDataDest->I0, augDataSource->I0, augDataSource->T*sizeof(int));
    memcpy(augDataDest->I1, augDataSource->I1, augDataSource->T*sizeof(int));
}







/***************** param ******************/
parameters *createParam(){
    parameters *param = (parameters *) malloc(sizeof(parameters));
    if(param == NULL){
	fprintf(stderr, "\n[in: alloc.c->createParam]\nNo memory left for creating parameters. Exiting.\n");
	exit(1);
    }

    param->beta = gsl_matrix_calloc(2,2);

    return param;
}





void freeParam(parameters *param){
    gsl_matrix_free(param->beta);
    free(param);
}





void copyParam(parameters * paramDest, parameters * paramSource){
    gsl_matrix_memcpy (paramDest->beta,paramSource->beta);
    paramDest->betaWardOut = paramSource->betaWardOut;
    paramDest->betaOutOut = paramSource->betaOutOut;

    /* paramDest->Sp = paramSource->Sp; */
    paramDest->Se = paramSource->Se;

    paramDest->Pi = paramSource->Pi;

    paramDest->mu = paramSource->mu;
    paramDest->sigma = paramSource->sigma;

    paramDest->nu1 = paramSource->nu1;
    paramDest->nu2 = paramSource->nu2;

    paramDest->tau = paramSource->tau;
    paramDest->alpha = paramSource->alpha;
    paramDest->weightNaGen = paramSource->weightNaGen;
}







/************ MCMC internals **************/
mcmcInternals *createMcmcInternals(){
    mcmcInternals *MCMCSettings = (mcmcInternals *) malloc(sizeof(mcmcInternals));
    if(MCMCSettings == NULL){
	fprintf(stderr, "\n[in: alloc.c->createMcmcInternals]\nNo memory left for creating MCMCSettings. Exiting.\n");
	exit(1);
    }

    MCMCSettings->Sigma_beta = gsl_matrix_calloc(2,2);

    return MCMCSettings;
}





void printStdProp(mcmcInternals *MCMCSettings){
    int i,j;

    for (i=0;i<2;i++){
	for (j=0;j<2;j++){
	    printf("Std proposal for beta_%d,%d: %lg\n",i,j,gsl_matrix_get(MCMCSettings->Sigma_beta,i,j));
	}
    }
    printf("Std proposal for betaWardOut: %lg\n",MCMCSettings->Sigma_betaWardOut);
    printf("Std proposal for betaOutOut: %lg\n",MCMCSettings->Sigma_betaOutOut);
    printf("Std proposal for mu: %lg\n",MCMCSettings->Sigma_mu);
    printf("Std proposal for sigma: %lg\n",MCMCSettings->Sigma_sigma);
    printf("Std proposal for nu1: %lg\n",MCMCSettings->Sigma_nu1);
    printf("Std proposal for nu2: %lg\n",MCMCSettings->Sigma_nu2);
    printf("Std proposal for tau: %lg\n",MCMCSettings->Sigma_tau);
    printf("Std proposal for alpha: %lg\n",MCMCSettings->Sigma_alpha);

    fflush(stdout);
}





void freeMcmcInternals(mcmcInternals *MCMCSettings){
    gsl_matrix_free(MCMCSettings->Sigma_beta);

    free(MCMCSettings);
}







/************ Acceptance ***************/
acceptance *createAcceptance(){
    acceptance *accept = (acceptance *) malloc(sizeof(acceptance));
    if(accept == NULL){
	fprintf(stderr, "\n[in: alloc.c->createAcceptance]\nNo memory left for creating accept. Exiting.\n");
	exit(1);
    }

    accept->PourcAcc_beta = gsl_matrix_calloc(2,2);
    accept->PourcAcc_betaWardOut=0;
    accept->PourcAcc_betaOutOut=0;
    accept->PourcAcc_mu=0;
    accept->PourcAcc_sigma=0;
    accept->PourcAcc_nu1=0;
    accept->PourcAcc_nu2=0;
    accept->PourcAcc_tau=0;
    accept->PourcAcc_alpha=0;

    return accept;
}





void reInitiateAcceptance(acceptance *accept){
    int i,j;
    for(i=0;i<2;i++){
	for(j=0;j<2;j++){
	    gsl_matrix_set(accept->PourcAcc_beta,i,j,0);
	}
    }

    accept->PourcAcc_betaWardOut=0;
    accept->PourcAcc_betaOutOut=0;
    accept->PourcAcc_mu=0;
    accept->PourcAcc_sigma=0;
    accept->PourcAcc_nu1=0;
    accept->PourcAcc_nu2=0;
    accept->PourcAcc_tau=0;
    accept->PourcAcc_alpha=0;

}





void printAcceptance(acceptance *accept, NbProposals *NbProp){
    int i,j;

    for(i=0;i<2;i++){
	for(j=0;j<2;j++){
	    printf("Prob accept beta_%d,%d\t%lg\n",i,j,gsl_matrix_get(accept->PourcAcc_beta,i,j)/gsl_matrix_get(NbProp->NbProp_beta,i,j));
	    fflush(stdout);
	}
    }

    printf("Prob accept betaWardOut\t%lg\n",accept->PourcAcc_betaWardOut/NbProp->NbProp_betaWardOut);
    fflush(stdout);
    printf("Prob accept betaOutOut\t%lg\n",accept->PourcAcc_betaOutOut/NbProp->NbProp_betaOutOut);
    fflush(stdout);
    printf("Prob accept mu\t%lg\n",accept->PourcAcc_mu/NbProp->NbProp_mu);
    fflush(stdout);
    printf("Prob accept sigma\t%lg\n",accept->PourcAcc_sigma/NbProp->NbProp_sigma);
    fflush(stdout);
    printf("Prob accept nu1\t%lg\n",accept->PourcAcc_nu1/NbProp->NbProp_nu1);
    fflush(stdout);
    printf("Prob accept nu2\t%lg\n",accept->PourcAcc_nu2/NbProp->NbProp_nu2);
    fflush(stdout);
    printf("Prob accept tau\t%lg\n",accept->PourcAcc_tau/NbProp->NbProp_tau);
    fflush(stdout);
    printf("Prob accept alpha\t%lg\n",accept->PourcAcc_alpha/NbProp->NbProp_alpha);
    fflush(stdout);

}





void freeAcceptance(acceptance *accept){
    gsl_matrix_free(accept->PourcAcc_beta);
    free(accept);
}






/************* Is Acceptance OK *************/
isAcceptOK *createIsAcceptOK(){
    isAcceptOK *acceptOK = (isAcceptOK *) malloc(sizeof(isAcceptOK));
    if(acceptOK == NULL){
	fprintf(stderr, "\n[in: alloc.c->createIsAcceptOK]\nNo memory left for creating acceptOK. Exiting.\n");
	exit(1);
    }

    acceptOK->IsAccOK_beta = gsl_matrix_calloc(2,2);

    acceptOK->IsAccOK_betaWardOut=0;
    acceptOK->IsAccOK_betaOutOut=0;
    acceptOK->IsAccOK_mu=0;
    acceptOK->IsAccOK_sigma=0;
    acceptOK->IsAccOK_nu1=0;
    acceptOK->IsAccOK_nu2=0;
    acceptOK->IsAccOK_tau=0;
    acceptOK->IsAccOK_alpha=0;

    return acceptOK;
}





void freeIsAcceptOK(isAcceptOK *acceptOK){
    gsl_matrix_free(acceptOK->IsAccOK_beta);

    free(acceptOK);
}








/************ NbProposals ***************/
NbProposals *createNbProposals(){
    NbProposals *NbProp = (NbProposals *) malloc(sizeof(NbProposals));
    if(NbProp == NULL){
	fprintf(stderr, "\n[in: alloc.c->createNbProposals]\nNo memory left for creating NbProp. Exiting.\n");
	exit(1);
    }

    NbProp->NbProp_beta = gsl_matrix_calloc(2,2);

    NbProp->NbProp_betaWardOut=0;
    NbProp->NbProp_betaOutOut=0;
    NbProp->NbProp_mu=0;
    NbProp->NbProp_sigma=0;
    NbProp->NbProp_nu1=0;
    NbProp->NbProp_nu2=0;
    NbProp->NbProp_tau=0;
    NbProp->NbProp_alpha=0;

    return NbProp;
}





void reInitiateNbProp(NbProposals * NbProp){
    int i,j;

    for(i=0 ; i<2 ; i++)
	{
	    for(j=0 ; j<2 ; j++)
		{
		    gsl_matrix_set(NbProp->NbProp_beta,i,j,0);
		}
	}

    NbProp->NbProp_betaWardOut=0;
    NbProp->NbProp_betaOutOut=0;
    NbProp->NbProp_mu=0;
    NbProp->NbProp_sigma=0;
    NbProp->NbProp_nu1=0;
    NbProp->NbProp_nu2=0;
    NbProp->NbProp_tau=0;
    NbProp->NbProp_alpha=0;

}





void freeNbProposals(NbProposals *NbProp){
    gsl_matrix_free(NbProp->NbProp_beta);

    free(NbProp);
}







/************** OUTPUT FILES **************/
output_files *createFILES(char *workspace){
    output_files *fich = (output_files *) malloc(sizeof(output_files));
    if(fich == NULL){
	fprintf(stderr, "\n[in: alloc.c->createFILES]\nNo memory left for creating fich. Exiting.\n");
	exit(1);
    }

    char fileName[300];

    strcpy(fileName, workspace);
    strcat(fileName,"LogL.txt");
    fich->LogL = fopen(fileName,"w");
    if ( fich->LogL == NULL )
	{
	    printf("A problem occurred while opening a file. Check whether a file exists and is already opened.\n");
	    exit(1);
	}

    strcpy(fileName, workspace);
    strcat(fileName,"ColonDates.txt");
    fich->ColonDates = fopen(fileName,"w");
    if ( fich->ColonDates == NULL )
	{
	    printf("A problem occurred while opening a file. Check whether a file exists and is already opened.\n");
	    exit(1);
	}

    strcpy(fileName, workspace);
    strcat(fileName,"EndColonDates.txt");
    fich->EndColonDates = fopen(fileName,"w");
    if ( fich->EndColonDates == NULL )
	{
	    printf("A problem occurred while opening a file. Check whether a file exists and is already opened.\n");
	    exit(1);
	}

    strcpy(fileName, workspace);
    strcat(fileName,"Parameters.txt");
    fich->Parameters = fopen(fileName,"w");
    if (fich->Parameters == NULL )
	{
	    printf("A problem occurred while opening a file. Check whether a file exists and is already opened.\n");
	    exit(1);
	}

    return fich;
}






void freeFILES(output_files *fich){
    fclose(fich->LogL);
    fclose(fich->ColonDates);
    fclose(fich->EndColonDates);
    fclose(fich->Parameters);

    free(fich);
}
