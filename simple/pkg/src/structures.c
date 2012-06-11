#include "common.h"
#include "matvec.h"
#include "genclasses.h"
#include "structures.h"


/*
  ======================
  >>>>> DATA <<<<<
  ======================
*/
data *alloc_data(int n, int length){
    /* allocate pointer */
    data *out = (data *) malloc(sizeof(data));
    if(out == NULL){
	fprintf(stderr, "\n[in: alloc.c->alloc_data]\nNo memory left for creating data. Exiting.\n");
	exit(1);
    }

    /* fill in integers */
    out->n = n;
    out->length = length;

    /* dates: collection times for each sequence */
    out->dates = alloc_vec_int(n);

    /* dna: list of DNA sequences */
    out->dna = alloc_list_dnaseq(n, length);

    /* RANDOM NUMBER GENERATOR */
    time_t t = time(NULL); /* time in seconds, used to change the seed of the random generator */
    /* time_t t = 1; /\* time in seconds, used to change the seed of the random generator *\/ */
    gsl_rng_env_setup();
    const gsl_rng_type *typ=gsl_rng_default;
    out->rng = gsl_rng_alloc(typ);
    gsl_rng_set(out->rng,t); /* changes the seed of the random generator */

    return out;
} /* end alloc_data */





void free_data(data *in){
    free_vec_int(in->dates);
    free_list_dnaseq(in->dna);
    gsl_rng_free(in->rng);
    free(in);
} /* end free_data*/




void print_data(data *in){
    printf("\n= Collection dates =\n");
    print_vec_int(in->dates);
    printf("\n= Sequences =");
    print_list_dnaseq(in->dna);
} /* end print_data*/







/* /\**************** augdata ****************\/ */
/* augdata *alloc_AugData(int NbPatients, int T, int NbSequences){ */
/*     augdata *augData = (augdata *) malloc(sizeof(augdata)); */
/*     if(augData == NULL){ */
/* 	fprintf(stderr, "\n[in: alloc.c->alloc_AugData]\nNo memory left for creating augData. Exiting.\n"); */
/* 	exit(1); */
/*     } */


/*     augData->C = (int *) calloc(NbPatients, sizeof(int)); */
/*     if(augData->C == NULL){ */
/* 	fprintf(stderr, "\n[in: alloc.c->alloc_AugData]\nNo memory left for creating augData. Exiting.\n"); */
/* 	exit(1); */
/*     } */

/*     augData->E = (int *) calloc(NbPatients, sizeof(int)); */
/*     if(augData->E == NULL){ */
/* 	fprintf(stderr, "\n[in: alloc.c->alloc_AugData]\nNo memory left for creating augData. Exiting.\n"); */
/* 	exit(1); */
/*     } */

/*     augData->I0 = (int *) calloc(T, sizeof(int)); */
/*     if(augData->I0 == NULL){ */
/* 	fprintf(stderr, "\n[in: alloc.c->alloc_AugData]\nNo memory left for creating augData. Exiting.\n"); */
/* 	exit(1); */
/*     } */

/*     augData->I1 = (int *) calloc(T, sizeof(int)); */
/*     if(augData->I1 == NULL){ */
/* 	fprintf(stderr, "\n[in: alloc.c->alloc_AugData]\nNo memory left for creating augData. Exiting.\n"); */
/* 	exit(1); */
/*     } */

/*     augData->alpha = gsl_matrix_calloc(NbSequences,NbSequences); */
/*     if(augData->alpha == NULL){ */
/* 	fprintf(stderr, "\n[in: alloc.c->alloc_AugData]\nNo memory left for creating augData. Exiting.\n"); */
/* 	exit(1); */
/*     } */

/*     augData->tau = gsl_matrix_calloc(NbSequences,NbSequences); */
/*     if(augData->tau == NULL){ */
/* 	fprintf(stderr, "\n[in: alloc.c->alloc_AugData]\nNo memory left for creating augData. Exiting.\n"); */
/* 	exit(1); */
/*     } */

/*     augData->NbPatients = NbPatients; */
/*     augData->T = T; */
/*     augData->NbSequences = NbSequences; */

/*     return augData; */
/* } */





/* void freeAugData(augdata *augData){ */
/*     free(augData->C); */
/*     free(augData->E); */
/*     free(augData->I0); */
/*     free(augData->I1); */
/*     gsl_matrix_free(augData->alpha); */
/*     gsl_matrix_free(augData->tau); */
/*     free(augData); */
/* } */




/* void print_augData(augdata *augData){ */
/*     int i; */
/*     printf("\nNb of patients: %d, time span 0-%d", augData->NbPatients, augData->T); */
/*     printf("\nColonisation times: \n"); */
/*     for(i=0;i<augData->NbPatients;i++) printf("%d ", augData->C[i]); */
/*     printf("\nClearance times: \n"); */
/*     for(i=0;i<augData->NbPatients;i++) printf("%d ", augData->E[i]); */
/*     printf("\nNb of colonised individuals at each time step, ward 0: \n"); */
/*     for(i=0;i<augData->T;i++) printf("%d ", augData->I0[i]); */
/*     printf("\nNb of colonised individuals at each time step, ward 1: \n"); */
/*     for(i=0;i<augData->T;i++) printf("%d ", augData->I1[i]); */
/*     fflush(stdout); */
/*     printf("\nMatrix of alpha_kq: \n"); */
/*     gsl_matrix_fprintf(stdout, augData->alpha, "%.3f"); */
/*     fflush(stdout); */
/*     printf("\nMatrix of tau_kq: \n"); */
/*     gsl_matrix_fprintf(stdout, augData->tau, "%.3f"); */
/*     fflush(stdout); */
/* } */




/* void copyAugData(augdata *augDataDest, augdata *augDataSource){ */
/*     augDataDest->NbPatients = augDataSource->NbPatients; */
/*     augDataDest->T = augDataSource->T; */
/*     memcpy(augDataDest->C, augDataSource->C, augDataDest->NbPatients*sizeof(int)); */
/*     memcpy(augDataDest->E, augDataSource->E, augDataDest->NbPatients*sizeof(int)); */
/*     memcpy(augDataDest->I0, augDataSource->I0, augDataSource->T*sizeof(int)); */
/*     memcpy(augDataDest->I1, augDataSource->I1, augDataSource->T*sizeof(int)); */
/*     gsl_matrix_memcpy(augDataDest->alpha, augDataSource->alpha); */
/*     gsl_matrix_memcpy(augDataDest->tau, augDataSource->tau); */
/* } */






/* /\***************** param ******************\/ */
/* parameters *alloc_Param(){ */
/*     parameters *param = (parameters *) malloc(sizeof(parameters)); */
/*     if(param == NULL){ */
/* 	fprintf(stderr, "\n[in: alloc.c->alloc_Param]\nNo memory left for creating parameters. Exiting.\n"); */
/* 	exit(1); */
/*     } */

/*     param->beta = gsl_matrix_calloc(2,2); */
/*     param->betaWardOut = 0.1; */
/*     param->betaOutOut = 0.1; */
/*     param->Se = 0.0; */
/*     param->Pi = 0.0; */
/*     param->mu = 0.0; */
/*     param->sigma = 0.0; */
/*     param->nu1 = 0.0; */
/*     param->kappa = 0.0; */
/*     param->weightNaGen = 0.0; */

/*     return param; */
/* } */





/* void freeParam(parameters *param){ */
/*     gsl_matrix_free(param->beta); */
/*     free(param); */
/* } */





/* void copyParam(parameters * paramDest, parameters * paramSource){ */
/*     gsl_matrix_memcpy (paramDest->beta,paramSource->beta); */
/*     paramDest->betaWardOut = paramSource->betaWardOut; */
/*     paramDest->betaOutOut = paramSource->betaOutOut; */

/*     /\* paramDest->Sp = paramSource->Sp; *\/ */
/*     paramDest->Se = paramSource->Se; */

/*     paramDest->Pi = paramSource->Pi; */

/*     paramDest->mu = paramSource->mu; */
/*     paramDest->sigma = paramSource->sigma; */

/*     paramDest->nu1 = paramSource->nu1; */
/*     paramDest->kappa = paramSource->kappa; */

/*     paramDest->weightNaGen = paramSource->weightNaGen; */
/* } */






/* void print_param(parameters *param){ */
/*     printf("\nBeta matrix:\n"); */
/*     gsl_matrix_fprintf(stdout, param->beta, "%.3f"); */
/*     printf("\nBeta ward<-out: %.3f", param->betaWardOut); */
/*     printf("\nBeta out<-out: %.3f", param->betaOutOut); */
/*     printf("\nsensibility of the test: %.3f", param->Se); */
/*     printf("\nprobability of being colonized at first admission: %.3f", param->Pi); */
/*     printf("\nmu/sigma - duration of colonisation: %.3f %.3f", param->mu, param->sigma); */
/*     printf("\nnu1: %.3f   kappa: %3f   nu2: %.3f", param->nu1, param->kappa, param->nu1*param->kappa); */
/*     printf("\nweightNaGen: %.3f", param->weightNaGen); */
/*     fflush(stdout); */
/* } */


/* void readParameters(char* workspace, parameters * param, hospDurationParam *paramHosp) */
/* { */
/*     FILE * paramInit; */
/*     char val[30]; */
/*     char file[200]; */

/*     strcpy(file, workspace); */
/*     strcat(file,"param.txt"); */
/*     if ((paramInit=fopen(file,"r"))==NULL) */
/* 	{ */
/* 	    printf("Cannot read param.txt"); */
/* 	    exit(2); */
/* 	} */

/*     fscanf(paramInit,"%s %lf",val,gsl_matrix_ptr (param->beta,0,0)); */
/*     printf("%s %g\n",val,gsl_matrix_get(param->beta,0,0)); */

/*     fscanf(paramInit,"%s %lf",val,gsl_matrix_ptr (param->beta,0,1)); */
/*     printf("%s %g\n",val,gsl_matrix_get(param->beta,0,1)); */

/*     fscanf(paramInit,"%s %lf",val,gsl_matrix_ptr (param->beta,1,0)); */
/*     printf("%s %g\n",val,gsl_matrix_get(param->beta,1,0)); */

/*     fscanf(paramInit,"%s %lf",val,gsl_matrix_ptr (param->beta,1,1)); */
/*     printf("%s %g\n",val,gsl_matrix_get(param->beta,1,1)); */

/*     fscanf(paramInit,"%s %lf",val,&param->betaWardOut); */
/*     printf("%s %g\n",val,param->betaWardOut); */

/*     fscanf(paramInit,"%s %lf",val,&param->betaOutOut); */
/*     printf("%s %g\n",val,param->betaOutOut); */

/*     /\*fscanf(paramInit,"%s %lf",val,&param->Sp); */
/*       printf("%s %g\n",val,param->Sp);*\/ */

/*     fscanf(paramInit,"%s %lf",val,&param->Se); */
/*     printf("%s %g\n",val,param->Se); */

/*     fscanf(paramInit,"%s %lf",val,&param->Pi); */
/*     printf("%s %g\n",val,param->Pi); */

/*     fscanf(paramInit,"%s %lf",val,&paramHosp->mu); */
/*     printf("%s %g\n",val,paramHosp->mu); */

/*     fscanf(paramInit,"%s %lf",val,&paramHosp->sigma); */
/*     printf("%s %g\n",val,paramHosp->sigma); */

/*     fscanf(paramInit,"%s %lf",val,&param->mu); */
/*     printf("%s %g\n",val,param->mu); */

/*     fscanf(paramInit,"%s %lf",val,&param->sigma); */
/*     printf("%s %g\n",val,param->sigma); */

/*     fscanf(paramInit,"%s %lf",val,&param->nu1); */
/*     printf("%s %g\n",val,param->nu1); */

/*     fscanf(paramInit,"%s %lf",val,&param->kappa); */
/*     printf("%s %g\n",val,param->kappa); */

/*     fclose(paramInit); */

/* } */

/* /\***************** hospDurationParam ******************\/ */
/* hospDurationParam *alloc_HospDurationParam(){ */
/*     hospDurationParam *param = (hospDurationParam *) malloc(sizeof(hospDurationParam)); */
/*     if(param == NULL){ */
/* 	fprintf(stderr, "\n[in: alloc.c->alloc_HospDurationParam]\nNo memory left for creating hospDurationParam. Exiting.\n"); */
/* 	exit(1); */
/*     } */

/*     param->mu = 0.0; */
/*     param->sigma = 0.0; */

/*     return param; */
/* } */


/* void print_HospDurationParam(hospDurationParam *param){ */
/*     printf("\nmu/sigma - duration of hospitalisation: %.3f %.3f", param->mu, param->sigma); */
/*     fflush(stdout); */
/* } */



/* void freeHospDurationParam(hospDurationParam *in){ */
/*     free(in); */
/* } */

/* /\************* Is Acceptance OK *************\/ */
/* isAcceptOK *alloc_IsAcceptOK(){ */
/*     isAcceptOK *acceptOK = (isAcceptOK *) malloc(sizeof(isAcceptOK)); */
/*     if(acceptOK == NULL){ */
/* 	fprintf(stderr, "\n[in: alloc.c->alloc_IsAcceptOK]\nNo memory left for creating acceptOK. Exiting.\n"); */
/* 	exit(1); */
/*     } */

/*     acceptOK->IsAccOK_beta = gsl_matrix_calloc(2,2); */

/*     acceptOK->IsAccOK_betaWardOut=0.0; */
/*     acceptOK->IsAccOK_betaOutOut=0.0; */
/*     acceptOK->IsAccOK_mu=0.0; */
/*     acceptOK->IsAccOK_sigma=0.0; */
/*     acceptOK->IsAccOK_nu1=0.0; */
/*     acceptOK->IsAccOK_kappa=0.0; */
   
/*     return acceptOK; */
/* } */





/* void freeIsAcceptOK(isAcceptOK *acceptOK){ */
/*     gsl_matrix_free(acceptOK->IsAccOK_beta); */

/*     free(acceptOK); */
/* } */








/* /\************ NbProposals ***************\/ */
/* NbProposals *alloc_NbProposals(){ */
/*     NbProposals *NbProp = (NbProposals *) malloc(sizeof(NbProposals)); */
/*     if(NbProp == NULL){ */
/* 	fprintf(stderr, "\n[in: alloc.c->alloc_NbProposals]\nNo memory left for creating NbProp. Exiting.\n"); */
/* 	exit(1); */
/*     } */

/*     NbProp->NbProp_beta = gsl_matrix_calloc(2,2); */

/*     NbProp->NbProp_betaWardOut=0.0; */
/*     NbProp->NbProp_betaOutOut=0.0; */
/*     NbProp->NbProp_mu=0.0; */
/*     NbProp->NbProp_sigma=0.0; */
/*     NbProp->NbProp_nu1=0.0; */
/*     NbProp->NbProp_kappa=0.0; */

/*     return NbProp; */
/* } */





/* void reInitiateNbProp(NbProposals * NbProp){ */
/*     int i,j; */

/*     for(i=0 ; i<2 ; i++) */
/* 	{ */
/* 	    for(j=0 ; j<2 ; j++) */
/* 		{ */
/* 		    gsl_matrix_set(NbProp->NbProp_beta,i,j,0.0); */
/* 		} */
/* 	} */

/*     NbProp->NbProp_betaWardOut=0.0; */
/*     NbProp->NbProp_betaOutOut=0.0; */
/*     NbProp->NbProp_mu=0.0; */
/*     NbProp->NbProp_sigma=0.0; */
/*     NbProp->NbProp_nu1=0.0; */
/*     NbProp->NbProp_kappa=0.0; */

/* } */





/* void freeNbProposals(NbProposals *NbProp){ */
/*     gsl_matrix_free(NbProp->NbProp_beta); */

/*     free(NbProp); */
/* } */







/* /\************** OUTPUT FILES **************\/ */
/* output_files *alloc_FILES(char *workspace){ */
/*     output_files *fich = (output_files *) malloc(sizeof(output_files)); */
/*     if(fich == NULL){ */
/* 	fprintf(stderr, "\n[in: alloc.c->alloc_FILES]\nNo memory left for creating fich. Exiting.\n"); */
/* 	exit(1); */
/*     } */

/*     char fileName[300]; */

/*     strcpy(fileName, workspace); */
/*     strcat(fileName,"LogL.txt"); */
/*     fich->LogL = fopen(fileName,"w"); */
/*     if ( fich->LogL == NULL ) */
/* 	{ */
/* 	    printf("A problem occurred while opening a file. Check whether a file exists and is already opened.\n"); */
/* 	    exit(1); */
/* 	} */

/*     strcpy(fileName, workspace); */
/*     strcat(fileName,"ColonDates.txt"); */
/*     fich->ColonDates = fopen(fileName,"w"); */
/*     if ( fich->ColonDates == NULL ) */
/* 	{ */
/* 	    printf("A problem occurred while opening a file. Check whether a file exists and is already opened.\n"); */
/* 	    exit(1); */
/* 	} */

/*     strcpy(fileName, workspace); */
/*     strcat(fileName,"EndColonDates.txt"); */
/*     fich->EndColonDates = fopen(fileName,"w"); */
/*     if ( fich->EndColonDates == NULL ) */
/* 	{ */
/* 	    printf("A problem occurred while opening a file. Check whether a file exists and is already opened.\n"); */
/* 	    exit(1); */
/* 	} */

/*     strcpy(fileName, workspace); */
/*     strcat(fileName,"Parameters.txt"); */
/*     fich->Parameters = fopen(fileName,"w"); */
/*     if (fich->Parameters == NULL ) */
/* 	{ */
/* 	    printf("A problem occurred while opening a file. Check whether a file exists and is already opened.\n"); */
/* 	    exit(1); */
/* 	} */

/*     return fich; */
/* } */






/* void freeFILES(output_files *fich){ */
/*     fclose(fich->LogL); */
/*     fclose(fich->ColonDates); */
/*     fclose(fich->EndColonDates); */
/*     fclose(fich->Parameters); */

/*     free(fich); */
/* } */



/* to check */



/* /\************ MCMC internals **************\/ */
/* mcmcInternals *alloc_McmcInternals(){ */
/*     mcmcInternals *MCMCSettings = (mcmcInternals *) malloc(sizeof(mcmcInternals)); */
/*     if(MCMCSettings == NULL){ */
/* 	fprintf(stderr, "\n[in: alloc.c->alloc_McmcInternals]\nNo memory left for creating MCMCSettings. Exiting.\n"); */
/* 	exit(1); */
/*     } */

/*     MCMCSettings->Sigma_beta = gsl_matrix_calloc(2,2); */
/*     if(MCMCSettings->Sigma_beta == NULL){ */
/* 	fprintf(stderr, "\n[in: alloc.c->alloc_McmcInternals]\nNo memory left for creating MCMCSettings. Exiting.\n"); */
/* 	exit(1); */
/*     } */

/*     MCMCSettings->NbSimul = 10000; */
/*     MCMCSettings->SubSample = 200; */
/*     MCMCSettings->BurnIn = 5000; */
/*     MCMCSettings->Sigma_betaWardOut = 0.1; */
/*     MCMCSettings->Sigma_betaOutOut = 0.1; */
/*     MCMCSettings->Sigma_mu = 0.1; */
/*     MCMCSettings->Sigma_sigma = 0.1; */
/*     MCMCSettings->Sigma_nu1 = 0.1; */
/*     MCMCSettings->Sigma_kappa = 0.1; */
/*     MCMCSettings->Sigma_tau = 0.1; */
/*     MCMCSettings->Sigma_alpha = 0.1; */

/*     return MCMCSettings; */
/* } */





/* void printStdProp(mcmcInternals *MCMCSettings){ */
/*     int i,j; */

/*     for (i=0;i<2;i++){ */
/* 	for (j=0;j<2;j++){ */
/* 	    printf("Std proposal for beta_%d,%d: %lg\n",i,j,gsl_matrix_get(MCMCSettings->Sigma_beta,i,j)); */
/* 	} */
/*     } */
/*     printf("Std proposal for betaWardOut: %lg\n",MCMCSettings->Sigma_betaWardOut); */
/*     printf("Std proposal for betaOutOut: %lg\n",MCMCSettings->Sigma_betaOutOut); */
/*     printf("Std proposal for mu: %lg\n",MCMCSettings->Sigma_mu); */
/*     printf("Std proposal for sigma: %lg\n",MCMCSettings->Sigma_sigma); */
/*     printf("Std proposal for nu1: %lg\n",MCMCSettings->Sigma_nu1); */
/*     printf("Std proposal for kappa: %lg\n",MCMCSettings->Sigma_kappa); */
/*     printf("Std proposal for tau: %lg\n",MCMCSettings->Sigma_tau); */
/*     printf("Std proposal for alpha: %lg\n",MCMCSettings->Sigma_alpha); */

/*     fflush(stdout); */
/* } */





/* void freeMcmcInternals(mcmcInternals *MCMCSettings){ */
/*     gsl_matrix_free(MCMCSettings->Sigma_beta); */

/*     free(MCMCSettings); */
/* } */







/* /\************ Acceptance ***************\/ */
/* acceptance *alloc_Acceptance(){ */
/*     acceptance *accept = (acceptance *) malloc(sizeof(acceptance)); */
/*     if(accept == NULL){ */
/* 	fprintf(stderr, "\n[in: alloc.c->alloc_Acceptance]\nNo memory left for creating accept. Exiting.\n"); */
/* 	exit(1); */
/*     } */

/*     accept->PourcAcc_beta = gsl_matrix_calloc(2,2); */
/*     accept->PourcAcc_betaWardOut=0.0; */
/*     accept->PourcAcc_betaOutOut=0.0; */
/*     accept->PourcAcc_mu=0.0; */
/*     accept->PourcAcc_sigma=0.0; */
/*     accept->PourcAcc_nu1=0.0; */
/*     accept->PourcAcc_kappa=0.0; */

/*     return accept; */
/* } */





/* void reInitiateAcceptance(acceptance *accept){ */
/*     int i,j; */
/*     for(i=0;i<2;i++){ */
/* 	for(j=0;j<2;j++){ */
/* 	    gsl_matrix_set(accept->PourcAcc_beta,i,j,0.0); */
/* 	} */
/*     } */

/*     accept->PourcAcc_betaWardOut=0.0; */
/*     accept->PourcAcc_betaOutOut=0.0; */
/*     accept->PourcAcc_mu=0.0; */
/*     accept->PourcAcc_sigma=0.0; */
/*     accept->PourcAcc_nu1=0.0; */
/*     accept->PourcAcc_kappa=0.0; */

/* } */





/* void printAcceptance(acceptance *accept, NbProposals *NbProp){ */
/*     int i,j; */

/*     for(i=0;i<2;i++){ */
/* 	for(j=0;j<2;j++){ */
/* 	    printf("Prob accept beta_%d,%d\t%lg\n",i,j,gsl_matrix_get(accept->PourcAcc_beta,i,j)/gsl_matrix_get(NbProp->NbProp_beta,i,j)); */
/* 	    fflush(stdout); */
/* 	} */
/*     } */

/*     printf("Prob accept betaWardOut\t%lg\n",accept->PourcAcc_betaWardOut/NbProp->NbProp_betaWardOut); */
/*     fflush(stdout); */
/*     printf("Prob accept betaOutOut\t%lg\n",accept->PourcAcc_betaOutOut/NbProp->NbProp_betaOutOut); */
/*     fflush(stdout); */
/*     printf("Prob accept mu\t%lg\n",accept->PourcAcc_mu/NbProp->NbProp_mu); */
/*     fflush(stdout); */
/*     printf("Prob accept sigma\t%lg\n",accept->PourcAcc_sigma/NbProp->NbProp_sigma); */
/*     fflush(stdout); */
/*     printf("Prob accept nu1\t%lg\n",accept->PourcAcc_nu1/NbProp->NbProp_nu1); */
/*     fflush(stdout); */
/*     printf("Prob accept kappa\t%lg\n",accept->PourcAcc_kappa/NbProp->NbProp_kappa); */
/*     fflush(stdout); */

/* } */





/* void freeAcceptance(acceptance *accept){ */
/*     gsl_matrix_free(accept->PourcAcc_beta); */
/*     free(accept); */
/* } */







/*
  ======================
  >>>>> TESTING <<<<<
  ======================
*/

int main(){
    data * dat = alloc_data(10,100);
    print_data(dat);
    free_data(dat);

    return 0;
}



/* 
   gcc instructions:

   gcc -o structures matvec.c genclasses.c structures.c -lgsl -lgslcblas -g

  ./structures

   valgrind --leak-check=full --track-origins=yes -v structures

*/
