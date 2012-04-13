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
/* Reading data                                                               */
/******************************************************************************/

/*
  Fill in a nb_data object using input form R.
*/
void importNbData(int *nbAdmVec, int *nbPosSwab, int *nbNegSwab, int *nbColPatients, int *nbPatients, int *duration, int *idxColPatients, int *nbSeqPat, int *nbSeq, nb_data *nb){
    int i;

    /* fill in simple values */
    nb->NbPatients = *nbPatients;
    nb->T = *duration;
    nb->NbColonisedPatients = *nbColPatients;
    nb->NbSequences = *nbSeq;

    /* fill in vectors values */
    for(i=0;i<*nbPatients;i++){
	nb->NbAdmissions[i] = nbAdmVec[i];
	nb->NbPosSwabs[i] = nbPosSwab[i];
	nb->NbNegSwabs[i] = nbNegSwab[i];
	nb->M[i] = nbSeqPat[i];
    }

    for(i=0;i<*nbColPatients;i++){
	nb->indexColonisedPatients[i] = idxColPatients[i];
    }

} /* end importNbData*/







/*
  Fill in a nb_data object using input form R.
  == legend (all lists in R are given patient-wise) ==
  - tAdmVec: result of unlisting of patients dates of admissions
  - tDisVec: result of unlisting of patients dates of discharge
  - tPosSwab: result of unlisting of patients dates of positive swabs
  - tNegSwab: result of unlisting of patients dates of positive swabs
  - tNegSwab: result of unlisting of hospital presence data
*/
void importRawData(int *wardVec, int *tAdmVec, int *tDisVec, int *tPosSwab, int *tNegSwab, int *hospPres, int *idxSeqVec, int *totNbSeq, double *tCollecVec, nb_data *nb, raw_data *data){
    int i, k, counter;

    /* FILL IN SIMPLE VALUES */
    data->NbPatients = nb->NbPatients;
    data->T = nb->T;

    /* FILL IN VECTORS / ARRAYS */
    /* wards */
    for(i=0;i<nb->NbPatients;i++){
	data->ward[i] = wardVec[i];
    }

    /* admission/discharge times */
    counter = 0;
    for(i=0;i<nb->NbPatients;i++){
	for(k=0;k<nb->NbAdmissions[i];k++){
	    /* data->A[i][k] = tAdmVec[counter]; */
	    /* data->D[i][k] = tDisVec[counter++]; */
	    gsl_vector_set(data->A[i], k, tAdmVec[counter]);
	    gsl_vector_set(data->D[i], k, tDisVec[counter++]);
	}
    }

    /* positive swab times */
    counter = 0;
    for(i=0;i<nb->NbPatients;i++){
	for(k=0;k<nb->NbPosSwabs[i];k++){
	    /* data->P[i][k] = tPosSwab[counter++]; */
	    gsl_vector_set(data->P[i], k, tPosSwab[counter++]);
	}
    }

    /* negative swab times */
    counter = 0;
    for(i=0;i<nb->NbPatients;i++){
	for(k=0;k<nb->NbNegSwabs[i];k++){
	    /* data->N[i][k] = tNegSwab[counter++]; */
	    gsl_vector_set(data->N[i], k, tNegSwab[counter++]);
	}
    }

    /* hospital presence across time */
    counter = 0;
    for(i=0;i<nb->NbPatients;i++){
	for(k=0;k<nb->T;k++){
	    /* data->IsInHosp[i][k] = hospPres[counter++]; */
	    gsl_vector_set(data->IsInHosp[i], k, hospPres[counter++]);
	}
    }
 
    /* indices of DNA sequences for each patient */
    counter = 0;
    for(i=0;i<nb->NbPatients;i++){
	for(k=0;k<nb->M[i];k++){
	    data->S[i][k] = idxSeqVec[counter++];
	}
    }

    /* collection dates of DNA sequences */
    for(i=0;i<*totNbSeq;i++){
	data->Tcollec[i] = tCollecVec[i];
    }

    /* nb of sequences per patient */
    for(i=0;i<nb->NbPatients;i++){
	data->M[i] = nb->M[i];
    }

} /* end importRawData*/















/*
  Old text-file based version.
*/
void readFakeNbData(nb_data *nb){
    int i, NbPatients = nb->NbPatients;
    FILE *fich;
    int V;

    for(i=0;i<NbPatients;i++){
	nb->NbAdmissions[i]=1;
    }

    fich = fopen("nbNeg.txt","r");
    if(fich == NULL){
	printf("A problem occurred while opening nbNeg.txt. Check that the file exists and is not opened.\n");
	fflush(stdout);
	exit(1);
    }

    for (i=0;i<NbPatients;i++){
	if( fscanf(fich,"%d",&V) != 1)
	    {
		printf("A problem occurred while reading the file\n");
		fflush(stdout);
		break;
	    }
	nb->NbNegSwabs[i]=V;
    }
    fclose(fich);

    fich = fopen("nbPos.txt","r");
    if ( fich == NULL ){
	printf("A problem occurred while opening nbPos.txt. Check that the file exists and is not opened.\n");
	fflush(stdout);
	exit(1);
    }

    for (i=0;i<NbPatients;i++){
	if( fscanf(fich,"%d",&V) != 1)
	    {
		printf("A problem occurred while reading the file\n");
		fflush(stdout);
		break;
	    }
	nb->NbPosSwabs[i]=V;
    }

    fclose(fich);
}







/*
  Old text-file based version.
*/
void readFakeData(nb_data *nb, raw_data *data){
    FILE *fich;
    int i, k, NbPatients = nb->NbPatients;
    int V;

    fich = fopen("Admission.txt","r");
    if ( fich == NULL ){
	printf("A problem occurred while opening Admission.txt. Check that the file exists and is not opened.\n");
	fflush(stdout);
	exit(1);
    }

    for (i=0;i<NbPatients;i++){
	for (k=0 ; k<nb->NbAdmissions[i] ; k++)
	    {
		if( fscanf(fich,"%d",&V) != 1)
		    {
			printf("A problem occurred while reading the file\n");
			fflush(stdout);
			break;
		    }
		gsl_vector_set(data->A[i],k,V);
	    }
    }
    fclose(fich);

    fich = fopen("Discharge.txt","r");
    if ( fich == NULL ){
	printf("A problem occurred while opening Discharge.txt. Check that the file exists and is not opened.\n");
	fflush(stdout);
	exit(1);
    }

    for (i=0;i<NbPatients;i++){
	for (k=0 ; k<nb->NbAdmissions[i] ; k++)
	    {
		if( fscanf(fich,"%d",&V) != 1)
		    {
			printf("A problem occurred while reading the file\n");
			fflush(stdout);
			break;
		    }
		gsl_vector_set(data->D[i],k,V);
	    }
    }
    fclose(fich);

    fich = fopen("PositiveSwabDates.txt","r");
    if ( fich == NULL ){
	printf("A problem occurred while opening PositiveSwabDates.txt. Check that the file exists and is not opened.\n");
	fflush(stdout);
	exit(1);
    }

    for (i=0;i<NbPatients;i++){
	for (k=0 ; k<nb->NbPosSwabs[i] ; k++)
	    {
		if( fscanf(fich,"%d",&V) != 1)
		    {
			printf("A problem occurred while reading the file\n");
			fflush(stdout);
			break;
		    }
		gsl_vector_set(data->P[i],k,V);
	    }
    }
    fclose(fich);

    fich = fopen("NegativeSwabDates.txt","r");
    if ( fich == NULL ){
	printf("A problem occurred while opening NegativeSwabDates.txt. Check that the file exists and is not opened.\n");
	fflush(stdout);
	exit(1);
    }

    for (i=0;i<NbPatients;i++){
	for (k=0 ; k<nb->NbNegSwabs[i] ; k++)
	    {
		if( fscanf(fich,"%d",&V) != 1)
		    {
			printf("A problem occurred while reading the file\n");
			fflush(stdout);
			break;
		    }
		gsl_vector_set(data->N[i],k,V);
	    }
    }
    fclose(fich);

    fich = fopen("Ward.txt","r");
    if ( fich == NULL ){
	printf("A problem occurred while opening Ward.txt. Check that the file exists and is not opened.\n");
	fflush(stdout);
	exit(1);
    }

    for (i=0;i<NbPatients;i++){
	if( fscanf(fich,"%d",&V) != 1)
	    {
		printf("A problem occurred while reading the file\n");
		fflush(stdout);
		break;
	    }
	data->ward[i]=V;
    }
    fclose(fich);

    for (i=0;i<NbPatients;i++){
	data->PatientIndex[i]=i;
    }

    CalculIsInHosp(nb, data);

}







/******************************************************************************/
/* Function to read initial parameters values                                 */
/******************************************************************************/

/*************************************************************************/
/* preparing output file                                                 */
/*************************************************************************/

void prepAllFiles(output_files * Files, int NbPatients){
    int i;

    fprintf(Files->LogL,"LogLikelihood\n");
    fflush(Files->LogL);

    fprintf(Files->Parameters,"beta[0,0]\t");
    fprintf(Files->Parameters,"beta[0,1]\t");
    fprintf(Files->Parameters,"beta[1,0]\t");
    fprintf(Files->Parameters,"beta[1,1]\t");
    fprintf(Files->Parameters,"betaWardOut\t");
    fprintf(Files->Parameters,"betaOutOut\t");
    /* fprintf(Files->Parameters,"Sp\t"); */
    fprintf(Files->Parameters,"Se\t");
    fprintf(Files->Parameters,"Pi\t");
    fprintf(Files->Parameters,"mu\t");
    fprintf(Files->Parameters,"sigma\t");
    fprintf(Files->Parameters,"nu1\t");
    fprintf(Files->Parameters,"nu2\t");
    fprintf(Files->Parameters,"tau\t");
    fprintf(Files->Parameters,"alpha\n");
    fflush(Files->Parameters);

    for(i=0;i<NbPatients;i++){
	fprintf(Files->ColonDates,"C[%d]\t",i);
	fprintf(Files->EndColonDates,"E[%d]\t",i);
    }
    fprintf(Files->ColonDates,"\n");
    fprintf(Files->EndColonDates,"\n");
    fflush(Files->ColonDates);
    fflush(Files->EndColonDates);

}







/*************************************************************************/
/* writing output files                                                  */
/*************************************************************************/

void writeAllFiles(output_files * Files, parameters * param, nb_data *nb, raw_data * data, aug_data *augData, dna_dist *dnainfo){
    /* Writing results in the output files */

    double L;
    int i;

    L = fullLoglikelihoodWithPrior(data, nb, augData, dnainfo, param);
    fprintf(Files->LogL,"%lf\n",L);
    fflush(Files->LogL);

    fprintf(Files->Parameters,"%lf\t",gsl_matrix_get(param->beta,0,0));
    fprintf(Files->Parameters,"%lf\t",gsl_matrix_get(param->beta,0,1));
    fprintf(Files->Parameters,"%lf\t",gsl_matrix_get(param->beta,1,0));
    fprintf(Files->Parameters,"%lf\t",gsl_matrix_get(param->beta,1,1));
    fprintf(Files->Parameters,"%lf\t",param->betaWardOut);
    fprintf(Files->Parameters,"%lf\t",param->betaOutOut);
    /* fprintf(Files->Parameters,"%lf\t",param->Sp); */
    fprintf(Files->Parameters,"%lf\t",param->Se);
    fprintf(Files->Parameters,"%lf\t",param->Pi);
    fprintf(Files->Parameters,"%lf\t",param->mu);
    fprintf(Files->Parameters,"%lf\t",param->sigma);
    fprintf(Files->Parameters,"%lf\t",param->nu1);
    fprintf(Files->Parameters,"%lf\t",param->nu2);
    fprintf(Files->Parameters,"%lf\t",param->tau);
    fprintf(Files->Parameters,"%lf\n",param->alpha);
    fflush(Files->Parameters);

    for(i=0 ; i<nb->NbPatients ; i++){
	if(nb->NbPosSwabs[i]>0)
	    {
		fprintf(Files->ColonDates,"%d\t",augData->C[i]);
		fprintf(Files->EndColonDates,"%d\t",augData->E[i]);
	    }else
	    {
		fprintf(Files->ColonDates,"NA\t");
		fprintf(Files->EndColonDates,"NA\t");
	    }
    }
    fprintf(Files->ColonDates,"\n");
    fprintf(Files->EndColonDates,"\n");
    fflush(Files->ColonDates);
    fflush(Files->EndColonDates);

}

