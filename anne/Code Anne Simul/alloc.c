#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_eigen.h>

#include <omp.h>

#include "common.h"
#include "alloc.h"
#include "output.h"

extern gsl_rng * rng;

/**************** nb_data ****************/

nb_data * createNbData()
{
	int i;
	nb_data *nb = (nb_data *) malloc(sizeof(nb_data));

	nb->NbAdmissions = (int *) calloc(NbPatientsMax, sizeof(int));
	nb->NbPosSwabs = (int *) calloc(NbPatientsMax, sizeof(int));
	nb->NbNegSwabs = (int *) calloc(NbPatientsMax, sizeof(int));

	for(i=0 ; i<NbPatientsMax ; i++)
	{
		nb->NbAdmissions[i]=0;
		nb->NbPosSwabs[i]=0;
		nb->NbNegSwabs[i]=0;
	}

	return nb;
}

void freeNbData(nb_data *nb)
{
	free(nb->NbAdmissions);
	free(nb->NbPosSwabs);
	free(nb->NbNegSwabs);

	free(nb);
}
/******************************************/

/**************** raw_data ****************/

raw_data *createRawData()
{
	raw_data *data = (raw_data *) malloc(sizeof(raw_data));

	int i,t;

	data->ward = (int *) calloc(NbPatientsMax, sizeof(int));
	data->timeSeq = (int *) calloc(NbPatientsMax, sizeof(int));

	for(i=0 ; i<NbPatientsMax ; i++)
	{
		data->A[i] = gsl_vector_calloc(Tmax);
		data->D[i] = gsl_vector_calloc(Tmax);
		data->P[i] = gsl_vector_calloc(Tmax);
		data->N[i] = gsl_vector_calloc(Tmax);
	}

	data->C = (int *) calloc(NbPatientsMax, sizeof(int));
	data->E = (int *) calloc(NbPatientsMax, sizeof(int));

	for(i=0 ; i<NbPatientsMax ; i++)
	{
		data->C[i]=100000;
		data->E[i]=100010;
	}

	data->I0 = (int *) calloc(Tmax, sizeof(int));
	data->I1 = (int *) calloc(Tmax, sizeof(int));

	for(t=0 ; t<Tmax ; t++)
	{
		data->I0[t]=0;
		data->I1[t]=0;
	}

	return data;
}

void freeRawData(raw_data *data)
{
	int i;

	free(data->ward);
	free(data->timeSeq);

	for(i=0 ; i<NbPatientsMax ; i++)
	{
		gsl_vector_free(data->A[i]);
		gsl_vector_free(data->D[i]);
		gsl_vector_free(data->P[i]);
		gsl_vector_free(data->N[i]);
	}

	free(data->C);
	free(data->E);

	free(data->I0);
	free(data->I1);

	free(data);
}
/******************************************/
/***************** param ******************/

parameters *createParam()
{
	parameters *param = (parameters *) malloc(sizeof(parameters));

	param->beta = gsl_matrix_calloc(2,2);

	return param;
}


void readParameters(char* workspace, parameters * param)
{
    FILE * paramInit;
	char val[30];
	char file[200];

	strcpy(file, workspace);
	strcat(file,"param.txt");
	if ((paramInit=fopen(file,"r"))==NULL)
	{
		printf("Cannot read param.txt");
		exit(2);
	}

	fscanf(paramInit,"%s %lf",val,gsl_matrix_ptr (param->beta,0,0));
	printf("%s %g\n",val,gsl_matrix_get(param->beta,0,0));

	fscanf(paramInit,"%s %lf",val,gsl_matrix_ptr (param->beta,0,1));
	printf("%s %g\n",val,gsl_matrix_get(param->beta,0,1));

	fscanf(paramInit,"%s %lf",val,gsl_matrix_ptr (param->beta,1,0));
	printf("%s %g\n",val,gsl_matrix_get(param->beta,1,0));

	fscanf(paramInit,"%s %lf",val,gsl_matrix_ptr (param->beta,1,1));
	printf("%s %g\n",val,gsl_matrix_get(param->beta,1,1));

	fscanf(paramInit,"%s %lf",val,&param->betaWardOut);
	printf("%s %g\n",val,param->betaWardOut);

	fscanf(paramInit,"%s %lf",val,&param->betaOutOut);
	printf("%s %g\n",val,param->betaOutOut);

	fscanf(paramInit,"%s %lf",val,&param->Sp);
	printf("%s %g\n",val,param->Sp);

	fscanf(paramInit,"%s %lf",val,&param->Se);
	printf("%s %g\n",val,param->Se);

	fscanf(paramInit,"%s %lf",val,&param->Pi);
	printf("%s %g\n",val,param->Pi);

	fscanf(paramInit,"%s %lf",val,&param->muHosp);
	printf("%s %g\n",val,param->muHosp);

	fscanf(paramInit,"%s %lf",val,&param->sigmaHosp);
	printf("%s %g\n",val,param->sigmaHosp);

	fscanf(paramInit,"%s %lf",val,&param->mu);
	printf("%s %g\n",val,param->mu);

	fscanf(paramInit,"%s %lf",val,&param->sigma);
	printf("%s %g\n",val,param->sigma);

	fscanf(paramInit,"%s %lf",val,&param->nu1);
	printf("%s %g\n",val,param->nu1);

	fscanf(paramInit,"%s %lf",val,&param->nu2);
	printf("%s %g\n",val,param->nu2);

	fscanf(paramInit,"%s %lf",val,&param->tau);
	printf("%s %g\n",val,param->tau);

	fscanf(paramInit,"%s %lf",val,&param->alpha);
	printf("%s %g\n",val,param->alpha);

    fclose(paramInit);

}

void freeParam(parameters *param)
{
	gsl_matrix_free(param->beta);

	free(param);
}

/******************************************/

