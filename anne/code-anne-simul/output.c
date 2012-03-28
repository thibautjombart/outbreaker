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

void writeData(int NbCases, char* workspace, nb_data *nbData, raw_data *data)
{
	char file[200];
	FILE *Ward, *Admission, *Discharge, *ColonDate, *EndColonDate, *PosSwabDates, *NegSwabDates, *NbAdmissions, *NbPosSwabs, *NbNegSwabs;
	int i,k;

	strcpy(file, workspace);
	strcat(file,"Ward.txt");
	if ((Ward=fopen(file,"w"))==NULL)
	{
		printf("Cannot open Ward.txt");
		exit(2);
	}
	for (i=0 ; i<NbCases ; i++)
	{
	  	fprintf(Ward,"%d\n",data->ward[i]);
	   	fflush(Ward);
	}

	strcpy(file, workspace);
	strcat(file,"Admission.txt");
	if ((Admission=fopen(file,"w"))==NULL)
	{
		printf("Cannot open Admission.txt");
		exit(2);
	}
	for (i=0 ; i<NbCases ; i++)
	{
		for(k=0 ; k<nbData->NbAdmissions[i] ; k++)
		{
			fprintf(Admission,"%d\t",(int)gsl_vector_get(data->A[i],k));
			fflush(Admission);
		}
		fprintf(Admission,"\n");
		fflush(Admission);
	}

	strcpy(file, workspace);
	strcat(file,"Discharge.txt");
	if ((Discharge=fopen(file,"w"))==NULL)
	{
		printf("Cannot open Discharge.txt");
		exit(2);
	}
	for (i=0 ; i<NbCases ; i++)
	{
	  	fprintf(Discharge,"%d\n",(int)gsl_vector_get(data->D[i],0));
	   	fflush(Discharge);
	}

	strcpy(file, workspace);
	strcat(file,"ColonDate.txt");
	if ((ColonDate=fopen(file,"w"))==NULL)
	{
		printf("Cannot open ColonDate.txt");
		exit(2);
	}
	for (i=0 ; i<NbCases ; i++)
	{
	  	fprintf(ColonDate,"%d\n",data->C[i]);
	   	fflush(ColonDate);
	}

	strcpy(file, workspace);
	strcat(file,"EndColonDate.txt");
	if ((EndColonDate=fopen(file,"w"))==NULL)
	{
		printf("Cannot open EndColonDate.txt");
		exit(2);
	}
	for (i=0 ; i<NbCases ; i++)
	{
	  	fprintf(EndColonDate,"%d\n",data->E[i]);
	   	fflush(EndColonDate);
	}


	strcpy(file, workspace);
	strcat(file,"NbAdmissions.txt");
	if ((NbAdmissions=fopen(file,"w"))==NULL)
	{
		printf("Cannot open NbAdmissions.txt");
		exit(2);
	}
	for (i=0 ; i<NbCases ; i++)
	{
	  	fprintf(NbAdmissions,"%d\n",nbData->NbAdmissions[i]);
	   	fflush(NbAdmissions);
	}

	strcpy(file, workspace);
	strcat(file,"nbPos.txt");
	if ((NbPosSwabs=fopen(file,"w"))==NULL)
	{
		printf("Cannot open nbPos.txt");
		exit(2);
	}
	for (i=0 ; i<NbCases ; i++)
	{
	  	fprintf(NbPosSwabs,"%d\n",nbData->NbPosSwabs[i]);
	   	fflush(NbPosSwabs);
	}

	strcpy(file, workspace);
	strcat(file,"nbNeg.txt");
	if ((NbNegSwabs=fopen(file,"w"))==NULL)
	{
		printf("Cannot open nbNeg.txt");
		exit(2);
	}
	for (i=0 ; i<NbCases ; i++)
	{
	  	fprintf(NbNegSwabs,"%d\n",nbData->NbNegSwabs[i]);
	   	fflush(NbNegSwabs);
	}

	strcpy(file, workspace);
	strcat(file,"PositiveSwabDates.txt");
	if ((PosSwabDates=fopen(file,"w"))==NULL)
	{
		printf("Cannot open PositiveSwabDates.txt");
		exit(2);
	}
	for (i=0 ; i<NbCases ; i++)
	{
		if(nbData->NbPosSwabs[i]>0)
		{
			for(k=0 ; k<nbData->NbPosSwabs[i] ; k++)
			{
				fprintf(PosSwabDates,"%d\t",(int)gsl_vector_get(data->P[i],k));
				fflush(PosSwabDates);
			}
		}
		fprintf(PosSwabDates,"\n");
		fflush(PosSwabDates);
	}

	strcpy(file, workspace);
	strcat(file,"NegativeSwabDates.txt");
	if ((NegSwabDates=fopen(file,"w"))==NULL)
	{
		printf("Cannot open NegativeSwabDates.txt");
		exit(2);
	}
	for (i=0 ; i<NbCases ; i++)
	{
		if(nbData->NbNegSwabs[i]>0)
		{
			for(k=0 ; k<nbData->NbNegSwabs[i] ; k++)
			{
				fprintf(NegSwabDates,"%d\t",(int)gsl_vector_get(data->N[i],k));
				fflush(NegSwabDates);
			}
		}
		fprintf(NegSwabDates,"\n");
		fflush(NegSwabDates);
	}

	fclose(Ward);
	fclose(Admission);
	fclose(Discharge);
	fclose(ColonDate);
	fclose(EndColonDate);
	fclose(PosSwabDates);
	fclose(NegSwabDates);
	fclose(NbAdmissions);
	fclose(NbPosSwabs);
	fclose(NbNegSwabs);
}

extern gsl_rng * rng;


