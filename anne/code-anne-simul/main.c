#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h> // for matrix & vectors operations
#include <gsl/gsl_permutation.h> //for LU decomposition
#include <gsl/gsl_sf_gamma.h> // for combinations calculation
#include <gsl/gsl_eigen.h>

#include <time.h>

#include "common.h"
#include "alloc.h"
#include "output.h"

/******************************************************************************/
/* DECLARATION OF GLOBAL VARIABLES                                            */
/******************************************************************************/

gsl_rng * rng; /* to be included as an external variable eslewhere */

/******************************************************************************/
/* MAIN                                                                       */
/******************************************************************************/

#define OPTIONS "awrh"

int SimulEpid(parameters *param, nb_data *nbData, raw_data *data)
{
	int t,k,i,u, prov;
	int NbCases=0;
	double rand;
	double lambda, prob;
	int S0, S1;
	int Incid0, Incid1;
	double probaInfectedPerWard0;
	double probaInfectedPerWard1;
	
	int *bedNb = (int *) calloc(NbPatientsMax, sizeof(int));

	int *indexInfector = (int *) calloc(NbPatientsMax, sizeof(int)); /* -1 if infected from outside ; -100 never infected */

	for(i=0 ; i<NbPatientsMax ; i++)
	{
		indexInfector[i]=-100;
	}

	int *indexI0 = calloc(SizeWard0, sizeof(int));
	int *indexS0 = calloc(SizeWard0, sizeof(int));
	int *indexI1 = calloc(SizeWard1, sizeof(int));
	int *indexS1 = calloc(SizeWard1, sizeof(int));
	int *indexIncid0 = calloc(SizeWard0, sizeof(int));
	int *indexIncid1 = calloc(SizeWard1, sizeof(int));

	gsl_matrix * NbDischarges0 = gsl_matrix_calloc(SizeWard0,Tmax);
	gsl_matrix * NbDischarges1 = gsl_matrix_calloc(SizeWard1,Tmax);

	int *TotalNbDischarges0 = calloc(Tmax, sizeof(int));
	int *TotalNbDischarges1 = calloc(Tmax, sizeof(int));

	for(t=0 ; t<Tmax ; t++)
	{
		TotalNbDischarges0[t]=0;
		TotalNbDischarges1[t]=0;
	}

	int *AlreadyColonised = calloc(NbPatientsMax, sizeof(int));

	gsl_vector * isColonised[NbPatientsMax];
	for(i=0 ; i<NbPatientsMax ; i++)
	{
		isColonised[i] = gsl_vector_calloc(Tmax);
	}

	gsl_vector * isInHosp[NbPatientsMax];
	for(i=0 ; i<NbPatientsMax ; i++)
	{
		isInHosp[i] = gsl_vector_calloc(Tmax);
	}

	gsl_matrix * indexOfPatientsInWard0 = gsl_matrix_calloc(SizeWard0,Tmax);
	gsl_matrix * indexOfPatientsInWard1 = gsl_matrix_calloc(SizeWard1,Tmax);

	/* cases present in the 2 wards on day 0*/

	for(i=0 ; i<SizeWard0 ; i++) /* patients initially present in ward 0 */
	{
		data->ward[NbCases]=0;
		bedNb[NbCases]=i;
		nbData->NbAdmissions[NbCases]=1;
		gsl_vector_set(data->A[NbCases],0,-ceil(gsl_ran_gamma (rng, (param->muHosp/2)*(param->muHosp/2)/(param->sigmaHosp*param->sigmaHosp), param->sigmaHosp*param->sigmaHosp/(param->muHosp/2))));
		do
		{
			gsl_vector_set(data->D[NbCases],0,ceil(gsl_vector_get(data->A[NbCases],0) + gsl_ran_gamma (rng, (param->muHosp)*(param->muHosp)/(param->sigmaHosp*param->sigmaHosp), param->sigmaHosp*param->sigmaHosp/(param->muHosp))));
		}while(gsl_vector_get(data->D[NbCases],0)<=0);
		if(gsl_vector_get(data->D[NbCases],0)<Tmax)
		{
			prov =gsl_matrix_get(NbDischarges0,bedNb[NbCases],gsl_vector_get(data->D[NbCases],0));
			gsl_matrix_set(NbDischarges0,bedNb[NbCases],gsl_vector_get(data->D[NbCases],0),prov+1);
			TotalNbDischarges0[(int)gsl_vector_get(data->D[NbCases],0)]++;
		}
		if(GSL_MAX(0,gsl_vector_get(data->A[NbCases],0))<GSL_MIN(Tmax,gsl_vector_get(data->D[NbCases],0)))
		{
			for(u=GSL_MAX(0,gsl_vector_get(data->A[NbCases],0)) ; u<GSL_MIN(Tmax,gsl_vector_get(data->D[NbCases],0)) ; u++)
			{
				gsl_vector_set(isInHosp[NbCases],u,1);
				gsl_matrix_set(indexOfPatientsInWard0,i,u,NbCases);
			}
		}
		rand=gsl_ran_flat (rng, 0.0, 1.0);
		if(rand<param->Pi) /* colonised before admission */
		{
			AlreadyColonised[NbCases]=1;
			indexInfector[NbCases]=-1;
			data->C[NbCases]=gsl_vector_get(data->A[NbCases],0)-1;
			data->E[NbCases]=data->C[NbCases]+ceil(gsl_ran_gamma (rng, (param->mu)*(param->mu)/(param->sigma*param->sigma), param->sigma*param->sigma/(param->mu)));
			if(GSL_MAX(0,data->C[NbCases])<GSL_MIN(Tmax,data->E[NbCases]))
			{
				for(u=GSL_MAX(0,data->C[NbCases]) ; u<GSL_MIN(Tmax,data->E[NbCases]) ; u++)
				{
					gsl_vector_set(isColonised[NbCases],u,1);
				}
			}
		}else
		{
			AlreadyColonised[NbCases]=0;
		}
		NbCases++;
		printf("NbCases\t %d\n",NbCases);
		fflush(stdout);
	}
	for(i=0 ; i<SizeWard1 ; i++) /* patients initially present in ward 1 */
	{
		data->ward[NbCases]=1;
		bedNb[NbCases]=i;
		nbData->NbAdmissions[NbCases]=1;
		gsl_vector_set(data->A[NbCases],0,-ceil(gsl_ran_gamma (rng, (param->muHosp/2)*(param->muHosp/2)/(param->sigmaHosp*param->sigmaHosp), param->sigmaHosp*param->sigmaHosp/(param->muHosp/2))));
		do
		{
			gsl_vector_set(data->D[NbCases],0,ceil(gsl_vector_get(data->A[NbCases],0) + gsl_ran_gamma (rng, (param->muHosp)*(param->muHosp)/(param->sigmaHosp*param->sigmaHosp), param->sigmaHosp*param->sigmaHosp/(param->muHosp))));
		}while(gsl_vector_get(data->D[NbCases],0)<=0);
		if(gsl_vector_get(data->D[NbCases],0)<Tmax)
		{
			prov =gsl_matrix_get(NbDischarges1,bedNb[NbCases],gsl_vector_get(data->D[NbCases],0));
			gsl_matrix_set(NbDischarges1,bedNb[NbCases],gsl_vector_get(data->D[NbCases],0),prov+1);
			TotalNbDischarges1[(int)gsl_vector_get(data->D[NbCases],0)]++;
		}
		if(GSL_MAX(0,gsl_vector_get(data->A[NbCases],0))<GSL_MIN(Tmax,gsl_vector_get(data->D[NbCases],0)))
		{
			for(u=GSL_MAX(0,gsl_vector_get(data->A[NbCases],0)) ; u<GSL_MIN(Tmax,gsl_vector_get(data->D[NbCases],0)) ; u++)
			{
				gsl_vector_set(isInHosp[NbCases],u,1);
				gsl_matrix_set(indexOfPatientsInWard1,i,u,NbCases);
			}
		}
		rand=gsl_ran_flat (rng, 0.0, 1.0);
		if(rand<param->Pi) /* colonised before admission */
		{
			AlreadyColonised[NbCases]=1;
			indexInfector[NbCases]=-1;
			data->C[NbCases]=gsl_vector_get(data->A[NbCases],0)-1;
			data->E[NbCases]=data->C[NbCases]+ceil(gsl_ran_gamma (rng, (param->mu)*(param->mu)/(param->sigma*param->sigma), param->sigma*param->sigma/(param->mu)));
			if(GSL_MAX(0,data->C[NbCases])<GSL_MIN(Tmax,data->E[NbCases]))
			{
				for(u=GSL_MAX(0,data->C[NbCases]) ; u<GSL_MIN(Tmax,data->E[NbCases]) ; u++)
				{
					gsl_vector_set(isColonised[NbCases],u,1);
				}
			}
		}else
		{
			AlreadyColonised[NbCases]=0;
		}
		NbCases++;
		printf("NbCases\t %d\n",NbCases);
		fflush(stdout);
	}

	/* simulate the admissions/discharges up to Tmax */

	for(t=0 ; t<Tmax ; t++)
	{
		for(i=0 ; i<SizeWard0 ; i++)
		{
			if(gsl_matrix_get(NbDischarges0,i,t)>0)
			{
				data->ward[NbCases]=0;
				bedNb[NbCases]=i;
				nbData->NbAdmissions[NbCases]=1;
				gsl_vector_set(data->A[NbCases],0,t);
				gsl_vector_set(data->D[NbCases],0,ceil(gsl_vector_get(data->A[NbCases],0) + gsl_ran_gamma (rng, (param->muHosp)*(param->muHosp)/(param->sigmaHosp*param->sigmaHosp), param->sigmaHosp*param->sigmaHosp/(param->muHosp))));
				if(gsl_vector_get(data->D[NbCases],0)<Tmax)
				{
					prov =gsl_matrix_get(NbDischarges0,bedNb[NbCases],gsl_vector_get(data->D[NbCases],0));
					gsl_matrix_set(NbDischarges0,bedNb[NbCases],gsl_vector_get(data->D[NbCases],0),prov+1);
					TotalNbDischarges0[(int)gsl_vector_get(data->D[NbCases],0)]++;
				}
				if(GSL_MAX(0,gsl_vector_get(data->A[NbCases],0))<GSL_MIN(Tmax,gsl_vector_get(data->D[NbCases],0)))
				{
					for(u=GSL_MAX(0,gsl_vector_get(data->A[NbCases],0)) ; u<GSL_MIN(Tmax,gsl_vector_get(data->D[NbCases],0)) ; u++)
					{
						gsl_vector_set(isInHosp[NbCases],u,1);
						gsl_matrix_set(indexOfPatientsInWard0,i,u,NbCases);
					}
				}
				rand=gsl_ran_flat (rng, 0.0, 1.0);
				if(rand<param->Pi) /* colonised before admission */
				{
					AlreadyColonised[NbCases]=1;
					data->C[NbCases]=gsl_vector_get(data->A[NbCases],0)-1;
					data->E[NbCases]=data->C[NbCases]+ceil(gsl_ran_gamma (rng, (param->mu)*(param->mu)/(param->sigma*param->sigma), param->sigma*param->sigma/(param->mu)));
					if(GSL_MAX(0,data->C[NbCases])<GSL_MIN(Tmax,data->E[NbCases]))
					{
						for(u=GSL_MAX(0,data->C[NbCases]) ; u<GSL_MIN(Tmax,data->E[NbCases]) ; u++)
						{
							gsl_vector_set(isColonised[NbCases],u,1);
						}
					}
				}else
				{
					AlreadyColonised[NbCases]=0;
				}
				NbCases++;
				printf("NbCases\t %d\n",NbCases);
				fflush(stdout);
			}
		}

		for(i=0 ; i<SizeWard1 ; i++)
		{
			if(gsl_matrix_get(NbDischarges1,i,t)>0)
			{
				data->ward[NbCases]=1;
				bedNb[NbCases]=i;
				nbData->NbAdmissions[NbCases]=1;
				gsl_vector_set(data->A[NbCases],0,t);
				gsl_vector_set(data->D[NbCases],0,ceil(gsl_vector_get(data->A[NbCases],0) + gsl_ran_gamma (rng, (param->muHosp)*(param->muHosp)/(param->sigmaHosp*param->sigmaHosp), param->sigmaHosp*param->sigmaHosp/(param->muHosp))));
				if(gsl_vector_get(data->D[NbCases],0)<Tmax)
				{
					prov =gsl_matrix_get(NbDischarges1,bedNb[NbCases],gsl_vector_get(data->D[NbCases],0));
					gsl_matrix_set(NbDischarges1,bedNb[NbCases],gsl_vector_get(data->D[NbCases],0),prov+1);
					TotalNbDischarges1[(int)gsl_vector_get(data->D[NbCases],0)]++;
				}
				if(GSL_MAX(0,gsl_vector_get(data->A[NbCases],0))<GSL_MIN(Tmax,gsl_vector_get(data->D[NbCases],0)))
				{
					for(u=GSL_MAX(0,gsl_vector_get(data->A[NbCases],0)) ; u<GSL_MIN(Tmax,gsl_vector_get(data->D[NbCases],0)) ; u++)
					{
						gsl_vector_set(isInHosp[NbCases],u,1);
						gsl_matrix_set(indexOfPatientsInWard1,i,u,NbCases);
					}
				}
				rand=gsl_ran_flat (rng, 0.0, 1.0);
				if(rand<param->Pi) /* colonised before admission */
				{
					AlreadyColonised[NbCases]=1;
					data->C[NbCases]=gsl_vector_get(data->A[NbCases],0)-1;
					data->E[NbCases]=data->C[NbCases]+ceil(gsl_ran_gamma (rng, (param->mu)*(param->mu)/(param->sigma*param->sigma), param->sigma*param->sigma/(param->mu)));
					if(GSL_MAX(0,data->C[NbCases])<GSL_MIN(Tmax,data->E[NbCases]))
					{
						for(u=GSL_MAX(0,data->C[NbCases]) ; u<GSL_MIN(Tmax,data->E[NbCases]) ; u++)
						{
							gsl_vector_set(isColonised[NbCases],u,1);
						}
					}
				}else
				{
					AlreadyColonised[NbCases]=0;
				}
				NbCases++;
				printf("NbCases\t %d\n",NbCases);
				fflush(stdout);
			}
		}
	}

	/* simulate the colonisation events up to the end (for those not colonised at admission) */

	for(t=1 ; t<Tmax ; t++)
	{
		S0=0;
		S1=0;

		for(k=0 ; k<SizeWard0 ; k++)
		{
			i=gsl_matrix_get(indexOfPatientsInWard0,k,t);
			if(gsl_vector_get(isColonised[i],t-1)==1)
			{
				indexI0[data->I0[t-1]]=i;
				data->I0[t-1]++;
			}else
			{
				if(AlreadyColonised[i]==0)
				{
					indexS0[S0]=i;
					S0++;
				}
			}
		}

		for(k=0 ; k<SizeWard1 ; k++)
		{
			i=gsl_matrix_get(indexOfPatientsInWard1,k,t);
			if(gsl_vector_get(isColonised[i],t-1)==1)
			{
				indexI1[data->I1[t-1]]=i;
				data->I1[t-1]++;
			}else
			{
				if(AlreadyColonised[i]==0)
				{
					indexS1[S1]=i;
					S1++;
				}
			}
		}

		/* Total number of individuals infected on day t-1 in ward 0 */

		lambda = gsl_matrix_get(param->beta,0,0)*data->I0[t-1] + gsl_matrix_get(param->beta,0,1)*data->I1[t-1] + param->betaWardOut;
		prob = 1-exp(-lambda);
		Incid0=gsl_ran_binomial (rng, prob, S0);

		if(Incid0>0)
		{
			/* Choose them */
			gsl_ran_choose (rng, indexIncid0, Incid0, indexS0, S0, sizeof(int));

			/* draw their infector */

			probaInfectedPerWard0 = gsl_matrix_get(param->beta,0,0)*data->I0[t-1]/lambda;
			probaInfectedPerWard1 = gsl_matrix_get(param->beta,0,1)*data->I1[t-1]/lambda;

			for(i=0 ; i<Incid0 ; i++)
			{
				rand = gsl_rng_uniform(rng);
				if(rand<probaInfectedPerWard0) /* infected by someone in ward 0 */
				{
					indexInfector[indexIncid0[i]]= indexI0[gsl_rng_uniform_int (rng, data->I0[t-1])];
				}else
				{
					if(rand<probaInfectedPerWard0+probaInfectedPerWard1) /* infected by someone in ward 1 */
					{
						indexInfector[indexIncid0[i]]= indexI1[gsl_rng_uniform_int (rng, data->I1[t-1])];
					}else /* infected by outside */
					{
						indexInfector[indexIncid0[i]]= -1;
					}
				}
				/**********************************************************/
				/**********************************************************/
				/**********************************************************/
				/**********************************************************/
				/********** HERE NEED TO SIMULATE THEIR SEQUENCE **********/
				/**********************************************************/
				/**********************************************************/
				/**********************************************************/
				/**********************************************************/
			}

			/* Draw there time of end of colonistaion */
			for(i=0 ; i<Incid0 ; i++)
			{
				data->C[indexIncid0[i]]=t-1;
				data->E[indexIncid0[i]]=data->C[indexIncid0[i]]+ceil(gsl_ran_gamma (rng, (param->mu)*(param->mu)/(param->sigma*param->sigma), param->sigma*param->sigma/(param->mu)));

				/* Change isColonised and AlreadyColonised for those individuals */

				AlreadyColonised[indexIncid0[i]]=1;
				if(GSL_MAX(0,data->C[indexIncid0[i]])<GSL_MIN(Tmax,data->E[indexIncid0[i]]))
				{
					for(u=GSL_MAX(0,data->C[indexIncid0[i]]) ; u<GSL_MIN(Tmax,data->E[indexIncid0[i]]) ; u++)
					{
						gsl_vector_set(isColonised[indexIncid0[i]],u,1);
					}
				}
			}

		}

		/* Total number of individuals infected on day t-1 in ward 1 */

		lambda = gsl_matrix_get(param->beta,1,0)*data->I0[t-1] + gsl_matrix_get(param->beta,1,1)*data->I1[t-1] + param->betaWardOut;
		prob = 1-exp(-lambda);
		Incid1=gsl_ran_binomial (rng, prob, S1);

		if(Incid1>0)
		{
			/* Choose them */
			gsl_ran_choose (rng, indexIncid1, Incid1, indexS1, S1, sizeof(int));

			/* draw their infector */
			probaInfectedPerWard0 = gsl_matrix_get(param->beta,1,0)*data->I0[t-1]/lambda;
			probaInfectedPerWard1 = gsl_matrix_get(param->beta,1,1)*data->I1[t-1]/lambda;

			for(i=0 ; i<Incid1 ; i++)
			{
				rand = gsl_rng_uniform(rng);
				if(rand<probaInfectedPerWard0) /* infected by someone in ward 0 */
				{
					indexInfector[indexIncid1[i]]= indexI0[gsl_rng_uniform_int (rng, data->I0[t-1])];
				}else
				{
					if(rand<probaInfectedPerWard0+probaInfectedPerWard1) /* infected by someone in ward 1 */
					{
						indexInfector[indexIncid1[i]]= indexI1[gsl_rng_uniform_int (rng, data->I1[t-1])];
					}else /* infected by outside */
					{
						indexInfector[indexIncid1[i]]= -1;
					}
				}
				/**********************************************************/
				/**********************************************************/
				/**********************************************************/
				/**********************************************************/
				/********** HERE NEED TO SIMULATE THEIR SEQUENCE **********/
				/**********************************************************/
				/**********************************************************/
				/**********************************************************/
				/**********************************************************/
			}

			/* Draw there time of end of colonistaion */
			for(i=0 ; i<Incid1 ; i++)
			{
				data->C[indexIncid1[i]]=t-1;
				data->E[indexIncid1[i]]=data->C[indexIncid1[i]]+ceil(gsl_ran_gamma (rng, (param->mu)*(param->mu)/(param->sigma*param->sigma), param->sigma*param->sigma/(param->mu)));
			}

			/* Change isColonised and AlreadyColonised for those individuals */

			AlreadyColonised[indexIncid1[i]]=1;
			if(GSL_MAX(0,data->C[indexIncid1[i]])<GSL_MIN(Tmax,data->E[indexIncid1[i]]))
			{
				for(u=GSL_MAX(0,data->C[indexIncid1[i]]) ; u<GSL_MIN(Tmax,data->E[indexIncid1[i]]) ; u++)
				{
					gsl_vector_set(isColonised[indexIncid1[i]],u,1);
				}
			}
		}

	}

	/* simulate the testing */
	/* here testing occurs at every time step */
	/* you can then subsample to simulate a less regular testing scheme */

	for(t=0 ; t<Tmax ; t++)
	{
		for(k=0 ; k<SizeWard0 ; k++)
		{
			i=gsl_matrix_get(indexOfPatientsInWard0,k,t);

			if(gsl_vector_get(isColonised[i],t)==1) /* colonised */
			{
				rand=gsl_ran_flat (rng, 0.0, 1.0);
				if(rand<param->Se) /* test positive */
				{
					gsl_vector_set(data->P[i],nbData->NbPosSwabs[i],t);
					nbData->NbPosSwabs[i]++;
				}else /* false negative test */
				{
					gsl_vector_set(data->N[i],nbData->NbNegSwabs[i],t);
					nbData->NbNegSwabs[i]++;
				}
			}else /* not colonised */
			{
				rand=gsl_ran_flat (rng, 0.0, 1.0);
				if(rand<param->Sp) /* test negative */
				{
					gsl_vector_set(data->N[i],nbData->NbNegSwabs[i],t);
					nbData->NbNegSwabs[i]++;
				}else /* false positive test */
				{
					gsl_vector_set(data->P[i],nbData->NbPosSwabs[i],t);
					nbData->NbPosSwabs[i]++;
				}
			}
		}

		for(k=0 ; k<SizeWard1 ; k++)
		{
			i=gsl_matrix_get(indexOfPatientsInWard1,k,t);
			if(gsl_vector_get(isColonised[i],t)==1) /* colonised */
			{
				rand=gsl_ran_flat (rng, 0.0, 1.0);
				if(rand<param->Se) /* test positive */
				{
					gsl_vector_set(data->P[i],nbData->NbPosSwabs[i],t);
					nbData->NbPosSwabs[i]++;
				}else /* false negative test */
				{
					gsl_vector_set(data->N[i],nbData->NbNegSwabs[i],t);
					nbData->NbNegSwabs[i]++;
				}
			}else /* not colonised */
			{
				rand=gsl_ran_flat (rng, 0.0, 1.0);
				if(rand<param->Sp) /* test negative */
				{
					gsl_vector_set(data->N[i],nbData->NbNegSwabs[i],t);
					nbData->NbNegSwabs[i]++;
				}else /* false positive test */
				{
					gsl_vector_set(data->P[i],nbData->NbPosSwabs[i],t);
					nbData->NbPosSwabs[i]++;
				}
			}
		}
	}

	free(bedNb);

	free(indexI0);
	free(indexI1);
	free(indexS0);
	free(indexS1);
	free(indexIncid0);
	free(indexIncid1);

	free(NbDischarges0);
	free(NbDischarges1);

	free(indexInfector);

	free(AlreadyColonised);

	for(i=0 ; i<NbPatientsMax ; i++)
	{
		gsl_vector_free(isColonised[i]);
	}

	for(i=0 ; i<NbPatientsMax ; i++)
	{
		gsl_vector_free(isInHosp[i]);
	}

	gsl_matrix_free(indexOfPatientsInWard0);
	gsl_matrix_free(indexOfPatientsInWard1);

	return NbCases;
			
}

int main(int argc, char *argv[]) 
{
	
	int where = 1;

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


	/*****************************************************/
	/***     DEFINITION OF A RANDOM GENERATOR          ***/
	/*****************************************************/
		
    time_t t;
    t = 1;//time(NULL); // time in seconds, used to change the seed of the random generator
    const gsl_rng_type *typ;
    gsl_rng_env_setup();    
    typ=gsl_rng_default;
    rng=gsl_rng_alloc(typ);
    gsl_rng_set(rng,t); // changes the seed of the random generator

    /*****************************************************/

    /* parameters */
   	parameters *param = createParam();
   	readParameters(workspace,param);

   	/* data */
   	nb_data *nbData = createNbData();
   	raw_data *data = createRawData();

    /* Simulation algorithm */
   	int NbCases = SimulEpid(param, nbData, data);
   	writeData(NbCases, workspace, nbData, data);

  	/* Closing files and freeing memory */
    freeParam(param);
    gsl_rng_free(rng);
    
	//getchar();
    return 0;

}
