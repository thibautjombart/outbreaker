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

void fillingBeds(parameters *param, hospDurationParam *paramHosp, nb_data *nbData, raw_data *data, aug_data *augData, int *indexInfector, int *NbCases, int *bedNb, gsl_matrix * NbDischarges0, gsl_matrix * NbDischarges1, int *TotalNbDischarges0, int *TotalNbDischarges1, gsl_vector ** isInHosp,gsl_matrix * indexOfPatientsInWard0, gsl_matrix * indexOfPatientsInWard1, int *AlreadyColonised, gsl_vector ** isColonised)
{
	int i, prov, u, t;
	double rand;
	/* cases present in the 2 wards on day 0*/

		for(i=0 ; i<SizeWard0 ; i++) /* patients initially present in ward 0 */
		{
			data->ward[*NbCases]=0;
			bedNb[*NbCases]=i;
			nbData->NbAdmissions[*NbCases]=1;
			gsl_vector_set(data->A[*NbCases],0,-ceil(gsl_ran_gamma (rng, (paramHosp->mu/2)*(paramHosp->mu/2)/(paramHosp->sigma*paramHosp->sigma), paramHosp->sigma*paramHosp->sigma/(paramHosp->mu/2))));
			do
			{
				gsl_vector_set(data->D[*NbCases],0,ceil(gsl_vector_get(data->A[*NbCases],0) + gsl_ran_gamma (rng, (paramHosp->mu)*(paramHosp->mu)/(paramHosp->sigma*paramHosp->sigma), paramHosp->sigma*paramHosp->sigma/(paramHosp->mu))));
			}while(gsl_vector_get(data->D[*NbCases],0)<=0);
			if(gsl_vector_get(data->D[*NbCases],0)<nbData->T)
			{
				prov =gsl_matrix_get(NbDischarges0,bedNb[*NbCases],gsl_vector_get(data->D[*NbCases],0));
				gsl_matrix_set(NbDischarges0,bedNb[*NbCases],gsl_vector_get(data->D[*NbCases],0),prov+1);
				TotalNbDischarges0[(int)gsl_vector_get(data->D[*NbCases],0)]++;
			}
			if(GSL_MAX(0,gsl_vector_get(data->A[*NbCases],0))<GSL_MIN(nbData->T,gsl_vector_get(data->D[*NbCases],0)))
			{
				for(u=GSL_MAX(0,gsl_vector_get(data->A[*NbCases],0)) ; u<GSL_MIN(nbData->T,gsl_vector_get(data->D[*NbCases],0)) ; u++)
				{
					gsl_vector_set(isInHosp[*NbCases],u,1);
					gsl_matrix_set(indexOfPatientsInWard0,i,u,*NbCases);
				}
			}
			rand=gsl_ran_flat (rng, 0.0, 1.0);
			if(rand<param->Pi) /* colonised before admission */
			{
				AlreadyColonised[*NbCases]=1;
				indexInfector[*NbCases]=-1;
				augData->C[*NbCases]=gsl_vector_get(data->A[*NbCases],0)-1;
				augData->E[*NbCases]=augData->C[*NbCases]+ceil(gsl_ran_gamma (rng, (param->mu)*(param->mu)/(param->sigma*param->sigma), param->sigma*param->sigma/(param->mu)));
				if(GSL_MAX(0,augData->C[*NbCases])<GSL_MIN(nbData->T,augData->E[*NbCases]))
				{
					for(u=GSL_MAX(0,augData->C[*NbCases]) ; u<GSL_MIN(nbData->T,augData->E[*NbCases]) ; u++)
					{
						gsl_vector_set(isColonised[*NbCases],u,1);
					}
				}
			}else
			{
				AlreadyColonised[*NbCases]=0;
			}
			(*NbCases)++;
			printf("NbCases\t %d\n",*NbCases);
			fflush(stdout);
		}
		for(i=0 ; i<SizeWard1 ; i++) /* patients initially present in ward 1 */
		{
			data->ward[*NbCases]=1;
			bedNb[*NbCases]=i;
			nbData->NbAdmissions[*NbCases]=1;
			gsl_vector_set(data->A[*NbCases],0,-ceil(gsl_ran_gamma (rng, (paramHosp->mu/2)*(paramHosp->mu/2)/(paramHosp->sigma*paramHosp->sigma), paramHosp->sigma*paramHosp->sigma/(paramHosp->mu/2))));
			do
			{
				gsl_vector_set(data->D[*NbCases],0,ceil(gsl_vector_get(data->A[*NbCases],0) + gsl_ran_gamma (rng, (paramHosp->mu)*(paramHosp->mu)/(paramHosp->sigma*paramHosp->sigma), paramHosp->sigma*paramHosp->sigma/(paramHosp->mu))));
			}while(gsl_vector_get(data->D[*NbCases],0)<=0);
			if(gsl_vector_get(data->D[*NbCases],0)<nbData->T)
			{
				prov =gsl_matrix_get(NbDischarges1,bedNb[*NbCases],gsl_vector_get(data->D[*NbCases],0));
				gsl_matrix_set(NbDischarges1,bedNb[*NbCases],gsl_vector_get(data->D[*NbCases],0),prov+1);
				TotalNbDischarges1[(int)gsl_vector_get(data->D[*NbCases],0)]++;
			}
			if(GSL_MAX(0,gsl_vector_get(data->A[*NbCases],0))<GSL_MIN(nbData->T,gsl_vector_get(data->D[*NbCases],0)))
			{
				for(u=GSL_MAX(0,gsl_vector_get(data->A[*NbCases],0)) ; u<GSL_MIN(nbData->T,gsl_vector_get(data->D[*NbCases],0)) ; u++)
				{
					gsl_vector_set(isInHosp[*NbCases],u,1);
					gsl_matrix_set(indexOfPatientsInWard1,i,u,*NbCases);
				}
			}
			rand=gsl_ran_flat (rng, 0.0, 1.0);
			if(rand<param->Pi) /* colonised before admission */
			{
				AlreadyColonised[*NbCases]=1;
				indexInfector[*NbCases]=-1;
				augData->C[*NbCases]=gsl_vector_get(data->A[*NbCases],0)-1;
				augData->E[*NbCases]=augData->C[*NbCases]+ceil(gsl_ran_gamma (rng, (param->mu)*(param->mu)/(param->sigma*param->sigma), param->sigma*param->sigma/(param->mu)));
				if(GSL_MAX(0,augData->C[*NbCases])<GSL_MIN(nbData->T,augData->E[*NbCases]))
				{
					for(u=GSL_MAX(0,augData->C[*NbCases]) ; u<GSL_MIN(nbData->T,augData->E[*NbCases]) ; u++)
					{
						gsl_vector_set(isColonised[*NbCases],u,1);
					}
				}
			}else
			{
				AlreadyColonised[*NbCases]=0;
			}
			(*NbCases)++;
			printf("NbCases\t %d\n",*NbCases);
			fflush(stdout);
		}

		/* simulate the admissions/discharges up to nbData->T */

		for(t=0 ; t<nbData->T ; t++)
		{
			for(i=0 ; i<SizeWard0 ; i++)
			{
				if(gsl_matrix_get(NbDischarges0,i,t)>0)
				{
					data->ward[*NbCases]=0;
					bedNb[*NbCases]=i;
					nbData->NbAdmissions[*NbCases]=1;
					gsl_vector_set(data->A[*NbCases],0,t);
					gsl_vector_set(data->D[*NbCases],0,ceil(gsl_vector_get(data->A[*NbCases],0) + gsl_ran_gamma (rng, (paramHosp->mu)*(paramHosp->mu)/(paramHosp->sigma*paramHosp->sigma), paramHosp->sigma*paramHosp->sigma/(paramHosp->mu))));
					if(gsl_vector_get(data->D[*NbCases],0)<nbData->T)
					{
						prov =gsl_matrix_get(NbDischarges0,bedNb[*NbCases],gsl_vector_get(data->D[*NbCases],0));
						gsl_matrix_set(NbDischarges0,bedNb[*NbCases],gsl_vector_get(data->D[*NbCases],0),prov+1);
						TotalNbDischarges0[(int)gsl_vector_get(data->D[*NbCases],0)]++;
					}
					if(GSL_MAX(0,gsl_vector_get(data->A[*NbCases],0))<GSL_MIN(nbData->T,gsl_vector_get(data->D[*NbCases],0)))
					{
						for(u=GSL_MAX(0,gsl_vector_get(data->A[*NbCases],0)) ; u<GSL_MIN(nbData->T,gsl_vector_get(data->D[*NbCases],0)) ; u++)
						{
							gsl_vector_set(isInHosp[*NbCases],u,1);
							gsl_matrix_set(indexOfPatientsInWard0,i,u,*NbCases);
						}
					}
					rand=gsl_ran_flat (rng, 0.0, 1.0);
					if(rand<param->Pi) /* colonised before admission */
					{
						AlreadyColonised[*NbCases]=1;
						augData->C[*NbCases]=gsl_vector_get(data->A[*NbCases],0)-1;
						augData->E[*NbCases]=augData->C[*NbCases]+ceil(gsl_ran_gamma (rng, (param->mu)*(param->mu)/(param->sigma*param->sigma), param->sigma*param->sigma/(param->mu)));
						if(GSL_MAX(0,augData->C[*NbCases])<GSL_MIN(nbData->T,augData->E[*NbCases]))
						{
							for(u=GSL_MAX(0,augData->C[*NbCases]) ; u<GSL_MIN(nbData->T,augData->E[*NbCases]) ; u++)
							{
								gsl_vector_set(isColonised[*NbCases],u,1);
							}
						}
					}else
					{
						AlreadyColonised[*NbCases]=0;
					}
					(*NbCases)++;
					printf("NbCases\t %d\n",*NbCases);
					fflush(stdout);
				}
			}

			for(i=0 ; i<SizeWard1 ; i++)
			{
				if(gsl_matrix_get(NbDischarges1,i,t)>0)
				{
					data->ward[*NbCases]=1;
					bedNb[*NbCases]=i;
					nbData->NbAdmissions[*NbCases]=1;
					gsl_vector_set(data->A[*NbCases],0,t);
					gsl_vector_set(data->D[*NbCases],0,ceil(gsl_vector_get(data->A[*NbCases],0) + gsl_ran_gamma (rng, (paramHosp->mu)*(paramHosp->mu)/(paramHosp->sigma*paramHosp->sigma), paramHosp->sigma*paramHosp->sigma/(paramHosp->mu))));
					if(gsl_vector_get(data->D[*NbCases],0)<nbData->T)
					{
						prov =gsl_matrix_get(NbDischarges1,bedNb[*NbCases],gsl_vector_get(data->D[*NbCases],0));
						gsl_matrix_set(NbDischarges1,bedNb[*NbCases],gsl_vector_get(data->D[*NbCases],0),prov+1);
						TotalNbDischarges1[(int)gsl_vector_get(data->D[*NbCases],0)]++;
					}
					if(GSL_MAX(0,gsl_vector_get(data->A[*NbCases],0))<GSL_MIN(nbData->T,gsl_vector_get(data->D[*NbCases],0)))
					{
						for(u=GSL_MAX(0,gsl_vector_get(data->A[*NbCases],0)) ; u<GSL_MIN(nbData->T,gsl_vector_get(data->D[*NbCases],0)) ; u++)
						{
							gsl_vector_set(isInHosp[*NbCases],u,1);
							gsl_matrix_set(indexOfPatientsInWard1,i,u,*NbCases);
						}
					}
					rand=gsl_ran_flat (rng, 0.0, 1.0);
					if(rand<param->Pi) /* colonised before admission */
					{
						AlreadyColonised[*NbCases]=1;
						augData->C[*NbCases]=gsl_vector_get(data->A[*NbCases],0)-1;
						augData->E[*NbCases]=augData->C[*NbCases]+ceil(gsl_ran_gamma (rng, (param->mu)*(param->mu)/(param->sigma*param->sigma), param->sigma*param->sigma/(param->mu)));
						if(GSL_MAX(0,augData->C[*NbCases])<GSL_MIN(nbData->T,augData->E[*NbCases]))
						{
							for(u=GSL_MAX(0,augData->C[*NbCases]) ; u<GSL_MIN(nbData->T,augData->E[*NbCases]) ; u++)
							{
								gsl_vector_set(isColonised[*NbCases],u,1);
							}
						}
					}else
					{
						AlreadyColonised[*NbCases]=0;
					}
					(*NbCases)++;
					printf("NbCases\t %d\n",*NbCases);
					fflush(stdout);
				}
			}
		}
}

void SimulColonisationInBeds(parameters *param, nb_data *nbData, aug_data *augData, int *indexInfector, gsl_matrix * indexOfPatientsInWard0, gsl_matrix * indexOfPatientsInWard1, int *AlreadyColonised, gsl_vector ** isColonised, int *indexI0, int *indexI1, int *indexS0, int *indexS1, 	int *indexIncid0, 	int *indexIncid1)
{
	int t,k,i,u;
	double rand;
	double lambda, prob;
	int S0, S1;
	int Incid0, Incid1;
	double probaInfectedPerWard0;
	double probaInfectedPerWard1;

	for(t=1 ; t<nbData->T ; t++)
		{
			S0=0;
			S1=0;

			for(k=0 ; k<SizeWard0 ; k++)
			{
				i=gsl_matrix_get(indexOfPatientsInWard0,k,t);
				if(gsl_vector_get(isColonised[i],t-1)==1)
				{
					indexI0[augData->I0[t-1]]=i;
					augData->I0[t-1]++;
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
					indexI1[augData->I1[t-1]]=i;
					augData->I1[t-1]++;
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

			lambda = gsl_matrix_get(param->beta,0,0)*augData->I0[t-1] + gsl_matrix_get(param->beta,0,1)*augData->I1[t-1] + param->betaWardOut;
			prob = 1-exp(-lambda);
			Incid0=gsl_ran_binomial (rng, prob, S0);

			if(Incid0>0)
			{
				/* Choose them */
				gsl_ran_choose (rng, indexIncid0, Incid0, indexS0, S0, sizeof(int));

				/* draw their infector */

				probaInfectedPerWard0 = gsl_matrix_get(param->beta,0,0)*augData->I0[t-1]/lambda;
				probaInfectedPerWard1 = gsl_matrix_get(param->beta,0,1)*augData->I1[t-1]/lambda;

				for(i=0 ; i<Incid0 ; i++)
				{
					rand = gsl_rng_uniform(rng);
					if(rand<probaInfectedPerWard0) /* infected by someone in ward 0 */
					{
						indexInfector[indexIncid0[i]]= indexI0[gsl_rng_uniform_int (rng, augData->I0[t-1])];
					}else
					{
						if(rand<probaInfectedPerWard0+probaInfectedPerWard1) /* infected by someone in ward 1 */
						{
							indexInfector[indexIncid0[i]]= indexI1[gsl_rng_uniform_int (rng, augData->I1[t-1])];
						}else /* infected by outside */
						{
							indexInfector[indexIncid0[i]]= -1;
						}
					}

				}

				/* Draw there time of end of colonistaion */
				for(i=0 ; i<Incid0 ; i++)
				{
					augData->C[indexIncid0[i]]=t-1;
					augData->E[indexIncid0[i]]=augData->C[indexIncid0[i]]+ceil(gsl_ran_gamma (rng, (param->mu)*(param->mu)/(param->sigma*param->sigma), param->sigma*param->sigma/(param->mu)));

					/* Change isColonised and AlreadyColonised for those individuals */

					AlreadyColonised[indexIncid0[i]]=1;
					if(GSL_MAX(0,augData->C[indexIncid0[i]])<GSL_MIN(nbData->T,augData->E[indexIncid0[i]]))
					{
						for(u=GSL_MAX(0,augData->C[indexIncid0[i]]) ; u<GSL_MIN(nbData->T,augData->E[indexIncid0[i]]) ; u++)
						{
							gsl_vector_set(isColonised[indexIncid0[i]],u,1);
						}
					}
				}

			}

			/* Total number of individuals infected on day t-1 in ward 1 */

			lambda = gsl_matrix_get(param->beta,1,0)*augData->I0[t-1] + gsl_matrix_get(param->beta,1,1)*augData->I1[t-1] + param->betaWardOut;
			prob = 1-exp(-lambda);
			Incid1=gsl_ran_binomial (rng, prob, S1);

			if(Incid1>0)
			{
				/* Choose them */
				gsl_ran_choose (rng, indexIncid1, Incid1, indexS1, S1, sizeof(int));

				/* draw their infector */
				probaInfectedPerWard0 = gsl_matrix_get(param->beta,1,0)*augData->I0[t-1]/lambda;
				probaInfectedPerWard1 = gsl_matrix_get(param->beta,1,1)*augData->I1[t-1]/lambda;

				for(i=0 ; i<Incid1 ; i++)
				{
					rand = gsl_rng_uniform(rng);
					if(rand<probaInfectedPerWard0) /* infected by someone in ward 0 */
					{
						indexInfector[indexIncid1[i]]= indexI0[gsl_rng_uniform_int (rng, augData->I0[t-1])];
					}else
					{
						if(rand<probaInfectedPerWard0+probaInfectedPerWard1) /* infected by someone in ward 1 */
						{
							indexInfector[indexIncid1[i]]= indexI1[gsl_rng_uniform_int (rng, augData->I1[t-1])];
						}else /* infected by outside */
						{
							indexInfector[indexIncid1[i]]= -1;
						}
					}

				}

				/* Draw there time of end of colonistaion */
				for(i=0 ; i<Incid1 ; i++)
				{
					augData->C[indexIncid1[i]]=t-1;
					augData->E[indexIncid1[i]]=augData->C[indexIncid1[i]]+ceil(gsl_ran_gamma (rng, (param->mu)*(param->mu)/(param->sigma*param->sigma), param->sigma*param->sigma/(param->mu)));
				}

				/* Change isColonised and AlreadyColonised for those individuals */

				AlreadyColonised[indexIncid1[i]]=1;
				if(GSL_MAX(0,augData->C[indexIncid1[i]])<GSL_MIN(nbData->T,augData->E[indexIncid1[i]]))
				{
					for(u=GSL_MAX(0,augData->C[indexIncid1[i]]) ; u<GSL_MIN(nbData->T,augData->E[indexIncid1[i]]) ; u++)
					{
						gsl_vector_set(isColonised[indexIncid1[i]],u,1);
					}
				}
			}

		}
}


int SimulEpid(parameters *param, hospDurationParam *paramHosp, nb_data *nbData, raw_data *data, aug_data *augData, int *indexInfector)
{
	int t,k,i;
	int NbCases=0;
	double rand;

	int *bedNb = (int *) calloc(nbData->NbPatients, sizeof(int));

	for(i=0 ; i<nbData->NbPatients ; i++)
	{
		indexInfector[i]=-100;
	}

	int *indexI0 = calloc(SizeWard0, sizeof(int));
	int *indexS0 = calloc(SizeWard0, sizeof(int));
	int *indexI1 = calloc(SizeWard1, sizeof(int));
	int *indexS1 = calloc(SizeWard1, sizeof(int));
	int *indexIncid0 = calloc(SizeWard0, sizeof(int));
	int *indexIncid1 = calloc(SizeWard1, sizeof(int));

	gsl_matrix * NbDischarges0 = gsl_matrix_calloc(SizeWard0,nbData->T);
	gsl_matrix * NbDischarges1 = gsl_matrix_calloc(SizeWard1,nbData->T);

	int *TotalNbDischarges0 = calloc(nbData->T, sizeof(int));
	int *TotalNbDischarges1 = calloc(nbData->T, sizeof(int));

	for(t=0 ; t<nbData->T ; t++)
	{
		TotalNbDischarges0[t]=0;
		TotalNbDischarges1[t]=0;
	}

	int *AlreadyColonised = calloc(nbData->NbPatients, sizeof(int));

	gsl_vector * isColonised[nbData->NbPatients];
	for(i=0 ; i<nbData->NbPatients ; i++)
	{
		isColonised[i] = gsl_vector_calloc(nbData->T);
	}

	gsl_vector * isInHosp[nbData->NbPatients];
	for(i=0 ; i<nbData->NbPatients ; i++)
	{
		isInHosp[i] = gsl_vector_calloc(nbData->T);
	}

	gsl_matrix * indexOfPatientsInWard0 = gsl_matrix_calloc(SizeWard0,nbData->T);
	gsl_matrix * indexOfPatientsInWard1 = gsl_matrix_calloc(SizeWard1,nbData->T);

	/* Filling the beds with patienst */

	fillingBeds(param, paramHosp, nbData, data, augData, indexInfector, &NbCases, bedNb, NbDischarges0, NbDischarges1, TotalNbDischarges0, TotalNbDischarges1, isInHosp,indexOfPatientsInWard0, indexOfPatientsInWard1, AlreadyColonised, isColonised);

	/* simulate the colonisation events up to the end (for those not colonised at admission) */

	SimulColonisationInBeds(param, nbData, augData, indexInfector, indexOfPatientsInWard0, indexOfPatientsInWard1, AlreadyColonised, isColonised, indexI0, indexI1, indexS0, indexS1, 	indexIncid0, 	indexIncid1);

	/* simulate the testing */
	/* here testing occurs at every time step */
	/* you can then subsample to simulate a less regular testing scheme */

	for(t=0 ; t<nbData->T ; t++)
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
				//rand=gsl_ran_flat (rng, 0.0, 1.0);
				/* here Sp = 100% */
				// if(rand<param->Sp) /* test negative */
				//{
					gsl_vector_set(data->N[i],nbData->NbNegSwabs[i],t);
					nbData->NbNegSwabs[i]++;
				/*}else // false positive test
				{
					gsl_vector_set(data->P[i],nbData->NbPosSwabs[i],t);
					nbData->NbPosSwabs[i]++;
				}*/
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
				//rand=gsl_ran_flat (rng, 0.0, 1.0);
				/* here Sp = 100% */
				//if(rand<param->Sp) /* test negative */
				//{
					gsl_vector_set(data->N[i],nbData->NbNegSwabs[i],t);
					nbData->NbNegSwabs[i]++;
				/*}else // false positive test
				{
					gsl_vector_set(data->P[i],nbData->NbPosSwabs[i],t);
					nbData->NbPosSwabs[i]++;
				}*/
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

	free(AlreadyColonised);

	for(i=0 ; i<nbData->NbPatients ; i++)
	{
		gsl_vector_free(isColonised[i]);
	}

	for(i=0 ; i<nbData->NbPatients ; i++)
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
   	hospDurationParam *paramHosp = createHospDurationParam();
   	readParameters(workspace,param,paramHosp);

   	/* data */
   	int NbPatientsMax = 10000;
   	int Tmax=100;
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

    /* Simulation algorithm */
   	int NbCases = SimulEpid(param, paramHosp, nbData, data, augData,indexInfector);
   	writeData(NbCases, workspace, nbData, data, augData,indexInfector);

  	/* Closing files and freeing memory */
    freeParam(param);
    freeAugData(augData);
    freeNbData(nbData);
    freeRawData(data);
    gsl_rng_free(rng);
    
	//getchar();
    return 0;

}
