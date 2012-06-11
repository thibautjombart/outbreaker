#if 0

#include "common.h"
#include "init.h"
#include "InputOutput.h"
#include "logL.h"
#include "mcmc.h"
#include "moves.h"
#include "prior.h"
#include "alloc.h"
#include "genclasses.h"
#include "distances.h"
#include "genlike.h"
#include "tuneVariances.h"


/******************************************************************************/
/* Loglikelihood                                                              */
/******************************************************************************/

/**************** Observation level ****************/

/* observation level for patient i */
double ObsLevelPerCase (int i, raw_data * data, nb_data *nb, aug_data *augData, parameters * param){
    int j,k;
    double logpTruePos=0;
    double logpTrueNeg=0;
    double logpFalsePos=0;
    double logpFalseNeg=0;
    double L;

    if(augData->C[i]>=augData->E[i]){
	L = -100000.00;
    } else {

	for(j=0;j<nb->NbPosSwabs[i];j++){/* positive swabs */

	    if(gsl_vector_get(data->P[i],j)<augData->C[i] || gsl_vector_get(data->P[i],j)>=augData->E[i]){
		/* if(1-param->Sp!=0) */
		/* { */
		/* 	logpFalsePos+=log(1-param->Sp); /\* false positives *\/ */
		/* } else */
		/* { */
		logpFalsePos-=100000;
		/* } */
	    } else {
		if(param->Se!=0){ 
		    logpTruePos+=log(param->Se); /* true positives */
		} else {
		    logpTruePos-=100000;
		}
	    }
	}
    }

    for(k=0;k<nb->NbNegSwabs[i];k++){ /* positive swabs */
	if(gsl_vector_get(data->N[i],k)<augData->C[i] || gsl_vector_get(data->N[i],k)>=augData->E[i]){
	    /* if(param->Sp!=0) */
	    /* { */
	    logpTrueNeg+=0;/*log(param->Sp); */ /* true negatives */
	    /* } else */
	    /* { */
	    /* 	logpTrueNeg-=100000; */
	    /* } */
	} else {
	    if(1-param->Se!=0){
		logpFalseNeg+=log(1-param->Se); /* false negatives */
	    } else {
		logpFalseNeg-=100000;
	    }
	}
    }

    L=logpFalsePos+logpTruePos+logpTrueNeg+logpFalseNeg;

    return L;
}




/* observation level for all patients */
double ObsLevel (raw_data * data, nb_data *nb, aug_data *augData, parameters * param){
    int i;
    double L = 0;

    for(i=0;i<data->NbPatients;i++){
	L+=ObsLevelPerCase (i, data, nb, augData, param);
    }

    return L;
}







/**************** Transmission level ****************/
double DurationColonPerCase (int i, aug_data *augData, parameters * param){
    double like,L;
    double a = param->mu*param->mu/(param->sigma*param->sigma);
    double b = param->sigma*param->sigma/param->mu;
    /* continuous time version: double L = log(gsl_ran_gamma_pdf (augData->E[i]-augData->C[i], a, b)); */
    like = gsl_cdf_gamma_P(augData->E[i]-augData->C[i]+0.5, a, b) - gsl_cdf_gamma_P(augData->E[i]-augData->C[i]-0.5, a, b);
    if(like>0){
	L = log(like);
    } else {
	L=-100000;
    }

    return L;
}





double DurationColon (aug_data *augData, parameters * param){
    int i;
    double L = 0;

    for(i=0;i<augData->NbPatients;i++){
	L+=DurationColonPerCase (i, augData, param);
    }

    return L;
}





double ColonPerCase (int i, raw_data *data, nb_data *nb, aug_data *augData, dna_dist *dnainfo, parameters * param){
    double L = 0;
    int l, s, t, j, r, NbPatients=data->NbPatients, T=data->T;
    double Num = 0;
    double Denom = 0;

    double fij;
    double fiWardOut;
    double fiOutOut;

    /* printf("\ndata->IsInHosp[i]:\n"); */
    /* print_gsl_vector(data->IsInHosp[i], "%.0f "); */
    /* fflush(stdout); */

    /* printf("\nColonPerCase i=%d - step 1", i); */
    if(augData->C[i] < gsl_vector_get(data->A[i],0)) /* colonised at first admission */
	{
	    L+=log(param->Pi);
	} else {
	L+=log(1-param->Pi);

	if(augData->C[i]>=0 && augData->C[i]<T){
	    /* printf("\ngsl_vector_get(data->IsInHosp[i],augData->C[i]): %.0f\n", gsl_vector_get(data->IsInHosp[i],augData->C[i])); */
	    if(augData->C[i]>=0 && augData->C[i]<T && gsl_vector_get(data->IsInHosp[i],augData->C[i])==1){ /* colonised whilst in hospital */
		/* printf("\nColonPerCase i=%d - step 2", i); */
		/* fflush(stdout); */
    
		/* finding l such that individual i is colonised during l^th hospital stay (numbered from zero) */
		l=-1;
		while(l<nb->NbAdmissions[i]-1 && augData->C[i]>=gsl_vector_get(data->A[i],l+1)){
		    l++;
		}
		/* printf("\nColonPerCase i=%d - step 3", i); */
		/* fflush(stdout); */

		if(l>0) /* need to escape transmission during at least 1 whole hospital stay and one whole period outside hospital */
		    {
			for(s=0;s<l-1;s++) /* for each hospital stay preceeding the one where colonisation takes place */
			    {
				/* escaping transmission in hospital: */
				for(t=GSL_MAX(0,gsl_vector_get(data->A[i],s));t<GSL_MIN(gsl_vector_get(data->D[i],s),T);t++){
				    L -= gsl_matrix_get(param->beta,data->ward[i],0)*augData->I0[t] + gsl_matrix_get(param->beta,data->ward[i],1)*augData->I1[t];
				}
				L -= param->betaWardOut*(gsl_vector_get(data->D[i],s)-gsl_vector_get(data->A[i],s));

				/* escaping transmission outside hospital: */
				L -= param->betaOutOut*(gsl_vector_get(data->A[i],s+1)-gsl_vector_get(data->D[i],s));
			    }
		    }
		/* printf("\nColonPerCase i=%d - step 4", i); */
		/* fflush(stdout); */

		/* escaping transmission in hospital during the l^th stay, before being colonised : */
		for(t=GSL_MAX(0,gsl_vector_get(data->A[i],l));t<GSL_MIN(T,augData->C[i]);t++){
		    L -= gsl_matrix_get(param->beta,data->ward[i],0)*augData->I0[t] + gsl_matrix_get(param->beta,data->ward[i],1)*augData->I1[t];
		}
		L -= param->betaWardOut*(augData->C[i]-gsl_vector_get(data->A[i],l));

		/* infection at time step C_i : */
		L += log(1-exp(-gsl_matrix_get(param->beta,data->ward[i],0)*augData->I0[augData->C[i]] - gsl_matrix_get(param->beta,data->ward[i],1)*augData->I1[augData->C[i]] - param->betaWardOut));

		/* printf("\nColonPerCase i=%d - step 5", i); */
		/* fflush(stdout); */

		/* relative weights of each route of transmission incorporating genetic data : */
		for(j=0;j<NbPatients;j++){
		    /* calculate fij */
		    fij = genlike_ij(i,j,data, augData, dnainfo, param);
		    /* printf("\nColonPerCase i=%d - step 6", i); */
		    /* fflush(stdout); */

		    /*****************/
		    if(augData->C[i]>=0 && augData->C[i]<T && augData->C[j]<=augData->C[i] && augData->C[i]<augData->E[j] && gsl_vector_get(data->IsInHosp[j],augData->C[i])) /* individual j is colonised and in hospital at time C_i */
			{
			    r=-1;
			    while(r<nb->NbAdmissions[j]-1 && augData->C[i]>=gsl_vector_get(data->A[j],r+1)){
				r++;
			    }
			    Num+=gsl_matrix_get(param->beta,data->ward[i],data->ward[j])*fij;
			    Denom+=gsl_matrix_get(param->beta,data->ward[i],data->ward[j]);
			}
		}

		/* printf("\nColonPerCase i=%d - step 7", i); */
		/* fflush(stdout); */

		/* calculate fiWardOut */
		/**************** Here need to incorporate Thibaut's genetic likelihood ****************/
		fiWardOut = 1;
		/*****************/
		Num+=param->betaWardOut*fiWardOut;
		Denom+=param->betaWardOut;
		L += log(Num) - log(Denom);

	    } else { /* colonised between hospital stays or not colonised over the study period */

		/* finding l such that individual i is colonised after l^th hospital stay (numbered from zero) */
		l=-1;
		/* printf("nb->NbAdmissions[i]-1: %d   \naugData->C[i]: %d  \naugData->E[i]: %d gsl_vector_get(data->A[i],l+1): %.0f gsl_vector_get(data->D[i],l+1): %.0f\n", nb->NbAdmissions[i]-1, augData->C[i], augData->E[i], gsl_vector_get(data->A[i],l+1), gsl_vector_get(data->D[i],l+1)); */
		/* fflush(stdout); */

		while(l<nb->NbAdmissions[i]-1 && augData->C[i]>=gsl_vector_get(data->D[i],l+1)){
		    l++;
		}
		/* printf("\ni = %d  l = %d\n", i, l); */
		/* printf("\nColonPerCase i=%d - step 8", i); */
		/* fflush(stdout); */

		if(l<nb->NbAdmissions[i]){ /* colonised before the end of the period study */
		    /* printf("\nColonPerCase i=%d - step 8bis", i); */
		    /* printf("\ni = %d   s = %d   l = %d\n", i, s, l); */
		    /* fflush(stdout); */
		    for(s=0;s<l;s++) /* for each previous hospital stay */
			{
			    /* escaping transmission in hospital: */
			    for(t=GSL_MAX(0,gsl_vector_get(data->A[i],s)) ; t<GSL_MIN(T,gsl_vector_get(data->D[i],s)) ; t++){
				L -= gsl_matrix_get(param->beta,data->ward[i],0)*augData->I0[t] + gsl_matrix_get(param->beta,data->ward[i],1)*augData->I1[t];
			    }
			    /* printf("\nColonPerCase i=%d - step 9", i); */
			    /* fflush(stdout); */

			    L -= param->betaWardOut*(gsl_vector_get(data->D[i],s)-gsl_vector_get(data->A[i],s));
			    /* printf("\nColonPerCase i=%d - step 10", i); */
			    /* fflush(stdout); */
			}
		    if(l>0) /* need to escape transmission during at least 1 whole period inbetween hospital stays */
			{
			    for(s=0;s<l-1;s++) /* for each hospital stay preceeding the one where colonisation takes place */
				{
				    /* escaping transmission outside hospital: */
				    L -= param->betaOutOut*(gsl_vector_get(data->A[i],s+1)-gsl_vector_get(data->D[i],s));
				}
			    /* printf("\nColonPerCase i=%d - step 11", i); */
			    /* fflush(stdout); */

			}

		    /* escaping transmission outside hospital after the l^th stay, before being colonised : */
		    L -= param->betaOutOut*(augData->C[i]-gsl_vector_get(data->D[i],l));

		    /* infection at time step C_i : */
		    L += log(1-exp(-param->betaOutOut));

		    /* relative weights of each route of transmission incorporating genetic data : */
		    /* calculate fiOutOut */
		    /**************** Here need to incorporate Thibaut's genetic likelihood ****************/
		    fiOutOut = 1;
		    /*****************/
		    L += log(fiOutOut);

		}
	    }
	}
    }
    /* printf("\nColonPerCase i=%d - step 10", i); */

    return L;
}





double Colon(raw_data *data, nb_data *nb, aug_data *augData, dna_dist *dnainfo, parameters * param){
    int i;
    double L = 0;

    for(i=0;i<data->NbPatients;i++){
	/* printf("\ni: %d\n",i); */
	/* fflush(stdout); */
	L+=ColonPerCase (i, data, nb, augData, dnainfo, param);
    }

    return L;
}







/**************** Total ****************/
double fullLoglikelihoodPerCase(int i, raw_data * data, nb_data *nb, aug_data *augData, dna_dist *dnainfo, parameters * param){
    double L = ObsLevelPerCase (i, data, nb, augData, param) + DurationColonPerCase (i, augData, param) + ColonPerCase (i, data, nb, augData, dnainfo, param);
    return(L);
}





double fullLoglikelihood(raw_data * data, nb_data *nb, aug_data *augData, dna_dist *dnainfo, parameters * param){
    double L = ObsLevel (data, nb, augData, param);
    L+= DurationColon(augData, param);
    L+= Colon(data, nb, augData, dnainfo, param);
    return(L);
}





double fullLoglikelihoodWithPrior(raw_data * data, nb_data *nb, aug_data *augData, dna_dist *dnainfo, parameters * param){
    double L = fullLoglikelihood(data, nb, augData, dnainfo, param) + logprior(param);
    return(L);
}



#endif
