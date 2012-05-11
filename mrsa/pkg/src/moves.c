#include "common.h"
#include "init.h"
#include "InputOutput.h"
#include "logL.h"
#include "mcmc.h"
#include "moves.h"
#include "prior.h"
#include "alloc.h"
#include "tuneVariances.h"






/*********************** Move parameters ***********************/
void moveBeta(int i, int j, mcmcInternals * MCMCSettings, parameters * curParam, parameters *newParam, raw_data * data, nb_data *nb, aug_data *augData, dna_dist *dnainfo, acceptance *accept, NbProposals *NbProp){
    /* updating beta_{i,j} with Metropolis algorithm */
    /* proposal distribution = lognormal */

    double * newVal, *curVal;
    double * nbAccept;
    double * nbPropos;
    double sigmaProp;
    double r,z;
    double pAccept = 0;
  
    /* printf("\nMove beta: i = %d, j = %d", i, j); */
    curVal = gsl_matrix_ptr(curParam->beta,i,j);
    newVal = gsl_matrix_ptr(newParam->beta,i,j);
    /* printf("\n - step i"); */

    sigmaProp = gsl_matrix_get(MCMCSettings->Sigma_beta,i,j);
    /* printf("\n - step ii"); */
    nbAccept = gsl_matrix_ptr(accept->PourcAcc_beta,i,j);
    /* printf("\n - step iii"); */
    nbPropos = gsl_matrix_ptr(NbProp->NbProp_beta,i,j);
    /* printf("\n - step iv"); */

    /* printf("\nValue of sigmaProp: %.5f", sigmaProp); */
    /* printf("\nValue of newVal: %.5f", *newVal); */
    /* printf("\nValue of curVal: %.5f", *curVal); */
    /* printf("\nRandom value: %.5f", gsl_ran_flat (data->rng, 0.0, 100.0)); */
    *newVal = *curVal * gsl_ran_lognormal(data->rng,0.0,sigmaProp);
    /* printf("\n - step v "); */
    /* fflush(stdout); */

    pAccept += Colon(data, nb, augData, dnainfo, newParam);
    pAccept -= Colon(data, nb, augData, dnainfo, curParam);

    /* printf("\n - step vi");fflush(stdout); */
    pAccept +=  logpriorBeta(i,j, newParam) - logpriorBeta(i,j,curParam);
    /* printf("\n - step viii");fflush(stdout); */

    pAccept +=  log(*(newVal)) - log(*(curVal)); /* correction for lognormal */

    if (pAccept>0) r=0; else r=pAccept;
    z=gsl_rng_uniform(data->rng);
    if (log(z)<=r) {
    	*curVal = *newVal;
    	*nbAccept +=1;
    }else{
	*newVal = *curVal;
    }
    *nbPropos+=1;

}





void moveAllBeta(mcmcInternals * MCMCSettings, parameters * curParam, parameters *newParam, raw_data * data, nb_data *nb, aug_data *augData, dna_dist *dnainfo, acceptance *accept, NbProposals *NbProp){
    int i,j;

    for (i=0;i<2;i++){
	for(j=0;j<2;j++){
	    moveBeta(i,j,MCMCSettings, curParam, newParam, data, nb, augData, dnainfo, accept, NbProp);
	}
    }

}





void moveBetaOut( char what, mcmcInternals * MCMCSettings, parameters * curParam, parameters *newParam, raw_data * data, nb_data *nb, aug_data *augData, dna_dist *dnainfo, acceptance *accept, NbProposals *NbProp){
    /* what is one of "w" or "o" */
    /* updating betaWardOut ("w") or betaOutOut ("o") with Metropolis algorithm */
    /* proposal distribution = lognormal */

    double *newVal, *curVal;
    double * nbAccept;
    double * nbPropos;
    double sigmaProp;
    double r,z;
    double pAccept = 0;
    double (* logPrior) (parameters *);

    switch(what) {
    case 'w':
	curVal = &curParam->betaWardOut;
	newVal = &newParam->betaWardOut;
	sigmaProp = MCMCSettings->Sigma_betaWardOut;
	nbAccept = &accept->PourcAcc_betaWardOut;
	nbPropos = &NbProp->NbProp_betaWardOut;
	logPrior=logpriorBetaWardOut;
	break;
    case 'o':
	curVal = &curParam->betaOutOut;
	newVal = &newParam->betaOutOut;
	sigmaProp = MCMCSettings->Sigma_betaOutOut;
	nbAccept = &accept->PourcAcc_betaOutOut;
	nbPropos = &NbProp->NbProp_betaOutOut;
	logPrior=logpriorBetaOutOut;
	break;
    }

    *(newVal) = *(curVal)*gsl_ran_lognormal(data->rng,0,sigmaProp);

    pAccept += Colon(data, nb, augData, dnainfo, newParam);
    pAccept -= Colon(data, nb, augData, dnainfo, curParam);

    pAccept +=  logPrior(newParam) - logPrior(curParam);

    pAccept +=  log(*(newVal)) - log(*(curVal)); /* correction for lognormal */

    if (pAccept>0) r=0; else r=pAccept;
    z=gsl_rng_uniform(data->rng);
    if (log(z)<=r) {
	*curVal = *newVal;
	*nbAccept +=1;
    }else{
	*newVal = *curVal;
    }
    *nbPropos+=1;

}





void moveSe(parameters * curParam, parameters * newParam, raw_data * data, nb_data *nb, aug_data *augData, acceptance *accept){
    /* updating Se with a Gibbs sampler */
    double NbFalseNeg=0;
    double NbTruePos=0;
    int i,j,k;

    for(i=0;i<data->NbPatients;i++){
	for(j=0;j<nb->NbPosSwabs[i];j++){ /* positive swabs */
	    if(augData->C[i]<=gsl_vector_get(data->P[i],j) && gsl_vector_get(data->P[i],j)<=augData->E[i]){
		NbTruePos++; /* true positives */
	    }
	}
	for(k=0;k<nb->NbNegSwabs[i];k++){ /* positive swabs */

	    if(augData->C[i]<=gsl_vector_get(data->N[i],k) && gsl_vector_get(data->N[i],k)<=augData->E[i]){
		NbFalseNeg++; /* false negatives */
	    }
	}
    }

    curParam->Se = gsl_ran_beta (data->rng, NbTruePos+1, NbFalseNeg+1);
}





void movePi(parameters * curParam, parameters * newParam, raw_data * data, aug_data *augData, acceptance *accept){
    /* updating Pi with a Gibbs sampler */
    double ColonisedAtFirstAdmission=0;
    double NonColonisedAtFirstAdmission=0;
    int i;

    for(i=0;i<data->NbPatients;i++){
	if(augData->C[i]<gsl_vector_get(data->A[i],0)){
	    ColonisedAtFirstAdmission++;
	} else {
	    NonColonisedAtFirstAdmission++;
	}
    }

    curParam->Pi = gsl_ran_beta (data->rng, ColonisedAtFirstAdmission+1, NonColonisedAtFirstAdmission+1);
}





void moveDurationColon(char what,mcmcInternals * MCMCSettings, parameters * curParam, parameters *newParam, raw_data *data, aug_data *augData, acceptance *accept, NbProposals *NbProp){
    /* what is one of "m" or "s" */
    /* updating mu ("m") or sigma ("s") with Metropolis algorithm */
    /* proposal distribution = lognormal */
    double *newVal, *curVal;
    double * nbAccept;
    double * nbPropos;
    double sigmaProp;
    double r,z;
    double pAccept = 0;
    double (* logPrior) (parameters *);

    switch(what) {
    case 'm':
	curVal = &curParam->mu;
	newVal = &newParam->mu;
	sigmaProp = MCMCSettings->Sigma_mu;
	nbAccept = &accept->PourcAcc_mu;
	nbPropos = &NbProp->NbProp_mu;
	logPrior = logpriorMu;
	break;
    case 's':
	curVal = &curParam->sigma;
	newVal = &newParam->sigma;
	sigmaProp = MCMCSettings->Sigma_sigma;
	nbAccept = &accept->PourcAcc_sigma;
	nbPropos = &NbProp->NbProp_sigma;
	logPrior = logpriorSigma;
	break;
    }

    *(newVal) = *(curVal)*gsl_ran_lognormal(data->rng,0,sigmaProp);

    pAccept += DurationColon(augData, newParam);
    pAccept -= DurationColon(augData, curParam);

    pAccept +=  logPrior(newParam) - logPrior(curParam);

    pAccept +=  log(*(newVal)) - log(*(curVal)); /* correction for lognormal */

    if (pAccept>0) r=0; else r=pAccept;
    z=gsl_rng_uniform(data->rng);
    if (log(z)<=r) {
	*curVal = *newVal;
	*nbAccept +=1;
    }else{
	*newVal = *curVal;
    }
    *nbPropos+=1;

}





void moveNu(mcmcInternals * MCMCSettings, parameters * curParam, parameters *newParam, raw_data * data, nb_data *nb, aug_data *augData, dna_dist *dnainfo, acceptance *accept, NbProposals *NbProp){
    /* updating nu1 with Metropolis algorithm */
    /* proposal distribution = lognormal */
    double *newVal, *curVal;
    double * nbAccept;
    double * nbPropos;
    double sigmaProp;
    double r,z;
    double pAccept = 0;
    /* double (* logPrior) (parameters *); */

    curVal = &curParam->nu1;
    newVal = &newParam->nu1;
    sigmaProp = MCMCSettings->Sigma_nu1;
    nbAccept = &accept->PourcAcc_nu1;
    nbPropos = &NbProp->NbProp_nu1;

    *(newVal) = *(curVal)*gsl_ran_lognormal(data->rng,0,sigmaProp);

    pAccept += Colon(data, nb, augData, dnainfo, newParam);
    pAccept -= Colon(data, nb, augData, dnainfo, curParam);

    pAccept +=  logpriorNu1(newParam) - logpriorNu1(curParam);

    pAccept +=  log(*(newVal)) - log(*(curVal)); /* correction for lognormal */

    if (pAccept>0) r=0; else r=pAccept; /* r=log(min(1,ratioOfMHalgo))*/
    z=gsl_rng_uniform(data->rng);
    if (log(z)<=r) {
	*curVal = *newVal;
	*nbAccept +=1;
    }else{
	*newVal = *curVal;
    }
    *nbPropos+=1;

}




void moveKappa(mcmcInternals * MCMCSettings, parameters * curParam, parameters *newParam, raw_data * data, nb_data *nb, aug_data *augData, dna_dist *dnainfo, acceptance *accept, NbProposals *NbProp){
    /* updating kappa with Metropolis algorithm */
    /* proposal distribution = lognormal */
    double *newVal, *curVal;
    double * nbAccept;
    double * nbPropos;
    double sigmaProp;
    double r,z;
    double pAccept = 0;
    /* double (* logPrior) (parameters *); */

    curVal = &curParam->kappa;
    newVal = &newParam->kappa;
    sigmaProp = MCMCSettings->Sigma_kappa;
    nbAccept = &accept->PourcAcc_kappa;
    nbPropos = &NbProp->NbProp_kappa;

    *(newVal) = *(curVal)*gsl_ran_lognormal(data->rng,0,sigmaProp);

    pAccept += Colon(data, nb, augData, dnainfo, newParam);
    pAccept -= Colon(data, nb, augData, dnainfo, curParam);

    pAccept +=  logpriorKappa(newParam) - logpriorKappa(curParam);

    pAccept +=  log(*(newVal)) - log(*(curVal)); /* correction for lognormal */

    if (pAccept>0) r=0; else r=pAccept; /* r=log(min(1,ratioOfMHalgo))*/
    z=gsl_rng_uniform(data->rng);
    if (log(z)<=r) {
	*curVal = *newVal;
	*nbAccept +=1;
    }else{
	*newVal = *curVal;
    }
    *nbPropos+=1;

}






/* !!! This needs to be checked. !!! */
void moveTau(mcmcInternals * MCMCSettings, parameters * curParam, parameters *newParam, raw_data * data, nb_data *nb, aug_data *augData, dna_dist *dnainfo, acceptance *accept, NbProposals *NbProp){

    /* /\* updating tau with Metropolis algorithm *\/ */
    /* /\* proposal distribution = lognormal *\/ */
    /* double *newVal, *curVal; */
    /* double * nbAccept; */
    /* double * nbPropos; */
    /* double sigmaProp; */
    /* double r,z; */
    /* double pAccept = 0; */
 
    /* curVal = &curParam->tau; */
    /* newVal = &newParam->tau; */

    /* sigmaProp = MCMCSettings->Sigma_tau; */
    /* nbAccept = &accept->PourcAcc_tau; */
    /* nbPropos = &NbProp->NbProp_tau; */

    /* *(newVal) = *(curVal)*gsl_ran_lognormal(data->rng,0,sigmaProp); */

    /* pAccept += Colon(data, nb, augData, dnainfo, newParam); */
    /* pAccept -= Colon(data, nb, augData, dnainfo, curParam); */

    /* pAccept +=  logpriorTau(newParam) - logpriorTau(curParam); */

    /* pAccept +=  log(*(newVal)) - log(*(curVal)); /\* correction for lognormal *\/ */

    /* if (pAccept>0) r=0; else r=pAccept; */
    /* z=gsl_rng_uniform(data->rng); */
    /* if (log(z)<=r) { */
    /* 	*curVal = *newVal; */
    /* 	*nbAccept +=1; */
    /* }else{
	*newVal = *curVal;
    } */
    /* *nbPropos+=1; */

}





/* !!! This needs to be checked. !!! */
void moveAlpha(mcmcInternals * MCMCSettings, parameters * curParam, parameters *newParam, raw_data * data, nb_data *nb, aug_data *augData, dna_dist *dnainfo, acceptance *accept, NbProposals *NbProp){
    /* /\* updating alpha with Metropolis algorithm *\/ */
    /* /\* proposal distribution = truncated lognormal (no values >1) *\/ */


    /* double QCur, QNew; */
    /* double *newVal, *curVal; */
    /* double * nbAccept; */
    /* double * nbPropos; */
    /* double sigmaProp; */
    /* double r,z; */
    /* double pAccept = 0; */
   
    /* curVal = &curParam->alpha; */
    /* newVal = &newParam->alpha; */

    /* sigmaProp = MCMCSettings->Sigma_alpha; */
    /* nbAccept = &accept->PourcAcc_alpha; */
    /* nbPropos = &NbProp->NbProp_alpha; */

    /* do */
    /* 	{ */
    /* 	    *(newVal) = *(curVal)*gsl_ran_lognormal(data->rng,0,sigmaProp); */
    /* 	}while(*(newVal)>1); */

    /* QCur = gsl_cdf_gaussian_P(-log(*(curVal)),sigmaProp); */
    /* QNew = gsl_cdf_gaussian_P(-log(*(newVal)),sigmaProp); */

    /* pAccept += Colon(data, nb, augData, dnainfo, newParam); */
    /* pAccept -= Colon(data, nb, augData, dnainfo, curParam); */

    /* pAccept +=  logpriorAlpha(newParam) - logpriorAlpha(curParam); */

    /* pAccept +=  log(*(newVal)) - log(*(curVal)); /\* correction for lognormal *\/ */
    /* pAccept +=   log(QCur) - log(QNew); /\* correction for truncation (no values >1) *\/ */

    /* if (pAccept>0) r=0; else r=pAccept; */
    /* z=gsl_rng_uniform(data->rng); */
    /* if (log(z)<=r) { */
    /* 	*curVal = *newVal; */
    /* 	*nbAccept +=1; */
    /* } else{
	*newVal = *curVal;
    }*/
    /* *nbPropos+=1; */

}






/*********************** Move augmented data ***********************/
void moveC(int i, mcmcInternals * MCMCSettings, parameters * param, raw_data * data, nb_data *nb, aug_data *curAugData, aug_data *newAugData, dna_dist *dnainfo){
    /* updating C[i] with Metropolis algorithm */
    /* proposal distribution = uniform between [C_i-l;C_i+l] */

    int l = 1, NbPatients = data->NbPatients, T=data->T; /* half-width of the proposal interval */
    int random;
    int j,t;
    int curMinTime, curMaxTime;
    int newMinTime, newMaxTime;
    int curNbPoss;
    int newNbPoss;
    double QCur, QNew;
    int *newVal, *curVal;
    double r,z;
    double pAccept = 0;
 
    curVal = &curAugData->C[i];
    newVal = &newAugData->C[i];

    curMinTime = *curVal - l;
    curMaxTime = GSL_MIN(*curVal + l,curAugData->E[i]-1);
    curNbPoss = curMaxTime-curMinTime;

    random=gsl_rng_uniform_int (data->rng, curNbPoss);

    if(random<l){
	*newVal = curMinTime + random;
    } else {
	*newVal = curMinTime + random + 1;
    }

    if(*newVal < *curVal){
	if(data->ward[i]==0){
	    if(GSL_MAX(*newVal,0)<GSL_MIN(*curVal,T)){
		for (t=GSL_MAX(*newVal,0);t<GSL_MIN(*curVal,T);t++){
		    if(gsl_vector_get(data->IsInHosp[i],t)==1){
			newAugData->I0[t]++;
		    }
		}
	    }
	} else {
	    if(GSL_MAX(*newVal,0)<GSL_MIN(*curVal,T)){
		for (t=GSL_MAX(*newVal,0);t<GSL_MIN(*curVal,T);t++){
		    if(gsl_vector_get(data->IsInHosp[i],t)==1){
			newAugData->I1[t]++;
		    }
		}
	    }
	}
    } else {
	if(data->ward[i]==0){
	    if(GSL_MAX(*curVal,0)<GSL_MIN(*newVal,T)){
		for (t=GSL_MAX(*curVal,0);t<GSL_MIN(*newVal,T);t++){
		    if(gsl_vector_get(data->IsInHosp[i],t)==1){
			newAugData->I0[t]--;
		    }
		}
	    }
	} else {
	    if(GSL_MAX(*curVal,0)<GSL_MIN(*newVal,T)){
		for (t=GSL_MAX(*curVal,0);t<GSL_MIN(*newVal,T);t++){
		    if(gsl_vector_get(data->IsInHosp[i],t)==1){
			newAugData->I1[t]--;
		    }
		}
	    }
	}
    }

    newMinTime = *newVal - l;
    newMaxTime = GSL_MIN(*newVal + l,curAugData->E[i]-1);
    newNbPoss = newMaxTime-newMinTime;

    QCur = 1.0/newNbPoss;
    QNew = 1.0/curNbPoss;

    pAccept += ObsLevelPerCase(i, data, nb, newAugData, param) - ObsLevelPerCase(i, data, nb, curAugData, param);

    pAccept += DurationColonPerCase (i, newAugData, param) - DurationColonPerCase (i, curAugData, param);

    for(j=0;j<NbPatients;j++) /* only patients who are colonised after C_i are affected by the change in C_i */
	{
	    if(curAugData->C[j]>=GSL_MIN(*newVal,*curVal)){
		pAccept += ColonPerCase(j,data, nb, newAugData, dnainfo, param) - ColonPerCase(j,data, nb, curAugData, dnainfo, param);
	    }
	} /* this loop should be equivalent to pAccept += Colon(data, nb, newAugData, param) - Colon (data, nb, curAugData, param); */

    pAccept +=   log(QCur) - log(QNew); /* correction for truncation (C_i has to be < E_i) */

    /* pAccept2 = fullLoglikelihoodWithPrior(data, nb, newAugData, param) - fullLoglikelihoodWithPrior(data, nb, curAugData, param); */
    /* pAccept2 +=   log(QCur) - log(QNew); /\* correction for truncation (C_i has to be < E_i) *\/ */

    /* printf("%d\t%lg\t%lg\n",i,pAccept,pAccept2); */
    /* fflush(stdout); */
    /*if(abs(pAccept-pAccept2)>0.01){
      printf("difference pAccept pAccept2\n");
      fflush(stdout);
      system("pause");
      }
      if(pAccept>0){
      printf("pAccept>0\n");
      fflush(stdout);
      system("pause");
      }*/


    if (pAccept>0) r=0; else r=pAccept;
    z=gsl_rng_uniform(data->rng);
    if (log(z)<=r) {
	copyAugData(curAugData,newAugData);
    } else {
	copyAugData(newAugData,curAugData);
    }

}






void moveE(int i, mcmcInternals * MCMCSettings, parameters * param, raw_data * data, nb_data *nb, aug_data *curAugData, aug_data *newAugData, dna_dist *dnainfo){
    /* updating E[i] with Metropolis algorithm */
    /* proposal distribution = uniform between [E_i-l;E_i+l] */

    int l = 1, NbPatients = data->NbPatients, T = data->T; /* half-width of the proposal interval */
    int random;
    int j,t;
    int curMinTime, curMaxTime;
    int newMinTime, newMaxTime;
    int curNbPoss;
    int newNbPoss;
    double QCur, QNew;
    int *newVal, *curVal;
    double r,z;
    double pAccept = 0;
    
    /* printf("\naaa\n");fflush(stdout); */

    curVal = &curAugData->E[i];
    newVal = &newAugData->E[i];

    curMinTime = GSL_MAX(*curVal - l,curAugData->C[i]+1);
    curMaxTime = *curVal + l;
    curNbPoss = curMaxTime-curMinTime;

    random=gsl_rng_uniform_int (data->rng, curNbPoss);
    /* printf("\nbbb\n");fflush(stdout); */

    if(random<l){
	*newVal = curMaxTime - random;
    } else {
	*newVal = curMaxTime - random - 1;
    }
    /* printf("\nbbba\n");fflush(stdout); */

    if(*newVal < *curVal){
	if(data->ward[i]==0){
	    if(GSL_MAX(*newVal,0)<GSL_MIN(*curVal,T)){
		for (t=GSL_MAX(*newVal,0);t<GSL_MIN(*curVal,T);t++){
		    if(gsl_vector_get(data->IsInHosp[i],t)==1){
			newAugData->I0[t] = newAugData->I0[t]-1;
		    }
		}
	    }   
	    /* printf("\nbbbb\n");fflush(stdout); */

	} else {
	    if(GSL_MAX(*newVal,0)<GSL_MIN(*curVal,T)){
		for (t=GSL_MAX(*newVal,0);t<GSL_MIN(*curVal,T);t++){
		    if(gsl_vector_get(data->IsInHosp[i],t)==1){
			newAugData->I1[t] = newAugData->I1[t] -1;
		    }
		}
	    }
	}
    } else {
	if(data->ward[i]==0){
	    /* printf("\nbbbc\n");fflush(stdout); */
	    if(GSL_MAX(*curVal,0)<GSL_MIN(*newVal,T)){
	    /* printf("\nbite\n");fflush(stdout); */

		for (t=GSL_MAX(*curVal,0);t<GSL_MIN(*newVal,T);t++){
	    /* printf("\npoil\n");fflush(stdout); */
	    /* printf("\nt:%d i=%d\n", t, i); */
		    if(gsl_vector_get(data->IsInHosp[i],t)==1){
	    /* printf("\nzob\n");fflush(stdout); */
	    /* 		printf("\nt:%d\n", t); */
			newAugData->I0[t] = newAugData->I0[t]+1;
		    }
		}
	    }
	} else {
	    /* printf("\nbbbd\n");fflush(stdout); */

	    if(GSL_MAX(*curVal,0)<GSL_MIN(*newVal,T)){
		for (t=GSL_MAX(*curVal,0);t<GSL_MIN(*newVal,T);t++){
		    if(gsl_vector_get(data->IsInHosp[i],t)==1){
			newAugData->I1[t] = newAugData->I1[t]+1;
		    }
		}
	    }
	}
    }
	    /* printf("\nbbbe\n");fflush(stdout); */

    newMinTime = GSL_MAX(*newVal - l,curAugData->C[i]+1);
    newMaxTime = *newVal + l;
    newNbPoss = newMaxTime-newMinTime;
    /* printf("\nccc\n");fflush(stdout); */

    QCur = 1.0/newNbPoss;
    QNew = 1.0/curNbPoss;

    pAccept += ObsLevelPerCase(i, data, nb, newAugData, param) - ObsLevelPerCase(i, data, nb, curAugData, param);

    pAccept += DurationColonPerCase (i, newAugData, param) - DurationColonPerCase (i, curAugData, param);
    /* printf("\nddd\n");fflush(stdout); */

    for(j=0;j<NbPatients;j++) /* only patients who are colonised after E_i are affected by the change in E_i */
	{
	    if(curAugData->C[j]>=GSL_MIN(*newVal,*curVal)){
		pAccept += ColonPerCase(j,data, nb, newAugData, dnainfo, param) - ColonPerCase(j,data, nb, curAugData, dnainfo, param);
	    }
	} /* this loop should be equivalent to pAccept += Colon(data, nb, newAugData, param) - Colon (data, nb, curAugData, param); */

    pAccept +=   log(QCur) - log(QNew); /* correction for truncation (E_i has to be > C_i) */

    if (pAccept>0) r=0; else r=pAccept;
    z=gsl_rng_uniform(data->rng);
    /* printf("\neee\n");fflush(stdout); */

    if (log(z)<=r) {
	copyAugData(curAugData,newAugData);
    } else {
	copyAugData(newAugData,curAugData);
    }
   
}






void moveCandE(int i, mcmcInternals * MCMCSettings, parameters * param, raw_data * data, nb_data *nb, aug_data *curAugData, aug_data *newAugData, dna_dist *dnainfo){
    /* updating C[i] and E[i] with Metropolis algorithm */
    /* proposal distribution = C_i uniform between [C_i-l;C_i+l] and E_i shifted so that E_i - C_i constant */

    int l = 1, NbPatients = data->NbPatients, T = data->T; /* half-width of the proposal interval */
    int random;
    int j,t;
    int *newValC, *curValC;
    int *newValE, *curValE;
    double r,z;
    double pAccept = 0;
   
    curValC = &curAugData->C[i];
    newValC = &newAugData->C[i];

    curValE = &curAugData->E[i];
    newValE = &newAugData->E[i];

    random=gsl_rng_uniform_int (data->rng, 2*l);

    if(random<l){
	*newValC = *curValC - l + random;
	*newValE = *curValE - l + random;
    } else {
	*newValC = *curValC - l + random + 1;
	*newValE = *curValE - l + random + 1;
    }

    if(*newValC < *curValC){
	if(data->ward[i]==0){
	    if(GSL_MAX(*newValC,0)<GSL_MIN(*curValC,T)){
		for (t=GSL_MAX(*newValC,0);t<GSL_MIN(*curValC,T);t++){
		    if(gsl_vector_get(data->IsInHosp[i],t)==1){
			newAugData->I0[t]++;
		    }
		}
	    }
	} else {
	    if(GSL_MAX(*newValC,0)<GSL_MIN(*curValC,T)){
		for (t=GSL_MAX(*newValC,0);t<GSL_MIN(*curValC,T);t++){
		    if(gsl_vector_get(data->IsInHosp[i],t)==1){
			newAugData->I1[t]++;
		    }
		}
	    }
	}
    } else {
	if(data->ward[i]==0){
	    if(GSL_MAX(*curValC,0)<GSL_MIN(*newValC,T)){
		for (t=GSL_MAX(*curValC,0);t<GSL_MIN(*newValC,T);t++){
		    if(gsl_vector_get(data->IsInHosp[i],t)==1){
			newAugData->I0[t]--;
		    }
		}
	    }
	} else {
	    if(GSL_MAX(*curValC,0)<GSL_MIN(*newValC,T)){
		for (t=GSL_MAX(*curValC,0);t<GSL_MIN(*newValC,T);t++){
		    if(gsl_vector_get(data->IsInHosp[i],t)==1){
			newAugData->I1[t]--;
		    }
		}
	    }
	}
    }

    if(*newValE < *curValE){
	if(data->ward[i]==0){
	    if(GSL_MAX(*newValE,0)<GSL_MIN(*curValE,T)){
		for (t=GSL_MAX(*newValE,0);t<GSL_MIN(*curValE,T);t++){
		    if(gsl_vector_get(data->IsInHosp[i],t)==1){
			newAugData->I0[t]--;
		    }
		}
	    }
	} else {
	    if(GSL_MAX(*newValE,0)<GSL_MIN(*curValE,T)){
		for (t=GSL_MAX(*newValE,0);t<GSL_MIN(*curValE,T);t++){
		    if(gsl_vector_get(data->IsInHosp[i],t)==1){
			newAugData->I1[t]--;
		    }
		}
	    }
	}
    } else {
	if(data->ward[i]==0){
	    if(GSL_MAX(*curValE,0)<GSL_MIN(*newValE,T)){
		for (t=GSL_MAX(*curValE,0);t<GSL_MIN(*newValE,T);t++){
		    if(gsl_vector_get(data->IsInHosp[i],t)==1){
			newAugData->I0[t]++;
		    }
		}
	    }
	} else {
	    if(GSL_MAX(*curValE,0)<GSL_MIN(*newValE,T)){
		for (t=GSL_MAX(*curValE,0);t<GSL_MIN(*newValE,T);t++){
		    if(gsl_vector_get(data->IsInHosp[i],t)==1){
			newAugData->I1[t]++;
		    }
		}
	    }
	}
    }

    pAccept += ObsLevelPerCase(i, data, nb, newAugData, param) - ObsLevelPerCase(i, data, nb, curAugData, param);
    for(j=0;j<NbPatients;j++){ /* only patients who are colonised after C_i are affected by the change in C_i */

	if(curAugData->C[j]>=GSL_MIN(*newValC,*curValC)){
	    pAccept += ColonPerCase(j,data, nb, newAugData, dnainfo, param) - ColonPerCase(j,data, nb, curAugData, dnainfo, param);
	}
    } /* this loop should be equivalent to pAccept += Colon(data, nb, newAugData, param) - Colon (data, nb, curAugData, param); */
    if (pAccept>0) r=0; else r=pAccept;
    z=gsl_rng_uniform(data->rng);
    if (log(z)<=r) {
	copyAugData(curAugData,newAugData);
    } else {
	copyAugData(newAugData,curAugData);
    }

}





