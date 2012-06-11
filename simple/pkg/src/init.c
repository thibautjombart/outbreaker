

#include "common.h"

/* #include "init.h" */
/* #include "InputOutput.h" */
/* #include "logL.h" */
/* #include "mcmc.h" */
/* #include "moves.h" */
/* #include "tuneVariances.h" */



gsl_rng * create_gsl_rng(){
    time_t t = time(NULL); /* time in seconds, used to change the seed of the random generator */
    gsl_rng_env_setup();
    gsl_rng *out = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(out,t); /* changes the seed of the random generator */
    return out;
}


# if 0


void CalculIsInHosp(nb_data * nbData, raw_data *data){
    int i,k,t, T=data->T;
    for(i=0 ; i<nbData->NbPatients ; i++){
	for(k = 0 ; k<nbData->NbAdmissions[i] ; k++){
	    if(GSL_MAX(0,gsl_vector_get(data->A[i],k))<GSL_MIN(T,gsl_vector_get(data->D[i],k))){
		for(t=GSL_MAX(0,gsl_vector_get(data->A[i],k)) ; t<GSL_MIN(T,gsl_vector_get(data->D[i],k)) ; t++){
		    gsl_vector_set(data->IsInHosp[i],t,1);
		}
	    }
	}
    }
}





void InitAugData(parameters *param, nb_data * nbData, raw_data *data, aug_data *augData){
    int i,t, NbPatients = nbData->NbPatients, T = nbData->T;
    double a = param->mu*param->mu/(param->sigma*param->sigma);
    double b = param->sigma*param->sigma/param->mu;
    FILE *fichC, *fichE;
    int V;

    for(t=0 ; t<T ; t++){
	augData->I0[t]=0;
	augData->I1[t]=0;
    }

    /*fichC = fopen("TrueC.txt","r");
      if ( fichC == NULL ){
      printf("A problem occurred while opening TrueC.txt. Check that the file exists and is not opened.\n");
      fflush(stdout);
      exit(1);
      }

      fichE = fopen("TrueE.txt","r");
      if ( fichE == NULL ){
      printf("A problem occurred while opening TrueE.txt. Check that the file exists and is not opened.\n");
      fflush(stdout);
      exit(1);
      }*/

    for(i=0 ; i<NbPatients ; i++){
	if(nbData->NbPosSwabs[i]>0){
	    augData->C[i] = gsl_vector_get(data->P[i],0);
	    /* augData->E[i] = GSL_MAX(augData->C[i]+ceil(gsl_ran_gamma (param->rng, a, b)),1+gsl_vector_get(data->P[i],nbData->NbPosSwabs[i]-1)); */
	    augData->E[i] = augData->C[i]+ceil(param->mu);
	}else
	    {
		augData->C[i] = -100;
		/* augData->E[i] = augData->C[i]+ceil(gsl_ran_gamma (param->rng, a, b)); */
		augData->E[i] = augData->C[i]+ceil(param->mu);
	    }

	/*if( fscanf(fichC,"%d",&V) != 1){
	  printf("A problem occurred while reading the file\n");
	  fflush(stdout);
	  break;
	  }
	  augData->C[i] = V;

	  if( fscanf(fichE,"%d",&V) != 1){
	  printf("A problem occurred while reading the file\n");
	  fflush(stdout);
	  break;
	  }
	  augData->E[i] = V;*/


	if(data->ward[i]==0){
	    if(GSL_MAX(0,augData->C[i])<GSL_MIN(T,augData->E[i])){
		for(t=GSL_MAX(0,augData->C[i]) ; t<GSL_MIN(T,augData->E[i]) ; t++){
		    if(gsl_vector_get(data->IsInHosp[i],t)==1){
			augData->I0[t]++;
		    }
		}
	    }
	}
	if(data->ward[i]==1){
	    if(GSL_MAX(0,augData->C[i])<GSL_MIN(T,augData->E[i])){
		for(t=GSL_MAX(0,augData->C[i]) ; t<GSL_MIN(T,augData->E[i]) ; t++){
		    if(gsl_vector_get(data->IsInHosp[i],t)==1){
			augData->I1[t]++;
		    }
		}
	    }
	}
    }

    /* fclose(fichC); */
    /* fclose(fichE); */
}





void InitMCMCSettings(mcmcInternals *MCMCSettings){
    MCMCSettings->NbSimul = 110000;
    MCMCSettings->SubSample = 10;
    MCMCSettings->BurnIn = 10000;

    gsl_matrix_set(MCMCSettings->Sigma_beta,0,0,0.1);
    gsl_matrix_set(MCMCSettings->Sigma_beta,0,1,0.1);
    gsl_matrix_set(MCMCSettings->Sigma_beta,1,0,0.1);
    gsl_matrix_set(MCMCSettings->Sigma_beta,1,1,0.1);

    MCMCSettings->Sigma_betaWardOut=0.1;
    MCMCSettings->Sigma_betaOutOut=0.1;
    MCMCSettings->Sigma_mu=5;
    MCMCSettings->Sigma_sigma=1;
    MCMCSettings->Sigma_nu1=0.005;
    MCMCSettings->Sigma_kappa=0.005;
    MCMCSettings->Sigma_tau=10;
    MCMCSettings->Sigma_alpha=0.1;
}





void InitParam(parameters *param){
    int i,j;

    for(i=0 ; i<2 ; i++){
	for(j=0 ; j<2 ; j++){
	    gsl_matrix_set(param->beta,i,j,0.1);
	}
    }

    param->betaWardOut=0.1;
    param->betaOutOut=0.1;

    /* param->Sp = 1.0; */
    param->Se = 0.9;

    param->Pi = 0.1;

    param->mu = 5;
    param->sigma=1;

    param->nu1=1e-6;
    param->kappa=1.0;

}



#endif
