#include "common.h"
#include "matvec.h"
#include "genclasses.h"
#include "structures.h"


/*
  ======
   DATA
  ======
*/
data *alloc_data(int n, int length){
    /* allocate pointer */
    data *out = (data *) malloc(sizeof(data));
    if(out == NULL){
	fprintf(stderr, "\n[in: structures.c->alloc_data]\nNo memory left for creating data. Exiting.\n");
	exit(1);
    }

    /* fill in integers */
    out->n = n;
    out->length = length;

    /* dates: collection times for each sequence */
    out->dates = alloc_vec_int(n);

    /* dna: list of DNA sequences */
    out->dna = alloc_list_dnaseq(n, length);

    return out;
} /* end alloc_data */





void free_data(data *in){
    free_vec_int(in->dates);
    free_list_dnaseq(in->dna);
    free(in);
} /* end free_data*/




void print_data(data *in){
    fflush(stdout);
    printf("\n= Collection dates =\n");
    print_vec_int(in->dates);
    printf("\n= Sequences =");
    print_list_dnaseq(in->dna);
    fflush(stdout);
} /* end print_data*/




/* Create a data object using inputs from R */
data * Rinput2data(unsigned char * DNAbinInput, int *Tcollec, int *n, int *length){
    int i, j, count=0;
    data * out = alloc_data(*n, *length);

    /* FILL IN DATES */
    for(i=0;i<*n;i++){
	out->dates->values[i] = Tcollec[i];
    }

    /* FILL IN DNA DATA */
    /* avoid using DNAbin2list_dnaseq here as it re-allocates memory */
    for(i=0;i<*n;i++){
	for(j=0;j<*length;j++){
	    out->dna->list[i]->seq[j] = DNAbin2char(DNAbinInput[count++]);
	}
    }

    /* RETURN */
    return out;
} /* end Rinput2data */










/*
  =======
  GENTIME
  =======
*/

/* 'trunc' is the time at which w is truncated to zero */
/* 'maxK' is the maximum number of unobserved generations, for which we need convolutions */
/* 'maxT' must be larger than the largest time difference that can be observed between two related cases */
gentime *alloc_gentime(int maxK, int trunc){
  /* allocate pointer */
    gentime *out = (gentime *) malloc(sizeof(gentime));
    if(out == NULL){
	fprintf(stderr, "\n[in: structures.c->alloc_gentime]\nNo memory left for creating gentime. Exiting.\n");
	exit(1);
    }

    out->type = 0;
    out->param1 = 0.0;
    out->param2 = 0.0;
    out->param3 = 0.0;
    out->trunc = trunc>0 ? trunc : 1; /* make sur that p(0) is not zero */
    out->maxK = maxK>0 ? maxK : 1;

    /* allocate vector of densities */
    out->dens = alloc_mat_double(out->maxK, out->trunc*out->maxK);

    /* return */
    return out;
} /* end alloc_gentime */




void free_gentime(gentime *in){
    free_mat_double(in->dens);
    free(in);
} /* end free_gentime*/




void print_gentime(gentime *in){
    fflush(stdout);
    printf("\n= Description of generation time function =\n");
    printf("\nType: %d", in->type);
    printf("\nparam1: %.5f \tparam2: %.5f \tparam3: %.5f", in->param1, in->param2, in->param3);
    printf("\n= Pre-computed density (truncated to 0 at %d)=\n",in->trunc);
    print_mat_double(in->dens);
    fflush(stdout);
} /* end print_gentime*/




/* get density from generation time funtion at time 't' with 'kappa_i' generations*/
double gentime_dens(gentime *in, int t, int kappa_i){
    /* error if requested kappa_i does not exist */
    if(kappa_i > in->maxK){
	fprintf(stderr, "\n[in: structures.c->gentime_dens]\nTrying to get density for %d generations (max: %d). Exiting.\n", kappa_i, in->maxK);
	fflush(stdout);
	exit(1);
    }

    /* error if requested time too large */
    if(t >= in->maxK*in->trunc){
	fprintf(stderr, "\n[in: structures.c->gentime_dens]\nTrying to get density for %d time units (max: %d). Exiting.\n", t, in->maxK*in->trunc);
	fflush(stdout);
	exit(1);
    }

    /* otherwise fetch density value */
    if(t >= in->trunc || t < 0) return 0.0;

    double out=mat_double_ij(in->dens, kappa_i, t);
    return out;
}



/*
 =======
  PARAM
 =======
*/

param *alloc_param(int n){
  /* allocate pointer */
    param *out = (param *) malloc(sizeof(param));
    if(out == NULL){
	fprintf(stderr, "\n[in: structures.c->alloc_param]\nNo memory left for creating param. Exiting.\n");
	exit(1);
    }

    /* fill in integers */
    out->n = n;

    /* allocates vectors of integers */
    out->Tinf = alloc_vec_int(n);
    out->alpha = alloc_vec_int(n);
    out->kappa = alloc_vec_int(n);

    /* fill in doubles */
    out->mu1 = 0.0001;
    out->gamma = 1.0;
    out->pi = 1.0;

    /* return */
    return out;
} /* end alloc_param */




void free_param(param *in){
    free_vec_int(in->Tinf);
    free_vec_int(in->alpha);
    free_vec_int(in->kappa);
    free(in);
} /* end free_param*/




void print_param(param *in){
    fflush(stdout);
    printf("\n= Tinf (infection dates) =\n");
    print_vec_int(in->Tinf);
    printf("\n= Alpha_i (ancestries) =\n");
    print_vec_int(in->alpha);
    printf("\n= Kappa_i (generations from nearest ancestor) =\n");
    print_vec_int(in->kappa);
    printf("\n= mu1, mu2, gamma (transi, transver, coef) =\n");
    printf("%.5f   %.5f   %.5f", in->mu1, in->gamma*in->mu1, in->gamma);
    printf("\n= pi (proportion of observed cases =\n");
    printf("%.5f", in->pi);
    fflush(stdout);
} /* end print_param*/




void copy_param(param *in, param *out){
    /* copy atomic values */
    out->n = in->n;
    out->mu1 = in->mu1;
    out->gamma = in->gamma;
    out->pi = in->pi;
    copy_vec_int(in->Tinf,out->Tinf);
    copy_vec_int(in->alpha,out->alpha);
    copy_vec_int(in->kappa,out->kappa);
} /* end copy_param */









/*
 ============
  MCMC_PARAM
 ============
*/

mcmc_param *alloc_mcmc_param(int n){
  /* allocate pointer */
    mcmc_param *out = (mcmc_param *) malloc(sizeof(mcmc_param));
    if(out == NULL){
	fprintf(stderr, "\n[in: structures.c->alloc_mcmc_param]\nNo memory left for creating mcmc_param. Exiting.\n");
	exit(1);
    }

    /* DETERMINE THE NUMBER OF Tinf */
    /* set to N/10, minimum 1 */
    out->n_move_Tinf = (int) n/10;
    out->n_move_Tinf = out->n_move_Tinf < 1 ? 1 : out->n_move_Tinf;


    /* DETERMINE THE NUMBER OF KAPPA AND ALPHA TO MOVE */
    /* set to N/10, minimum 1 */
    out->n_move_alpha = (int) n/10;
    out->n_move_alpha = out->n_move_alpha < 1 ? 1 : out->n_move_alpha;
    out->n_move_kappa = out->n_move_alpha;

    /* ALLOCATE VECTORS */
    out->idx_move_Tinf = alloc_vec_int(out->n_move_Tinf);
    out->idx_move_alpha = alloc_vec_int(out->n_move_alpha);
    out->idx_move_kappa = alloc_vec_int(out->n_move_kappa);
    out->proba_vec = alloc_vec_double(n);
    out->all_idx = alloc_vec_int(n);


    /* FILL OUT INTEGERS */
    out->n_reject = 0;
    out->n_accept_mu1 = 0;
    out->n_reject_mu1 = 0;
    out->n_accept_gamma = 0;
    out->n_reject_gamma = 0;
    out->n_accept_Tinf = 0;
    out->n_reject_Tinf = 0;
    out->n_accept_alpha = 0;
    out->n_reject_alpha = 0;
    out->n_accept_kappa = 0;
    out->n_reject_kappa = 0;
    out->n_accept = 0;


    /* FILL IN DOUBLES */
    out->sigma_mu1 = 0.0;
    out->sigma_gamma = 0.0;


    /* RETURN */
    return out;
} /* end alloc_mcmc_param */




void free_mcmc_param(mcmc_param *in){
    free_vec_int(in->idx_move_Tinf);
    free_vec_int(in->idx_move_alpha);
    free_vec_int(in->idx_move_kappa);
    free_vec_int(in->all_idx);
    free_vec_int(in->candid_ances);
    free(in);
} /* end free_mcmc_param*/




void print_mcmc_param(mcmc_param *in){
    fflush(stdout);
    printf("\nsigma for mu1: %.10f",in->sigma_mu1);
    printf("\nsigma for gamma: %.10f",in->sigma_gamma);
    printf("\nnb moves for Tinf: %d",in->n_move_Tinf);
    printf("\nnb moves for alpha: %d",in->n_move_alpha);
    printf("\nnb moves for kappa: %d",in->n_move_kappa);

    printf("\nglobal nb. accepted: %d   nb. rejected: %d   (acc/rej ratio:%.3f)", in->n_accept, in->n_reject, (double) in->n_accept / in->n_reject);

    printf("\nmu1: nb. accepted: %d   nb. rejected: %d   (acc/rej ratio:%.3f)", in->n_accept_mu1, in->n_reject_mu1, (double) in->n_accept_mu1 / in->n_reject_mu1);

    printf("\ngamma: nb. accepted: %d   nb. rejected: %d   (acc/rej ratio:%.3f)", in->n_accept_gamma, in->n_reject_gamma, (double) in->n_accept_gamma / in->n_reject_gamma);

    printf("\nTinf: nb. accepted: %d   nb. rejected: %d   (acc/rej ratio:%.3f)", in->n_accept_Tinf, in->n_reject_Tinf, (double) in->n_accept_Tinf / in->n_reject_Tinf);

    printf("\nalpha: nb. accepted: %d   nb. rejected: %d   (acc/rej ratio:%.3f)", in->n_accept_alpha, in->n_reject_alpha, (double) in->n_accept_alpha / in->n_reject_alpha);

    printf("\nkappa: nb. accepted: %d   nb. rejected: %d   (acc/rej ratio:%.3f)", in->n_accept_kappa, in->n_reject_kappa, (double) in->n_accept_kappa / in->n_reject_kappa);

    printf("\nIndices of Tinf_i to move:\n");
    print_vec_int(in->idx_move_Tinf);
    printf("\nIndices of alpha_i to move:\n");
    print_vec_int(in->idx_move_alpha);
    printf("\nIndices of kappa_i to move:\n");
    print_vec_int(in->idx_move_kappa);

    printf("\nVector of all indices (0:(n-1)):\n");
    print_vec_int(in->all_idx);

    printf("\nVector of candidate ancestors:\n");
    print_vec_int(in->proba_vec);

    fflush(stdout);
} /* end print_mcmc_param*/






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



/* int main(){ */
/*     int i; */

/*     /\* data *\/ */
/*     data * dat = alloc_data(10,100); */
/*     printf("\nData\n"); */
/*     print_data(dat); */
/*     free_data(dat); */

/*     /\* gentime *\/ */
/*     gentime * gen = alloc_gentime(5, 20); */
/*     printf("\nGentime\n"); */
/*     print_gentime(gen); */
/*     printf("\nDensity for kappa=1, first 20 values:\n"); */
/*     for(i=0;i<20;i++) printf("%.6f ", gentime_dens(gen, i, 1)); */
/*     printf("\nDensity for kappa=2, first 20 values:\n"); */
/*     for(i=0;i<20;i++) printf("%.6f ", gentime_dens(gen, i, 2)); */
/*     printf("\nDensity for kappa=3, first 20 values:\n"); */
/*     for(i=0;i<20;i++) printf("%.6f ", gentime_dens(gen, i, 3)); */
/*     free_gentime(gen); */

/*     /\* param *\/ */
/*     param * par = alloc_param(10); */
/*     printf("\nParam\n"); */
/*     print_param(par); */
/*     free_param(par); */

/*     /\* /\\* mcmcParam *\\/ *\/ */
/*     mcmc_param * mcmcPar = alloc_mcmc_param(10); */
/*     printf("\nMcmcParam\n"); */
/*     print_mcmc_param(mcmcPar); */
/*     free_mcmc_param(mcmcPar); */


/*     printf("\n\n"); */
/*     return 0; */
/* } */



/*
   gcc instructions:

   gcc -o structures matvec.c genclasses.c structures.c -lgsl -lgslcblas -g -Wall

  ./structures

   valgrind --leak-check=full --track-origins=yes -v structures

*/
