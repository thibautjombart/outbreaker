
#include "common.h"
#include "structures.h"
#include "matvec.h"
#include "genclasses.h"
#include "distances.h"
#include "init.h"
#include "prior.h"
#include "likelihood.h"
#include "moves.h"
#include "mcmc.h"



/*
  ===================
  AUXILIARY FUNCTIONS
  ===================
*/

/* parameters are output in the following order:
   chain-number, posterior, likelihood, prior, mu1, mu2, gamma, pi, phi, Tinf_1, ..., Tinf_n, alpha_1, ..., alpha_n, kappa_1, ..., kappa_n

   notes:
   - the output text file ("output.txt") is tab-delimited
   - indices are provided from 1 to n, i.e. not as C indices (from 0 to n-1)
*/

void fprint_chains(FILE *file, data *dat, dna_dist *dnainfo, gentime *gen, param *par, int step, gsl_rng *rng, bool quiet){
    int i;
    double like, prior;

    /* OUTPUT TO FILE */
    /* chain number */
    fprintf(file,"\n%d", step);

    /* posterior, likelihood, prior */
    like = loglikelihood_all(dat, dnainfo, gen, par, rng);
    prior = logprior_all(par);
    fprintf(file,"\t%.15f", like + prior);
    fprintf(file,"\t%.15f", like);
    fprintf(file,"\t%.15f", prior);

    /* parameters */
    fprintf(file,"\t%.15f", par->mu1);
    fprintf(file,"\t%.15f", par->mu1 * par->gamma);
    fprintf(file,"\t%.15f", par->gamma);
    fprintf(file,"\t%.15f", par->pi);
    /* fprintf(file,"\t%.15f", par->phi); */
    for(i=0;i<par->Tinf->length;i++){
	fprintf(file, "\t%d", vec_int_i(par->Tinf, i));
    }
    for(i=0;i<par->Tinf->length;i++){
	fprintf(file, "\t%d", vec_int_i(par->alpha, i)+1);
    }
    for(i=0;i<par->Tinf->length;i++){
	fprintf(file, "\t%d", vec_int_i(par->kappa, i));
    }

    /* OUTPUT TO SCREEN */
    if(!quiet){
	printf("\n%d\t", step);
	fprintf(file,"\t%.15f", like*prior);
	fprintf(file,"\t%.15f", like);
	fprintf(file,"\t%.15f", prior);
	printf("\t%.15f", par->mu1);
	printf("\t%.15f", par->mu1 * par->gamma);
	printf("\t%.15f", par->gamma);
	printf("\t%.15f", par->pi);
	/* printf("\t%.15f", par->phi); */
	for(i=0;i<par->Tinf->length;i++){
	    printf("\t%d", vec_int_i(par->Tinf, i));
	}
	for(i=0;i<par->Tinf->length;i++){
	    printf("\t%d", vec_int_i(par->alpha, i)+1);
	}
	for(i=0;i<par->Tinf->length;i++){
	    printf("\t%d", vec_int_i(par->kappa, i));
	}
    }
} /* end fprint_chains */






/* print mcmc parameter (e.g. acceptance/rejection) to file 
   order is as follows:
   step | global prop accept | accept_mu1 | sigma_mu1 | sigma_gamma
*/
void fprint_mcmc_param(FILE *file, mcmc_param *mcmcPar, int step){
    double temp=0.0;
    /* OUTPUT TO FILE */
    fprintf(file,"\n%d\t", step);
    temp = (double) mcmcPar->n_accept_mu1 / (double) (mcmcPar->n_accept_mu1 + mcmcPar->n_reject_mu1);
    fprintf(file,"\t%.5f", temp);
    temp = (double) mcmcPar->n_accept_gamma / (double) (mcmcPar->n_accept_gamma + mcmcPar->n_reject_gamma);
    fprintf(file,"\t%.5f", temp);
    temp = (double) mcmcPar->n_accept_pi / (double) (mcmcPar->n_accept_pi + mcmcPar->n_reject_pi);
    fprintf(file,"\t%.5f", temp);
    /* temp = (double) mcmcPar->n_accept_phi / (double) (mcmcPar->n_accept_phi + mcmcPar->n_reject_phi); */
    /* fprintf(file,"\t%.5f", temp); */
    temp = (double) mcmcPar->n_accept_Tinf / (double) (mcmcPar->n_accept_Tinf + mcmcPar->n_reject_Tinf);
    fprintf(file,"\t%.5f", temp);
    fprintf(file,"\t%.15f", mcmcPar->sigma_mu1);
    fprintf(file,"\t%.15f", mcmcPar->sigma_gamma);
    fprintf(file,"\t%.15f", mcmcPar->sigma_pi);
    /* fprintf(file,"\t%.15f", mcmcPar->sigma_phi); */
    fprintf(file,"\t%d", mcmcPar->n_like_zero);
}









/*
   ================
   TUNING FUNCTIONS
   ================
*/

/*
   AIM: get ~40% acceptance for univariate param, ~20% for multivariate param
*/
void tune_mu1(mcmc_param * in, gsl_rng *rng){
    /* get acceptance proportion */
    double paccept = (double) in->n_accept_mu1 / (double) (in->n_accept_mu1 + in->n_reject_mu1);

    /* acceptable zone: 35-45% acceptance */
    if(paccept<0.25) {
	in->sigma_mu1 /= 1.5;
	in->n_accept_mu1 = 0;
	in->n_reject_mu1 = 0;
    } else if (paccept>0.50) {
	in->sigma_mu1 *= 1.5;
	in->n_accept_mu1 = 0;
	in->n_reject_mu1 = 0;
	/* do not allow sigma to be > 1 (for lognormal not to go crazy) */
	if(in->sigma_mu1>1.0){
	    in->sigma_mu1 = 1.0;
	    in->tune_mu1 = FALSE;
	}
    } else {
	in->tune_mu1 = FALSE;
    }
}





void tune_gamma(mcmc_param * in, gsl_rng *rng){
    /* get acceptance proportion */
    double paccept = (double) in->n_accept_gamma / (double) (in->n_accept_gamma + in->n_reject_gamma);

    /* acceptable zone: 35-45% acceptance */
    if(paccept<0.25) {
	in->sigma_gamma /= 1.5;
	in->n_accept_gamma = 0;
	in->n_reject_gamma = 0;
    } else if (paccept>0.50) {
	in->sigma_gamma *= 1.5;
	in->n_accept_gamma = 0;
	in->n_reject_gamma = 0;
	/* do not allow sigma to be > 1 (for lognormal not to go crazy) */
	if(in->sigma_gamma>1.0){
	    in->sigma_gamma = 1.0;
	    in->tune_gamma = FALSE;
	}
    } else {
	in->tune_gamma = FALSE;
    }
}





void tune_pi(mcmc_param * in, gsl_rng *rng){
    /* get acceptance proportion */
    double paccept = (double) in->n_accept_pi / (double) (in->n_accept_pi + in->n_reject_pi);

    /* acceptable zone: 35-45% acceptance */
    if(paccept<0.25) {
	in->sigma_pi /= 1.5;
	in->n_accept_pi = 0;
	in->n_reject_pi = 0;
    } else if (paccept>0.50) {
	in->sigma_pi *= 1.5;
	in->n_accept_pi = 0;
	in->n_reject_pi = 0;
	/* do not allow sigma to be > 1 (for lognormal not to go crazy) */
	if(in->sigma_pi>1.0){
	    in->sigma_pi = 1.0;
	    in->tune_pi = FALSE;
	}
    } else {
	in->tune_pi = FALSE;
    }
}






/* void tune_phi(mcmc_param * in, gsl_rng *rng){ */
/*     /\* get acceptance proportion *\/ */
/*     double paccept = (double) in->n_accept_phi / (double) (in->n_accept_phi + in->n_reject_phi); */

/*     /\* acceptable zone: 35-45% acceptance *\/ */
/*     if(paccept<0.25) { */
/* 	in->sigma_phi /= 1.5; */
/* 	in->n_accept_phi = 0; */
/* 	in->n_reject_phi = 0; */
/*     } else if (paccept>0.50) { */
/* 	in->sigma_phi *= 1.5; */
/* 	in->n_accept_phi = 0; */
/* 	in->n_reject_phi = 0; */
/* 	/\* do not allow sigma to be > 1 (for lognormal not to go crazy) *\/ */
/* 	if(in->sigma_phi>1.0){ */
/* 	    in->sigma_phi = 1.0; */
/* 	    in->tune_phi = FALSE; */
/* 	} */
/*     } else { */
/* 	in->tune_phi = FALSE; */
/*     } */
/* } */




/* void tune_Tinf(mcmc_param * in, gsl_rng *rng){ */
/*     /\* get acceptance proportion *\/ */
/*     double paccept = (double) in->n_accept_Tinf / (double) (in->n_accept_Tinf + in->n_reject_Tinf); */

/*     /\* Note: Tinf treated as univariate as each value is accepted/rejected independently *\/ */
/*     /\* acceptable zone: 35-45% acceptance *\/ */
/*     if(paccept<0.35) { */
/* 	in->lambda_Tinf /= 1.5; */
/*     } else if (paccept>0.45) in->lambda_Tinf *= 1.5; */
/* } */






/*
  ===============================================
  METROPOLIS-HASTING ALGORITHM FOR ALL PARAMETERS
  ===============================================
*/

/* PRELIM MCMC FOR FINDING OUTLIERS */
void mcmc_find_import(vec_int *areOutliers, int outEvery, int tuneEvery, bool quiet, param *par, 
		      data *dat, dna_dist *dnainfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng){

    int i, j, nbTermsLike = 0;
    double medLogLike = 0.0;

    /* OUTPUT TO SCREEN - HEADER */
    if(!quiet){
	printf("step\tpost\tlike\tprior\tmu1\tmu2\tgamma\tpi");
	for(i=0;i<dat->n;i++){
	    printf("\tTinf_%d", i+1);
	}
	for(i=0;i<dat->n;i++){
	    printf("\talpha_%d", i+1);
	}
	for(i=0;i<dat->n;i++){
	    printf("\tkappa_%d", i+1);
	}
    }


    /* CREATE TEMPORARY PARAMETERS */
    /* ! do not alter 'par' or mcmcPar !*/
    param *localPar = alloc_param(dat->n), *tempPar = alloc_param(dat->n);
    copy_param(par,localPar);
    copy_param(par,tempPar);

    mcmc_param *localMcmcPar = alloc_mcmc_param(dat->n);
    copy_mcmc_param(mcmcPar, localMcmcPar);

    /* CREATE TEMPORARY VECTOR STORING INDIVIDUAL LIKELIHOODS */
    vec_double *indivLogLike = alloc_vec_double(dat->n);

    /* RUN MCMC */
    for(i=2;i<=localMcmcPar->find_import_at;i++){
	/* printf("i: %d ",i);fflush(stdout); */
	/* COLLECT INFORMATION ABOUT INDIVIDUAL LIKELIHOODS */
	if(i>=localMcmcPar->burnin && i % outEvery == 0){
	    /* printf("\ni=%d - computing individual likelihoods\n",i);fflush(stdout); */
	    for(j=0;j<dat->n;j++){
		indivLogLike->values[j] += loglikelihood_gen_i(j,dat, dnainfo, localPar, rng);
	    }
	    /* printf("\nlikelihood vector:\n");fflush(stdout); */
	    /* print_vec_double(indivLogLike); */
	    nbTermsLike++;
	}

	/* TUNING */
	if(i % tuneEvery == 0 && localMcmcPar->tune_all){
	    tune_mu1(localMcmcPar,rng);
	    tune_gamma(localMcmcPar,rng);
	    tune_pi(localMcmcPar,rng);
	    /* tune_phi(localMcmcPar,rng); */
	    /* localMcmcPar->tune_all = localMcmcPar->tune_mu1 || localMcmcPar->tune_gamma || localMcmcPar->tune_pi || localMcmcPar->tune_phi; */
	    localMcmcPar->tune_all = localMcmcPar->tune_mu1 || localMcmcPar->tune_gamma || localMcmcPar->tune_pi;
	}

	/* MOVEMENTS */
	/* move mutation rates */
	if(localMcmcPar->move_mut){
	    /* move mu1 */
	    move_mu1(localPar, tempPar, dat, dnainfo, localMcmcPar, rng);

	    /* move gamma */
	    move_gamma(localPar, tempPar, dat, dnainfo, localMcmcPar, rng);
	}

	/* move pi */
	if(localMcmcPar->move_pi) move_pi(localPar, tempPar, dat, localMcmcPar, rng);

	/* /\* move phi *\/ */
	/* if(localMcmcPar->move_phi) move_phi(localPar, tempPar, dat, localMcmcPar, rng); */

	/* move Tinf */
	if(localMcmcPar->move_Tinf) move_Tinf(localPar, tempPar, dat, dnainfo, gen, localMcmcPar, rng);

	/* /\* move alpha_i*\/ */
	/* move_alpha(localPar, tempPar, dat, dnainfo, gen, localMcmcPar, rng); */

	/* /\* move kappa_i*\/ */
	/* move_kappa(localPar, tempPar, dat, dnainfo, gen, localMcmcPar, rng); */

	/* move alpha_i and kappa_i*/
	move_alpha_kappa(localPar, tempPar, dat, dnainfo, gen, localMcmcPar, rng);

    } /* end of MCMC */


    /* FIND IMPORTED CASES */
    /* compute individual average log-like */
    for(j=0;j<dat->n;j++){
	indivLogLike->values[j] = vec_double_i(indivLogLike,j)/((double) nbTermsLike);
    }

    /* compute general average log-like */
    medLogLike = mean_vec_double(indivLogLike);
    /* printf("\nAverage loglike: %f\n", medLogLike);fflush(stdout); */
    /* printf("\nIndividual loglike:\n");fflush(stdout); */
    /* print_vec_double(indivLogLike); */

    /* browse each likelihood, define outliers */
    /* printf("\n\nLooking for outliers...\n"); */
    for(j=0;j<dat->n;j++){
	/* outliers = likelihood 1000 times lower than the median */
	/* printf("\nIndiv %d: loglike difference= %.5f", j+1, medLogLike - vec_double_i(indivLogLike,j));fflush(stdout); */
	if((medLogLike - vec_double_i(indivLogLike,j)) > log(1000)){
	    areOutliers->values[j] = 1;
	    printf("\nIndividual %d identified as imported case\n",j+1);fflush(stdout);
	} else {
	    areOutliers->values[j] = 0;
	}
    } /* end setting outliers */


    /* FREE TEMPORARY PARAMETERS */
    free_param(localPar);
    free_param(tempPar);
    free_mcmc_param(localMcmcPar);
    free_vec_double(indivLogLike);
} /* end mcmc_find_import */








/*
   ==============
   MAIN MCMC ALGO
   ==============
*/
void mcmc(int nIter, int outEvery, char outputFile[256], char mcmcOutputFile[256], int tuneEvery, 
	  bool quiet, param *par, data *dat, dna_dist *dnainfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng){

    int i, j, nbTermsLike = 0;
    vec_int *areOutliers = alloc_vec_int(dat->n);

    /* OPEN OUTPUT FILE */
    FILE *file = fopen(outputFile,"w"), *mcmcFile = fopen(mcmcOutputFile,"w");
    if(file==NULL){
	fprintf(stderr, "\n[in: mcmc.c->mcmc]\nCannot open output file %s.\n", outputFile);
	exit(1);
    }
    if(mcmcFile==NULL){
	fprintf(stderr, "\n[in: mcmc.c->mcmc]\nCannot open output file %s.\n", mcmcOutputFile);
	exit(1);
    }


    /* OUTPUT TO OUTFILE - HEADER */
    /* fprintf(file, "step\tpost\tlike\tprior\tmu1\tmu2\tgamma\tpi\tphi"); */
    fprintf(file, "step\tpost\tlike\tprior\tmu1\tmu2\tgamma\tpi");
    for(i=0;i<dat->n;i++){
	fprintf(file, "\tTinf_%d", i+1);
    }
    for(i=0;i<dat->n;i++){
	fprintf(file, "\talpha_%d", i+1);
    }
    for(i=0;i<dat->n;i++){
	fprintf(file, "\tkappa_%d", i+1);
    }

    /* OUTPUT TO MCMCOUTFILE - HEADER */
    /* fprintf(mcmcFile, "step\tp_accept_mu1\tp_accept_gamma\tp_accept_pi\tp_accept_phi\tp_accept_Tinf"); */
    /* fprintf(mcmcFile, "\tsigma_mu1\tsigma_gamma\tsigma_pi\tsigma_phi\tn_like_zero"); */
    fprintf(mcmcFile, "step\tp_accept_mu1\tp_accept_gamma\tp_accept_pi\tp_accept_Tinf");
    fprintf(mcmcFile, "\tsigma_mu1\tsigma_gamma\tsigma_pi\tn_like_zero");


    /* OUTPUT TO SCREEN - HEADER */
    if(!quiet){
	printf("step\tpost\tlike\tprior\tmu1\tmu2\tgamma\tpi");
	for(i=0;i<dat->n;i++){
	    printf("\tTinf_%d", i+1);
	}
	for(i=0;i<dat->n;i++){
	    printf("\talpha_%d", i+1);
	}
	for(i=0;i<dat->n;i++){
	    printf("\tkappa_%d", i+1);
	}
    }

    fprint_chains(file, dat, dnainfo, gen, par, 1, rng, quiet);
    fprint_mcmc_param(mcmcFile, mcmcPar, 1);

    mcmcPar->step_notune = nIter;


    /* PRELIM STEP - FINDING OUTLIERS */
    if(mcmcPar->find_import){
	mcmc_find_import(areOutliers, outEvery, tuneEvery, quiet, par, dat, dnainfo, gen, mcmcPar, rng);

	/* RESTORE INITIAL TUNING SETTINGS AND PARAM */
	/* mcmcPar->tune_all = TRUE; */
	/* copy_param(par,tempPar); */
	/* mcmcPar->step_notune = nIter; */

	/* printf("\nContent of 'areOutliers:\n");fflush(stdout); */
	/* print_vec_int(areOutliers); */
	for(i=0;i<dat->n;i++){
	    if(vec_int_i(areOutliers,i)==1){
		par->alpha->values[i] = -1;
		par->kappa->values[i] = 1;
		mcmcPar->move_alpha->values[i] = 0.0;
	    }
	}
    } /* END PRELIM MCMC FOR FINDING OUTLIERS */


    /* CREATE TEMPORARY PARAMETERS */
    param *tempPar = alloc_param(dat->n);
    copy_param(par,tempPar);

     /* RUN MAIN MCMC */
    for(i=2;i<=nIter;i++){
	/* /\* debugging *\/ */
	/* printf("\n\n = MCMC iteration %d =\n",i); */
	/* fflush(stdout); */

	/* OUTPUT TO FILES */
	if(i % outEvery == 0){
	    fprint_chains(file, dat, dnainfo, gen, par, i, rng, quiet);
	    fprint_mcmc_param(mcmcFile, mcmcPar, i);
	}

	/* TUNING */
	if(i % tuneEvery == 0 && mcmcPar->tune_all){
	    tune_mu1(mcmcPar,rng);
	    tune_gamma(mcmcPar,rng);
	    tune_pi(mcmcPar,rng);
	    /* tune_phi(mcmcPar,rng); */
	    /* mcmcPar->tune_all = mcmcPar->tune_mu1 || mcmcPar->tune_gamma || mcmcPar->tune_pi || mcmcPar->tune_phi; */
	    mcmcPar->tune_all = mcmcPar->tune_mu1 || mcmcPar->tune_gamma || mcmcPar->tune_pi;
	    if(!mcmcPar->tune_all) {
		mcmcPar->step_notune = i;
		/* printf("\nStopped tuning at chain %d\n",i);fflush(stdout); */
	    }
	}

	/* /\* debugging *\/ */
	/* double logLike = loglikelihood_all(dat, dnainfo, gen, par); */
	/* printf("\n\n = Initial Log-likelihood value (in mcmc, before movement): %f\n", logLike); */
	/* fflush(stdout); */

	/* check_loglikelihood_all(dat, dnainfo, gen, par); */

	/* MOVEMENTS */
	if(mcmcPar->move_mut){/* move mu1 */
	    move_mu1(par, tempPar, dat, dnainfo, mcmcPar, rng);

	    /* move gamma */
	    move_gamma(par, tempPar, dat, dnainfo, mcmcPar, rng);
	}

	/* move pi */
	if(mcmcPar->move_pi) move_pi(par, tempPar, dat, mcmcPar, rng);

	/* /\* move phi *\/ */
	/* if(mcmcPar->move_phi) move_phi(par, tempPar, dat, mcmcPar, rng); */

	/* move Tinf */
	/* printf("\nTinf:"); */
	/* print_vec_int(par->Tinf); */
	/* fflush(stdout); */
	if(mcmcPar->move_Tinf) move_Tinf(par, tempPar, dat, dnainfo, gen, mcmcPar, rng);

	/* /\* move alpha_i*\/ */
	/* move_alpha(par, tempPar, dat, dnainfo, gen, mcmcPar, rng); */

	/* /\* move kappa_i*\/ */
	/* move_kappa(par, tempPar, dat, dnainfo, gen, mcmcPar, rng); */

	/* move alpha_i and kappa_i*/
	move_alpha_kappa(par, tempPar, dat, dnainfo, gen, mcmcPar, rng);

    } /* end of mcmc */


    /* CLOSE OUTPUT OUTFILE */
    fclose(file);
    fclose(mcmcFile);

    /* FREE TEMPORARY PARAMETERS */
    free_param(tempPar);
    free_vec_int(areOutliers);
} /* end mcmc */









/*
>>>> TESTING <<<<
*/

/* int main(){ */
/*   /\* DECLARATIONS *\/ */
/*     int TIMESPAN, i; */
/*     data *dat; */
/*     gentime *gen; */
/*     param *par; */
/*     dna_dist * dnainfo; */
/*     mcmc_param * mcmcPar; */

/*     double logPrior, logLike, logPost; */

/*     /\* INITIALIZE RNG *\/ */
/*     gsl_rng *rng = create_gsl_rng(time(NULL)); */


/*     /\* CONVERT DATA *\/ */
/*     dat = alloc_data(3,10); */
/*     dat->dates->values[0] = 0; */
/*     dat->dates->values[1] = 2; */
/*     dat->dates->values[1] = 3; */
/*     dat->dna->list[0]->seq[0] = 'a'; */
/*     dat->dna->list[1]->seq[0] = 'a'; */
/*     dat->dna->list[2]->seq[0] = 't'; */
/*     printf("\n>>> Data <<<\n"); */
/*     print_data(dat); */


/*     /\* GET TIME SPAN *\/ */
/*     TIMESPAN = max_vec_int(dat->dates) - min_vec_int(dat->dates); */
/*     printf("\nTimespan is %d\n",TIMESPAN); */


/*     /\* CREATE AND INIT GENERATION TIME *\/ */
/*     gen = alloc_gentime(TIMESPAN, 5); */
/*     init_gentime(gen, 1, 1.0, 0.0, 0.0); */
/*     printf("\n>>> gentime info <<<\n"); */
/*     print_gentime(gen); */
/*     printf("sizes of rows in gen: "); */
/*     for(i=0;i<gen->dens->n;i++) printf("%d ", gen->dens->rows[i]->length); */


/*      /\* CREATE AND INIT PARAMETERS *\/ */
/*     par = alloc_param(3); */
/*     par->alpha->values[0] = -1; */
/*     par->alpha->values[1] = 0; */
/*     par->alpha->values[2] = 0; */
/*     par->kappa->values[0] = 1; */
/*     par->kappa->values[1] = 1; */
/*     par->kappa->values[2] = 1; */
/*     par->mu1 = 0.0001; */
/*     par->gamma = 1.0; */
/*     par->pi = 0.5; */
/*     printf("\nParameters (par)\n"); */
/*     print_param(par); */


/*     /\* ALLOCATE MCMCPAR *\/ */
/*     mcmcPar = alloc_mcmc_param(3); */
/*     init_mcmc_param(mcmcPar, dat); */
/*     printf("\nMCMC parameters (mcmcPar)\n"); */
/*     print_mcmc_param(mcmcPar); */


/*     /\* COMPUTE GENETIC DISTANCES *\/ */
/*     dnainfo = compute_dna_distances(dat->dna); */
/*     printf("\n>>> DNA info <<<\n"); */
/*     print_dna_dist(dnainfo); */


/*     /\* COMPUTE PRIORS *\/ */
/*     logPrior = logprior_all(par); */
/*     printf("\nPrior value (log): %.10f\n", logPrior); */

/*    /\* COMPUTE LIKELIHOOD *\/ */
/*     logLike = loglikelihood_all(dat, dnainfo, gen, par); */
/*     printf("\nLog-likelihood value: %.10f\n", logLike); */

/*     /\* COMPUTE POSTERIOR *\/ */
/*     logPost = logposterior_all(dat, dnainfo, gen, par); */
/*     printf("\nLog-posterior value: %.10f\n", logPost); */


/*     /\* PROCEED TO MCMC *\/ */
/*     int nIter=10000, outEvery=100; */
/*     char outFile[256] = "output.txt"; */

/*     mcmc(nIter, outEvery, outFile, FALSE, par, dat, dnainfo, gen, mcmcPar, rng); */

/*     printf("\n\n");fflush(stdout); */

/*     /\* /\\* RUNTIME TEST *\\/ *\/ */
/*     /\* int ITER=10e6, i; *\/ */
/*     /\* time_t t1, t2; *\/ */
/*     /\* time(&t1); *\/ */
/*     /\* printf("\nRuntime (%d computations of posterior): \n", ITER); *\/ */
/*     /\* for(i=0;i<ITER;i++){ *\/ */
/*     /\* 	logPost = logposterior_all(dat, dnainfo, gen, par); *\/ */
/*     /\* } *\/ */
/*     /\* time(&t2); *\/ */
/*     /\* printf("\nellapsed time: %d seconds\n", (int) t2 - (int) t1); *\/ */


/*     /\* FREE / RETURN *\/ */
/*     gsl_rng_free(rng); */
/*     free_data(dat); */
/*     free_gentime(gen); */
/*     free_dna_dist(dnainfo); */
/*     free_param(par); */
/*     free_mcmc_param(mcmcPar); */

/*     return 0; */
/* } */





/*
  gcc instructions

  gcc -o mcmc matvec.c genclasses.c structures.c init.c distances.c prior.c likelihood.c moves.c mcmc.c -lgsl -lgslcblas -Wall -g


  gcc -o mcmc matvec.c genclasses.c structures.c init.c distances.c prior.c likelihood.c moves.c mcmc.c -lgsl -lgslcblas -Wall -O3

 ./mcmc

  valgrind --leak-check=full -v mcmc

*/

