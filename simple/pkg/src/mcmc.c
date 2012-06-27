
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
   iteration-number, mu1, mu2, gamma, Tinf_1, ..., Tinf_n, alpha_1, ..., alpha_n, kappa_1, ..., kappa_n

   notes:
   - the output text file ("output.txt") is tab-delimited
   - indices are provided from 1 to n, i.e. not as C indices (from 0 to n-1)
*/

void fprint_param(FILE *file, param *par, int step, bool quiet){
    int i;

    /* OUTPUT TO FILE */
    fprintf(file,"\n%d\t", step);
    fprintf(file,"\t%.15f", par->mu1);
    fprintf(file,"\t%.15f", par->mu1 * par->gamma);
    fprintf(file,"\t%.15f", par->gamma);
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
	printf("\t%.15f", par->mu1);
	printf("\t%.15f", par->mu1 * par->gamma);
	printf("\t%.15f", par->gamma);
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
} /* end fprint_param */






/* print mcmc parameter (e.g. acceptance/rejection) to file 
   order is as follows:
   step | global prop accept | accept_mu1 | sigma_mu1 | sigma_gamma
*/
void fprint_mcmc_param(FILE *file, mcmc_param *mcmcPar, int step){
    double temp=0.0;
    /* OUTPUT TO FILE */
    fprintf(file,"\n%d\t", step);
    fprintf(file,"\t%.lf", update_get_accept_rate(mcmcPar));
    temp = (double) mcmcPar->n_accept_mu1 / (double) (mcmcPar->n_accept_mu1 + mcmcPar->n_reject_mu1);
    fprintf(file,"\t%.5f", temp);
    temp = (double) mcmcPar->n_accept_gamma / (double) (mcmcPar->n_accept_gamma + mcmcPar->n_reject_gamma);
    fprintf(file,"\t%.5f", temp);
    fprintf(file,"\t%.15f", mcmcPar->sigma_mu1);
    fprintf(file,"\t%.15f", mcmcPar->sigma_gamma);
    fprintf(file,"\t%d", mcmcPar->n_like_zero);
}






/*
   UPDATE GLOBAL ACCEPTANCE RATE
*/
double update_get_accept_rate(mcmc_param *in){
    in->n_accept = in->n_accept_mu1 + in->n_accept_gamma + in->n_accept_Tinf + in->n_accept_alpha + in->n_accept_kappa;
    in->n_reject = in->n_reject_mu1 + in->n_reject_gamma + in->n_reject_Tinf + in->n_reject_alpha + in->n_reject_kappa;
    return (double) in->n_accept / (double) (in->n_accept+in->n_reject);
}










/*
   ================
   TUNING FUNCTIONS
   ================
*/

void tune_mu1(mcmc_param * in, gsl_rng *rng){
    /* get acceptance proportion */
    double paccept = (double) in->n_accept_mu1 / (double) (in->n_accept_mu1 + in->n_reject_mu1);

    /* acceptable zone: 35-45% acceptance */
    if(paccept<0.35) {
	in->sigma_mu1 /= 1.5;
    } else if (paccept>0.45) in->sigma_mu1 *= 1.5;
}




void tune_gamma(mcmc_param * in, gsl_rng *rng){
    /* get acceptance proportion */
    double paccept = (double) in->n_accept_gamma / (double) (in->n_accept_gamma + in->n_reject_gamma);

    /* acceptable zone: 35-45% acceptance */
    if(paccept<0.35) {
	in->sigma_gamma /= 1.5;
    } else if (paccept>0.45) in->sigma_gamma *= 1.5;
}







/*
   ===============================================
   METROPOLIS-HASTING ALGORITHM FOR ALL PARAMETERS
   ===============================================
*/
void mcmc(int nIter, int outEvery, char outputFile[256], char mcmcOutputFile[256], int tuneEvery, bool quiet, param *par, data *dat, dna_dist *dnainfo, gentime *gen, mcmc_param *mcmcPar, gsl_rng *rng){

    int i;

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
    fprintf(file, "step\tmu1\tmu2\tgamma");
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
    fprintf(mcmcFile, "step\tp_accept\tp_accept_mu1\tp_accept_gamma\tsigma_mu1\tsigma_gamma\tn_like_zero");


    /* OUTPUT TO SCREEN - HEADER */
    if(!quiet){
	printf("step\tmu1\tmu2\tgamma");
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

    fprint_param(file, par, 1, quiet);
    fprint_mcmc_param(mcmcFile, mcmcPar, 1);


    /* CREATE TEMPORARY PARAMETERS */
    param *tempPar = alloc_param(dat->n);
    copy_param(par,tempPar);


    /* RUN NITER CHAINS */
    for(i=2;i<=nIter;i++){
	if(i % outEvery == 0){
	    fprint_param(file, par, i, quiet);
	    fprint_mcmc_param(mcmcFile, mcmcPar, i);
	}

	/* TUNING */
	if(i % tuneEvery == 0){
	    tune_mu1(mcmcPar,rng);
	    tune_gamma(mcmcPar,rng);
	}

	/* MOVEMENTS */
	/* move mu1 */
	move_mu1(par, tempPar, dat, dnainfo, mcmcPar, rng);

	/* move gamma */
	move_gamma(par, tempPar, dat, dnainfo, mcmcPar, rng);

	/* move Tinf */
	printf("\nTinf:");
	print_vec_int(par->Tinf);
	fflush(stdout);
	move_Tinf(par, tempPar, dat, dnainfo, gen, mcmcPar, rng);

	/* move alpha_i*/
	move_alpha(par, tempPar, dat, dnainfo, gen, mcmcPar, rng);

	/* move kappa_i*/
	/* move_kappa(par, tempPar, dat, dnainfo, gen, mcmcPar, rng); */

	/* /\* move alpha and kappa *\/ */
	/* move_alpha_kappa(par, tempPar, dat, dnainfo, gen, mcmcPar, rng); */
    }


    /* CLOSE OUTPUT OUTFILE */
    fclose(file);
    fclose(mcmcFile);

    /* FREE TEMPORARY PARAMETERS */
    free_param(tempPar);
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

