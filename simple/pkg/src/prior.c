#if 0

#include "common.h"
#include "init.h"
#include "InputOutput.h"
#include "logL.h"
#include "mcmc.h"
#include "moves.h"
#include "prior.h"
#include "alloc.h"
#include "tuneVariances.h"

/* extern gsl_rng * rng; */







/*****************************************************************************/
/*                        PRIORS                                             */
/*****************************************************************************/

/******************************************************************************/
/* Transmission rates                                                         */
/******************************************************************************/

/* flat exponential priors */
/* Note : here we exclude the case were any transmission rate is 0 */

double logpriorBeta(int i, int j, parameters * param){
    double h = log(gsl_ran_exponential_pdf (gsl_matrix_get(param->beta,i,j), 1000));
    return h;
}





double logpriorBetaWardOut(parameters * param){
    double h = log(gsl_ran_exponential_pdf(param->betaWardOut, 1000));
    return h;
}





double logpriorBetaOutOut(parameters * param){
    double h = log(gsl_ran_exponential_pdf(param->betaOutOut, 1000));
    return h;
}

/******************************************************************************/
/* Probability of being infected at first admission                           */
/******************************************************************************/

double logpriorPi(parameters * param){
    /* uniform prior */
    double h = 0;
    return h;
}

/******************************************************************************/
/* Specificity and sensitivity of testing                                     */
/******************************************************************************/

/* uniform priors */
double logpriorSp(parameters * param){

    double h = 0;
    return h;
}

double logpriorSe(parameters * param){

    double h = 0;
    return h;
}







/******************************************************************************/
/* Mean and std of duration of colonisation period                            */
/******************************************************************************/
double logpriorMu(parameters * param){

    double h = log(gsl_ran_exponential_pdf(param->mu, 1000));
    return h;
}





double logpriorSigma(parameters * param){

    double h = log(gsl_ran_exponential_pdf(param->sigma, 1000));
    return h;
}







/******************************************************************************/
/* Mutation rates                                                             */
/******************************************************************************/
double logpriorNu1(parameters * param){

    double h = log(gsl_ran_exponential_pdf(param->nu1, 1000));
    return h;
}





double logpriorKappa(parameters * param){

    double h = log(gsl_ran_exponential_pdf(param->kappa, 1000));
    return h;
}








/******************************************************************************/
/* Tau                                                                        */
/******************************************************************************/
/* double logpriorTau(parameters * param){ */

/*     double h = log(gsl_ran_exponential_pdf(param->tau, 1000)); */
/*     return h; */
/* } */







/******************************************************************************/
/* Alpha                                                                        */
/******************************************************************************/
/* uniform prior */
/* double logpriorAlpha(parameters * param){ */

/*     double h = 0; */
/*     return h; */
/* } */







/******************************************************************************/
/* Total                                                                      */
/******************************************************************************/
double logprior (parameters * param){
    int i,j;
    double h=0;

    for(i=0 ; i<2 ; i++)
	{
	    for(j=0 ; j<2 ; j++)
		{
		    h+=logpriorBeta(i,j,param);
		}
	}

    h+=logpriorBetaWardOut(param)+logpriorBetaOutOut(param);

    /* h+=logpriorSp(param)+logpriorSe(param); */
    h+=logpriorSe(param);

    h+=logpriorMu(param)+logpriorSigma(param);

    h+=logpriorNu1(param)+logpriorKappa(param);

    /* h+=logpriorTau(param)+logpriorAlpha(param); */

    return h;
}


#endif
