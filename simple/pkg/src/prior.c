
#include "common.h"
#include "prior.h"


/* p(alpha_i): = 1/(n-1) if alpha_i \neq i, and zero otherwise */
double logprior_alpha_i(int i, param *par){
    double out = (i==vec_int_i(par->alpha,i)) ? 0.0 : (double) 1/(par->n-1);
    return log(out);
}



/* p(kappa_i) = NB(1, 1-pi) with pi: prop obs cases */
double logprior_kappa_i(int i, param *par){
    double out = gsl_ran_negativeBINOMIALshit_pdf(vec_int_i(par->kappa,i)-1, 1, par->pi); /* TO BE REPLACED WITH ADEQUATE FUNCTION */
    return log(out);
}




/* p(mu_1) = Unif(0,1) = 1*/
double logprior_mu_1(){
    return 0.0; /* log(1) -> 0 */
}



/* p(gamma) = logN(1,1.25) */
double logprior_gamma(param *par){
    double out = gsl_ran_somethingforLogNormal(param->gamma, 1, 1.25); /* TO BE REPLACED WITH ADEQUATE FUNCTION CALL */
    return log(out);
}





double logprior_all(param *par){
    int i;
    double out=0.0;

    /* result is the sum of priors over all parameters and individuals */
    for(i=0;i<par->n;i++){
	out += logprior_alpha_i(i,par);
	out += logprior_kappa_i(i,par);
    }

    out += logprior_mu1();
    out += logprior_gamma(par);

    return(out);
}





/* /\******************************************************************************\/ */
/* /\* Total                                                                      *\/ */
/* /\******************************************************************************\/ */
/* double logprior (parameters * param){ */
/*     int i,j; */
/*     double h=0; */

/*     for(i=0 ; i<2 ; i++) */
/* 	{ */
/* 	    for(j=0 ; j<2 ; j++) */
/* 		{ */
/* 		    h+=logpriorBeta(i,j,param); */
/* 		} */
/* 	} */

/*     h+=logpriorBetaWardOut(param)+logpriorBetaOutOut(param); */

/*     /\* h+=logpriorSp(param)+logpriorSe(param); *\/ */
/*     h+=logpriorSe(param); */

/*     h+=logpriorMu(param)+logpriorSigma(param); */

/*     h+=logpriorNu1(param)+logpriorKappa(param); */

/*     /\* h+=logpriorTau(param)+logpriorAlpha(param); *\/ */

/*     return h; */
/* } */

