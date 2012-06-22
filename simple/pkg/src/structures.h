
/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk) and Anne Cori (a.cori@imperial.ac.uk), 2012.
  Licence: GPL >=2.
*/


/*
  This file contains the defition of the data structures used in this program, ad elementary methods (allocation, de-allocation, printing, accessing/setting content).
*/


#ifndef __STRUCTURES_H
#define __STRUCTURES_H

#include "common.h"
#include "matvec.h"
#include "genclasses.h"


/*
  =========
   CLASSES
  =========
*/

/* data: contains only observed data */
 typedef struct{
     int n, length; /* n: number of observations; length: sequence length */
     vec_int * dates; /* collection dates*/
     list_dnaseq * dna; /* sequences */
 } data;




/* descriptors/parameters of the generation time function 'w' */
typedef struct{
    int type; /* type of distribution for generation time */
    double param1, param2, param3; /* parameters for the generation time distribution 'w' */
    int trunc; /* value of truncation; p(x>=truc)=0 */
    int maxK; /* maximum value of kappa_i (i.e. max nb of generations between two cases) */
    mat_double *dens; /* pre-computed values of density: row 'i' gives densities for kappa=i at t=0, 1, ..., trunc-1*/
} gentime;




/* param: contains augmented data and parameters */
typedef struct{
    int n; /* number of observations, length of the vec_int objects */
    vec_int *Tinf; /* times of infection */
    vec_int *alpha; /* idx of closest ancestor for each case; -1 = unknown */
    vec_int *kappa; /* number of generations before a case and its closest ancestor */
    double mu1; /* rate of transitions */
    double gamma; /* so that rate of transversions mu2 = gamma x mu1 */
    double pi; /* proportion of observed cases */
} param;




/* mcmc_param: contains parameters of mcmc */
typedef struct{
    int n_accept, n_reject; /* global accept/reject*/
    int n_accept_mu1, n_reject_mu1; /* accept/reject for mu1 */
    int n_accept_gamma, n_reject_gamma; /* accept/reject for gamma */
    double sigma_mu1; /* sigma for normal proposal for mu1 */
    double sigma_gamma; /* sigma for normal proposal for gamma */
} mcmc_param;









/*
  =========
   METHODS
  =========
*/

/*
  ======
   DATA
  ======
*/
data *alloc_data(int n, int length);

void free_data(data *in);

void print_data(data *in);

data * Rinput2data(unsigned char * DNAbinInput, int *Tcollec, int *n, int *length);





/*
  =======
  GENTIME
  =======
*/

gentime *alloc_gentime(int maxK, int trunc);

void free_gentime(gentime *in);

void print_gentime(gentime *in);

double gentime_dens(gentime *in, int t, int kappa_i);


/*
 =======
  PARAM
 =======
*/

param *alloc_param(int n);

void free_param(param *in);

void print_param(param *in);

void copy_param(param *in, param *out);




/*
 ============
  MCMC_PARAM
 ============
*/

mcmc_param *alloc_mcmc_param(int n);

void free_mcmc_param(mcmc_param *in);

void print_mcmc_param(mcmc_param *in);









/* typedef struct{ */
/*     /\* to be estimated *\/ */
/*     gsl_matrix *beta; /\* person to person transmission rates between and within wards *\/ */
/*     double betaWardOut; /\* force of transmission from outside applied to patients in the wards *\/ */
/*     double betaOutOut; /\* force of transmission applied to individuals outside the wards *\/ */

/*     /\* double Sp; /\\* specificity of the test *\\/ assumed = 100% *\/ */
/*     double Se; /\* sensibility of the test *\/ */

/*     double Pi; /\* probability of being colonized at first admission *\/ */

/*     double mu; /\* mean duration of colonisation *\/ */
/*     double sigma; /\* std of the duration of colonisation *\/ */
/*     /\* *************** MAYBE NEED TO REPARAMETERIZE MU, SIGMA INTO MU, CV *************** *\/ */

/*     double nu1; /\* rate of transitions *\/ */
/*     /\* double nu2; /\\* rate of tranversions *\\/ *\/ */
/*     double kappa; /\* nu2 = kappa*nu1 *\/ */

/*     /\* double tau; /\\* time to the most recent common ancestor *\\/ *\/ */
/*     /\* double alpha; /\\* probability that two sampled isolates belong to the same lineage *\\/ *\/ */
/*     double weightNaGen; /\* weight used to replace missing genetic likelihood *\/ */
/* } parameters; */


/* typedef struct{ */
/*     double mu;  //mean duration of hospitalisation */
/*     double sigma;  //std of the duration of hospitalisation */
/* } hospDurationParam; */





/* typedef struct{ */
/*     long NbSimul; /\* nb of iterations of Metropolis-Hastings *\/ */
/*     int SubSample; /\* results recorded every SubSample iterations *\/ */
/*     int BurnIn; /\* nb of iterations considered as the Burn in period *\/ */

/*     /\* standard deviation for the proposition laws *\/ */
/*     gsl_matrix *Sigma_beta; */
/*     double Sigma_betaWardOut; */
/*     double Sigma_betaOutOut; */
/*     double Sigma_mu; */
/*     double Sigma_sigma; */
/*     double Sigma_nu1; */
/*     double Sigma_kappa; */
/*     double Sigma_tau; */
/*     double Sigma_alpha; */
/* } mcmcInternals; */





/* typedef struct{ */
/*     /\* average probabilities of acceptance *\/ */
/*     gsl_matrix *PourcAcc_beta; */
/*     double PourcAcc_betaWardOut; */
/*     double PourcAcc_betaOutOut; */
/*     double PourcAcc_mu; */
/*     double PourcAcc_sigma; */
/*     double PourcAcc_nu1; */
/*     double PourcAcc_kappa; */
/*     double PourcAcc_tau; */
/*     double PourcAcc_alpha; */

/* } acceptance; */





/* typedef struct{ */
/*     gsl_matrix *IsAccOK_beta; */
/*     double IsAccOK_betaWardOut; */
/*     double IsAccOK_betaOutOut; */
/*     double IsAccOK_mu; */
/*     double IsAccOK_sigma; */
/*     double IsAccOK_nu1; */
/*     double IsAccOK_kappa; */
/*     double IsAccOK_tau; */
/*     double IsAccOK_alpha; */
/* } isAcceptOK; */





/* typedef struct{ */
/*     gsl_matrix *NbProp_beta; */
/*     double NbProp_betaWardOut; */
/*     double NbProp_betaOutOut; */
/*     double NbProp_mu; */
/*     double NbProp_sigma; */
/*     double NbProp_nu1; */
/*     double NbProp_kappa; */
/*     double NbProp_tau; */
/*     double NbProp_alpha; */
/* } NbProposals; */





/* typedef struct{ */
/*     FILE *LogL; */
/*     FILE *Parameters; */
/*     FILE *ColonDates; */
/*     FILE *EndColonDates; */
/* } output_files; */



#endif
