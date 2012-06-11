
/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk) and Anne Cori (a.cori@imperial.ac.uk), 2012.
  Licence: GPL >=2.
*/


/*
  This file contains the defition of the data structures used in this program.
*/


#ifndef __STRUCTURES_H
#define __STRUCTURES_H

#include "common.h"
#include "matvec.h"
#include "genclasses.h"

 typedef struct{
     int n, length; /* n: number of observations; length: sequence length */
     vec_int * dates; /* collection dates*/
     list_dnaseq * dna; /* sequences */
     gsl_rng * rng;
 } data;



/* typedef struct{ */
/*     int n; /\* number of infections *\/ */
/*     vec_int * T_inf; /\* infection dates *\/ */
/* } aug_data; */








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
