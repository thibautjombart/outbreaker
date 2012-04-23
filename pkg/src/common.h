/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk) and Anne Cori (a.cori@imperial.ac.uk), 2012.
  Licence: GPL >=2.

*/

#ifndef __COMMON_H
#define __COMMON_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

/* Calls to GNU Scientific Library */
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>

#define NEARZERO 0.0000000001
#define TRUE 1
#define FALSE 0

#define SizeWard0 8
#define SizeWard1 12
/* #define NbPatients 142 */
/* #define T 100 /\* number of time steps corresponding to the study period *\/ */

typedef short int bool;



typedef struct{
    int NbPatients, T, NbSequences; /* T: number of time steps corresponding to the study period */
    int * NbAdmissions; /* vector of size NbPatients */
    int * NbPosSwabs; /* vector of size NbPatients */
    int * NbNegSwabs; /* vector of size NbPatients */
    int NbColonisedPatients; /* those who have at least one positive swab */
    int * indexColonisedPatients; /* those who have at least one positive swab */
    int * M; /* number of sequences in each patient */
} nb_data;

typedef struct{
    /* EPI DATA */
    int NbPatients, T, NbSequences; /* T: number of time steps corresponding to the study period */
    int * ward; /* vector of length NbPatients, 0 for adult, 1 for pediatric */
    /* int * timeSeq; /\* time of sequence for each patient *\/ */
    /* int * PatientIndex; /\* vector from 0 to (NbPatients-1) *\/ */
    gsl_vector ** A; /* admission times. A[i] is a vector of size nb_data->NbAdmissions[i] */
    gsl_vector ** D; /* discharge times. D[i] is a vector of size nb_data->NbAdmissions[i] */
    gsl_vector ** P; /* times of positive swabs. P[i] is a vector of size nb_data->NbPosSwabs[i] */
    gsl_vector ** N; /* times of negative swabs. N[i] is a vector of size nb_data->NbNegSwabs[i] */
    /* gsl_vector *IsInHosp[NbPatients]; /\* 0 when patient is outside hospital, 1 when inside. IsInHosp[i] is a vector of size T *\/ */
    gsl_vector ** IsInHosp; /* 0 when patient is outside hospital, 1 when inside. IsInHosp[i] is a vector of size T */

    /* GENETIC DATA */
    int ** S; /* indices of sequences collected for each patient; S[i] has length M[i] */
    int * Tcollec; /* collection times of the sequences; length: NbSequences*/
    int * M; /* number of sequences collected for each patient */

    /* RANDOM NUMBER GENERATOR */
    gsl_rng * rng; /* nice here since raw data don't change */
} raw_data;




typedef struct{
    int NbPatients, T;
    int *C; /* colonisation times */
    int *E; /* times of end of colonisation */
    int *I0; /* number of colonised individuals in ward 0 at each time step */
    int *I1; /* number of colonised individuals in ward 1 at each time step */
} aug_data;



typedef struct{
    /* to be estimated */
    gsl_matrix *beta; /* person to person transmission rates between and within wards */
    double betaWardOut; /* force of transmission from outside applied to patients in the wards */
    double betaOutOut; /* force of transmission applied to individuals outside the wards */

    /* double Sp; /\* specificity of the test *\/ assumed = 100% */
    double Se; /* sensibility of the test */

    double Pi; /* probability of being colonized at first admission */

    double mu; /* mean duration of colonisation */
    double sigma; /* std of the duration of colonisation */
    /* *************** MAYBE NEED TO REPARAMETERIZE MU, SIGMA INTO MU, CV *************** */

    double nu1; /* rate of transitions */
    double nu2; /* rate of tranversions */
    /*************** MAYBE NEED TO REPARAMETERIZE NU1, NU2 INTO NU1, coeff de proportionalite *************/

    double tau; /* time to the most recent common ancestor */
    double alpha; /* probability that two sampled isolates belong to the same lineage */
    double weightNaGen; /* weight used to replace missing genetic likelihood */
} parameters;


typedef struct{
    double mu;  //mean duration of hospitalisation
    double sigma;  //std of the duration of hospitalisation
} hospDurationParam;





typedef struct{
    long NbSimul; /* nb of iterations of Metropolis-Hastings */
    int SubSample; /* results recorded every SubSample iterations */
    int BurnIn; /* nb of iterations considered as the Burn in period */

    /* standard deviation for the proposition laws */
    gsl_matrix *Sigma_beta;
    double Sigma_betaWardOut;
    double Sigma_betaOutOut;
    double Sigma_mu;
    double Sigma_sigma;
    double Sigma_nu1;
    double Sigma_nu2;
    double Sigma_tau;
    double Sigma_alpha;
} mcmcInternals;





typedef struct{
    /* average probabilities of acceptance */
    gsl_matrix *PourcAcc_beta;
    double PourcAcc_betaWardOut;
    double PourcAcc_betaOutOut;
    double PourcAcc_mu;
    double PourcAcc_sigma;
    double PourcAcc_nu1;
    double PourcAcc_nu2;
    double PourcAcc_tau;
    double PourcAcc_alpha;

} acceptance;





typedef struct{
    gsl_matrix *IsAccOK_beta;
    double IsAccOK_betaWardOut;
    double IsAccOK_betaOutOut;
    double IsAccOK_mu;
    double IsAccOK_sigma;
    double IsAccOK_nu1;
    double IsAccOK_nu2;
    double IsAccOK_tau;
    double IsAccOK_alpha;
} isAcceptOK;





typedef struct{
    gsl_matrix *NbProp_beta;
    double NbProp_betaWardOut;
    double NbProp_betaOutOut;
    double NbProp_mu;
    double NbProp_sigma;
    double NbProp_nu1;
    double NbProp_nu2;
    double NbProp_tau;
    double NbProp_alpha;
} NbProposals;





typedef struct{
    FILE *LogL;
    FILE *Parameters;
    FILE *ColonDates;
    FILE *EndColonDates;
} output_files;



#endif
