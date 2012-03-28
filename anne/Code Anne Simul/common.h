#ifndef __COMMON_H
#define __COMMON_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#define Tmax 100 // duration of the study
#define NbPatientsMax 2000
#define SizeWard0 8
#define SizeWard1 12

typedef struct
{
	int *NbAdmissions;
	int *NbPosSwabs;
	int *NbNegSwabs;
} nb_data;

typedef struct
{
	int *ward; /* 0 for adult, 1 for pediatric */
	gsl_vector * A[NbPatientsMax]; /* admission times. A[i] is a vector of size nb_data->NbAdmissions[i] */
	gsl_vector * D[NbPatientsMax]; /* discharge times. D[i] is a vector of size nb_data->NbAdmissions[i] */
	gsl_vector *P[NbPatientsMax]; /* times of positive swabs. P[i] is a vector of size nb_data->NbPosSwabs[i] */
	gsl_vector *N[NbPatientsMax]; /* times of positive swabs. N[i] is a vector of size nb_data->NbNegSwabs[i] */
	int *C; /* colonisation times */
	int *E; /* times of end of colonisation */
	int *I0; /* number of colonised individuals in ward 0 at each time step */
	int *I1; /* number of colonised individuals in ward 1 at each time step */
	int *timeSeq; /* time of sequence for each patient */
} raw_data;

typedef struct
{
    /* to be estimated */
	gsl_matrix *beta; /* person to person transmission rates between and within wards */
	double betaWardOut; /* force of transmission from outside applied to patients in the wards */
	double betaOutOut; /* force of transmission applied to individuals outside the wards */

	double Sp; /* specificity of the test */
	double Se; /* sensibility of the test */

    double Pi; /* probability of being colonized at first admission */

    double muHosp; /* mean duration of hospitalisation */
    double sigmaHosp; /* std of the duration of hospitalisation */

    double mu; /* mean duration of colonisation */
    double sigma; /* std of the duration of colonisation */
    /* *************** MAYBE NEED TO REPARAMETERIZE MU, SIGMA INTO MU, CV *************** */

    double nu1; /* rate of transitions */
    double nu2; /* rate of tranversions */
    /*************** MAYBE NEED TO REPARAMETERIZE NU1, NU2 INTO NU1, coeff de proportionalite *************/

    double tau; /* time to the most recent common ancestor */
    double alpha; /* probability that two sampled isolates belong to the same lineage */

} parameters;

#endif
