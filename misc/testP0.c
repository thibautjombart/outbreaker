/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), January 2012.
  Licence: GPL >=2.

*/


#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
/* Calls to GNU Scientific Library */
#include <gsl/gsl_rng.h> /* random nb generators */
#include <gsl/gsl_randist.h> /* rng with specific distributions */



int main(){
	int i;

	/* INIT GSL RNG*/
	
	time_t t = time(NULL); // time in seconds, used to change the seed of the random generator
	gsl_rng * rng;
 	const gsl_rng_type *typ;
	gsl_rng_env_setup();
	typ=gsl_rng_default;
	rng=gsl_rng_alloc(typ);
	gsl_rng_set(rng,t); // changes the seed of the random generator

	printf("\n==Test Poisson==");fflush(stdout);
	for(i=0;i<20;i++){
	    printf("\nP(x=%d|lambda=0):",i);fflush(stdout);
	    printf(" %f",gsl_ran_poisson_pdf(i, 0.0));fflush(stdout);
	}

	printf("\n\n==Test exponential==");fflush(stdout);
	for(i=0;i<20;i++){
	    printf("\nP(x=%d|lambda=0):",i);fflush(stdout);
	    printf(" %f",gsl_ran_exponential_pdf(i, 1.0));fflush(stdout);
	}

	printf("\n\n==Test lognormal==");fflush(stdout);
	for(i=0;i<20;i++){
	    printf("\nP(x=%d|lambda=0):",i);fflush(stdout);
	    printf(" %f",gsl_ran_lognormal_pdf(i, 1.0, 1.25));fflush(stdout);
	}

}


/* gcc line

   gcc -o testP0 testP0.c -lgsl -lgslcblas
 
*/
