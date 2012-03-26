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
	const int M = 5e7;
	time_t t1, t2, t;
	double mu= 1e-7;
	int L=3e6, T, N;

	/* INIT GSL RNG*/
	t = time(NULL); // time in seconds, used to change the seed of the random generator
	gsl_rng * rng;
 	const gsl_rng_type *typ;
	gsl_rng_env_setup();
	typ=gsl_rng_default;
	rng=gsl_rng_alloc(typ);
	gsl_rng_set(rng,t); // changes the seed of the random generator


	printf("\n - performing %d computations of Bernoulli prob mass fct - ", M);
	time(&t1);
	for(i=0;i<M;i++){
		T=gsl_rng_uniform_int(rng,100); /* time */
		N=gsl_rng_uniform_int(rng, 20); /* nb of mutations*/
		gsl_ran_binomial_pdf((unsigned int) N, mu, T*L);
	}
	time(&t2);
	printf("\nTime ellapsed: %d ", (int) (t2-t1));


	printf("\n - performing %d computations of Poisson prob mass fct- ", M);
	time(&t1);
	for(i=0;i<M;i++){
		T=gsl_rng_uniform_int(rng,100); /* time */
		N=gsl_rng_uniform_int(rng, 20); /* nb of mutations*/
		gsl_ran_poisson_pdf((unsigned int) N, mu*T*L);
	}
	time(&t2);
	printf("\nTime ellapsed: %d \n", (int) (t2-t1));

}


/* gcc line 

   gcc -o bernouliVsPoisson bernouliVsPoisson.c -lgsl -lgslcblas
 
*/
