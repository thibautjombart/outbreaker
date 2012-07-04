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
	/* INIT GSL RNG*/
	time_t t = time(NULL); // time in seconds, used to change the seed of the random generator
	gsl_rng * rng;
 	const gsl_rng_type *typ;
	gsl_rng_env_setup();
	typ=gsl_rng_default;
	rng=gsl_rng_alloc(typ);
	gsl_rng_set(rng,t); // changes the seed of the random generator

	/* OPEN OUTPUT FILE */
	FILE *file = fopen("test.csv","w");
	fprintf(file, "lambda,p(1|lambda)");
	printf("\nmet1\tmet2");

	/* VARIABLES */
	int i, N=10000;
	double mu = 8, temp1, temp2, sigma=0.5;
	

	for(i=0;i<N;i++){
	    /* values centred on mu - method 1 */
	    temp1 = gsl_ran_lognormal(rng,log(mu),sigma);
	    
	    /* Poisson density there */
	    temp2 =  gsl_ran_poisson_pdf((unsigned int) 1, temp1);
	    
	    /* print to screen/file */
	    printf("\n%.15f\t%.15f", temp1,temp2);
	    fprintf(file, "\n%.15f,%.15f", temp1,temp2);
	
	}

	printf("\nPoisson(0|0)=%.15f\n",gsl_ran_poisson_pdf((unsigned int) 0, 0));
	
	printf("\n\n");fflush(stdout);
	fclose(file);
}


/* gcc line 

   gcc -o testRandist testRandist.c -lgsl -lgslcblas
 
*/
