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
	const int M = 1e9;
	time_t t1, t2, t;
	double mu= 1e-7;
	int L=3e6, T, toto;

	/* INIT GSL RNG*/
	t = time(NULL); // time in seconds, used to change the seed of the random generator
	gsl_rng * rng;
 	const gsl_rng_type *typ;
	gsl_rng_env_setup();
	typ=gsl_rng_default;
	rng=gsl_rng_alloc(typ);
	gsl_rng_set(rng,t); // changes the seed of the random generator


	printf("\n - task using successive if - ", M);
	time(&t1);
	for(i=0;i<M;i++){
		T=gsl_rng_uniform_int(rng,5);		
		if(T==0) toto=0;
		else if(T==1) toto=0;
		else if(T==2) toto=0;
		else if(T==3) toto=0;
		else if(T==4) toto=0;
		else toto=0;

	}
	time(&t2);
	printf("\nTime ellapsed: %d ", (int) (t2-t1));


	printf("\n - using switch - ", M);
	time(&t1);
	for(i=0;i<M;i++){
		T=gsl_rng_uniform_int(rng,5);		
		switch(T){
		case 0:
			toto=0;
			break;
		case 1:
			toto=0;
			break;	
		case 2:
			toto=0;
			break;	
		case 3:
			toto=0;
			break;	
		case 4:
			toto=0;
			break;	
		default:
			toto=0;
		}
	}
	time(&t2);
	printf("\nTime ellapsed: %d \n", (int) (t2-t1));
		
}


/* gcc line 
   
   gcc -o switch switch.c -lgsl -lgslcblas && switch

   gcc -o switch switch.c -lgsl -lgslcblas -O2 && switch
   
*/
