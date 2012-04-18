/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), January 2012.
  Licence: GPL >=2.

*/

#include "common.h"
#include "genclasses.h"
#include "simgen.h"




/*
  ===================
  AUXILIARY FUNCTIONS
  ===================
*/

char transi(char in){
    switch(in){
    case 'a':
	return 'g';
    case 'g':
	return 'a';
    case 'c':
	return 't';
    case 't':
	return 'c';
    default:
	fprintf(stderr, "\n[in: simgen.c->transi]\nUnknown character %c.\n", in);
	exit(1);
    }
}




char transv(char in, gsl_rng *rng){
    double x = gsl_ran_flat(rng, 0.0, 1.0);
    switch(in){
    case 'a':
	return x<0.5 ? 't': 'c';
    case 'g':
	return x<0.5 ? 't': 'c';
    case 'c':
    	return x<0.5 ? 'a': 'g';
    case 't':
    	return x<0.5 ? 'a': 'g';
    default:
	fprintf(stderr, "\n[in: simgen.c->transv]\nUnknown character %c.\n", in);
	exit(1);
    }
}





dnaseq *create_haplo(int length, gsl_rng *rng){
    int i;
    double x;
    dnaseq *out = create_dnaseq(length);

    for(i=0;i<length;i++){
	x = gsl_ran_flat(rng, 0.0, 4.0);
	if(x<1.0) {
	    out->seq[i] = 'a';
	} else if(x<2.0){
	    out->seq[i] = 'g';
	} else if(x<3.0){
	    out->seq[i] = 't';
	} else {
	    out->seq[i] = 'c';
	}
    }

    return out;
}







/*
  =======
  TESTING
  =======
*/


int main(){
    time_t t = time(NULL); /* time in seconds, used to change the seed of the random generator */
    const gsl_rng_type *typ;
    gsl_rng_env_setup();
    typ=gsl_rng_default;
    gsl_rng * rng=gsl_rng_alloc(typ);
    gsl_rng_set(rng,t); /* changes the seed of the random generator */

    int i;

    dnaseq *seq;

    /* transitions */
    printf("\n== Transitions ==");
    printf("\na:%c", transi('a'));
    printf("\ng:%c", transi('g'));
    printf("\nt:%c", transi('t'));
    printf("\nc:%c", transi('c'));

    printf("\n\n== Transversions ==");
    for(i=0;i<5;i++) printf("\na:%c", transv('a',rng));
    printf("\n");
    for(i=0;i<5;i++) printf("\ng:%c", transv('g',rng));
    printf("\n");
    for(i=0;i<5;i++) printf("\nt:%c", transv('t',rng));
    printf("\n");
     for(i=0;i<5;i++) printf("\nc:%c", transv('c',rng));
    printf("\n");


    printf("\n== Haplotype creation ==\n");
    seq = create_haplo(30, rng);
    print_dnaseq(seq);


    free_dnaseq(seq);
    gsl_rng_free(rng);
    return 0;
}


/* gcc instructions:

gcc -o simgen genclasses.c simgen.c -lgsl -lgslcblas && ./simgen


valgrind --leak-check=full simgen


*/



