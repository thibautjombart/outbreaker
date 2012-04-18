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




/* Make an existing haplotype evolve using :
   - nu1: rate of transitions 
   - nu2: rate of transversions 
   - double t: time of evolution between in and out
*/
void evolve_haplo(dnaseq *in, double nu1, double nu2, double t, gsl_rng *rng){
    int i, nbtransi, nbtransver, posi=0;

    /* handle transitions */
    nbtransi = gsl_ran_poisson(rng, nu1 * t * (double) in->length);
    for(i=0;i<nbtransi;i++){
	posi = gsl_rng_uniform_int(rng, in->length);
	in->seq[posi] = transi(in->seq[posi]);
    }

    /* handle transversions */
    nbtransver = gsl_ran_poisson(rng, nu2 * t * (double) in->length);
    for(i=0;i<nbtransver;i++){
	posi = gsl_rng_uniform_int(rng, in->length);
	in->seq[posi] = transv(in->seq[posi], rng);
    }
} /* end evolve_haplo */






/* Replicate an haplotype, creating a new sequence, using:
   - nu1: rate of transitions 
   - nu2: rate of transversions 
   - double t: time of evolution between in and out
*/
dnaseq *replicate_haplo(dnaseq *in, double nu1, double nu2, double t, gsl_rng *rng){
    dnaseq *out = create_dnaseq(in->length);

    /* copy haplotype */
    copy_dnaseq(in, out);

    /* evolve haplotype */
    evolve_haplo(out, nu1, nu2, t, rng);

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

    dnaseq *seq1, *seq2, *seq3;

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
    seq1 = create_haplo(30, rng);
    print_dnaseq(seq1);
    
    printf("\n== Haplotype copy ==\n");
    seq2 = create_dnaseq(30);
    copy_dnaseq(seq1, seq2);
    print_dnaseq(seq2);

    printf("\n== Haplotype replication ==\n");
    printf("\nref:");
    seq3 = replicate_haplo(seq1, 0.1, 0.2, 1.0, rng);
    printf("\nref:"); 
    print_dnaseq(seq2);
    printf("\nnew:"); 
    print_dnaseq(seq3);

    printf("\n== Evolution across several time steps ==\n");
    copy_dnaseq(seq1, seq3);
    for(i=0;i<20;i++){
	evolve_haplo(seq3, 0.05, 0.1, 1.0, rng);
	print_dnaseq(seq3);
    }

    free_dnaseq(seq1);
    free_dnaseq(seq2);
    free_dnaseq(seq3);
    gsl_rng_free(rng);
    return 0;
}


/* gcc instructions:

gcc -o simgen genclasses.c simgen.c -lgsl -lgslcblas && ./simgen


valgrind --leak-check=full simgen


*/



