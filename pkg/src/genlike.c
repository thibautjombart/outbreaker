/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), January 2012.
  Licence: GPL >=2.

*/

#include "common.h"
#include "matvec.h"
#include "genclasses.h"
#include "distances.h"
#include "param.h"
#include "genlike.h"







/*
   =================
   === AUXILIARY ===
   =================
*/



/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

double compute_genlike(int i, int j, double ti, double tj, double nu1, double nu2, double alpha, double tau, struct dna_dist *dnainfo, struct param *par){
	double out, T = ti>tj ? ti-tj : tj-ti;
	int temp;

	/* SORT I/J TO OPTIMIZE CACHE USAGE */
	if(i>j){
		temp=i;
		i=j;
		j=temp;
	}

	/* GET PSEUDO-LIKELIHOOD */
	if(dnainfo->nbcommon->rows[i]->values[j] < 1) { /* if no genetic data */
		out = par->weightNaGen;
	} else { /* if genetic info availavable */
		out = alpha * ( gsl_ran_poisson_pdf((unsigned int) dnainfo->transi->rows[i]->values[j], nu1*T*dnainfo->nbcommon->rows[i]->values[j]) 
				+ gsl_ran_poisson_pdf((unsigned int) dnainfo->transv->rows[i]->values[j], nu2*T*dnainfo->nbcommon->rows[i]->values[j]) ) 
			+ (1.0-alpha) * ( gsl_ran_poisson_pdf((unsigned int) dnainfo->transi->rows[i]->values[j], nu1*(T+2.0*tau)*dnainfo->nbcommon->rows[i]->values[j]) 
				       + gsl_ran_poisson_pdf((unsigned int) dnainfo->transv->rows[i]->values[j], nu2*(T+2.0*tau)*dnainfo->nbcommon->rows[i]->values[j]) ) ;
	}

	/* RETURN */
	return out;
} /* end compute_genlike */




/*
   =========================
   === TESTING FUNCTIONS ===
   =========================
*/


int main(){
	const int N=5, L=10;
	int i,j;
	double alpha=0.5, tau=2.0;

	/* create a list of sequences */
	struct list_dnaseq * test = create_list_dnaseq(N, L);

	for(i=0;i<N;i++){
		for(j=0;j<L;j++){
			if(i*j % 5 ==0) test->list[i]->seq[j] = 'a';
			else if(i*j % 3 ==0)test->list[i]->seq[j] = 't';
			else if(i*j % 2 ==0)test->list[i]->seq[j] = 'g';
			else test->list[i]->seq[j] = 'c';
		}
	}

	for(i=5;i<L;i++)
		test->list[0]->seq[i] = '-';
	for(i=0;i<5;i++)
		test->list[N-1]->seq[i] = '-';

	print_list_dnaseq(test);

	/* COMPUTE DISTANCES */
	struct dna_dist *distinfo = compute_dna_distances(test);
	print_dna_dist(distinfo);


	/* COMPUTE LIKELIHOOD */
	struct param *par = create_param();
	double out=0, nu1 = 0.01, nu2=0.02, t_vec[5]={0.0, 10.0, 12.0, 50.0, 100.0}, T;
	int count=0;

	par->weightNaGen = 0.001; /* near zero if no data */
	for(i=0;i<N-1;i++){
		for(j=i+1;j<N;j++){
			T = t_vec[i]>t_vec[j] ? t_vec[i] - t_vec[j] : t_vec[j]-t_vec[i];
			printf("\npair %d <-> %d:",i+1,j+1);
			printf("\nnb transi %d:", get_transi(distinfo,i,j));
			printf(" (esperance: %.1f)", T*nu1*get_nbcommon(distinfo,i,j));
			printf("\nnb transv %d:", get_transv(distinfo,i,j));
			printf(" (esperance: %.1f)", T*nu2*get_nbcommon(distinfo,i,j));
			out = compute_genlike(i, j, t_vec[i], t_vec[j], nu1, nu2, alpha, tau, distinfo, par);
			printf("\npseudo-likelihood %d <-> %d: %.5f\n",i+1,j+1, out);
		}
	}
	printf("\n");

	/* free and return */
	free_list_dnaseq(test);
	free_dna_dist(distinfo);
	free_param(par);

	return 0;
}



/*
  gcc instructions

  gcc -o genlike matvec.c genclasses.c distances.c param.c genlike.c -lgsl -lgslcblas && ./genlike

  valgrind --leak-check=full genlike

*/
