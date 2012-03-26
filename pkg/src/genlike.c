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
/* 
   computes the likelihood of a sequences "S_i" given that "j infected i" and "S_j", with:
   - s_i: indices of sequences in i (target patient)
   - s_j: indices of sequences in j (infector)
   - t_i, t_j: collection times for sets of sequences S_i and S_j
   - m_i, m_j: number of sequences for patients i and j; length of vectors s_i,s_j,t_i,t_j
   - dnainfo: the pre-computed pairwise comparisons of DNA sequences
   
*/
double genlike_ij(int *s_i, int *s_j, double *t_i, double *t_j, int m_i, int m_j, double nu1, double nu2, double alpha, double tau, struct dna_dist *dnainfo, struct param *par){
	/* double out, T = ti>tj ? ti-tj : 0.0, Tabs = ti>tj ? ti-tj : tj-ti; */
	double out, T, Xi1, Xi2, Xi3, Xi4, Pk;
	int k, q, r, transi, tranv, common, nb_comp;


	/* Compute Pk for each sequence 'k' in S_i */
	/* Pk = p(s_i^k| s_i^1, ..., s_i^{k-1}, S_j, i<-j) */

	nb_comp = 0; /* count the number of effective sequence comparisons used in Pk */

	if((m_i > 0 && m_j > 0) || m_i>1){ /* likelihood tractable if at least a pair is available */
		Pk = 0.0;
		for(k=0;k<m_i;k++){
			/* for a given 'k' */
			/* ancestries from S_j */
			Xi1=0.0;
			Xi3=0.0;
			for(q=0;q<m_j;q++){
				/* stuff used to compute Poisson mass function */
				T = t_i[k]>t_j[q] ? t_i[k]-t_j[q] : t_j[q]-t_i[k]; /* time difference (absolute value) */
				tansi = matint_ij(dnainfo->transi,s_i[k],s_j[q]); /* transitions */
				tansv = matint_ij(dnainfo->transv,s_i[k],s_j[q]); /* transversions */
				common = matint_ij(dnainfo->nbcommon,s_i[k],s_j[q]); /* nb of nucleotides compared */

				if(common > 0){
					nb_comp++;

					/* direct ancestries */
					if(t_i[k] > t_j[q]){
						/* transitions */
						Xi1 += gsl_ran_poisson_pdf((unsigned int) matint_ij(dnainfo->transi,s_i[k],s_j[q]), nu1*T*common);
						/* transversions */
						Xi1 += gsl_ran_poisson_pdf((unsigned int) matint_ij(dnainfo->transv,s_i[k],s_j[q]), nu2*T*common);
					}

					/* indirect ancestries */
					/* transitions */
					Xi3 += gsl_ran_poisson_pdf((unsigned int) matint_ij(dnainfo->transi,s_i[k],s_j[q]), nu1*(T+2.0*tau)*common);
					/* transversions */
					Xi3 += gsl_ran_poisson_pdf((unsigned int) matint_ij(dnainfo->transv,s_i[k],s_j[q]), nu2*(T+2.0*tau)*common);
				}
			}

			/* ancestries from S_i */
			Xi2=0.0;
			Xi4=0.0;
			for(r=0;r<(k-1);r++){
				/* stuff used to compute Poisson mass function */
				T = t_i[k]>t_i[r] ? t_i[k]-t_i[r] : t_i[r]-t_i[k]; /* time difference (absolute value) */
				tansi = matint_ij(dnainfo->transi,s_i[k],s_i[r]); /* transitions */
				tansv = matint_ij(dnainfo->transv,s_i[k],s_i[r]); /* transversions */
				common = matint_ij(dnainfo->nbcommon,s_i[k],s_i[r]); /* nb of nucleotides compared */

				if(common > 0){
					nb_comp++;

					/* direct ancestries */
					if(t_i[k] > t_i[r]){
						/* transitions */
						Xi2 += gsl_ran_poisson_pdf((unsigned int) matint_ij(dnainfo->transi,s_i[k],s_i[r]), nu1*T*common);
						/* transversions */
						Xi2 += gsl_ran_poisson_pdf((unsigned int) matint_ij(dnainfo->transv,s_i[k],s_i[r]), nu2*T*common);
					}

					/* indirect ancestries */
					/* transitions */
					Xi4 += gsl_ran_poisson_pdf((unsigned int) matint_ij(dnainfo->transi,s_i[k],s_i[r]), nu1*(T+2.0*tau)*common);
					/* transversions */
					Xi4 += gsl_ran_poisson_pdf((unsigned int) matint_ij(dnainfo->transv,s_i[k],s_i[r]), nu2*(T+2.0*tau)*common);
				}
			}

			/* likelihood (p_s_i^k | s_i^1, ..., s_i^{k-1}, S_j, i<-j)*/
			Pk = alpha * (Xi1 + Xi2) + (1.0-alpha) * (Xi3 + Xi4);

			/* update general likelihood */
			out += log(Pk);
		}
	}

	/* Two cases when likelihood can't be computed:
	   - no pair of sequences to be compared
	   - no compared sequence had nucleotide in common */
	if(nb_comp==0){
		out = log(par->weightNaGen);
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
