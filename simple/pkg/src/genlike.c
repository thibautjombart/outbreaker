
/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), January 2012.
  Licence: GPL >=2.

*/
#if 0

#include "common.h"
#include "alloc.h"
#include "matvec.h"
#include "genclasses.h"
#include "distances.h"
#include "genlike.h"





/*
  ===============================
  === MAIN EXTERNAL FUNCTIONS ===
  ===============================
*/
/*
   computes the likelihood of a sequences "S_i" given that "j infected i" and "S_j", with:
   == in 'data' ==
   - s_i: indices of sequences in i (target patient)
   - s_j: indices of sequences in j (infector)
   - t_i, t_j: collection times for sets of sequences S_i and S_j
   - m_i, m_j: number of sequences for patients i and j; length of vectors s_i,s_j,t_i,t_j

   == in 'dnainfo' ==
   - the pre-computed pairwise comparisons of DNA sequences

   == in 'param' ==
   - the various parameters needed
*/
double genlike_ij(int i, int j, raw_data *data, aug_data *augData, dna_dist *dnainfo, parameters *param){

    /* extract variables from input objects */
    int *s_i=data->S[i], *s_j=data->S[j], m_i=data->M[i], m_j=data->M[j];
    double t_i[m_i], t_j[m_j];
    double nu1=param->nu1, nu2=param->nu1*param->kappa, alpha, tau;

    /* get alpha and tau */
    alpha = gsl_matrix_get(augData->alpha,i,j);
    tau = gsl_matrix_get(augData->tau,i,j);

    /* variables used in computations */
    double out, Tabs, Xi1, Xi2, Xi3, Xi4, Pk;
    int k, q, r, transi, transv, common, nb_comp, nb_comp_k;

    /* fill in vectors of collection dates */
    for(k=0;k<m_i;k++){
	t_i[k] = data->Tcollec[s_i[k]];
    }
    for(k=0;k<m_j;k++){
	t_j[k] = data->Tcollec[s_j[k]];
    }


    /* Compute Pk for each sequence 'k' in S_i */
    /* Pk = p(s_i^k| s_i^1, ..., s_i^{k-1}, S_j, i<-j) */

    nb_comp = 0; /* count the number of effective sequence comparisons used overall */
    out=0.0; /* important initialization here */

    if((m_i > 0 && m_j > 0) || m_i>1){ /* likelihood tractable if at least a pair is available */
	for(k=0;k<m_i;k++){
	    /* initialize k-specific variables */
	    Pk = 0.0;
	    nb_comp_k = 0;

	    /* for a given 'k' */
	    /* ancestries from S_j */
	    Xi1=0.0;
	    Xi3=0.0;
	    for(q=0;q<m_j;q++){
		/* stuff used to compute Poisson mass function */
		Tabs = t_i[k]>t_j[q] ? t_i[k]-t_j[q] : t_j[q]-t_i[k]; /* time difference (absolute value) */
		transi = matint_ij(dnainfo->transi,s_i[k],s_j[q]); /* transitions */
		transv = matint_ij(dnainfo->transv,s_i[k],s_j[q]); /* transversions */
		common = matint_ij(dnainfo->nbcommon,s_i[k],s_j[q]); /* nb of nucleotides compared */

		if(common > 0){
		    nb_comp_k++;

		    /* direct ancestries */
		    if(t_i[k] > t_j[q]){
			/* transitions */
			Xi1 += gsl_ran_poisson_pdf((unsigned int) matint_ij(dnainfo->transi,s_i[k],s_j[q]), nu1*Tabs*common);
			/* transversions */
			Xi1 += gsl_ran_poisson_pdf((unsigned int) matint_ij(dnainfo->transv,s_i[k],s_j[q]), nu2*Tabs*common);
		    }

		    /* indirect ancestries */
		    /* transitions */
		    Xi3 += gsl_ran_poisson_pdf((unsigned int) matint_ij(dnainfo->transi,s_i[k],s_j[q]), nu1*(Tabs+2.0*tau)*common);
		    /* transversions */
		    Xi3 += gsl_ran_poisson_pdf((unsigned int) matint_ij(dnainfo->transv,s_i[k],s_j[q]), nu2*(Tabs+2.0*tau)*common);
		}
	    }

	    /* ancestries from S_i */
	    Xi2=0.0;
	    Xi4=0.0;
	    for(r=0;r<k;r++){
		/* stuff used to compute Poisson mass function */
		Tabs = t_i[k]>t_i[r] ? t_i[k]-t_i[r] : t_i[r]-t_i[k]; /* time difference (absolute value) */
		transi = matint_ij(dnainfo->transi,s_i[k],s_i[r]); /* transitions */
		transv = matint_ij(dnainfo->transv,s_i[k],s_i[r]); /* transversions */
		common = matint_ij(dnainfo->nbcommon,s_i[k],s_i[r]); /* nb of nucleotides compared */

		if(common > 0){
		    nb_comp_k++;

		    /* direct ancestries */
		    if(t_i[k] > t_i[r]){
			/* transitions */
			Xi2 += gsl_ran_poisson_pdf((unsigned int) matint_ij(dnainfo->transi,s_i[k],s_i[r]), nu1*Tabs*common);
			/* transversions */
			Xi2 += gsl_ran_poisson_pdf((unsigned int) matint_ij(dnainfo->transv,s_i[k],s_i[r]), nu2*Tabs*common);
		    }

		    /* indirect ancestries */
		    /* transitions */
		    Xi4 += gsl_ran_poisson_pdf((unsigned int) matint_ij(dnainfo->transi,s_i[k],s_i[r]), nu1*(Tabs+2.0*tau)*common);
		    /* transversions */
		    Xi4 += gsl_ran_poisson_pdf((unsigned int) matint_ij(dnainfo->transv,s_i[k],s_i[r]), nu2*(Tabs+2.0*tau)*common);
		}
	    }

	    /* compute likelihood if sequences were compared */
	    if(nb_comp_k>0){
		nb_comp++;

		/* likelihood (p_s_i^k | s_i^1, ..., s_i^{k-1}, S_j, i<-j)*/
		Pk = alpha * (Xi1 + Xi2) + (1.0-alpha) * (Xi3 + Xi4);
		/* printf("\nvalue of pk: %f \n",Pk); */

		/* update general likelihood */
		out += log(Pk);
	    }
	}
    }

    /* Two cases when likelihood can't be computed:
       - no pair of sequences to be compared
       - no compared sequence had nucleotide in common */
    /* printf("\nnb comparisons: %d", nb_comp); */
    if(nb_comp==0){
	out = log(param->weightNaGen);
    } else {
	/* log-like are averaged, not summed */
	out = out/((double) nb_comp);
    }


    /* RETURN */
    return out;
} /* end compute_genlike */








/*
  =========================
  === TESTING FUNCTIONS ===
  =========================
*/

//* * TESTING WITH R INTERFACE *\/ */
/* void test_genlike(unsigned char *DNAbinInput, int *n, int *length, int *s_i, int *s_j, double *t_i, double *t_j, int *m_i, int *m_j, double *nu1, double *nu2, double *alpha, double *tau, double *out){ */


/*     /\* CHECK THAT ARGUMENT ARE PASSED ALL RIGHT *\/ */
/*     /\* printf("\ns_i:"); *\/ */
/*     /\* for(i=0;i<*m_i;i++){ *\/ */
/*     /\* 	printf("%d ",s_i[i]); *\/ */
/*     /\* } *\/ */
/*     /\* printf("\ns_j:"); *\/ */
/*     /\* for(i=0;i<*m_j;i++){ *\/ */
/*     /\* 	printf("%d ",s_j[i]); *\/ */
/*     /\* } *\/ */
/*     /\* printf("\nt_i:"); *\/ */
/*     /\* for(i=0;i<*m_i;i++){ *\/ */
/*     /\* 	printf("%.4f ",t_i[i]); *\/ */
/*     /\* } *\/ */
/*     /\* printf("\nt_j:"); *\/ */
/*     /\* for(i=0;i<*m_j;i++){ *\/ */
/*     /\* 	printf("%.4f ",t_j[i]); *\/ */
/*     /\* } *\/ */
/*     /\* printf("\nnu1: %.4f nu2: %.4f alpha:%.4f tau: %.4f out:%.4f\n",*nu1,*nu2, *alpha, *tau, *out); *\/ */


/*     /\* MAKE DNAINFO AND PARAM *\/ */
/*     parameters *param = createParam(); */
/*     list_dnaseq * dna = DNAbin2list_dnaseq(DNAbinInput, n, length); */
/*     /\* print_list_dnaseq(dna); *\/ */
/*     dna_dist *distinfo = compute_dna_distances(dna); */

/*     param->weightNaGen = 0.0000001; /\* near zero if no data *\/ */

/*     /\* COMPUTE LIKELIHOOD *\/ */
/*     *out = genlike_ij(s_i, s_j, t_i, t_j, *m_i, *m_j, *nu1, *nu2, *alpha, *tau, distinfo, param); */

/*     /\* FREE STUFF *\/ */
/*     free_list_dnaseq(dna); */
/*     free_dna_dist(distinfo); */
/*     freeParam(param); */

/* } */







/* /\* PURE C TESTING *\/ */
/* int main(){ */
/*     int N=5, L=10; */
/*     int i,j, s_i[2] = {0,1}, s_j[3] = {2,3,4}; */
/*     double alpha=0.5, tau=2.0, t_i[2] = {0.0, 10.0}, t_j[3] = {12.0, 50.0, 100.0}, out, nu1 = 0.01, nu2=0.02; */
/*     parameters *param = createParam(); */
/*     param->weightNaGen = 0.001; /\* near zero if no data *\/ */

/*     /\* create a list of sequences *\/ */
/*     list_dnaseq * test = create_list_dnaseq(N, L); */

/*     for(i=0;i<N;i++){ */
/* 	for(j=0;j<L;j++){ */
/* 	    if(i*j % 5 ==0) test->list[i]->seq[j] = 'a'; */
/* 	    else if(i*j % 3 ==0)test->list[i]->seq[j] = 't'; */
/* 	    else if(i*j % 2 ==0)test->list[i]->seq[j] = 'g'; */
/* 	    else test->list[i]->seq[j] = 'c'; */
/* 	} */
/*     } */

/*     for(i=5;i<L;i++) */
/* 	test->list[0]->seq[i] = '-'; */
/*     for(i=0;i<5;i++) */
/* 	test->list[N-1]->seq[i] = '-'; */

/*     print_list_dnaseq(test); */

/*     /\* COMPUTE DISTANCES *\/ */
/*     dna_dist *distinfo = compute_dna_distances(test); */
/*     print_dna_dist(distinfo); */


/*     /\* COMPUTE LIKELIHOOD *\/ */
/*     /\* likelihood, sequences available *\/ */
/*     printf("\nsi: %d %d\n",s_i[0], s_i[1]); */
/*     printf("\nsj: %d %d %d\n",s_j[0], s_j[1], s_j[2]); */
/*     printf("\nti: %f %f\n",t_i[0], t_i[1]); */
/*     printf("\ntj: %f %f %f\n",t_j[0], t_j[1], t_j[2]); */
/*     printf("\nnu1: %f nu2: %f alpha:%f tau: %f\n",nu1,nu2, alpha, tau); */

/*     out = genlike_ij(s_i, s_j, t_i, t_j, 2, 3, nu1, nu2, alpha, tau, distinfo, param); */
/*     log(out); */
/*     //printf("\npseudo log-likelihood (i: 0,1; j: 2,3,4): %.5f\n", out); */

/*     /\* likelihood, only within-host sequences available *\/ */
/*     /\* int *nothing=NULL; *\/ */
/*     /\* out = genlike_ij(s_i, nothing, t_i, t_j, 2, 0, nu1, nu2, alpha, tau, distinfo, param); *\/ */
/*     //printf("\npseudo log-likelihood (i: 0,1; j:NA): %.5f\n", out); */

/*     /\* likelihood, no sequence available *\/ */
/*     /\* out = genlike_ij(nothing, nothing, t_i, t_j, 0, 0, nu1, nu2, alpha, tau, distinfo, param); *\/ */
/*     /\* printf("\npseudo log-likelihood (i: NA, j:NA): %.5f\n", out); *\/ */


/*     printf("\n"); */

/*     /\* free and return *\/ */
/*     free_list_dnaseq(test); */
/*     free_dna_dist(distinfo); */
/*     freeParam(param); */

/*     return 0; */
/* } */



/*
  gcc instructions

  gcc -o genlike -g -Wall alloc.c matvec.c genclasses.c distances.c param.c genlike.c -lgsl -lgslcblas && ./genlike

  valgrind --leak-check=full genlike

*/



# endif
