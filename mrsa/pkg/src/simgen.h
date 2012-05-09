/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

*/

#ifndef __SIMGEN_H
#define __SIMGEN_H


/*
   This is the structure representing genetic sequences during epidemics. 
   Each patient 'i' hosts a nb of lineages nbLineages[i], stored as a 
   list_dnaseq object in dna[i].
*/
typedef struct{
    int nbPatients; /* nb of patients */
    int length; /* haplotype length */
    int maxNbLineages;
    /* int *nbLineages; /\* nb of lineages in each patient; vector of size nPatients*\/ */
    list_dnaseq **dna; /* haplotypes in each patient; vector of size nPatients*/
} epid_dna;





/*
  ===================
  AUXILIARY FUNCTIONS
  ===================
*/

char transi(char in);

char transv(char in, gsl_rng *rng);




/*
  ============
  CONSTRUCTORS
  ============
*/

epid_dna * create_epid_dna(int nPatients, int maxNlineages, int haploLength);



/*
  ===========
  DESTRUCTORS
  ===========
*/

void free_epid_dna(epid_dna *in);





/*
   ==========================
   === EXTERNAL FUNCTIONS ===
   ==========================
*/

void print_epid_dna(epid_dna *in);

dnaseq *create_haplo(int length, gsl_rng *rng);

void evolve_haplo(dnaseq *in, double nu1, double nu2, double t, gsl_rng *rng);

void replicate_haplo(dnaseq *in, dnaseq *out,double nu1, double nu2, double t, gsl_rng *rng);

void make_distant_lineage(dnaseq *in, dnaseq *out, int dist, double nu1, double nu2, gsl_rng *rng);

void evolve_epid_dna(epid_dna *in, int *ances, double mu_dist, double sigma_dist, double lambda_nlin, double nu1, double nu2, int *dates, gsl_rng *rng);


list_dnaseq *swab_dna_patient(epid_dna *in, int patient, int N, double nu1, double nu2, int colDate, int *swabDates, gsl_rng *rng);


list_dnaseq *sample_epid_dna(epid_dna *in, nb_data *nb_data, raw_data *data, aug_data *augData, double lambda_nseq, double nu1, double nu2, gsl_rng *rng);


#endif
