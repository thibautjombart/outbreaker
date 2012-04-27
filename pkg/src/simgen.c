/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), January 2012.
  Licence: GPL >=2.

*/

#include "common.h"
#include "alloc.h"
#include "genclasses.h"
#include "simgen.h"
#include "matvec.h"



/*
  ============
  CONSTRUCTORS
  ============
*/

epid_dna * create_epid_dna(int nbPatients, int maxNlineages, int haploLength){
    int i;

    /* ALLOCATE OUTPUT */
    epid_dna *out = (epid_dna *) malloc(sizeof(epid_dna));
    if(out==NULL){
	fprintf(stderr, "\n[in: simgen.c->create_epid_dna]\nNo memory left for creating list of DNA sequences. Exiting.\n");
	exit(1);
    }

    /* FILL/ALLOCATE CONTENT */
    /* out->nbLineages = (int *) calloc(nbPatients, sizeof(int)); */
    /* if(out->nbLineages==NULL){ */
    /* 	fprintf(stderr, "\n[in: simgen.c->create_epid_dna]\nNo memory left for creating list of DNA sequences. Exiting.\n"); */
    /* 	exit(1); */
    /* } */

    out->dna = (list_dnaseq **) calloc(nbPatients, sizeof(list_dnaseq *));
    for(i=0;i<nbPatients;i++){
	out->dna[i] = create_list_dnaseq(maxNlineages, haploLength);
    }

    out->nbPatients = nbPatients;
    out->length = haploLength;
    out->maxNbLineages = maxNlineages;
    return out;
}





/*
  ===========
  DESTRUCTORS
  ===========
*/

void free_epid_dna(epid_dna *in){
    int i;
    /* free(in->nbLineages); */
    for(i=0;i<in->nbPatients;i++){
	free_list_dnaseq(in->dna[i]);
    }
    free(in->dna);
    free(in);
}





/*
  ===================
  AUXILIARY FUNCTIONS
  ===================
*/

/* find one transition */
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
	fprintf(stderr, "\n[in: simgen.c->transi]\nUnknown character '%c'.\n", in);
	exit(1);
    }
}



/* find one transversion */
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
	fprintf(stderr, "\n[in: simgen.c->transv]\nUnknown character '%c'.\n", in);
	exit(1);
    }
}


/* PRINT */
void print_epid_dna(epid_dna *in){
    int i;
    printf("\n== Epidemic DNA info==\n");
    printf("\n%d patients   haplotype length:%d\n", in->nbPatients, in->length);
    printf("\nnb of lineages per patient:\n");
    for(i=0;i<in->nbPatients;i++) printf("%d ", in->dna[i]->n);
    printf("\n");

    for(i=0;i<in->nbPatients;i++){
	printf("\nPatient %d:", i);
	print_list_dnaseq(in->dna[i]);
    }

    printf("\n");
    fflush(stdout);
}



/*
  ==================
  EXTERNAL FUNCTIONS
  ==================
*/

/* create a new haplotype of a given length */
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
   - nu2: rate of transversions; nu2=kappa*nu1
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






/* Replicate an haplotype using:
   - nu1: rate of transitions 
   - nu2: rate of transversions; nu2=kappa*nu1
   - double t: time of evolution between in and out
*/
void replicate_haplo(dnaseq *in, dnaseq *out, double nu1, double nu2, double t, gsl_rng *rng){
    /* copy haplotype */
    copy_dnaseq(in, out);

    /* evolve haplotype */
    evolve_haplo(out, nu1, nu2, t, rng);
}





/* Create a new lineage at a given distance from an existing haplotype, using:
   - in: a reference haplotype
   - dist: the distance in nb of mutations from the reference
   - nu1: rate of transitions
   - nu2: rate of transversions; nu2=kappa*nu1
   - double t: time of evolution between in and out
*/
void make_distant_lineage(dnaseq *in, dnaseq *out, int dist, double nu1, double nu2, gsl_rng *rng){
    int i, posi, nbTransiTransver[2];
    double p[2];
    p[0] = nu1;
    p[1] = nu2;

    /* COPY HAPLOTYPE */
    copy_dnaseq(in, out);

    /* EVOLVE HAPLOTYPE */
    /* find nb of transitions / transversions */
    gsl_ran_multinomial (rng, 2, dist, p, (unsigned int *) nbTransiTransver);

  
    /* handle transitions */
    for(i=0;i<nbTransiTransver[0];i++){
	posi = gsl_rng_uniform_int(rng, out->length);
	out->seq[posi] = transi(out->seq[posi]);
    }

    /* handle transversions */
    for(i=0;i<nbTransiTransver[1];i++){
	posi = gsl_rng_uniform_int(rng, out->length);
	out->seq[posi] = transv(out->seq[posi], rng);
    }

} /* end make_distant_lineage */





/*
  Simulate the evolution of pathogens in a set of patients, given:
  - ances: a vector of ancestries (values indicate 'infecting' patient; -1=external)
  - mu_dist,sigma_dist: the mean/sd for the distance between lineages (lognormal)
  - lambda_nlin: lambda for the number of lineages in a external infection (nlin ~ pois(lambda)+1)
  - nu1: rate of transitions
  - nu2: rate of transversions; nu2=kappa*nu1
  - dates: dates at which patients are infected; ordered by increasing dates (forward-time)

  Sensible values:
  - mu_dist,sigma_dist: ([3,4] , 0.1)
  - lambda_nlin: [0,1]

*/

void evolve_epid_dna(epid_dna *in, int *ances, double mu_dist, double sigma_dist, double lambda_nlin, double nu1, double nu2, int *dates, gsl_rng *rng){
    int i, j, k, N = in->nbPatients, L=in->length, nlin, deltaT, dist, seqIdx;

    /* REFERENCE HAPLOTYPE */
    dnaseq *ref = create_haplo(L, rng);

    /* FIND ORDER IN THE COLONIZATION DATES */
    /* define auxiliary vectors */
    vec_int *inDates, *sortedDates, *datesIdx;
    inDates = create_vec_int(N);
    sortedDates = create_vec_int(N);
    datesIdx = create_vec_int(N);
    for(i=0;i<N;i++){
	inDates->values[i] = dates[i];
    }
    sort_vec_int(inDates, sortedDates, datesIdx);


    /* OUTSIDE TRANSMISSIONS FIRST */
    for(k=0;k<N;k++){
	i = vecint_i(datesIdx, k);
	/* printf("\nSimulating outside haplotype %d - date:%d\n",k, dates[i]); */
	if(ances[i]<0){ /* i.e. colonization from outside */
			/* determine the number of lineages */
	    nlin = 1 + gsl_ran_poisson(rng, lambda_nlin);

	    /* make sure not to exceed the max number of sequences per patient */
	    if(nlin>in->dna[i]->n) nlin = in->dna[i]->n;

	    /* draw haplotypes */
	    for(j=0;j<nlin;j++){
		dist = (int) gsl_ran_lognormal(rng, mu_dist, sigma_dist);
		make_distant_lineage(ref, in->dna[i]->list[j], dist, nu1, nu2, rng);
	    }

	    /* update the number of sequences in patient */
	    in->dna[i]->n = nlin;

	    /* free unused sequences */
	    for(j=nlin;j<in->maxNbLineages;j++){
		free_dnaseq(in->dna[i]->list[j]);
	    }
	}
    }


    /* PATIENT->PATIENT TRANSMISSIONS */
    for(k=0;k<N;k++){
	i = vecint_i(datesIdx, k);
	/* printf("\nSimulating transmitted haplotype %d - date:%d\n",k, dates[i]); */
	/* make sure the ancestor is known */
	if(ances[i]>=N){
	    fprintf(stderr, "\n[in: simgen.c->evolve_epid_dna]\nUnknown ancestor index %d (%d patients)\n", ances[i], N);
	    exit(1);
	}
	if(ances[i]>=0){ /* i.e. infection from a known patient */
	    /* determine the number of lineages */
	    nlin = 1 + gsl_ran_poisson(rng, lambda_nlin);

	    /* make sure not to exceed the max number of sequences */
	    if(nlin>in->maxNbLineages) nlin = in->maxNbLineages;
	    if(nlin==0){
		printf("\nLikely issue: patient %d infects patient %d but has no know pathogen sequence.\n", ances[i], i);
	    }

	    /* replicate haplotypes */
	    deltaT = dates[i] - dates[ances[i]];
	    for(j=0;j<nlin;j++){
		seqIdx = gsl_rng_uniform_int(rng, in->dna[ances[i]]->n);
		replicate_haplo(in->dna[ances[i]]->list[seqIdx], in->dna[i]->list[j], nu1, nu2, deltaT, rng);
	    }

	    /* update the number of sequences in patient i */
	    in->dna[i]->n = nlin;

	    /* free unused sequences */
	    for(j=nlin;j<in->maxNbLineages;j++){
		free_dnaseq(in->dna[i]->list[j]);
	    }
	}
    }

    /* free local pointers */
    free_dnaseq(ref);
    free_vec_int(sortedDates);
    free_vec_int(datesIdx);
    free_vec_int(inDates);
} /* end evolve_epid_dna */





/*
  get N sequences from a patient 'patient' at times collecDates
  - evolution of the sequences occurs between colDate and collecDates
  - 'patient': index of the patient (from 0 to NbPatients-1)
  - 'N': nb of sequences sampled; length of vector 'collecDates'
*/
list_dnaseq *swab_dna_patient(epid_dna *in, int patient, int N, double nu1, double nu2, int colonDate, int *collecDates, gsl_rng *rng){
    int i, ancesSeqId, deltaT;

    /* CHECKS */
    if(patient>=in->nbPatients || patient<0){
	fprintf(stderr, "\n[in: simgen.c->get_swab_patient]\nSwabbed patient index %d unknown (there are %d patients).\n", patient, in->nbPatients);
	exit(1);
    }


    /* CREATE OUTPUT AND FILL IT IN  */
    list_dnaseq *out = create_list_dnaseq(N, in->length);

    for(i=0;i<N;i++){
	if(collecDates[i]<colonDate){
	    fprintf(stderr, "\n[in: simgen.c->get_swab_patient]\nSwab date (%d) before colonisation data (%d).\n", collecDates[i], colonDate);
	    exit(1);
	}

	/* find deltaT */
	deltaT = collecDates[i] - colonDate;

	/* select ancestral strain at random */
	ancesSeqId = gsl_rng_uniform_int(rng, in->dna[patient]->n);

	/* replicate it */
	replicate_haplo(in->dna[patient]->list[ancesSeqId], out->list[i], nu1, nu2, deltaT, rng);
    }

    /* return */
    return out;
} /* end swab_dna_patient */






/*
  Sample genetic data from outbreak and fill in epid and dna data
  - 'dna_data': list_dnaseq used to output the dna sequences; must not have been allocated before!
  - lambda_nseq: parameter of the Poisson distribution giving the number of sequences for a given positive swab
  - nu1, nu2: transitions, transversions
  - colonDates: vector of colonisation dates for each patient, length nbData->NbPatients

*/
list_dnaseq *sample_epid_dna(epid_dna *in, nb_data *nb_data, raw_data *data, aug_data *augData, double lambda_nseq, double nu1, double nu2, gsl_rng *rng){
    int i, j, counter=0;
    int *colonDates=augData->C;
    int nbPatients=data->NbPatients;
    int *swabDates; /* temporary vector storing positive swab dates */
    list_dnaseq *out;

    int **listCollecDates = (int **) calloc(nbPatients, sizeof(int *)); /* temporary list of vectors storing collection dates */

    if(listCollecDates == NULL){
	fprintf(stderr, "\n[in: simgen.c->sample_epid_dna]\nNo memory left for creating listCollecDates. Exiting.\n");
	exit(1);
    }

    /* create temporary list storing swab results */
    list_dnaseq **listSwabSeq = (list_dnaseq **) malloc(nbPatients*sizeof(list_dnaseq *));
    if(listSwabSeq == NULL){
	fprintf(stderr, "\n[in: simgen.c->sample_epid_dna]\nNo memory left for creating listSwabSeq. Exiting.\n");
	exit(1);
    }


    /* INITIALIZE NB OF SEQUENCES PER PATIENTS */
    for(i=0;i<nbPatients;i++) data->M[i] = 0;
    data->NbSequences = 0;


    /* FOR EACH PATIENT, GET SEQUENCES FROM SWABS */
    for(i=0;i<nbPatients;i++){
	/* create temporary vector storing positive swab dates for i */
	/* (sampling procedure does not handle gsl_vectors) */
	swabDates = calloc(nb_data->NbPosSwabs[i], sizeof(int));
	if(swabDates == NULL){
	    fprintf(stderr, "\n[in: simgen.c->sample_epid_dna]\nNo memory left for creating swabDates. Exiting.\n");
	    exit(1);
	}

	/* fill swab dates in temporary vector */
	for(j=0;j<nb_data->NbPosSwabs[i];j++){
	    swabDates[j] = (int) gsl_vector_get(data->P[i], j);
	}

	/* determine the total number of sequences for patient i */
	data->M[i] = gsl_ran_poisson(rng, (double) nb_data->NbPosSwabs[i]*lambda_nseq);
	data->NbSequences = data->NbSequences + data->M[i];

	/* determine the collection dates for these sequences  */
	listCollecDates[i] = calloc(data->M[i], sizeof(int));
	for(j=0;j<data->M[i];j++){
	    /* vector used locally */
	    listCollecDates[i][j] = swabDates[gsl_rng_uniform_int(rng, nb_data->NbPosSwabs[i])];
	}


	/* get DNA sequences from swabs */
	listSwabSeq[i] = swab_dna_patient(in, i, data->M[i], nu1, nu2, colonDates[i], listCollecDates[i], rng);

	/* free local (patient-specific) pointers */
	free(swabDates);
    } /* end for patient i */


    /* BUILD FINAL LIST_DNASEQ FROM THE LIST OF SAMPLES */

    /* allocate output for dna sequences */
    /* dna_data must not have been allocated before */
    out = create_list_dnaseq(data->NbSequences, in->length);

    /* reallocate vector of collection times */
    /* msg = realloc(data->Tcollec, data->NbSequences*sizeof(int)); */
    free(data->Tcollec);
    data->Tcollec = (int *) calloc(data->NbSequences,sizeof(int)) ;
    if(data->Tcollec == NULL){
	fprintf(stderr, "\n[in: simgen.c->sample_epid_dna]\nNo memory left for creating swabDates. Exiting.\n");
	exit(1);
    }


    counter = 0;
    for(i=0;i<nbPatients;i++){
	/* realloc S vectors (S[i]: indices of sequences sampled in patient i */
	/* msg = realloc(data->S[i], data->M[i]*sizeof(int)); */
	free(data->S[i]);
	data->S[i] = (int *) calloc(data->M[i],sizeof(int));
	if(data->S[i] == NULL){
	    fprintf(stderr, "\n[in: simgen.c->sample_epid_dna]\nNo memory left for creating swabDates. Exiting.\n");
	    exit(1);
	}

	/* fill in data */
	for(j=0;j<data->M[i];j++){
	    /* fill in sequence index data */
	    data->S[i][j] = counter;

	    /* copy DNA sequences */
	    copy_dnaseq(listSwabSeq[i]->list[j],out->list[counter]);

	    /* copy collection times */
	    data->Tcollec[counter++] = listCollecDates[i][j];
	}
    }

    /* REALLOCATE MATRICES OF LAMBDA AND TAU IN AUGDATA */
    gsl_matrix_free(augData->alpha);
    gsl_matrix_free(augData->tau);
    augData->alpha = gsl_matrix_calloc(data->NbSequences, data->NbSequences);
    if(augData->alpha == NULL){
	fprintf(stderr, "\n[in: simgen.c->sample_epid_dna]\nNo memory left for creating matrix of alpha_kq in augData. Exiting.\n");
	exit(1);
    }

    augData->tau = gsl_matrix_calloc(data->NbSequences, data->NbSequences);
    if(augData->tau == NULL){
	fprintf(stderr, "\n[in: simgen.c->sample_epid_dna]\nNo memory left for creating matrix of tau_kq in augData. Exiting.\n");
	exit(1);
    }


    /* FREE LOCAL VARIABLES AND RETURN */
    for(i=0;i<nbPatients;i++) {
    	free_list_dnaseq(listSwabSeq[i]);
    	free(listCollecDates[i]);
    }
    free(listSwabSeq);
    free(listCollecDates);

    return out;
} /* end sample_epid_dna */





/*
  =======
  TESTING
  =======
*/

/* /\* TESTS OF BASIC ROUTINES *\/ */
/* int main(){ */
/*     time_t t = time(NULL); /\* time in seconds, used to change the seed of the random generator *\/ */
/*     const gsl_rng_type *typ; */
/*     gsl_rng_env_setup(); */
/*     typ=gsl_rng_default; */
/*     gsl_rng * rng=gsl_rng_alloc(typ); */
/*     gsl_rng_set(rng,t); /\* changes the seed of the random generator *\/ */

/*     int i; */

/*     dnaseq *seq1, *seq2, *seq3; */

/*     /\* transitions *\/ */
/*     printf("\n== Transitions =="); */
/*     printf("\na:%c", transi('a')); */
/*     printf("\ng:%c", transi('g')); */
/*     printf("\nt:%c", transi('t')); */
/*     printf("\nc:%c", transi('c')); */

/*     printf("\n\n== Transversions =="); */
/*     for(i=0;i<5;i++) printf("\na:%c", transv('a',rng)); */
/*     printf("\n"); */
/*     for(i=0;i<5;i++) printf("\ng:%c", transv('g',rng)); */
/*     printf("\n"); */
/*     for(i=0;i<5;i++) printf("\nt:%c", transv('t',rng)); */
/*     printf("\n"); */
/*      for(i=0;i<5;i++) printf("\nc:%c", transv('c',rng)); */
/*     printf("\n"); */


/*     printf("\n== Haplotype creation ==\n"); */
/*     seq1 = create_haplo(30, rng); */
/*     print_dnaseq(seq1); */
    
/*     printf("\n== Haplotype copy ==\n"); */
/*     seq2 = create_dnaseq(30); */
/*     copy_dnaseq(seq1, seq2); */
/*     print_dnaseq(seq2); */

/*     printf("\n== Haplotype replication ==\n"); */
/*     printf("\nref:"); */
/*     seq3 = create_dnaseq(30); */
/*     replicate_haplo(seq1, seq3, 0.1, 0.2, 1.0, rng); */
/*     printf("\nref:");  */
/*     print_dnaseq(seq2); */
/*     printf("\nnew:");  */
/*     print_dnaseq(seq3); */

/*     printf("\n== Evolution across several time steps ==\n"); */
/*     copy_dnaseq(seq1, seq3); */
/*     for(i=0;i<20;i++){ */
/* 	evolve_haplo(seq3, 0.05, 0.1, 1.0, rng); */
/* 	print_dnaseq(seq3); */
/*     } */

/*     free_dnaseq(seq1); */
/*     free_dnaseq(seq2); */
/*     free_dnaseq(seq3); */
/*     gsl_rng_free(rng); */
/*     return 0; */
/* } */




/* /\* TESTS OF EVOLVE_EPID_DNA *\/ */
/* int main(){ */
/*     int i,j; */

/*     /\* time_t t = time(NULL); /\\* time in seconds, used to change the seed of the random generator *\\/ *\/ */
/*     time_t t = 1; /\* time in seconds, used to change the seed of the random generator *\/ */
/*     const gsl_rng_type *typ; */
/*     gsl_rng_env_setup(); */
/*     typ=gsl_rng_default; */
/*     gsl_rng * rng=gsl_rng_alloc(typ); */
/*     gsl_rng_set(rng,t); /\* changes the seed of the random generator *\/ */

/*     epid_dna *out = create_epid_dna(10, 5, 50); /\* nb patients, max nb lineages, haplo length *\/ */

/*     int ances[10] = {-1, -1, 1, 1, 2, 3, 3, 2, 6, 6}; /\* vector indicating ancestries *\/ */
/*     int dates[10] = {0, 0, 1, 1, 1, 2, 2, 2, 4, 50}; /\* collonisation dates *\/ */
/*     double mu_dist=3.0, sigma_dist=0.1, lambda_nlin=2; */
/*     double nu1=0.02, nu2=0.05; */
/*     double lambdaNseq=1.0; */

/*     /\* OUTBREAK SIMULATION *\/ */
/*     evolve_epid_dna(out, ances, mu_dist, sigma_dist, lambda_nlin, nu1, nu2, dates, rng); */

/*     print_epid_dna(out); */


/*     /\* SAMPLING PROCEDURE *\/ */

/*     /\* make nb_data and raw_data *\/ */
/*     /\* 10 indiv in outbreak, 5 with positive swabs *\/ */
/*     nb_data *nb = createNbData(10, 100, 0, 5); */
/*     for(i=0;i<10;i++) nb->NbAdmissions[i] = 1; */

/*     nb->NbColonisedPatients = 5; */
 
/*     /\* get nb of positive swabs per patients *\/ */
/*     for(i=0;i<10;i++){ */
/* 	nb->NbPosSwabs[i] = i < 5 ? i+1: 0; */
/*     } */

/*     raw_data *data = createRawData(nb); */

/*     for(i=0;i<10;i++){ */
/* 	/\* resize the vector of positive swab dates *\/ */
/* 	if(data->P[i] != NULL) { */
/* 	    gsl_vector_free(data->P[i]); */
/* 	    data->P[i] = gsl_vector_alloc(nb->NbPosSwabs[i]); */
/* 	} */

/* 	/\* fill in positive swab dates *\/ */
/* 	if(i<5){ */
/* 	    for(j=0;j<nb->NbPosSwabs[i];j++){ */
/* 		gsl_vector_set(data->P[i], j, gsl_ran_flat(rng, 5.0, 20.0)); */
/* 		nb->indexColonisedPatients[i]=i; */
/* 	    } */
/* 	} */
/*     } */

/*     /\* printf("\nInitial nb data: \n"); *\/ */
/*     /\* print_nbData(nb); *\/ */
/*     /\* fflush(stdout); *\/ */

/*     /\* printf("\n\nInitial raw data: \n"); *\/ */
/*     /\* print_rawData(data); *\/ */
/*     /\* fflush(stdout); *\/ */

/*     list_dnaseq *dna; */

/*     dna = sample_epid_dna(out, nb, data, lambdaNseq, nu1, nu2, dates, rng); */

/*     printf("\nSampling done. \n"); */
/*     fflush(stdout); */

/*     printf("\n\n= Final data = \n"); */
/*     print_rawData(data); */

/*     printf("\n\n= Final sample of DNA sequences = \n"); */
/*     print_list_dnaseq(dna); */

/*     free_epid_dna(out); */
/*     freeNbData(nb); */
/*     freeRawData(data); */
/*     free_list_dnaseq(dna); */
/*     gsl_rng_free(rng); */
/*     return 0; */
/* } */






/* 

   gcc instructions:

   gcc -o simgen alloc.c matvec.c genclasses.c simgen.c -lgsl -lgslcblas -g -Wall && ./simgen

   valgrind --leak-check=full simgen


*/



