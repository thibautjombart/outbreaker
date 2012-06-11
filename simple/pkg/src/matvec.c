/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), January 2012.
  Licence: GPL >=2.

*/

#include "common.h"
#include "matvec.h"



/*
  ====================
  === CONSTRUCTORS ===
  ====================
*/

/* ALLOC A VECTOR OF INTEGERS OF SIZE N */
vec_int * alloc_vec_int(int n){
    vec_int *out = (vec_int *) malloc(sizeof(vec_int));
    if(out == NULL){
	fprintf(stderr, "\n[in: matvec.c->alloc_vec_int]\nNo memory left for creating vector of integers. Exiting.\n");
	exit(1);
    }

    /* NOTE out->values is not allocated when n=0 */
    if(n>0){
	out->values = (int *) calloc(n, sizeof(int));
	if(out->values == NULL){
	    fprintf(stderr, "\n[in: matvec.c->alloc_vec_int]\nNo memory left for creating vector of integers. Exiting.\n");
	    exit(1);
	}
    }

    out->length = n;

    return(out);
}





/* ALLOC A VECTOR OF DOUBLEEGERS OF SIZE N */
vec_double * alloc_vec_double(int n){
    vec_double *out = (vec_double *) malloc(sizeof(vec_double));
    if(out == NULL){
	fprintf(stderr, "\n[in: matvec.c->alloc_vec_double]\nNo memory left for creating vector of double. Exiting.\n");
	exit(1);
    }

    /* NOTE out->values is not allocated when n=0 */
    if(n>0){
	out->values = (double *) calloc(n, sizeof(double));
	if(out->values == NULL){
	    fprintf(stderr, "\n[in: matvec.c->alloc_vec_double]\nNo memory left for creating vector of double. Exiting.\n");
	    exit(1);
	}
    }

    out->length = n;

    return(out);
}




/* /\* ALLOC A VECTOR OF INTEGERS OF SIZE N INITIALIZED TO ZERO *\/ */
/* vec_int * alloc_vec_int_zero(int n){ */
/* 	vec_int *out = (vec_int *) malloc(sizeof(vec_int)); */
/* 	if(out == NULL){ */
/* 		fprintf(stderr, "\n[in: matvec.c->alloc_vec_int]\nNo memory left for creating vector of integers. Exiting.\n"); */
/* 		exit(1); */
/* 	} */

/* 	out->values = (int *) calloc(n, sizeof(int)); */
/* 	if(out->values == NULL){ */
/* 		fprintf(stderr, "\n[in: matvec.c->alloc_vec_int]\nNo memory left for creating vector of integers. Exiting.\n"); */
/* 		exit(1); */
/* 	} */

/* 	out->length = n; */

/* 	return(out); */
/* } */






/* ALLOC EMPTY MAT_INT BETWEEN N OBJECTS */
/* (values initialized to 0) */
mat_int * alloc_mat_int(int n){
    int i;
    mat_int *out;

    /* allocate output */
    out = (mat_int *) malloc(sizeof(mat_int));
    if(out == NULL){
	fprintf(stderr, "\n[in: matvec.c->alloc_mat_int]\nNo memory left for creating distance matrix. Exiting.\n");
	exit(1);
    }

    /* fill in content */
    out->rows = (vec_int **) calloc(n, sizeof(vec_int *));
    if(out->rows == NULL){
	fprintf(stderr, "\n[in: matvec.c->alloc_mat_int]\nNo memory left for creating distance matrix. Exiting.\n");
	exit(1);
    }

    for(i=0;i<n;i++){
	out->rows[i] = alloc_vec_int(n);
    }

    out->n = n;

    /* return */
    return out;
}





/* ALLOC EMPTY MAT_DOUBLE BETWEEN N OBJECTS */
/* (values initialized to 0) */
mat_double * alloc_mat_double(int n){
    int i;
    mat_double *out;

    /* allocate output */
    out = (mat_double *) malloc(sizeof(mat_double));
    if(out == NULL){
	fprintf(stderr, "\n[in: matvec.c->alloc_mat_double]\nNo memory left for creating distance matrix. Exiting.\n");
	exit(1);
    }

    /* fill in content */
    out->rows = (vec_double **) calloc(n, sizeof(vec_double *));
    if(out->rows == NULL){
	fprintf(stderr, "\n[in: matvec.c->alloc_mat_double]\nNo memory left for creating distance matrix. Exiting.\n");
	exit(1);
    }

    for(i=0;i<n;i++){
	out->rows[i] = alloc_vec_double(n);
    }

    out->n = n;

    /* return */
    return out;
}




/*
  ===================
  === DESTRUCTORS ===
  ===================
*/


void free_vec_int(vec_int *in){
    if(in->length > 0) free(in->values);
    free(in);
}


void free_mat_int(mat_int *in){
    int i;
    if(in->n > 0) {
	for(i=0;i<in->n;i++)
	    free_vec_int(in->rows[i]);
    }
    free(in->rows);
    free(in);
}




void free_vec_double(vec_double *in){
    if(in->length > 0) free(in->values);
    free(in);
}


void free_mat_double(mat_double *in){
    int i;
    if(in->n > 0) {
	for(i=0;i<in->n;i++)
	    free_vec_double(in->rows[i]);
    }
    free(in->rows);
    free(in);
}




/*
  ===============================
  === MAIN EXTERNAL FUNCTIONS ===
  ===============================
*/

int vec_int_i(vec_int *in, int i){
    if(i >= in->length) {
	fprintf(stderr, "\nTrying to access value %d in a vector of size %d\n",i,in->length);
	exit(1);
    }
    return in->values[i];
}




int mat_int_ij(mat_int *in, int i, int j){
    if(i >= in->n) {
	fprintf(stderr, "\nTrying to access item %d in a list of size %d\n",i,in->n);
	exit(1);
    }
    return vec_int_i(in->rows[i], j);
}





double vec_double_i(vec_double *in, int i){
    if(i >= in->length) {
	fprintf(stderr, "\nTrying to access value %d in a vector of size %d\n",i,in->length);
	exit(1);
    }
    return in->values[i];
}




double mat_double_ij(mat_double *in, int i, int j){
    if(i >= in->n) {
	fprintf(stderr, "\nTrying to access item %d in a list of size %d\n",i,in->n);
	exit(1);
    }
    return vec_double_i(in->rows[i], j);
}






/* print method */
void print_vec_int(vec_int *in){
    int i;
    printf("\nVector of %d values: ", in->length);
    /* for(i=0;i<in->length;i++) printf("%d ", in->values[i]); */
    for(i=0;i<in->length;i++) printf("%d ", vec_int_i(in,i));
    printf("\n");
}




/* print method */
void print_mat_int(mat_int *in){
    int i,j;

    for(i=0;i<in->n;i++){
	printf("\n");
	for(j=0;j<in->n;j++)
	    /* printf("%d ", in->rows[i]->values[j]); */
	    printf("%d ", mat_int_ij(in,i,j));
    }
    printf("\n");
}






/* print method */
void print_vec_double(vec_double *in){
    int i;
    printf("\nVector of %d values: ", in->length);
    /* for(i=0;i<in->length;i++) printf("%d ", in->values[i]); */
    for(i=0;i<in->length;i++) printf("%.3f ", vec_double_i(in,i));
    printf("\n");
}





/* print method */
void print_mat_double(mat_double *in){
    int i,j;

    for(i=0;i<in->n;i++){
	printf("\n");
	for(j=0;j<in->n;j++)
	    /* printf("%d ", in->rows[i]->values[j]); */
	    printf("%.3f ", mat_double_ij(in,i,j));
    }
    printf("\n");
}






/* alternative print method for gsl vectors */
void print_gsl_vector(gsl_vector *in, char format[256]){
    int i;
    for(i=0;i<in->size;i++){
	printf(format, in->data[i]);
    }
    printf("\n");
    fflush(stdout);
}



/* check if an integer 'x' is in a vector of integers, and returns the matching position */
int in_vec_int(int x, vec_int *vec){
    int i=0;
    while(i<vec->length && x!=vec_int_i(vec, i)) i++; /* note: condition needs to be in this order */
    if(i==vec->length || vec->length<1) return -1; /* -1 will mean: no match*/
    return i;
}



/* find max value in a vector of integers */
int max_vec_int(vec_int *vec){
    if(vec->length<1) return (int) NAN;
    int i, out=vec_int_i(vec,0);
    for(i=0;i<vec->length;i++) if(out<vec_int_i(vec,i)) out=vec_int_i(vec,i);
    return out;
}



/* find min value in a vector of integers */
int min_vec_int(vec_int *vec){
    if(vec->length<1) return (int) NAN;
    int i, out=vec_int_i(vec,0);
    for(i=0;i<vec->length;i++) if(out>vec_int_i(vec,i)) out=vec_int_i(vec,i);
    return out;
}



/* permut the values of a vector of integers */
void permut_vec_int(vec_int *in, gsl_rng * rng){
    if(in->length<1) return;
    /* if(in->length != out->length){ */
    /* 	fprintf(stderr, "\n[in: matvec.c->permut_vec_int]\nInconsistent vector sizes: in = %d, out = %d",in->length,out->length); */
    /* 	exit(1); */
    /* } */

    gsl_ran_shuffle(rng, in->values, in->length, sizeof (int));
}




/* sample values of a vector of integers with/without replacement */
void sample_vec_int(vec_int *in, vec_int *out, bool replace, gsl_rng * rng){
    if(in->length<1 || out->length<1) return;
    if(out->length > in->length && !replace){
  	fprintf(stderr, "\n[in: matvec.c->sample_vec_int]\nReplace is FALSE but sample size (%d) is bigger than input vector (%d)",out->length,in->length);
    	exit(1);
    }

    if(replace){
	gsl_ran_sample(rng, out->values, out->length, in->values, in->length, sizeof (int));
    } else {
	gsl_ran_choose(rng, out->values, out->length, in->values, in->length, sizeof (int));
    }
}






/* sort a vector of integers (ascending order)
   - in: input vector
   - out: vector of sorted values
   - idx: vector of indices (using R notation: out = in[idx])
*/
void sort_vec_int(vec_int *in, vec_int *out, vec_int *idx){
    if(out->length > in->length){
	fprintf(stderr, "\n[in: matvec.c->sort_vec_int_index]\nInput and output vectors have different lengths (in:%d out:%d)",in->length, out->length);
    	exit(1);
    }

    int i, j, curMin, curMinIdx;
    idx->length=0;

    for(j=0;j<in->length;j++){
	/* printf("\n- sorting value %d\n",j); */
	/* find minimal value and its index, discarding already sorted indices */
	curMin=max_vec_int(in);
	curMinIdx=0;

	for(i=0;i<in->length;i++){
	    /* if(in_vec_int(i, idx)>=0){ */
	    /* 	printf("\nindex %d found in array idx (position:%d)", i, in_vec_int(i, idx)); */
	    /* } */
	    if(in_vec_int(i, idx)<0 && curMin>=vec_int_i(in,i)) {
		/* printf("\nentering the loop, i=%d\n",i); */
		/* printf("\nidx:"); print_vec_int(idx); */
		/* printf("\nmatch i in idx: %d\n", in_vec_int(i, idx)); */
		curMin=vec_int_i(in,i);
		curMinIdx = i;
	    }
	}

	/* printf("\n- minimum %d found at index %d",curMin,curMinIdx); */

	/* update vectors */
	out->values[j] = curMin;
	idx->values[j] = curMinIdx;
	idx->length = idx->length + 1;
    }

} /* end sort_vec_int */




/* 
   =========
   COPYING
   =========
*/

void copy_vec_int(vec_int *in, vec_int *out){
    int i;
    if(in->length != out->length){
	fprintf(stderr, "\n[in: matvec.c->copy_vec_int]\nInput and output vectors have different lengths (in:%d out:%d)",in->length, out->length);
    	exit(1);
    }

    for(i=0;i<in->length;i++){
	out->values[i] = in->values[i];
    }
}



void copy_vec_double(vec_double *in, vec_double *out){
    int i;
    if(in->length != out->length){
	fprintf(stderr, "\n[in: matvec.c->copy_vec_double]\nInput and output vectors have different lengths (in:%d out:%d)",in->length, out->length);
    	exit(1);
    }

    for(i=0;i<in->length;i++){
	out->values[i] = in->values[i];
    }
}




void copy_mat_int(mat_int *in, mat_int *out){
    int i;
    if(in->n != out->n){
	fprintf(stderr, "\n[in: matvec.c->copy_mat_int]\nInput and output matrices have different numbers of rows (in:%d out:%d)",in->n, out->n);
    	exit(1);
    }

    for(i=0;i<in->n;i++){
	copy_vec_int(in->rows[i],out->rows[i]);
    }
}




void copy_mat_double(mat_double *in, mat_double *out){
    int i;
    if(in->n != out->n){
	fprintf(stderr, "\n[in: matvec.c->copy_mat_double]\nInput and output matrices have different numbers of rows (in:%d out:%d)",in->n, out->n);
    	exit(1);
    }

    for(i=0;i<in->n;i++){
	copy_vec_double(in->rows[i],out->rows[i]);
    }
}








/*
  =========================
  === TESTING FUNCTIONS ===
  =========================
*/


int main(){
    /* RANDOM NUMBER GENERATOR */
    time_t t = time(NULL); /* time in seconds, used to change the seed of the random generator */
    const gsl_rng_type *typ;
    gsl_rng_env_setup();
    typ=gsl_rng_default;
    gsl_rng * rng=gsl_rng_alloc(typ);
    gsl_rng_set(rng,t); /* changes the seed of the random generator */

    int i, N = 10;
    mat_int * test = alloc_mat_int(N);

    print_mat_int (test);
    free_mat_int(test);

    vec_int *myVec = alloc_vec_int(30), *toto;
    for(i=0;i<30;i++){
	myVec->values[i] = 30-i;
    }
    printf("\nVector\n");
    print_vec_int(myVec);
    
    printf("\nMin/Max: %d, %d\n", min_vec_int(myVec), max_vec_int(myVec));

    toto = alloc_vec_int(15);
    sample_vec_int(myVec, toto, 1, rng);
    printf("\n15 sampled values - with replacement \n");
    print_vec_int(toto);

    sample_vec_int(myVec, toto, 1, rng);
    printf("\nanother 15 sampled values - with replacement \n");
    print_vec_int(toto);

    sample_vec_int(myVec, toto, 0, rng);
    printf("\n15 sampled values - without replacement \n");
    print_vec_int(toto);

    sample_vec_int(myVec, toto, 0, rng);
    printf("\nanother 15 sampled values - without replacement \n");
    print_vec_int(toto);

    permut_vec_int(myVec,rng);
    printf("\npermut the vector myVec \n");
    print_vec_int(myVec);

    permut_vec_int(myVec,rng);
    printf("\nanother permutation of the vector myVec \n");
    print_vec_int(myVec);
    
    printf("\n== sorting a vector ==\n");
    vec_int *idx, *sortedVec;
    idx = alloc_vec_int(30);
    sortedVec = alloc_vec_int(30);
    sort_vec_int(myVec, sortedVec, idx);
    printf("\nvector to sort:");
    print_vec_int(myVec);
    printf("\nsorted vector:");
    print_vec_int(sortedVec);
    printf("\nindices:");
    print_vec_int(idx);

    /* copies */
    printf("\nCopy of sorted vector\n");
    vec_int *copyVec = alloc_vec_int(sortedVec->length);
    copy_vec_int(sortedVec,copyVec);
    print_vec_int(copyVec);


    printf("\nCopy of a matrix\n");
    mat_double *mat = alloc_mat_double(2);
    mat->rows[0]->values[0] = 1.1;
    mat->rows[0]->values[1] = 2.1;
    mat->rows[1]->values[1] = 666.0;
    printf("\nmat\n");
    print_mat_double(mat);
    mat_double *mat2 = alloc_mat_double(2);
    printf("\nmat2 before copy\n");
    print_mat_double(mat2);
    copy_mat_double(mat,mat2);
    printf("\nmat2 after copy\n");
    print_mat_double(mat);

    /* vec_int *a = alloc_vec_int(10); */
    /* for(i=0;i<10;i++){ */
    /* 	a->values[i] = i; */
    /* } */
    /* printf("\nvector a: \n"); */
    /* print_vec_int(a); */
    /* for(i=0;i<10;i++){ */
    /* 	printf("\n%d matches in a at position %d", i, in_vec_int(i,a)); */
    /* } */

    free_vec_int(toto);
    free_vec_int(myVec);
    free_vec_int(idx);
    free_vec_int(sortedVec);
    free_vec_int(copyVec);
    free_mat_double(mat);
    free_mat_double(mat2);
    gsl_rng_free(rng);
    return 0;
}



/*
  gcc instructions

  gcc -o matvec matvec.c -lgsl -lgslcblas && ./matvec

  valgrind --leak-check=full matvec

*/
