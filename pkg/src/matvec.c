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

/* CREATE A VECTOR OF INTEGERS OF SIZE N */
struct vec_int * create_vec_int(int n){
	struct vec_int *out = (struct vec_int *) malloc(sizeof(struct vec_int));
	if(out == NULL){
		fprintf(stderr, "\n[in: distances.c->create_vec_int]\nNo memory left for creating vector of integers. Exiting.\n");
		exit(1);
	}

	/* NOTE out->values is not allocated when n=0 */
	if(n>0){
		out->values = (int *) calloc(n, sizeof(int));
		if(out->values == NULL){
			fprintf(stderr, "\n[in: distances.c->create_vec_int]\nNo memory left for creating vector of integers. Exiting.\n");
			exit(1);
		}
	}

	out->length = n;

	return(out);
}





/* /\* CREATE A VECTOR OF INTEGERS OF SIZE N INITIALIZED TO ZERO *\/ */
/* struct vec_int * create_vec_int_zero(int n){ */
/* 	struct vec_int *out = (struct vec_int *) malloc(sizeof(struct vec_int)); */
/* 	if(out == NULL){ */
/* 		fprintf(stderr, "\n[in: distances.c->create_vec_int]\nNo memory left for creating vector of integers. Exiting.\n"); */
/* 		exit(1); */
/* 	} */

/* 	out->values = (int *) calloc(n, sizeof(int)); */
/* 	if(out->values == NULL){ */
/* 		fprintf(stderr, "\n[in: distances.c->create_vec_int]\nNo memory left for creating vector of integers. Exiting.\n"); */
/* 		exit(1); */
/* 	} */

/* 	out->length = n; */

/* 	return(out); */
/* } */






/* CREATE EMPTY MAT_INT BETWEEN N OBJECTS */
/* (values initialized to 0) */
struct mat_int * create_mat_int(int n){
	int i;
	struct mat_int *out;

	/* allocate output */
	out = (struct mat_int *) malloc(sizeof(struct mat_int));
	if(out == NULL){
		fprintf(stderr, "\n[in: distances.c->create_mat_int]\nNo memory left for creating distance matrix. Exiting.\n");
		exit(1);
	}

	/* fill in content */
	out->rows = (struct vec_int **) calloc(n, sizeof(struct vec_int *));
	if(out->rows == NULL){
		fprintf(stderr, "\n[in: distances.c->create_mat_int]\nNo memory left for creating distance matrix. Exiting.\n");
		exit(1);
	}

	for(i=0;i<n;i++){
		out->rows[i] = create_vec_int(n);
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


void free_vec_int(struct vec_int *in){
	if(in->length > 0) free(in->values);
	free(in);
}


void free_mat_int(struct mat_int *in){
	int i;
	if(in->n > 0) {
		for(i=0;i<in->n;i++)
			free_vec_int(in->rows[i]);
	}
	free(in->rows);
	free(in);
}






/*
   ===============================
   === MAIN EXTERNAL FUNCTIONS ===
   ===============================
*/

int vecint_i(struct vec_int *in, int i){
	if(i >= in->length) {
		fprintf(stderr, "\nTrying to access value %d in a vector of size %d\n",i,in->length);
		exit(1);
	}
	return in->values[i];
}




int matint_ij(struct mat_int *in, int i, int j){
	if(i >= in->n) {
		fprintf(stderr, "\nTrying to access item %d in a list of size %d\n",i,in->n);
		exit(1);
	}
	return vecint_i(in->rows[i], j);
}




void print_vec_int(struct vec_int *in){
	int i;
	printf("\nVector of %d values: ", in->length);
	/* for(i=0;i<in->length;i++) printf("%d ", in->values[i]); */
	for(i=0;i<in->length;i++) printf("%d ", vecint_i(in,i));
	printf("\n");
}




void print_mat_int(struct mat_int *in){
	int i,j;

	for(i=0;i<in->n;i++){
	printf("\n");
		for(j=0;j<in->n;j++)
			/* printf("%d ", in->rows[i]->values[j]); */
			printf("%d ", matint_ij(in,i,j));
	}
	printf("\n");
}






/*
   =========================
   === TESTING FUNCTIONS ===
   =========================
*/


/* int main(){ */
/* 	const int N = 10; */
/* 	struct mat_int * test = create_mat_int(N); */

/* 	print_mat_int (test); */

/* 	free_mat_int(test); */
/* 	return 0; */
/* } */



/*
  gcc instructions

  gcc -o matvec matvec.c && ./matvec

  valgrind --leak-check=full matvec

*/
