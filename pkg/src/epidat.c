/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), January 2012.
  Licence: GPL >=2.

*/

#include "common.h"
#include "epidat.h"




/*
  ============
  CONSTRUCTORS
  ============
*/

/* CREATE EPIDAT OBJECT - ONE DNA SEQUENCE */
struct epidat * create_epidat(int length){
	/* ALLOCATE OUTPUT */
	struct epidat *out = (struct epidat *) malloc(sizeof(struct epidat));
	if(out==NULL){
		fprintf(stderr, "\n[in: classes.c->create_epidat]\nNo memory left for creating DNA sequence. Exiting.\n");
		exit(1);
	}

	/* FILL/ALLOCATE CONTENT */
	out->seq =  (char*) malloc(length*sizeof(char));
	out->length = length;
	return out;
}




/* CREATE LIST_EPIDAT OBJECT - A LIST OF ALIGNED DNA SEQUENCES */
struct list_epidat * create_list_epidat(int n, int length){
	int i;

	/* ALLOCATE OUTPUT */
	struct list_epidat *out = (struct list_epidat *) malloc(sizeof(struct list_epidat));
	if(out==NULL){
		fprintf(stderr, "\n[in: classes.c->create_list_epidat]\nNo memory left for creating list of DNA sequences. Exiting.\n");
		exit(1);
	}

	/* FILL/ALLOCATE CONTENT */
	out->list = (struct epidat **) malloc(n*sizeof(struct epidat *));
	if(out->list==NULL){
		fprintf(stderr, "\n[in: classes.c->create_list_epidat]\nNo memory left for creating list of DNA sequences. Exiting.\n");
		exit(1);
	}

	for(i=0;i<n;i++){
		out->list[i] = create_epidat(length);
	}
	out->n = n;
	out->length = length;
	return out;
}








/*
  ===========
  DESTRUCTORS
  ===========
*/

void free_epidat(struct epidat *in){
	if(in!=NULL){
		free(in->seq);
		free(in);
	}
}


void free_list_epidat(struct list_epidat *in){
	int i;
	if(in!=NULL){
		for(i=0;i<in->n;i++){
			free_epidat(in->list[i]);
		}
		free(in->list);
		free(in);
	}
}



/*
  ===================
  AUXILIARY FUNCTIONS
  ===================
*/

void print_epidat(struct epidat *in){
	int i;
	for(i=0;i<in->length;i++){
		printf("%c", in->seq[i]);
	}
	printf("\n");
}





void print_list_epidat(struct list_epidat *in){
	int i;
	printf("\nList of DNA %d sequences (size: %d)\n", in->n, in->length);
	for(i=0;i<in->n;i++){
		printf("%d: ", i+1);
		print_epidat(in->list[i]);
	}
	printf("\n");
}



/*
  ==================
  EXTERNAL FUNCTIONS
  ==================
*/


/* convert DNAbin class to list_epidat */
struct list_epidat * DNAbin2list_epidat(unsigned char * in, int *n, int *length){
	int i, j, count=0;

	/* CREATE OUTPUT */
	struct list_epidat *out = create_list_epidat(*n,*length);

	/* FILL IN THE OUTPUT */
	for(i=0;i<*n;i++){
		for(j=0;j<*length;j++){
			out->list[i]->seq[j] = DNAbin2char(in[count++]);
		}
	}

	printf("\nlist_epidat in C:\n");
	print_list_epidat(out);
	printf("\n");

	/* RETURN */
	return out;
}





/*
  =======
  TESTING
  =======
*/


/* int main(){ */
/* 	const int N=5, L=30; */
/* 	int i,j; */

/* 	/\* create a list of sequences *\/ */
/* 	struct list_epidat * test = create_list_epidat(N, L); */

/* 	for(i=0;i<N;i++){ */
/* 		for(j=0;j<L;j++){ */
/* 			if(i*j % 5 ==0) test->list[i]->seq[j] = 'a'; */
/* 			else if(i*j % 3 ==0)test->list[i]->seq[j] = 't'; */
/* 			else if(i*j % 2 ==0)test->list[i]->seq[j] = 'g'; */
/* 			else test->list[i]->seq[j] = 'c'; */
/* 		} */
/* 	} */

/* 	print_list_epidat(test); */

/* 	/\* free and return *\/ */
/* 	free_list_epidat(test); */
/* 	return 0; */
/* } */


/* gcc instructions:

gcc -o genclasses genclasses.c && ./genclasses


valgrind --leak-check=full genclasses


*/



