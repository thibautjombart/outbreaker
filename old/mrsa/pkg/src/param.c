/* /\* */
/*   Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), January 2012. */
/*   Licence: GPL >=2. */

/* *\/ */

/* #include "common.h" */
/* #include "matvec.h" */
/* #include "param.h" */



/* /\* */
/*    ==================== */
/*    === CONSTRUCTORS === */
/*    ==================== */
/* *\/ */


/* struct param * create_param(){ */
/* 	struct param * out = (struct param *) malloc(sizeof(struct param)); */
/* 	if(out==NULL){ */
/* 		fprintf(stderr, "\n[in: param.c->create_param]\nNo memory left for creating list of parameters. Exiting.\n"); */
/* 		exit(1); */
/* 	} */

/* 	return out; */
/* } */





/* /\* */
/*    =================== */
/*    === DESTRUCTORS === */
/*    =================== */
/* *\/ */


/* void free_param(struct param * in){ */
/* 	free(in); */
/* } */







/* /\* */
/*    =============================== */
/*    === MAIN EXTERNAL FUNCTIONS === */
/*    =============================== */
/* *\/ */





/* /\* */
/*    ========================= */
/*    === TESTING FUNCTIONS === */
/*    ========================= */
/* *\/ */


/* /\* int main(){ *\/ */
/* /\* 	struct mat_int * test = create_mat_int(N); *\/ */

/* /\* 	print_mat_int (test); *\/ */

/* /\* 	free_mat_int(test); *\/ */
/* /\* 	return 0; *\/ */
/* /\* } *\/ */



/* /\* */
/*   gcc instructions */

/*   gcc -o param param.c && ./param */

/*   valgrind --leak-check=full param */

/* *\/ */
