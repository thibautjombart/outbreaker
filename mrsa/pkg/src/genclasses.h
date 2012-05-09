/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), January 2012.
  Licence: GPL >=2.

*/

#ifndef __GENCLASSES_H
#define __GENCLASSES_H


/*
  ===============
  DATA STRUCTURES
  ===============
*/
typedef struct{
	char *seq;
	int length;
} dnaseq;



typedef struct{
	dnaseq **list;
	int n, length; /* n: nb of sequences; length: length of sequences */
} list_dnaseq;





/*
  ============
  CONSTRUCTORS
  ============
*/
 dnaseq * create_dnaseq(int length);

 list_dnaseq * create_list_dnaseq(int n, int length);




/*
  ===========
  DESTRUCTORS
  ===========
*/

void free_dnaseq(dnaseq *in);

void free_list_dnaseq(list_dnaseq *in);




/*
  ===================
  AUXILIARY FUNCTIONS
  ===================
*/

void print_dnaseq(dnaseq *in);

void print_list_dnaseq(list_dnaseq *in);

char DNAbin2char(unsigned char in);





/*
  ==================
  EXTERNAL FUNCTIONS
  ==================
*/

list_dnaseq * DNAbin2list_dnaseq(unsigned char * in, int *n, int *length);

void copy_dnaseq(dnaseq *in, dnaseq *out);


#endif
