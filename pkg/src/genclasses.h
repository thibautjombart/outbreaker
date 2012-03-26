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
struct dnaseq{
	char *seq;
	int length;
};



struct list_dnaseq{
	struct dnaseq **list;
	int n, length; /* n: nb of sequences; length: length of sequences */
};





/*
  ============
  CONSTRUCTORS
  ============
*/
struct dnaseq * create_dnaseq(int length);

struct list_dnaseq * create_list_dnaseq(int n, int length);




/*
  ===========
  DESTRUCTORS
  ===========
*/

void free_dnaseq(struct dnaseq *in);

void free_list_dnaseq(struct list_dnaseq *in);




/*
  ===================
  AUXILIARY FUNCTIONS
  ===================
*/

void print_dnaseq(struct dnaseq *in);

void print_list_dnaseq(struct list_dnaseq *in);

char DNAbin2char(unsigned char in);





/*
  ==================
  EXTERNAL FUNCTIONS
  ==================
*/

struct list_dnaseq * DNAbin2list_dnaseq(unsigned char * in, int *n, int *length);


#endif
