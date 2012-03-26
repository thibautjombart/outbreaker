/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), January 2012.
  Licence: GPL >=2.

*/

#ifndef __EPIDAT_H
#define  __EPIDAT_H



/*
  ===============
  DATA STRUCTURES
  ===============
*/
struct epidat{
	char *seq;
	int length;
};



struct list_epidat{
	struct epidat **list;
	int n, length; /* n: nb of sequences; length: length of sequences */
};





/*
  ============
  CONSTRUCTORS
  ============
*/
struct epidat * create_epidat(int length);

struct list_epidat * create_list_epidat(int n, int length);




/*
  ===========
  DESTRUCTORS
  ===========
*/

void free_epidat(struct epidat *in);

void free_list_epidat(struct list_epidat *in);




/*
  ===================
  AUXILIARY FUNCTIONS
  ===================
*/

void print_epidat(struct epidat *in);

void print_list_epidat(struct list_epidat *in);





/*
  ==================
  EXTERNAL FUNCTIONS
  ==================
*/

struct list_epidat * DNAbin2list_epidat(unsigned char * in, int *n, int *length);


#endif
