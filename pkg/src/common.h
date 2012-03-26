/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), January 2012.
  Licence: GPL >=2.

*/


#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
/* Calls to GNU Scientific Library */
#include <gsl/gsl_rng.h> /* random nb generators */
#include <gsl/gsl_randist.h> /* rng with specific distributions */

#define NEARZERO 0.0000000001
#define TRUE 1
#define FALSE 0

typedef short int bool;
