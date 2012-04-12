#ifndef __COMMON_H
#include "common.h"
#endif

#ifndef __DISTANCES_H
#include "distances.h"
#endif

#ifndef __INPUTOUTPUT_H
#define __INPUTOUTPUT_H


/* import functions R -> C */
void importNbData(int *nbAdmVec, int *nbPosSwab, int *nbNegSwab, int *nbColPatients, int *nbPatients, int *duration, int *idxColPatients, int *nbSeqPat, nb_data *nb);

void importRawData(int *wardVec, int *tAdmVec, int *tDisVec, int *tPosSwab, int *tNegSwab, int *hospPres, int *idxSeqVec, int *totNbSeq, double *tCollecVec, nb_data *nb, raw_data *data);


/* output functions */
void prepAllFiles(output_files *, int NbPatients);
void writeAllFiles(output_files *, parameters *, nb_data *, raw_data *, aug_data *, dna_dist *);



/* deprecated */
void readFakeNbData(nb_data *);
void readFakeData(nb_data *, raw_data *);

#endif


