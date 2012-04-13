#ifndef __COMMON_H

#include "common.h"
#endif

#ifndef __ALLOC_H
#define __ALLOC_H

nb_data *createNbData(int NbPatients, int T, int NbSequences);
void freeNbData(nb_data *);
void print_nbData(nb_data *);

raw_data *createRawData(nb_data *);
void freeRawData(raw_data *);
void print_rawData(raw_data *);

aug_data *createAugData(int NbPatients, int T);
void freeAugData(aug_data *);
void copyAugData(aug_data *, aug_data *);
void print_augData(aug_data *);

parameters *createParam();
void freeParam(parameters *);
void copyParam(parameters * , parameters * );
void print_param(parameters *);

mcmcInternals *createMcmcInternals();
void printStdProp(mcmcInternals *);
void freeMcmcInternals(mcmcInternals*);

acceptance *createAcceptance();
void reInitiateAcceptance(acceptance *);
void printAcceptance(acceptance *, NbProposals *);
void freeAcceptance(acceptance*);

isAcceptOK *createIsAcceptOK();
void freeIsAcceptOK(isAcceptOK *);

NbProposals *createNbProposals();
void reInitiateNbProp(NbProposals *);
void freeNbProposals(NbProposals *);

output_files *createFILES(char*);
void freeFILES(output_files *);

#endif

