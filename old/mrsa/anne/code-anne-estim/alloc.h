#ifndef __COMMON_H
#include "common.h"
#endif



#ifndef __ALLOC_H
#define __ALLOC_H

nb_data *createNbData();
void freeNbData(nb_data *);

raw_data *createRawData(nb_data *);
void freeRawData(raw_data *);

aug_data *createAugData();
void freeAugData(aug_data *);
void copyAugData(aug_data *, aug_data *);

parameters *createParam();
void freeParam(parameters *);
void copyParam(parameters * , parameters * );

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

