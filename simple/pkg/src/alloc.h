#ifndef __COMMON_H

#include "common.h"
#endif

#ifndef __ALLOC_H
#define __ALLOC_H



aug_data *alloc_AugData(int NbPatients, int T, int NbSequences);
void freeAugData(aug_data *);
void copyAugData(aug_data *, aug_data *);
void print_augData(aug_data *);

parameters *alloc_Param();
void freeParam(parameters *);
void copyParam(parameters * , parameters * );
void print_param(parameters *);

mcmcInternals *alloc_McmcInternals();
void printStdProp(mcmcInternals *);
void freeMcmcInternals(mcmcInternals*);

acceptance *alloc_Acceptance();
void reInitiateAcceptance(acceptance *);
void printAcceptance(acceptance *, NbProposals *);
void freeAcceptance(acceptance*);

isAcceptOK *alloc_IsAcceptOK();
void freeIsAcceptOK(isAcceptOK *);

NbProposals *alloc_NbProposals();
void reInitiateNbProp(NbProposals *);
void freeNbProposals(NbProposals *);

output_files *alloc_FILES(char*);
void freeFILES(output_files *);

#endif

