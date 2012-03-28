#ifndef __COMMON_H

#include "common.h"
#endif

#ifndef __ALLOC_H
#define __ALLOC_H

nb_data *createNbData();
void freeNbData(nb_data *);

raw_data *createRawData();
void freeRawData(raw_data *);

parameters *createParam();
void readParameters(char* , parameters * );
void freeParam(parameters *);

#endif

