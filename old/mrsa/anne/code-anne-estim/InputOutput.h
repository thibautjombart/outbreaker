//Input / Output
#ifndef __COMMON_H
#include "common.h"
#endif

#ifndef __INPUTOUTPUT_H
#define __INPUTOUTPUT_H

void readFakeNbData(nb_data *);
void readFakeData(nb_data *, raw_data *);

void prepAllFiles(output_files *);

void writeAllFiles(output_files * ,parameters * , nb_data *, raw_data * , aug_data *);

#endif

