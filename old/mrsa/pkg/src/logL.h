#ifndef __COMMON_H
#include "common.h"
#endif

#ifndef __DISTANCES_H
#include "distances.h"
#endif

#ifndef __LOGVRAIS_H
#define __LOGVRAIS_H

double ObsLevelPerCase(int , raw_data * , nb_data *, aug_data *, parameters *);
double ObsLevel(raw_data * , nb_data *, aug_data *, parameters *);

double DurationColonPerCase(int , aug_data *, parameters *);

double DurationColon(aug_data *, parameters *);

double ColonPerCase(int, raw_data *, nb_data *, aug_data *, dna_dist *, parameters *);

double Colon(raw_data *, nb_data *, aug_data *, dna_dist *, parameters *);

double fullLoglikelihoodPerCase(int , raw_data * , nb_data *, aug_data *, dna_dist *, parameters *);

double fullLoglikelihood(raw_data * , nb_data *, aug_data *, dna_dist *, parameters *);

double fullLoglikelihoodWithPrior(raw_data * , nb_data *, aug_data *, dna_dist *, parameters *);

#endif








