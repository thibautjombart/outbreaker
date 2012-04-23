#ifndef __COMMON_H
#include "common.h"
#endif

#ifndef __MOVES_H
#define __MOVES_H

void moveBeta(int, int, mcmcInternals *, parameters *, raw_data *, nb_data *, aug_data *, dna_dist *, acceptance *, NbProposals *);

void moveAllBeta(mcmcInternals * , parameters * , raw_data * , nb_data *, aug_data *, dna_dist *, acceptance *, NbProposals *);

void moveBetaOut(char, mcmcInternals * , parameters * , raw_data * , nb_data *, aug_data *, dna_dist *, acceptance *, NbProposals *);

void moveSp(parameters * , raw_data * , nb_data *, aug_data *, acceptance *);

void moveSe(parameters * , raw_data * , nb_data *, aug_data *, acceptance *);

void movePi(parameters * , raw_data * , aug_data *, acceptance *);

void moveDurationColon(char, mcmcInternals * , parameters * , raw_data * , aug_data *, acceptance *, NbProposals *);

void moveMutationRate(int , mcmcInternals * , parameters * , raw_data * , nb_data *, aug_data *, dna_dist *, acceptance *, NbProposals *);

void moveTau(mcmcInternals * , parameters * , raw_data * , nb_data *, aug_data *, dna_dist *, acceptance *, NbProposals *);

void moveAlpha(mcmcInternals * , parameters * , raw_data * , nb_data *, aug_data *, dna_dist *, acceptance *, NbProposals *);

void moveC(int, mcmcInternals * , parameters * , raw_data * , nb_data *, aug_data *, dna_dist *);

void moveE(int, mcmcInternals * , parameters * , raw_data * , nb_data *, aug_data *, dna_dist *);

void moveCandE(int, mcmcInternals * , parameters * , raw_data * , nb_data *, aug_data *, dna_dist *);


#endif
