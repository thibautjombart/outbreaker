/*
  Coded by Thibaut Jombart (t.jombart@imperial.ac.uk), September 2011.
  Distributed with the epidemics package for the R software.
  Licence: GPL >=2.

*/

#ifndef __SIMEPID_H
#define __SIMEPID_H

void fillingBeds(parameters *param, hospDurationParam *paramHosp, nb_data *nbData, raw_data *data, aug_data *augData, int *indexInfector, int *NbCases, int *bedNb, gsl_matrix * NbDischarges0, gsl_matrix * NbDischarges1, int *TotalNbDischarges0, int *TotalNbDischarges1, gsl_vector ** isInHosp,gsl_matrix * indexOfPatientsInWard0, gsl_matrix * indexOfPatientsInWard1, int *AlreadyColonised, gsl_vector ** isColonised);

void SimulColonisationInBeds(parameters *param, nb_data *nbData, raw_data *data, aug_data *augData, int *indexInfector, gsl_matrix * indexOfPatientsInWard0, gsl_matrix * indexOfPatientsInWard1, int *AlreadyColonised, gsl_vector ** isColonised, int *indexI0, int *indexI1, int *indexS0, int *indexS1, 	int *indexIncid0, 	int *indexIncid1);

void SimulTesting(parameters *param, nb_data *nbData, raw_data *data, gsl_matrix * indexOfPatientsInWard0, gsl_matrix * indexOfPatientsInWard1, gsl_vector ** isColonised);

int SimulEpid(parameters *param, hospDurationParam *paramHosp, nb_data *nbData, raw_data *data, aug_data *augData, int *indexInfector);

#endif
