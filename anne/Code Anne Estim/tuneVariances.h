#ifndef __COMMON_H

#include "common.h"
#endif

#ifndef __TUNEVARIANCES_H
#define __TUNEVARIANCES_H

int IsAcceptOKBeta(int, int, isAcceptOK * );
int IsAcceptOKBetaWardOut(isAcceptOK * );
int IsAcceptOKBetaOutOut(isAcceptOK * );
int IsAcceptOKMu(isAcceptOK * );
int IsAcceptOKSigma(isAcceptOK * );
int IsAcceptOKNu1(isAcceptOK * );
int IsAcceptOKNu2(isAcceptOK * );
int IsAcceptOKTau(isAcceptOK * );
int IsAcceptOKAlpha(isAcceptOK * );

void updateAccept(NbProposals * , acceptance * );

void updateMCMCSettingsBeta(int , int , NbProposals * ,acceptance * ,isAcceptOK *, mcmcInternals *);
void updateMCMCSettingsBetaWardOut(NbProposals * ,acceptance * ,isAcceptOK *, mcmcInternals *);
void updateMCMCSettingsBetaOutOut(NbProposals * ,acceptance * ,isAcceptOK *, mcmcInternals *);
void updateMCMCSettingsMu(NbProposals * ,acceptance * ,isAcceptOK *, mcmcInternals *);
void updateMCMCSettingsSigma(NbProposals * ,acceptance * ,isAcceptOK *, mcmcInternals *);
void updateMCMCSettingsNu1(NbProposals * ,acceptance * ,isAcceptOK *, mcmcInternals *);
void updateMCMCSettingsNu2(NbProposals * ,acceptance * ,isAcceptOK *, mcmcInternals *);
void updateMCMCSettingsTau(NbProposals * ,acceptance * ,isAcceptOK *, mcmcInternals *);
void updateMCMCSettingsAlpha(NbProposals * ,acceptance * ,isAcceptOK *, mcmcInternals *);

#endif

