#ifndef __COMMON_H
#include "common.h"
#endif

#ifndef __PRIOR_H
#define __PRIOR_H

double logpriorBeta(int , int , parameters * );

double logpriorBetaWardOut(parameters * );

double logpriorBetaOutOut(parameters * );

double logpriorPi(parameters * );

double logpriorSp(parameters * );

double logpriorSe(parameters * );

double logpriorMu(parameters * );

double logpriorSigma(parameters * );

double logpriorNu1(parameters * param);

double logpriorNu2(parameters * param);

/* double logpriorAlpha(parameters * param); */

/* double logpriorTau(parameters * param); */

double logprior (parameters * );

#endif
