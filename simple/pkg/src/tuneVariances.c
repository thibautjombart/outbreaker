#if 0


#include "common.h"
#include "init.h"
#include "InputOutput.h"
#include "logL.h"
#include "mcmc.h"
#include "moves.h"
#include "prior.h"
#include "alloc.h"
#include "tuneVariances.h"

/* extern gsl_rng * rng; */


int IsAcceptOKBeta(int i, int j, isAcceptOK * acceptOK){
    int res=1;
    if(gsl_matrix_get(acceptOK->IsAccOK_beta,i,j) < 3){res=0;}
    return res;
}





int IsAcceptOKBetaWardOut(isAcceptOK * acceptOK){
    int res=1;
    if(acceptOK->IsAccOK_betaWardOut < 3){res=0;}
    return res;
}






int IsAcceptOKBetaOutOut(isAcceptOK * acceptOK){
    int res=1;
    if(acceptOK->IsAccOK_betaOutOut < 3){res=0;}
    return res;
}






int IsAcceptOKMu(isAcceptOK * acceptOK){
    int res=1;
    if(acceptOK->IsAccOK_mu < 3){res=0;}
    return res;
}






int IsAcceptOKSigma(isAcceptOK * acceptOK){
    int res=1;
    if(acceptOK->IsAccOK_sigma < 3){res=0;}
    return res;
}






int IsAcceptOKNu1(isAcceptOK * acceptOK){
    int res=1;
    if(acceptOK->IsAccOK_nu1 < 3){res=0;}
    return res;
}






int IsAcceptOKKappa(isAcceptOK * acceptOK){
    int res=1;
    if(acceptOK->IsAccOK_kappa < 3){res=0;}
    return res;
}






/* int IsAcceptOKTau(isAcceptOK * acceptOK){ */
/*     int res=1; */
/*     if(acceptOK->IsAccOK_tau < 3){res=0;} */
/*     return res; */
/* } */






/* int IsAcceptOKAlpha(isAcceptOK * acceptOK){ */
/*     int res=1; */
/*     if(acceptOK->IsAccOK_alpha < 3){res=0;} */
/*     return res; */
/* } */






void updateAccept(NbProposals * nbProp, acceptance * accept){
    int i,j;
    double prov;

    for(i=0 ; i<2 ; i++)
	{
	    for(j=0 ; j<2 ; j++)
		{
		    prov=gsl_matrix_get(accept->PourcAcc_beta,i,j);
		    gsl_matrix_set(accept->PourcAcc_beta,i,j,prov/gsl_matrix_get(nbProp->NbProp_beta,i,j));
		}
	}

    accept->PourcAcc_betaWardOut/=nbProp->NbProp_betaWardOut;
    accept->PourcAcc_betaOutOut/=nbProp->NbProp_betaOutOut;

    accept->PourcAcc_mu/=nbProp->NbProp_mu;
    accept->PourcAcc_sigma/=nbProp->NbProp_sigma;

    accept->PourcAcc_nu1/=nbProp->NbProp_nu1;
    accept->PourcAcc_kappa/=nbProp->NbProp_kappa;
}






void updateMCMCSettingsBeta(int i, int j, NbProposals * nbProp,acceptance * accept,isAcceptOK *acceptOK, mcmcInternals *MCMCSettings){
    double prov,perc;

    perc=gsl_matrix_get(accept->PourcAcc_beta,i,j)/gsl_matrix_get(nbProp->NbProp_beta,i,j);
    if(perc<0.1 || perc>0.4)
	{
	    prov=gsl_matrix_get(MCMCSettings->Sigma_beta,i,j);
	    gsl_matrix_set(MCMCSettings->Sigma_beta,i,j,prov*perc/0.25);
	    gsl_matrix_set(acceptOK->IsAccOK_beta,i,j,0);
	}else
	{
	    prov=gsl_matrix_get(acceptOK->IsAccOK_beta,i,j);
	    gsl_matrix_set(acceptOK->IsAccOK_beta,i,j,prov+1);
	}
}






void updateMCMCSettingsBetaWardOut(NbProposals * nbProp,acceptance * accept,isAcceptOK *acceptOK, mcmcInternals *MCMCSettings){
    double perc;

    perc=accept->PourcAcc_betaWardOut/nbProp->NbProp_betaWardOut;
    if(perc<0.1 || perc>0.4)
	{
	    MCMCSettings->Sigma_betaWardOut*=perc/0.25;
	    acceptOK->IsAccOK_betaWardOut=0;
	}else
	{
	    acceptOK->IsAccOK_betaWardOut+=1;
	}
}






void updateMCMCSettingsBetaOutOut(NbProposals * nbProp,acceptance * accept,isAcceptOK *acceptOK, mcmcInternals *MCMCSettings){
    double perc;

    perc=accept->PourcAcc_betaOutOut/nbProp->NbProp_betaOutOut;
    if(perc<0.1 || perc>0.4)
	{
	    MCMCSettings->Sigma_betaOutOut*=perc/0.25;
	    acceptOK->IsAccOK_betaOutOut=0;
	}else
	{
	    acceptOK->IsAccOK_betaOutOut+=1;
	}
}






void updateMCMCSettingsMu(NbProposals * nbProp,acceptance * accept,isAcceptOK *acceptOK, mcmcInternals *MCMCSettings){
    double perc;

    perc=accept->PourcAcc_mu/nbProp->NbProp_mu;
    if(perc<0.1 || perc>0.4)
	{
	    MCMCSettings->Sigma_mu*=perc/0.25;
	    acceptOK->IsAccOK_mu=0;
	}else
	{
	    acceptOK->IsAccOK_mu+=1;
	}
}






void updateMCMCSettingsSigma(NbProposals * nbProp,acceptance * accept,isAcceptOK *acceptOK, mcmcInternals *MCMCSettings){
    double perc;

    perc=accept->PourcAcc_sigma/nbProp->NbProp_sigma;
    if(perc<0.1 || perc>0.4)
	{
	    MCMCSettings->Sigma_sigma*=perc/0.25;
	    acceptOK->IsAccOK_sigma=0;
	}else
	{
	    acceptOK->IsAccOK_sigma+=1;
	}
}






void updateMCMCSettingsNu1(NbProposals * nbProp,acceptance * accept,isAcceptOK *acceptOK, mcmcInternals *MCMCSettings){
    double perc;

    perc=accept->PourcAcc_nu1/nbProp->NbProp_nu1;
    if(perc<0.1 || perc>0.4)
	{
	    MCMCSettings->Sigma_nu1*=perc/0.25;
	    acceptOK->IsAccOK_nu1=0;
	}else
	{
	    acceptOK->IsAccOK_nu1+=1;
	}
}





void updateMCMCSettingsKappa(NbProposals * nbProp,acceptance * accept,isAcceptOK *acceptOK, mcmcInternals *MCMCSettings){
    double perc;

    perc=accept->PourcAcc_kappa/nbProp->NbProp_kappa;
    if(perc<0.1 || perc>0.4)
	{
	    MCMCSettings->Sigma_kappa*=perc/0.25;
	    acceptOK->IsAccOK_kappa=0;
	}else
	{
	    acceptOK->IsAccOK_kappa+=1;
	}
}






/* void updateMCMCSettingsTau(NbProposals * nbProp,acceptance * accept,isAcceptOK *acceptOK, mcmcInternals *MCMCSettings){ */
/*     double perc; */

/*     perc=accept->PourcAcc_tau/nbProp->NbProp_tau; */
/*     if(perc<0.1 || perc>0.4) */
/* 	{ */
/* 	    MCMCSettings->Sigma_tau*=perc/0.25; */
/* 	    acceptOK->IsAccOK_tau=0; */
/* 	}else */
/* 	{ */
/* 	    acceptOK->IsAccOK_tau+=1; */
/* 	} */
/* } */





/* void updateMCMCSettingsAlpha(NbProposals * nbProp,acceptance * accept,isAcceptOK *acceptOK, mcmcInternals *MCMCSettings){ */
/*     double perc; */

/*     perc=accept->PourcAcc_alpha/nbProp->NbProp_alpha; */
/*     if(perc<0.1 || perc>0.4) */
/* 	{ */
/* 	    MCMCSettings->Sigma_alpha*=perc/0.25; */
/* 	    acceptOK->IsAccOK_alpha=0; */
/* 	}else */
/* 	{ */
/* 	    acceptOK->IsAccOK_alpha+=1; */
/* 	} */
/* } */



#endif
