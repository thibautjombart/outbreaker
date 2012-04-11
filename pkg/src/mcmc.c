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

/******************************************************************************/
/* Metropolis-Hastings                                                        */
/******************************************************************************/


void metro (mcmcInternals * MCMCSettings, parameters * param, raw_data * data, nb_data *nb, aug_data *augData, acceptance *accept, isAcceptOK *acceptOK, NbProposals *nbProp, output_files * Files){
    int i,j,l;
    int NbChangeAugData = 10; /* TO BE TUNED  */
    int *AugDataToMove = (int *) calloc(NbChangeAugData,sizeof(int));
    if(AugDataToMove == NULL){
	fprintf(stderr, "\n[in: mcmc.c->metro]\nNo memory left for creating AugDataToMove. Exiting.\n");
	exit(1);
    }

    int nbStepsVar=1000; /* step at which variances of proposal laws are updated */
    int VarPropOK = 0;
    gsl_matrix *isAcceptOKBeta = gsl_matrix_calloc(2,2);
    double isAcceptOKBetaWardOut;
    double isAcceptOKBetaOutOut;
    double isAcceptOKMu;
    double isAcceptOKSigma;
    double isAcceptOKNu1;
    double isAcceptOKNu2;
    double isAcceptOKAlpha;
    double isAcceptOKTau;
    int EndBurnIn=0;
    int NbParamToBeTuned = 12;

    /* printf("Starting metro (p=%d)\n",p); */
    /* fflush(stdout); */

    /* writeAllFiles(Files, param, nb, data, augData); */

    for(l=0 ; l < MCMCSettings->NbSimul ; l++)
	{
	    /*printf("%d\n",l);
	      fflush(stdout);*/

	    /* Updates the std of proposal distributions if needed */

	    if(l>0 && (l%nbStepsVar)==0 && EndBurnIn==0)
		{
		    VarPropOK = 0;
		    for (i=0 ; i<2 ; i++)
			{
			    for(j=0 ; j<2 ; j++)
				{
				    if(IsAcceptOKBeta(i,j,acceptOK)==0)
					{
					    updateMCMCSettingsBeta(i, j, nbProp,accept,acceptOK, MCMCSettings);
					}else
					{
					    gsl_matrix_set(isAcceptOKBeta,i,j,1);
					    VarPropOK ++;
					}
				}
			}
		    if(IsAcceptOKBetaWardOut(acceptOK)==0)
			{
			    updateMCMCSettingsBetaWardOut(nbProp,accept,acceptOK, MCMCSettings);
			}else
			{
			    isAcceptOKBetaWardOut=1;
			    VarPropOK ++;
			}
		    if(IsAcceptOKBetaOutOut(acceptOK)==0)
			{
			    updateMCMCSettingsBetaOutOut(nbProp,accept,acceptOK, MCMCSettings);
			}else
			{
			    isAcceptOKBetaOutOut=1;
			    VarPropOK ++;
			}

		    if(IsAcceptOKMu(acceptOK)==0)
			{
			    updateMCMCSettingsMu(nbProp,accept,acceptOK, MCMCSettings);
			}else
			{
			    isAcceptOKMu=1;
			    VarPropOK ++;
			}
		    if(IsAcceptOKSigma(acceptOK)==0)
			{
			    updateMCMCSettingsSigma(nbProp,accept,acceptOK, MCMCSettings);
			}else
			{
			    isAcceptOKSigma=1;
			    VarPropOK ++;
			}

		    if(IsAcceptOKNu1(acceptOK)==0)
			{
			    updateMCMCSettingsNu1(nbProp,accept,acceptOK, MCMCSettings);
			}else
			{
			    isAcceptOKNu1=1;
			    VarPropOK ++;
			}
		    if(IsAcceptOKNu2(acceptOK)==0)
			{
			    updateMCMCSettingsNu2(nbProp,accept,acceptOK, MCMCSettings);
			}else
			{
			    isAcceptOKNu2=1;
			    VarPropOK ++;
			}

		    if(IsAcceptOKTau(acceptOK)==0)
			{
			    updateMCMCSettingsTau(nbProp,accept,acceptOK, MCMCSettings);
			}else
			{
			    isAcceptOKTau=1;
			    VarPropOK ++;
		    	}
		    if(IsAcceptOKAlpha(acceptOK)==0)
			{
			    updateMCMCSettingsAlpha(nbProp,accept,acceptOK, MCMCSettings);
			}else
			{
			    isAcceptOKAlpha=1;
			    VarPropOK ++;
			}

		    if(VarPropOK<NbParamToBeTuned)
			{
			    reInitiateAcceptance(accept);
			    reInitiateNbProp(nbProp);
			    printf("starting preliminary simulation %u \n",l);
			    fflush(stdout);
			}else
			{
			    EndBurnIn++;
			    l=0;
			    printStdProp(MCMCSettings);
			}
		}

	    if(EndBurnIn==1 && (l%nbStepsVar)==0)
		{
		    printf("starting simulation %u \n",l);
		    fflush(stdout);
		}

	    /*** Moves ***/

	    moveAllBeta(MCMCSettings , param, data, nb, augData, accept,nbProp);
	    /* moveBetaOut('w', MCMCSettings, param, data, nb, augData, accept,nbProp); */
	    /* moveBetaOut('o', MCMCSettings, param, data, nb, augData, accept,nbProp); */

	    /* moveSp(param, data, nb, augData, accept); /\*  to remove if we assume Sp = 100% *\/ */
	    /* moveSe(param, data, nb, augData, accept); */

	    /* movePi(param, data, augData, accept); */

	    /* moveDurationColon('m',MCMCSettings, param, augData, accept,nbProp); */
	    /* moveDurationColon('s',MCMCSettings, param, augData, accept,nbProp); */

	    /* moveMutationRate(1, MCMCSettings, param, data, nb, augData, accept,nbProp);*/
	    /* moveMutationRate(2, MCMCSettings, param, data, nb, augData, accept,nbProp); */

	    /* moveTau(MCMCSettings, param, data, nb, augData, accept,nbProp); */
	    /* moveAlpha(MCMCSettings, param, data, nb, augData, accept,nbProp); */

	    gsl_ran_choose (param->rng, AugDataToMove, NbChangeAugData, nb->indexColonisedPatients, nb->NbColonisedPatients, sizeof(int));
	    for(i=0 ; i<NbChangeAugData ; i++)
		{
		    moveC(AugDataToMove[i], MCMCSettings, param, data,nb, augData) ;
		}

	    gsl_ran_choose (param->rng, AugDataToMove, NbChangeAugData, nb->indexColonisedPatients, nb->NbColonisedPatients, sizeof(int));
	    for(i=0 ; i<NbChangeAugData ; i++)
		{
		    moveE(AugDataToMove[i], MCMCSettings, param, data,nb, augData) ;
		}

	    gsl_ran_choose (param->rng, AugDataToMove, NbChangeAugData, nb->indexColonisedPatients, nb->NbColonisedPatients, sizeof(int));
	    for(i=0 ; i<NbChangeAugData ; i++)
		{
		    moveCandE(AugDataToMove[i], MCMCSettings, param, data,nb, augData) ;
		}

	    /* Writing results in an output file */
	    if (l>=MCMCSettings->BurnIn && l%MCMCSettings->SubSample==0 && EndBurnIn==1)
	 	/* if (EndBurnIn==1) */
	 	{
		    writeAllFiles(Files, param, nb, data, augData);
		    /* printf("L = %lg\n",fullLoglikelihoodWithPrior(data, nb, augData, param)); */
	 	}
	}

    /* PRINT ACCESPTANCE RATE */
    printAcceptance(accept,nbProp);

    /* FREE ALLOCATED MEMORY */
    free(AugDataToMove);
    gsl_matrix_free(isAcceptOKBeta);

}
