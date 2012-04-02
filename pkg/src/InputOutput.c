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
/* Reading data                                                               */
/******************************************************************************/


void readFakeNbData(nb_data *nb){
    int i;
    FILE *fich;
    int V;

    for(i=0 ; i<NbPatients ; i++)
	{
	    nb->NbAdmissions[i]=1;
	}

    fich = fopen("nbNeg.txt","r");
    if ( fich == NULL )
	{
	    printf("A problem occurred while opening nbNeg.txt. Check that the file exists and is not opened.\n");
	    fflush(stdout);
	    exit(1);
	}

    for (i=0 ; i<NbPatients ; i++)
	{
	    if( fscanf(fich,"%d",&V) != 1)
		{
		    printf("A problem occurred while reading the file\n");
		    fflush(stdout);
		    break;
		}
	    nb->NbNegSwabs[i]=V;
	}

    fich = fopen("nbPos.txt","r");
    if ( fich == NULL )
	{
	    printf("A problem occurred while opening nbPos.txt. Check that the file exists and is not opened.\n");
	    fflush(stdout);
	    exit(1);
	}

    for (i=0 ; i<NbPatients ; i++)
	{
	    if( fscanf(fich,"%d",&V) != 1)
		{
		    printf("A problem occurred while reading the file\n");
		    fflush(stdout);
		    break;
		}
	    nb->NbPosSwabs[i]=V;
	}

}





void readFakeData(nb_data *nb, raw_data *data){
    FILE *fich;
    int i,k;
    int V;

    fich = fopen("Admission.txt","r");
    if ( fich == NULL )
	{
	    printf("A problem occurred while opening Admission.txt. Check that the file exists and is not opened.\n");
	    fflush(stdout);
	    exit(1);
	}

    for (i=0 ; i<NbPatients ; i++)
	{
	    for (k=0 ; k<nb->NbAdmissions[i] ; k++)
		{
		    if( fscanf(fich,"%d",&V) != 1)
			{
			    printf("A problem occurred while reading the file\n");
			    fflush(stdout);
			    break;
			}
		    gsl_vector_set(data->A[i],k,V);
		}
	}
    fclose(fich);

    fich = fopen("Discharge.txt","r");
    if ( fich == NULL )
	{
	    printf("A problem occurred while opening Discharge.txt. Check that the file exists and is not opened.\n");
	    fflush(stdout);
	    exit(1);
	}

    for (i=0 ; i<NbPatients ; i++)
	{
	    for (k=0 ; k<nb->NbAdmissions[i] ; k++)
		{
		    if( fscanf(fich,"%d",&V) != 1)
			{
			    printf("A problem occurred while reading the file\n");
			    fflush(stdout);
			    break;
			}
		    gsl_vector_set(data->D[i],k,V);
		}
	}
    fclose(fich);

    fich = fopen("PositiveSwabDates.txt","r");
    if ( fich == NULL )
	{
	    printf("A problem occurred while opening PositiveSwabDates.txt. Check that the file exists and is not opened.\n");
	    fflush(stdout);
	    exit(1);
	}

    for (i=0 ; i<NbPatients ; i++)
	{
	    for (k=0 ; k<nb->NbPosSwabs[i] ; k++)
		{
		    if( fscanf(fich,"%d",&V) != 1)
			{
			    printf("A problem occurred while reading the file\n");
			    fflush(stdout);
			    break;
			}
		    gsl_vector_set(data->P[i],k,V);
		}
	}
    fclose(fich);

    fich = fopen("NegativeSwabDates.txt","r");
    if ( fich == NULL )
	{
	    printf("A problem occurred while opening NegativeSwabDates.txt. Check that the file exists and is not opened.\n");
	    fflush(stdout);
	    exit(1);
	}

    for (i=0 ; i<NbPatients ; i++)
	{
	    for (k=0 ; k<nb->NbNegSwabs[i] ; k++)
		{
		    if( fscanf(fich,"%d",&V) != 1)
			{
			    printf("A problem occurred while reading the file\n");
			    fflush(stdout);
			    break;
			}
		    gsl_vector_set(data->N[i],k,V);
		}
	}
    fclose(fich);

    fich = fopen("Ward.txt","r");
    if ( fich == NULL )
	{
	    printf("A problem occurred while opening Ward.txt. Check that the file exists and is not opened.\n");
	    fflush(stdout);
	    exit(1);
	}

    for (i=0 ; i<NbPatients ; i++)
	{
	    if( fscanf(fich,"%d",&V) != 1)
		{
		    printf("A problem occurred while reading the file\n");
		    fflush(stdout);
		    break;
		}
	    data->ward[i]=V;
	}
    fclose(fich);

    for (i=0 ; i<NbPatients ; i++)
	{
	    data->PatientIndex[i]=i;
	}

    CalculIsInHosp(nb, data);

}







/******************************************************************************/
/* Function to read initial parameters values                                 */
/******************************************************************************/

/*************************************************************************/
/* preparing output file                                                 */
/*************************************************************************/

void prepAllFiles(output_files * Files){
    int i;
	
    fprintf(Files->LogL,"LogLikelihood\n");
    fflush(Files->LogL);

    fprintf(Files->Parameters,"beta[0,0]\t");
    fprintf(Files->Parameters,"beta[0,1]\t");
    fprintf(Files->Parameters,"beta[1,0]\t");
    fprintf(Files->Parameters,"beta[1,1]\t");
    fprintf(Files->Parameters,"betaWardOut\t");
    fprintf(Files->Parameters,"betaOutOut\t");
    /* fprintf(Files->Parameters,"Sp\t"); */
    fprintf(Files->Parameters,"Se\t");
    fprintf(Files->Parameters,"Pi\t");
    fprintf(Files->Parameters,"mu\t");
    fprintf(Files->Parameters,"sigma\t");
    fprintf(Files->Parameters,"nu1\t");
    fprintf(Files->Parameters,"nu2\t");
    fprintf(Files->Parameters,"tau\t");
    fprintf(Files->Parameters,"alpha\n");
    fflush(Files->Parameters);

    for(i=0 ; i<NbPatients ; i++)
	{
	    fprintf(Files->ColonDates,"C[%d]\t",i);
	    fprintf(Files->EndColonDates,"E[%d]\t",i);
	}
    fprintf(Files->ColonDates,"\n");
    fprintf(Files->EndColonDates,"\n");
    fflush(Files->ColonDates);
    fflush(Files->EndColonDates);

}







/*************************************************************************/
/* writing output files                                                  */
/*************************************************************************/

void writeAllFiles(output_files * Files, parameters * param, nb_data *nb, raw_data * data, aug_data *augData){
    /* Writing results in the output files */

    double L;
    int i;

    L = fullLoglikelihoodWithPrior(data, nb, augData, param);
    fprintf(Files->LogL,"%lf\n",L);
    fflush(Files->LogL);

    fprintf(Files->Parameters,"%lf\t",gsl_matrix_get(param->beta,0,0));
    fprintf(Files->Parameters,"%lf\t",gsl_matrix_get(param->beta,0,1));
    fprintf(Files->Parameters,"%lf\t",gsl_matrix_get(param->beta,1,0));
    fprintf(Files->Parameters,"%lf\t",gsl_matrix_get(param->beta,1,1));
    fprintf(Files->Parameters,"%lf\t",param->betaWardOut);
    fprintf(Files->Parameters,"%lf\t",param->betaOutOut);
    /* fprintf(Files->Parameters,"%lf\t",param->Sp); */
    fprintf(Files->Parameters,"%lf\t",param->Se);
    fprintf(Files->Parameters,"%lf\t",param->Pi);
    fprintf(Files->Parameters,"%lf\t",param->mu);
    fprintf(Files->Parameters,"%lf\t",param->sigma);
    fprintf(Files->Parameters,"%lf\t",param->nu1);
    fprintf(Files->Parameters,"%lf\t",param->nu2);
    fprintf(Files->Parameters,"%lf\t",param->tau);
    fprintf(Files->Parameters,"%lf\n",param->alpha);
    fflush(Files->Parameters);

    for(i=0 ; i<NbPatients ; i++)
	{
	    if(nb->NbPosSwabs[i]>0)
		{
		    fprintf(Files->ColonDates,"%d\t",augData->C[i]);
		    fprintf(Files->EndColonDates,"%d\t",augData->E[i]);
		}else
		{
		    fprintf(Files->ColonDates,"NA\t");
		    fprintf(Files->EndColonDates,"NA\t");
		}
	}
    fprintf(Files->ColonDates,"\n");
    fprintf(Files->EndColonDates,"\n");
    fflush(Files->ColonDates);
    fflush(Files->EndColonDates);
	
}

