/**
 * @file essIO.c
 * @author Amir M. Abdol
 * @contact a.m.abdol@uva.nl
 * @date Created on Sep 2015
 *
 * @brief Implementation of IO routines
 */

#include "ess.h"
#include "fly_io.h"
#include <unistd.h>
#include <ctype.h>

void read_cli_params(eSSType *eSSParams, int argc, char **argv){

   // char *cvalue = NULL;
   // int index;
   int c;
   printf("Reading the command line paremters...\n");
		while ((c = getopt (argc, argv, "m:drwuslo:")) != -1)
	     switch (c)
	       {
	       case 'm':		// maxiter
	         eSSParams->maxiter = atoi(optarg);
	         break;
	       case 'd':
	         eSSParams->debug = 1;
	         break;
	       case 'u':
	         eSSParams->user_guesses = 1;
	         break;
	       case 'w':
	         eSSParams->warmStart = 1;
	         break;
	       case 's':
	         eSSParams->collectStats = 1;
	         break;
	       case 'r':
	         eSSParams->saveOutput = 1;
	         break;
	       case 'l':
	         eSSParams->perform_LocalSearch = 1;
	         break;
	       case 'o':
	         	eSSParams->local_method = optarg[0];
	         break;
	       case '?':
	       	  break;
	       default:
	         abort ();
	       }

}

void print_Set(eSSType *eSSParams, Set *set, int onlyCost){
	printf("-----------------------------------\n");
	for (int i = 0; i < set->size; ++i)
	{
		printf("%d: ", i); print_Ind(eSSParams, &(set->members[i]), onlyCost);
	}
	printf("\n");

}

void print_Ind(eSSType *eSSParams, individual *ind, int onlyCost){

	if (!onlyCost){
		for (int i = 0; i < eSSParams->n_Params; ++i)
		{
			printf("%.5lf, ", ind->params[i]);
		}
	}
	// printf("\t (cost: %lf)\t(dist: %lf)\n", ind->cost, ind->dist);
	printf("\t (cost: %lf)\n", ind->cost/*, ind->dist*/);

}


void write_Set(eSSType *eSSParams, Set *set, FILE *fpt, int iter){

	for (int i = 0; i < set->size; ++i)
	{
		write_Ind(eSSParams, &(set->members[i]), fpt, iter);
	}
}

void write_Ind(eSSType *eSSParams, individual *ind, FILE *fpt, int iter){

	if (iter != -1)
		fprintf(fpt, "%d\t", iter);

	for (int i = 0; i < eSSParams->n_Params; ++i)
	{
		fprintf(fpt, "%.5lf\t", ind->params[i]);
	}
	
	fprintf(fpt, "%lf\n", ind->cost);
}


void print_Stats(eSSType *eSSParams){

	printf("%s\n", KGRN);
	printf("Overall Statistics:\n");
	printf("\tn_iter: %d\n", eSSParams->iter);
	printf("\tn_successful_goBeyond: %d\n", eSSParams->stats->n_successful_goBeyond);
	printf("\tn_local_search_performed: %d\n", eSSParams->stats->n_local_search_performed);
	printf("\tn_successful_localSearch: %d\n", eSSParams->stats->n_successful_localSearch);
	printf("\tn_local_search_iterations: %d \t (avg: %d)\n", eSSParams->stats->n_local_search_iterations, eSSParams->stats->n_local_search_iterations / (eSSParams->stats->n_successful_localSearch + 1) );
	printf("\tn_refSet_randomized: %d\n", eSSParams->stats->n_refSet_randomized);
	printf("\tn_Stuck: %d\n", eSSParams->stats->n_Stuck);
	printf("\tn_successful_recombination: %d\n", eSSParams->stats->n_successful_recombination);
	if(eSSParams->compute_Set_Stats){
		printf("\tRefSet Mean Cost: %lf+/-%lf\n", eSSParams->refSet->mean_cost, eSSParams->refSet->std_cost);
	}
	printf("-----------------------------------------\n");
	printf("%s\n", KNRM);

}

void write_Stats(eSSType *eSSParams, FILE *fpt){

	fprintf(fpt, "%d\t", eSSParams->iter);
	fprintf(fpt, "%d\t", eSSParams->stats->n_successful_goBeyond);
	fprintf(fpt, "%d\t", eSSParams->stats->n_local_search_performed);
	fprintf(fpt, "%d\t", eSSParams->stats->n_successful_localSearch);
	fprintf(fpt, "%d\t", eSSParams->stats->n_local_search_iterations);
	fprintf(fpt, "%d\t", eSSParams->stats->n_Stuck);
	fprintf(fpt, "%d\t", eSSParams->stats->n_successful_recombination);
	fprintf(fpt, "\n");

}


void parse_double_row(eSSType *eSSParams, char *line, double *row){

    int i = 0;
    const char* tok;
    for (tok = strtok(line, "\t"); tok && *tok; i++, tok = strtok(NULL, "\t\n"))
    {
        row[i] = atof(tok);
    }
}


void parse_int_row(eSSType *eSSParams, char *line, int *row){

    int i = 0;
    const char* tok;
    for (tok = strtok(line, "\t"); tok && *tok; i++, tok = strtok(NULL, "\t\n"))
    {
        row[i] = atoi(tok);
    }
}

void print_Inputs(eSSType *eSSParams){
	printf("%s\n", KCYN);
	printf("\nInput Parameters:\n");
	printf("\tMaximum Iterations: % d\n", eSSParams->maxiter);
	printf("\tDebug: % d\n", eSSParams->debug);
	printf("\tWarm Start: % d\n", eSSParams->warmStart);
	printf("\t# of Sub Regions: % d\n", eSSParams->n_subRegions);
	printf("\t# of Parameters: % d\n", eSSParams->n_Params);
	printf("\tReference Set Size: % d\n", eSSParams->n_refSet);
	printf("\tCandidate Set Size: % d\n", eSSParams->n_candidateSet);
	printf("\tChildren Set Size: % d\n", eSSParams->n_childsSet);
	printf("\tStuck Tolerance: % d\n", eSSParams->maxStuck);
	printf("\tLocal Search Activated: %s\n", eSSParams->perform_LocalSearch == 1 ? "Yes" : "NO");
	printf("\tLocal Search Method: %s\n", eSSParams->local_method == 'l' ? "Levenberg-Marquardt" : "Nelder-Mead");
	printf("\tLocal Search Tolerance: %e\n", eSSParams->local_Tol);
	printf("\tLocal Search Max Iters: %d\n", eSSParams->local_maxIter);
	printf("\tLocal Search only on Best Sol: %s\n", eSSParams->local_onBest_Only == 1 ? "True" : "False");
	printf("------------------------------------------\n");
	printf("%s\n", KNRM);
}


void write_params_to_fly_output_standard(eSSType *eSSParams, Input *inp, char *inname){
    // int i, j;
    ParamList *ptab;
    
    printf("**********************\n%s\n", inname);
    
    printf("\nCreating output files for each parameters set.\n");
    
    ptab = inp->tra.array;      // Get the pointer to the 'inp' array.
    
    // printf("%d\n", eSSParams->i_archivesize);
    // printf("%d\n", eSSParams->i_hardl);
    for(int i=0;i<eSSParams->n_refSet;i++)
    {
        // printf("%d\n", eSSParams->i_totalno_var);
        for(int h=0;h<eSSParams->n_Params;h++){
        	*( ptab[h].param  ) = eSSParams->refSet->members[i].params[h];
            // printf("%lf,", eSSParams->d_archive[i][h]);
         //    printf("---%lf", eSSParams->d_func_archive[i][h]);
        	// printf("---%lf\t", eSSParams->d_solution[i][h]);

            // fprintf(fp, "%f\t", eSSParams->d_func_archive[i][h]);
        }
        // printf("\n");
        WriteParameters(inname, &(inp->zyg.parm), "eqparms", 9, inp->zyg.defs);
        
    }//End of for loop
    
    
    // for (i = 0; i < eSSParams->popsize; ++i){
    //     /* code */
    //     for (j = 0; j < inp->tra.size; ++j){
    //         /* Modifying the local parameters of 'inp' struct */
    //         *( ptab[j].param  ) = pop->ind[i].xreal[j];
    //     }
    //     // Each run will create an output file in form of 'inname_parm_XXXXXX' in which 'XXXXXXX'
    //     // will replace by random string.
    //     WriteParameters(inname, &(inp->zyg.parm), "eqparms", 9, inp->zyg.defs);
    // }
    printf("Done.\n");
}


