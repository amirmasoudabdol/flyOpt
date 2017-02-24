/**
 * @file ss.c
 * @author Amir M. Abdol
 * @contact a.m.abdol@uva.nl
 * @date Created on May 2015
 *
 * @brief Main optimization routines are implemented in this file.
 */

#include "ss.h"

FILE *ref_set_history_file;
FILE *best_sols_history_file;
FILE *freqs_matrix_file;
FILE *freq_mat_final_file;
FILE *prob_mat_final_file;
FILE *ref_set_final_file;
FILE *stats_file;

ScoreOutput out;


/**
 * @brief      Initialize the Scatter Search by calling several functions for allocation,
 * randomization and reporting. Check the code for more detail.
 *
 * @param      inp       A pointer to Input struct for linking Scatter Search routines and fly simulator.
 * @param      ssParams  A pointer to SSType struct used to store all necessary variables needed by 
 * Scatter Search.
 * @param      inname    The name of the input file.
 */
void InitSS(Input *inp, SSType *ssParams, Files *files){


	// Initializing the ScoreOutput struct.
    out.score          = 1e38;           // start with a very large number
    out.penalty        = 0;
    out.size_resid_arr = 0;
    out.jacobian       = NULL;
    out.residuals      = NULL;

	printf("\nInitializing Scatter Search...\n");

	// Initializing Mersenne Twister random number generator based on the 
	// seed value provided.
	InitRand(ssParams->seed);

	// Initializing the history files.
	init_report_files(ssParams, files);

	// Allocate memory of ssParams variables, and initialize some parameters
	init_ssParams(ssParams);

	if ( !ssParams->perform_warm_start ){

		init_scatter_set(ssParams, ssParams->scatter_set);
		evaluate_set(ssParams, ssParams->scatter_set, ssParams->scatter_set_size, inp, &out);

		init_ref_set(ssParams);
		quick_sort_set(ssParams, ssParams->ref_set, ssParams->ref_set_size);

		// Initialize the best solution
		ssParams->best = (individual *)malloc(sizeof(individual));

		// Assigning the best solutions
		ssParams->best = &(ssParams->ref_set->members[0]);
	}
	else{
		// Assume that 'final_ref_set_filename' point to a TSV file.
		// Rows count corresponds to the ref_set_size and columns count
		// refers to number of parameters to be tweaked.
		init_warm_start(ssParams);
	}	

#ifdef DEBUG
	write_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, ssParams->nreal, ref_set_history_file, 0, 'w');
	write_ind(ssParams, ssParams->best, ssParams->nreal, best_sols_history_file, 0, 'w');
#endif
}


/**
 * @brief      Main loop of Scatter Search. Basically, running the scatter search, 
 * computing stats, and writing the results into files.
 *
 */
void RunSS(Input *inp, SSType *ssParams, Files *files){

	int n_ref_set_update    = 0;
	int n_refinement        = 0;
	int n_duplicates        = 0;
	int n_function_evals    = 0;
	int n_flatzone_detected = 0;
	// bool wasChanged = false;

	printf("Starting the optimization procedure...\n");
	for (ssParams->n_iter = 1; ssParams->n_iter < ssParams->max_iter; 
		++ssParams->n_iter)
	{
		// Selecting the SubSets List
		select_subsets_list(ssParams, ssParams->ref_set, ssParams->ref_set_size);

		// Generate new candidates
		generate_candiates(ssParams);
		evaluate_set(ssParams, ssParams->candidates_set, ssParams->candidates_set_size, inp, &out);

		// Update refSet by replacing new candidates
		update_ref_set(ssParams);

		// Perform the local_search
		if (ssParams->perform_local_search && (ssParams->n_iter % ssParams->local_search_freq == 0)  ) {
			/* Nelder Mead by default */
			refine_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, 'n', inp, &out);
		}
		quick_sort_set(ssParams, ssParams->ref_set, ssParams->ref_set_size);
		
#ifdef DEBUG
		// Append the ref_set to the file
		write_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, 
			ssParams->nreal, ref_set_history_file, ssParams->n_iter, 'w');
		fflush(ref_set_history_file);

		// Append the best solution to the sol history file
		write_ind(ssParams, ssParams->best, ssParams->nreal, 
			best_sols_history_file, ssParams->n_iter, 'w');
		fflush(best_sols_history_file);
#endif

		// It basically check if all refSet members are in the same flat zone.
		if (ssParams->perform_stop_criteria)
			if (fabs(ssParams->ref_set->members[0].cost - ssParams->ref_set->members[ssParams->ref_set_size - 1].cost) < ssParams->stop_criteria){
				printf("\n%s   Stop by difference criteria!\n   The difference between the best and worst members of refSet is smaller than %lf\n\n%s", KRED, ssParams->stop_criteria, KNRM);
				break;
			}

		/*
		 * If user wants regeneration of the refSet and we're not any more in
		 * the first iteration, check if the number of duplicates is too much
		 * or it was time to do a periodic regeneration.
		 */
		if ((ssParams->perform_ref_set_regen && ssParams->n_iter != 1) &&
			((double)(ssParams->n_duplicates - n_duplicates) / 
				(double)ssParams->candidates_set_size > 0.7	|| 
				ssParams->n_iter % ssParams->ref_set_regen_freq == 0))
			 // || (double)(ssParams->n_flatzone_detected - n_flatzone_detected) / (double)ssParams->candidates_set_size > 0.7
		{
			re_gen_ref_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, 'n', inp, &out);
			ssParams->n_regen++;

			quick_sort_set(ssParams, ssParams->ref_set, ssParams->ref_set_size);
		}

#ifdef STATS
		write_int_matrix(ssParams, ssParams->freqs_matrix, ssParams->nreal, ssParams->p, freqs_matrix_file, ssParams->n_iter, 'w');
#endif		

		/* output some basic statistics to log file */
		if (ssParams->n_iter % 10 == 0) {
			fprintf(stats_file, "%d\t%d\t%g\t%g\t%g\t%d\t%d\t%d\t%d\t%d\n", 
				ssParams->n_iter, 
				ssParams->n_function_evals/* -  n_function_evals*/,
				ssParams->best->cost,
				average_cost_refset(ssParams, ssParams->ref_set, ssParams->ref_set_size),
				var_cost_refset(ssParams, ssParams->ref_set, ssParams->ref_set_size),
				ssParams->n_ref_set_update -  n_ref_set_update, 
				ssParams->n_refinement -  n_refinement, 
				ssParams->n_flatzone_detected -  n_flatzone_detected, 
				ssParams->n_duplicates -  n_duplicates, 
				ssParams->candidates_set_size);

#ifdef DEBUG
			/* in debug mode we want it all immediately on disk */
			fflush(stats_file);

			/* and write it to terminal */
			printf("\nStats - (%d):\n", ssParams->n_iter);
			// print_ind(ssParams, ssParams->best, ssParams->nreal);	
			printf("\t\tBest cost: %lf\n", ssParams->best->cost);
			printf("\t\t   as RMS: %lf\n", sqrt(ssParams->best->cost / ( double ) inp->zyg.ndp));
			printf("\t\t# Replacement in Reference Set: %d\n", ssParams->n_ref_set_update - n_ref_set_update);
			printf("\t\t# of Local Search Performed: %d\n", ssParams->n_refinement - n_refinement);
			printf("\t\t# Duplicates: %d\n", ssParams->n_duplicates - n_duplicates);
			printf("\t\t# Flatzone: %d\n", ssParams->n_flatzone_detected - n_flatzone_detected);
			printf("\t\t================= candidateSetSize: %d\n", ssParams->candidates_set_size);
#endif

			/* Reset delta counters after we've reported their values */
			n_ref_set_update    = ssParams->n_ref_set_update;
			n_refinement        = ssParams->n_refinement;
			n_duplicates        = ssParams->n_duplicates;
			n_function_evals    = ssParams->n_function_evals;
			n_flatzone_detected = ssParams->n_flatzone_detected;
		}
	} // End of the main loop

	/* AC: We do a last local search on the refset */
	refine_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, 'n', inp, &out);
	quick_sort_set(ssParams, ssParams->ref_set, ssParams->ref_set_size);

	/* Final output of basic statistics to log file */
	fprintf(stats_file, "%d\t%d\t%g\t%g\t%g\t%d\t%d\t%d\t%d\t%d\n",
		ssParams->n_iter, 
		ssParams->n_function_evals/* -  n_function_evals*/,
		ssParams->best->cost,
		average_cost_refset(ssParams, ssParams->ref_set, ssParams->ref_set_size),
		var_cost_refset(ssParams, ssParams->ref_set, ssParams->ref_set_size),
		ssParams->n_ref_set_update -  n_ref_set_update, 
		ssParams->n_refinement -  n_refinement, 
		ssParams->n_flatzone_detected -  n_flatzone_detected, 
		ssParams->n_duplicates -  n_duplicates, 
		ssParams->candidates_set_size);
	/* Mark end-of-file */
	fprintf(stats_file, "#eof\n");

	/* Generating final output to terminal */
	printf("\nReference Set:\n");
	print_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, ssParams->nreal);
	printf("\n====================================\n");
	printf("Best Solution:\n");
	print_ind(ssParams, ssParams->best, ssParams->nreal);
	printf("====================================\n");

	// // Writing stats and final refSet to file for possible warm_start
	// char *str_filename;
	// str_filename = (char *) calloc(MAX_RECORD + 1, sizeof(char));

	// sprintf(str_filename, "%s/%s", ssParams->output_path, ssParams->freq_mat_final_filename);
	// freq_mat_final_file = fopen(str_filename, "w");

	// sprintf(str_filename, "%s/%s", ssParams->output_path, ssParams->prob_mat_final_filename);
	// prob_mat_final_file = fopen(str_filename, "w");

	// sprintf(str_filename, "%s/%s", ssParams->output_path, ssParams->ref_set_final_filename);
	// ref_set_final_file = fopen(str_filename, "w");
	// free(str_filename);

	// printf("\nExporting %s, %s and %s for warm start...\n", ssParams->freq_mat_final_filename, ssParams->prob_mat_final_filename, ssParams->ref_set_final_filename);
	// write_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, ssParams->nreal, ref_set_final_file, -1, 'w');
	// write_int_matrix(ssParams, ssParams->freqs_matrix, ssParams->nreal, ssParams->p, freq_mat_final_file, -1, 'w');
	// write_double_matrix(ssParams, ssParams->probs_matrix, ssParams->nreal, ssParams->p, prob_mat_final_file, -1, 'w');

	/* Write refset as output configuration files */
	write_refset_eqparms(ssParams, files, inp);
	deallocate_ssParam(ssParams);

#ifdef DEBUG
	fclose(ref_set_history_file);
	fclose(best_sols_history_file);
#endif
#ifdef STATS
	fclose(freqs_matrix_file);
#endif
	fclose(stats_file);
}
