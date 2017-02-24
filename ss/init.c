/**
 * @file init.c
 * @author Amir M. Abdol
 * @contact a.m.abdol@uva.nl
 * @date Created on May 2015
 *
 * @brief Implementation for initialization of different Sets used by Scatter Search algorithm
 */

#include <string.h>

#include "ss.h"
#include "error.h"


/**
 * @brief      This function should call before almost everything. All the memory
 * allocations for different data structures happen here. Also, many default values
 * is being calculated and set. Procedure:
 * 1. Allocate memory for Reference Set
 * 2. Allocate memory for Scatter Set
 * 3. Allocate memory for Sub Sets
 * 4. Allocate and initialize statistics
 * 5. And more...
 *
 * @param      ssParams  The ss parameters
 */
void init_ssParams(SSType *ssParams){


	ssParams->n_refinement        = 0;
	ssParams->n_ref_set_update    = 0;
	ssParams->n_duplicates        = 0;
	ssParams->n_flatzone_detected = 0;
	ssParams->n_function_evals    = 0;
	ssParams->n_regen			  = 0;
	ssParams->n_duplicate_replaced = 0;
	ssParams->n_iter = 0;

	// Initialize the Reference Set
	ssParams->ref_set             = (Set *)malloc(sizeof(Set));
	allocate_set_memory(ssParams, ssParams->ref_set, ssParams->ref_set_size, ssParams->nreal);
	
	// Initialize the Scatter Set
	ssParams->scatter_set         = (Set *)malloc(sizeof(Set));
	allocate_set_memory(ssParams, ssParams->scatter_set, ssParams->scatter_set_size, ssParams->nreal);
	
	// Estimating the Subsets list size and allocating memory for it
	ssParams->subsets_list_size   = (ssParams->ref_set_size * ssParams->ref_set_size);		// Maximum possible size of subsets_list
	ssParams->subsets_list        = (Set *)malloc( ssParams->subsets_list_size * sizeof(Set) );
	for (int i = 0; i < ssParams->subsets_list_size; ++i){
		/* 	allocating memory for each subset, with the size equal to `pair_size` and lenth of `nreal`. */
		allocate_set_memory(ssParams, &(ssParams->subsets_list[i]), ssParams->pair_size, ssParams->nreal);
	}

	// Declare the candidate set with the maximum size possible!
	ssParams->candidates_set = (Set *)malloc(sizeof(Set));
	allocate_set_memory(ssParams, ssParams->candidates_set, ssParams->ref_set_size * ssParams->ref_set_size * 6, ssParams->nreal);

	ssParams->freqs_matrix        = (int **)malloc(ssParams->nreal * sizeof(int *));
	ssParams->probs_matrix        = (double **)malloc(ssParams->nreal * sizeof(double *));
	
	// Generating the sub regions matrixes
	ssParams->min_boundary_matrix = (double **)malloc( ssParams->nreal * sizeof(double *));
	ssParams->max_boundary_matrix = (double **)malloc( ssParams->nreal * sizeof(double *));

	// Calculating and Initializing sub set regions
	for (int i = 0; i < ssParams->nreal; ++i)
	{
		ssParams->min_boundary_matrix[i] = (double *)malloc((ssParams->p) * sizeof(double));
		ssParams->max_boundary_matrix[i] = (double *)malloc((ssParams->p) * sizeof(double));
		
		ssParams->freqs_matrix[i]        = (int *)malloc(ssParams->p * sizeof(int));
		ssParams->probs_matrix[i]        = (double *)malloc(ssParams->p * sizeof(double));

		for (int j = 1; j <= ssParams->p; ++j)
		{
			/* Building the sub region matrixes */	// One matrix would be enough but for clarity I use two.
			ssParams->min_boundary_matrix[i][j - 1] = ssParams->min_real_var[i] + ((ssParams->max_real_var[i] - ssParams->min_real_var[i]) / ssParams->p) * (j - 1);
			ssParams->max_boundary_matrix[i][j - 1] = ssParams->min_real_var[i] + ((ssParams->max_real_var[i] - ssParams->min_real_var[i]) / ssParams->p) * j; 
			
			/* Building the freqs_matrix: frequently_matrix */
			ssParams->freqs_matrix[i][j - 1]        = 1;
			ssParams->probs_matrix[i][j - 1]        = 1./(ssParams->p);		// Since all values are 1.
		}
	}

#ifdef STATS
	write_int_matrix(ssParams, ssParams->freqs_matrix, ssParams->nreal, ssParams->p, freqs_matrix_file, 0, 'w');
#endif

	// print_double_matrix(ssParams, ssParams->min_boundary_matrix, ssParams->nreal, ssParams->p);
	// print_double_matrix(ssParams, ssParams->max_boundary_matrix, ssParams->nreal, ssParams->p);

}

/**
 * @brief      Initialize report files. These files are being used by
 * during the optimization to keep the history of ref_set, best_sols in the
 * ref_set and frequency of results in each sub-regions of the search space
 * after each iteration.
 *
 * @param      ssParams  The ss parameters
 */
void init_report_files(SSType *ssParams, Files *files) {

#ifdef DEBUG
	printf("Initializing output files\n");
#endif
	char *temp_fname = (char *) calloc(MAX_RECORD + 1, sizeof(char));
	char *mode = "w";
	if (ssParams->perform_warm_start)
		mode = "a+";

#ifdef DEBUG
	sprintf(temp_fname, "%s_ref_history", files->outputfile);
	ref_set_history_file = fopen(temp_fname, mode);
    if( !ref_set_history_file )
        file_error( "fly_X error opening refset history file" );

	sprintf(temp_fname, "%s_best_history", files->outputfile);
	best_sols_history_file = fopen(temp_fname, mode);
    if( !best_sols_history_file )
        file_error( "fly_X error opening best history file" );
#endif

#ifdef STATS
	sprintf(temp_fname, "%s_freqs_history", files->outputfile);
	freqs_matrix_file = fopen(temp_fname, mode);
    if( !freqs_matrix_file )
        file_error( "fly_X error opening freq history file" );
#endif

	sprintf(temp_fname, "%s.log", files->outputfile);
	stats_file = fopen(temp_fname, mode);
    if( !stats_file )
        file_error( "fly_X error opening statistics log file" );
	write_stats_header(stats_file);

	free(temp_fname);
}	

/**
 * @brief      Initialize Scatter Set by randomly generating individuals based on 
 * the prefered discretization parameters `p` provided by users. 
 */
void init_scatter_set(SSType *ssParams, Set *set){

	// The first `p` members are uniformy selected from all sub-regions
	for (int k = 0; k < ssParams->p; ++k)
	{
		for (int i = 0; i < ssParams->nreal; ++i)
		{
			set->members[k].params[i] = rndreal(ssParams->min_boundary_matrix[i][k], ssParams->max_boundary_matrix[i][k]);
		}
	}

	// Indicate the index of selected sub-regions that the new value should be generated from.
	int a            = 0;
	double probs_sum = 0;
	double rnd       = 0;

	/* Generate new item to add to scatter_set
	   set->members[ssParams->p + k] is going to be generated */
	for (int k = 0; k < ssParams->scatter_set_size - ssParams->p; ++k)
	{		
		for (int i = 0; i < ssParams->nreal; ++i)
		{	
			rnd = rndreal(0, 1);
			for (int j = 0; j < ssParams->p; ++j)
			{	
				/* Compute `probs` and find the index */

				probs_sum += ssParams->probs_matrix[i][j];

				if ( rnd <= probs_sum ){
					a = j;
					
					// Updating the f matrix, and prob_matrix
					ssParams->freqs_matrix[i][a]++;
					
					// Updating the ssParams->probs_matrix
					// IMPROVE: This could be rewritten much nicer
					double f_col_prob = 0;
					for (int t = 0; t < ssParams->p; ++t){
						f_col_prob += (1. / ssParams->freqs_matrix[i][t]);
					}
					ssParams->probs_matrix[i][a] = (1. / ssParams->freqs_matrix[i][a] ) / (f_col_prob);

					probs_sum = 0;
					break;
				}
			}

			// sub-region index is selected now: `a`
			set->members[ssParams->p + k].params[i] = rndreal(ssParams->min_boundary_matrix[i][a], ssParams->max_boundary_matrix[i][a] );
		}

	}	// Scatter Set is now extented...

	// #ifdef STATS
	// 	write_int_matrix(ssParams, ssParams->freqs_matrix, ssParams->nreal, ssParams->p, freqs_matrix_file, 0, 'w');
	// #endif
}

/**
 * @brief      Generate the Reference Set using the procedure describes in Banga paper.
 * It assumes that the Reference Set is initialize and the proper memory allocated for it.
 * The procedure start by selected the best `h` members of currectly generated Scatter Set.
 * Then iterate over the rest of the Scatter Set to select the farthest member to the already
 * selected members in each iteration.
 *
 */
void init_ref_set(SSType *ssParams){

	printf("Forming the refSet...\n");

	int m = ssParams->scatter_set_size;
	int b = ssParams->ref_set_size;
	int h = b / 2;	// max_elite

	// Sorting the scatter_set
	quick_sort_set(ssParams, ssParams->scatter_set, ssParams->scatter_set_size);

	/* Copying first `h` best items of the scatter_set with respect to the cost */
	for (int i = 0; i < h; ++i){
		copy_ind(ssParams, &(ssParams->ref_set->members[i]), &(ssParams->scatter_set->members[i]));
	}

	/************************************************************/
	/* Expanding the ref_set by applying the diversity distance */
	/************************************************************/
	
	// dist_matrix will store the distances between the rest of scatter_set and the ref_set
	double **dist_matrix = (double **)malloc( (m - h) * sizeof(double *));;
	for (int i = 0; i < (m - h); ++i){
		dist_matrix[i] = (double *)malloc( (b) * sizeof(double));
	}

	double *min_dists = (double *)malloc((m - h) * sizeof(double));

	int min_index = 0;
	int max_index = 0;

	int i = 0;
	for (int k = h; k < b; ++k, --m)		// Iterate through the ref_set for expanding
	{										// --m reduces the depth that the second `for` should go through it
											// since the array is being shifted to the left after each iteration

		for (i = h; i < m; ++i)				// Iterate through the rest of scatter_set
		{							 

			for (int j = 0; j < k; ++j)		// Iterate trhough the selected members of ref_set
			{								
				/* Computing the distance */
				dist_matrix[i - (b / 2)][j] = euclidean_distance(ssParams, &(ssParams->scatter_set->members[i]), &(ssParams->ref_set->members[j]));
			}

			// Find the minimum of each row into a new variable, then finding it's max
			min_dists[i - (b / 2)] = min(dist_matrix[i - (b / 2)], k, &min_index);
		}
		
		// Now, the min_dists of vectors to all the members of R is computed,
		// we could compute the maximum of them and also it's index
		
		// Finding the place of the max/min distance
		max(min_dists, i - (b / 2), &max_index); 

		// Adding the vector in max_index to the ref_set
		copy_ind(ssParams, &(ssParams->ref_set->members[k]), &(ssParams->scatter_set->members[h + max_index]));

		// After assigment, the rest of should be shift into the left form the position that the 
		// vector is selected
		
		// It could be done more complicated!
		delete_and_shift(ssParams, ssParams->scatter_set, ssParams->scatter_set_size, (b / 2) + max_index);

	}

	#ifdef STATS
		write_int_matrix(ssParams, ssParams->freqs_matrix, ssParams->nreal, ssParams->p, freqs_matrix_file, 0, 'w');
	#endif

	free(min_dists);
	for (int i = 0; i < (ssParams->scatter_set_size - b/2); ++i){
		free(dist_matrix[i]);
	}
	free(dist_matrix);
}

/**
 * @brief      Performing the warmStart by initializing the Reference Set from `ssParams->freq_mat_final_filename`, 
 * `ssParams->prob_mat_final_filename`, `ssParams->ref_set_final_filename`.
 * 
 * @note 	   Files should exist and have the same number of rows (ref_set_size) and columns (nreal).
 *
 */
void init_warm_start(SSType *ssParams){
	int i;
    char line[4098];
    printf("Loading the data to perform warm start...\n");
	
	// Read refSet
    i = 0;
    FILE* refSetStream = fopen(ssParams->ref_set_final_filename, "r");
    while (fgets(line, 4098, refSetStream) && (i < ssParams->ref_set_size))
    {
    	double row[ssParams->nreal + 1];
        char* tmp = strdup(line);
        parse_double_row(ssParams, tmp, row);
        for (int j = 0; j < ssParams->nreal; ++j){
        	ssParams->ref_set->members[i].params[j] = row[j];
        }
        ssParams->ref_set->members[i].cost = row[ssParams->nreal];
        free(tmp);
        i++;
    }
    ssParams->best = (individual *)malloc(sizeof(individual));
	ssParams->best = &(ssParams->ref_set->members[0]);				// The first members of ref_set is always the best
    print_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, ssParams->nreal);
    print_ind(ssParams, ssParams->best, ssParams->nreal);

	// Read freqMat
    i = 0;
    FILE* freqMatStream = fopen(ssParams->freq_mat_final_filename, "r");
    while (fgets(line, 4098, freqMatStream) && (i < ssParams->nreal ))
    {
    	int row[ssParams->p];
        char* tmp = strdup(line);
        parse_int_row(ssParams, tmp, row);
        memcpy(ssParams->freqs_matrix[i], (int*)row, ssParams->p * sizeof(int));
        free(tmp);
        i++;
    }
    // print_int_matrix(ssParams, ssParams->freqs_matrix, ssParams->nreal, ssParams->p);

	// Read probMat
    i = 0;
    FILE* probMatStream = fopen(ssParams->prob_mat_final_filename, "r");
    while (fgets(line, 4098, probMatStream))
    {
    	double row[ssParams->p];
        char* tmp = strdup(line);
        parse_double_row(ssParams, tmp, row);
        memcpy(ssParams->probs_matrix[i], row, ssParams->p * sizeof(double));
        free(tmp);
        i++;
    }
    // print_double_matrix(ssParams, ssParams->probs_matrix, ssParams->nreal, ssParams->p);

}
