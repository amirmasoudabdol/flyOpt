/**
 * @file 
 * @author Amir M. Abdol
 * @contact a.m.abdol@uva.nl
 * @date Created on May 2015
 *
 * @brief Implementation of update and regeneration of Reference Set.
 */

#include "ss.h"
#include <string.h>

/**
 * @brief      The refinement procedure, based on the routine described in 
 * Egea et al. 2006. It only replaces the refSet members with better 
 * candidates by considering the fitness and the diversity (equality) to 
 * avoid duplicates.
 */
void update_ref_set(SSType *ssParams){

	/* Sort the candidate set */
	quick_sort_set(ssParams, ssParams->candidates_set, 
		ssParams->candidates_set_size);
	
	/* Index to members in the candidate and reference set */
	int i = 0;
	
	/* Replace the best ref set member if the best candidate is better */
	if (ssParams->candidates_set->members[i].cost < ssParams->ref_set->members[i].cost ){
		copy_ind(ssParams, &(ssParams->ref_set->members[i]), &(ssParams->candidates_set->members[i]));
		i++;
	}
	
	/* Continue until there is no better candidate. In each iteration, 
	 * compare the best candidate with the worst member of Reference Set. If 
	 * the candidate is better, then replace the worst member of Reference Set 
	 * with the candidate. 
     *
     * After each replacement, an Insertion Sort on the Reference Set keeps it 
     * ordered by cost.
     */
	while (ssParams->candidates_set->members[i].cost < ssParams->ref_set->members[ ssParams->ref_set_size - 1].cost )
	{
		/* Check if candidate is different enough from ref set members */
		int duplicate_index = is_exist(ssParams, ssParams->ref_set, ssParams->ref_set_size, &(ssParams->candidates_set->members[i]));
		if (duplicate_index == -1)
		{
			/* 
			 * Check if ref set members with different parameters have similar 
			 * cost.
			 */
			if ( ssParams->perform_flatzone_detection )
			{
			   	if ( !is_in_flatzone(ssParams, ssParams->ref_set , ssParams->ref_set_size, &(ssParams->candidates_set->members[i])) )
			   	{
					replace(ssParams, &(ssParams->ref_set->members[ssParams->ref_set_size - 1]), &(ssParams->candidates_set->members[i]));
			   	}
			}
			else
			{
				/* We don't care about flatzones, just replace */
				replace(ssParams, &(ssParams->ref_set->members[ssParams->ref_set_size - 1]), &(ssParams->candidates_set->members[i]));
			}
		}
		else
		{	
			/* 
			 * Candidate has a duplicate in the ref set, only replace it if 
			 * the candidate performs better.
			 */
			ssParams->n_duplicates++;
			if (ssParams->candidates_set->members[i].cost < ssParams->ref_set->members[duplicate_index].cost){
				replace(ssParams, &(ssParams->ref_set->members[duplicate_index]), &(ssParams->candidates_set->members[i]));
				ssParams->n_duplicate_replaced++;
			}
		}
		i++;
	}
}

/**
 * @brief      Replace `dest` with `src` in the Reference Set
 *
 * @param      ssParams         The ss parameters
 * @param      dest             The destination
 * @param      src              The source
 */
void replace(SSType *ssParams, individual *dest, individual *src){

	copy_ind(ssParams, dest, src);
	insertion_sort(ssParams, ssParams->ref_set, ssParams->ref_set_size);
	ssParams->n_ref_set_update++;

#ifdef STATS
	update_frequency_matrix(ssParams, dest);
#endif
}

/**
 * @brief      Regenerate the Reference Set. It keeps the first `g` members of the Reference Set (elite Members)
 * and it regenerate `b - g` new members.
 * 
 * @note 	   For more detail, check page 16 of [Novel metaheuristic for parameter estimation in nonlinear dynamic biological systems](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1654195/)
 *
 */
void re_gen_ref_set(SSType *ssParams, Set *set, int set_size, char type, Input *inp, ScoreOutput *out){
	int n = ssParams->nreal;			// Number of paramaters
	int b = ssParams->ref_set_size;		// Reference Set Size
	int g = ssParams->max_elite;		// Number of elite members in the Reference Set
	int h = b - g;						// Number of members to regenerated

	// refSet:
	//////////// h /////////////
	////////////////////////////
	//	h = b - g 	|		g //
	////////////////////////////

	// Generate a new Scatter Search
	int m = ssParams->scatter_set_size;
	init_scatter_set(ssParams, ssParams->scatter_set);
	
	double **M = (double **)malloc(n * sizeof(double*));
	for (int i = 0; i < n; ++i){
		M[i] = (double *)malloc( (b) * sizeof(double));
	}

	double **P = (double **)malloc(m * sizeof(double *));
	for (int i = 0; i < m; ++i){
		P[i] = (double *)malloc( b * sizeof(double *) );
	}

	double **tmp_row = (double **)malloc(1 * sizeof(double *));
	tmp_row[0] = (double *)malloc(n * sizeof(double));

	double *msp = (double *)malloc(m * sizeof(double));


	int max_index;
	int min_index;
	for (int k = h; k < b; ++k, --m)
	{
		compute_Mt(ssParams, ssParams->ref_set, ssParams->ref_set_size, M, n, k);

		for (int i = 0; i < m; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				tmp_row[0][j] = ssParams->ref_set->members[0].params[j] - ssParams->scatter_set->members[i].params[j];
			}
			matrix_product(ssParams, tmp_row, 1, n, M, n, k, &(P[i]), 1, k);
			msp[i] = max(P[i], k, &max_index);
		}
		min(msp, m, &min_index);

		evaluate_ind(ssParams, &(ssParams->scatter_set->members[min_index]), inp, out);

		copy_ind(ssParams, &(ssParams->ref_set->members[k]), &(ssParams->scatter_set->members[min_index]));
		delete_and_shift(ssParams, ssParams->scatter_set, ssParams->scatter_set_size, min_index);

		#ifdef STATS
			update_frequency_matrix(ssParams, &(ssParams->ref_set->members[k]));
		#endif
	}

	// @todo: Free arrays, M, P, msp

}

/**
 * @brief      Calculates the Mt. The matrix of difference between each paramter of each individual 
 * and the best memeber of the Reference Set.
 *
 */
void compute_Mt(SSType *ssParams, Set *set, int set_size, double **M, int m_row, int m_col){
	for (int i = 0; i < m_row; ++i){		
		for (int j = 0; j < m_col; ++j)
		{
			M[i][j] = ssParams->ref_set->members[0].params[i] - ssParams->ref_set->members[j + 1].params[i];
		}
	}
}

/**
 * @brief      Generate array of subsets for recombination process using the members of Reference Set
 */
void select_subsets_list(SSType *ssParams, Set *set, int set_size){

	int k = 0;
	for (int i = 0; i < ssParams->ref_set_size; ++i)
	{
		for (int j = i + 1; j < ssParams->ref_set_size; ++j)
		{

			if( !is_equal(ssParams, &(set->members[i]), &(set->members[j])) 
				 && !is_exist_in_subsets_list(ssParams, &(set->members[j]), &(set->members[i])) ){
				// @note: The candidateSet generation could be performed here to save some computation time.
				// 		 but the code will lose its clarity.
				
				copy_ind(ssParams, &(ssParams->subsets_list[k].members[0]), &(set->members[i]));

				copy_ind(ssParams, &(ssParams->subsets_list[k].members[1]), &(set->members[j]));

				k++;
				ssParams->subsets_list_size = k;
			}
		}
	}
}

