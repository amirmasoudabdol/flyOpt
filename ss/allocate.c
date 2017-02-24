/**
 * @file allocate.c
 * @author Amir M. Abdol
 * @contact a.m.abdol@uva.nl
 * @date Created on May 2015
 *
 * @brief Implementation of memory allocation and deallocation routines
 */


#include "ss.h"

/**
 * Allocate memory for an individual type
 */
/**
 * @brief      Allocate memory for an individual type
 *
 * @param      ssParams       
 * @param      ind            
 * @param[in]  member_length  Length of array to be allocated
 */
void allocate_ind_memory(SSType *ssParams, individual *ind, int member_length){

	ind->params = (double *)malloc(member_length * sizeof(double));
	ind->cost = 0;
	// ind->distance = 0;
}

/**
 * @brief      Allocate memory for subset of a population.
 * @note       It's not being used at the moment. It declared with the idea in mind
 * that the linear combination can be applied on more than two ::individual
 *
 * @param      ssParams  
 * @param      pair      Pointer to an array of ::individual
 */
void allocate_subset_memory(SSType *ssParams, individual *pair){

	
}

/*
	Allocate for set of individuals
 */
/**
 * @brief      Allocate memory for a set
 *
 * @param      ssParams       
 * @param      set            
 * @param[in]  set_size       The size of a set
 * @param[in]  member_length  The length of each ::individual in the set
 */
void allocate_set_memory(SSType *ssParams, Set *set, int set_size, int member_length){

	set->members = (individual *)malloc( set_size * sizeof(individual) );
	for (int i = 0; i < set_size; ++i)
	{
		allocate_ind_memory(ssParams, &(set->members[i]), member_length);
	}

}

/**
 * @brief      Free up the memory allocated by an ::individual
 *
 * @param      ssParams  
 * @param      ind       
 */
void deallocate_ind_memory(SSType *ssParams, individual *ind){

	free(ind->params);
}

/*
	Free memory of set of individuals
 */
/**
 * @brief      Free up the memory allocated by a ::Set
 *
 * @param      ssParams  
 * @param      set       
 * @param[in]  set_size  Size of set
 */
void deallocate_set_memory(SSType *ssParams, Set *set, int set_size){

	for (int i = 0; i < set_size; ++i)
	{
		deallocate_ind_memory(ssParams, &(set->members[i]));
	}
	free(set->members);

}

/**
 * @brief      Free up the memory allocated by a sub_set
 *
 * @param      ssParams  
 */
void deallocate_subsets_list_memory(SSType *ssParams){
						// max size of sub_sets_list
	for (int i = 0; i < ssParams->ref_set_size * ssParams->ref_set_size; ++i)
	{
		deallocate_set_memory(ssParams, &(ssParams->subsets_list[i]), ssParams->pair_size);
	}
}


/**
 * @brief      Free up the memory allocated by ssParams
 *
 * @param      ssParams  
 */
void deallocate_ssParam(SSType *ssParams){

	free(ssParams->min_real_var);
	free(ssParams->max_real_var);

	for (int i = 0; i < ssParams->nreal; ++i){
		free(ssParams->freqs_matrix[i]);
		free(ssParams->probs_matrix[i]);
		free(ssParams->min_boundary_matrix[i]);
		free(ssParams->max_boundary_matrix[i]);
	}
	free(ssParams->freqs_matrix);
	free(ssParams->probs_matrix);
	free(ssParams->min_boundary_matrix);
	free(ssParams->max_boundary_matrix);

	deallocate_set_memory(ssParams, ssParams->ref_set, ssParams->ref_set_size);
	free(ssParams->ref_set);
	deallocate_set_memory(ssParams, ssParams->candidates_set, ssParams->ref_set_size * ssParams->ref_set_size * 6);
	free(ssParams->candidates_set);
	if ( !ssParams->perform_warm_start )
	{
		deallocate_set_memory(ssParams, ssParams->scatter_set, ssParams->scatter_set_size);
		free(ssParams->scatter_set);
	}
	deallocate_subsets_list_memory(ssParams);
	free(ssParams->subsets_list);

	free(ssParams->ref_set_final_filename);
	free(ssParams->freq_mat_final_filename);
	free(ssParams->prob_mat_final_filename);
}
