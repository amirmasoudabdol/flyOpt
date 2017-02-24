/**
 * @file refine.c
 * @author Amir M. Abdol, 
 * @contact a.m.abdol@uva.nl
 * @date Created on May 2015
 *
 * @brief Implementation of Refining routines. Handing the link between refine.c and local_search.c
 */

#include "ss.h"
#include <string.h>

/** Helper functions to improve readability of if-statements in refine_set() */
inline bool
is_good_enough( SSType *ssParams, individual *ind ) {
	/*
	 * `sol' is often set to zero, implying the test is a simple cost smaller
	 * than a user-given number.
	 */
	return fabs(ind->cost - ssParams->sol) < ssParams->good_enough_score_diff;
}

inline bool
has_params_different_enough( SSType *ssParams, individual *ind, individual *closest ) {

	return euclidean_distance( ssParams, ind, closest ) > ssParams->different_enough_param_dist;
}

inline bool
has_cost_different_enough( SSType *ssParams, individual *ind, individual *closest ) {

	return ind->cost > closest->cost + (ssParams->different_cost_margin * closest->cost ) ||
		ind->cost < closest->cost - (ssParams->different_cost_margin * closest->cost );
}


/**
* @brief      Refine a set by applying local search on its individuals.
* - **Filter 1,** decide whether local search will be applied on individual or not. If the difference
* between the current individual and `sol` is smaller than `good_enough_score_filter`. Then local
* search will be performed. _This is to avoid performing local search to an individual with a very bad 
* cost for two reasons: 1. keeping it in the Reference Set provides more diversity, 2. It will avoid
* unnecessary local searches._
* - **Filter 2,** first finds the closest vector to the current individual, then decide if their costs is larger
* than 'different_cost_margin' or not. If yes, then it applies the local search. _This is to avoid
* performing local search on two similar individuals in the Reference Set._
* 
* @note To perform local search on every individual, choose a very large value for 
* `good_enough_score_diff` and only activate `filter_good_enough`.
*
* @note `different_enough_param_dist` indicates the distance between two individuals to
* be considered as the same in parameter space. You can use 'dist_epsilon' which is a value 
* close to zero.
*
* @param[in]  method    Local Search Method. 'n': Nelder-Mead. 't': Stochastic Hill Climbing.
*/
void 
refine_set(SSType *ssParams, Set *set, int set_size, char method, Input *inp, 
	ScoreOutput *out) {

	/* for each member of the reference set */
	int closest_member_index;
	int different_flag, cost_flag;

	/*
	   First determine if we filter individuals, then loop over the set. By 
	   doing the loop later, we allow the compiler to optimize with some
	   vectorization (if it can...)
	*/

	/* Both filters have to be applied! */
	if (ssParams->filter_good_enough && ssParams->filter_different_enough)
	{
		for (int i = 0; i < set_size; ++i)
		{
			if (is_good_enough( ssParams, &(set->members[i])))
			{
				// find closest member for comparison of parameters and score
				closest_member_index = closest_member(ssParams, set, set_size, &(set->members[i]), i);
				different_flag = has_params_different_enough(ssParams, &(set->members[i]), &(set->members[ closest_member_index ]) );
				cost_flag = has_cost_different_enough(ssParams, &(set->members[i]), &(set->members[ closest_member_index ]) );

				if (different_flag && cost_flag)
				{
#ifdef DEBUG					
					printf("\tPassed both filters, doing local_search\n");
#endif
					refine_individual(ssParams, set, set_size, &(set->members[i]), method, inp, out);
				}				
			}
		}

	} 
	/* Only apply the "well-scoring individual" filter */
	else if (ssParams->filter_good_enough) 
	{
		for (int i = 0; i < set_size; ++i)
		{
			if (is_good_enough( ssParams, &(set->members[i])))
			{
				refine_individual(ssParams, set, set_size, &(set->members[i]), method, inp, out);
			}
		}
	}
	/* Only apply the "really different parameter set" filter */
	else if (ssParams->filter_different_enough)
	{
		for (int i = 0; i < set_size; ++i)
		{
			// find closest member for comparison of parameters and score
			closest_member_index = closest_member(ssParams, set, set_size, &(set->members[i]), i);
			different_flag = has_params_different_enough(ssParams, &(set->members[i]), &(set->members[ closest_member_index ]) );
			cost_flag = has_cost_different_enough(ssParams, &(set->members[i]), &(set->members[ closest_member_index ]) );

			if (different_flag && cost_flag)
			{						
				refine_individual(ssParams, set, set_size, &(set->members[i]), method, inp, out);
			}				
		}
	}
}

/**
 * @brief      Apply local search on an individual
 *
 * @param[in]  method    Local Search Method. 'n': Nelder-Mead. 't': Stochastic Hill Climbing.
 */
void refine_individual(SSType *ssParams, Set *set, int set_size, 
	individual *ind, char method, Input *inp, ScoreOutput *out) {

	individual new_candidate;
	double *new_params = (double *)malloc( ssParams->nreal * sizeof(double));

    switch(method) {
	case 'n':
		nelder_mead( ssParams, ind, inp, out );

	case 't':
		// Stochastic Hill Climbing 
		for (int i = 0; i < ssParams->max_no_improve; ++i)
		{
			take_step(ssParams, ind->params, new_params);
			new_candidate.params = new_params;
			evaluate_ind(ssParams, &( new_candidate), inp, out);
	
			if (new_candidate.cost < ind->cost)
			{
				/* Replace ind->params with newly generated params */
				copy_ind(ssParams, ind, &( new_candidate));
			}

		}
		break;
	}

	/* Track stats on how often we refine an individual */
	ssParams->n_refinement++;
	free(new_params);
}
