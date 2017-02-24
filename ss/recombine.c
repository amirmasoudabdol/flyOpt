/**
 * @file recombine.c
 * @author Amir M. Abdol
 * @contact a.m.abdol@uva.nl
 * @date Created on May 2015
 *
 * @brief Implementation of linear combination operator between two individuals
 */

#include "ss.h"

/**
 * @brief      Use the sub_sets_list members to produce candidates, 
 * the recombined_set will be generated. 
 * Defining `x1 := member[0]` and `x2 := member[1]`, 4 types of 
 * candidates will be produced:
 * - Type 0: `x1 - d` or `x1 + d` (decided by a random number)
 * - Type 1: `x1 - d`
 * - Type 2: `x1 + d`
 * - Type 3: `x2 + d`
 * where `d = r.(x2 - x1)/2`, `r` being a vector or random numbers in `[0, 1]`.
 * 
 * Based on `x1` and `x2` being elite members of Reference Set or not, different numbers
 * of candidate is being generated in each case.
 */
void generate_candiates(SSType *ssParams){

	int candidates_count = 0;
	double *dists        = (double *)malloc( ssParams->nreal * sizeof(double));
	double mid_cost      = ssParams->ref_set->members[ssParams->max_elite].cost;
	double diff;

	for (int i = 0; i < ssParams->subsets_list_size; ++i)
	{
		for (int j = 0; j < ssParams->nreal; ++j)
		{
			/* Step sizes in every direction */
			diff     = ( ssParams->subsets_list[i].members[1].params[j] -  ssParams->subsets_list[i].members[0].params[j] );
			// dists[j] = sqrt(diff * diff) / 2;
			dists[j] = diff / 2;
		}

		/* Generating new candidates */
		if ( ssParams->subsets_list[i].members[0].cost < mid_cost && ssParams->subsets_list[i].members[1].cost < mid_cost )		// type 1
		{
			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '0');

			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[1]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '0');
			// Type 1
			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '1');
			// Type 2
			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '2');
			// Type 3
			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[1]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '3');
			// Type 3
			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[1]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '3');

		}else
		if ( ssParams->subsets_list[i].members[0].cost < mid_cost && ssParams->subsets_list[i].members[1].cost >= mid_cost )	// type 2
		{
			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '0');
			// Type 1 
			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '1');
			// Type 2
			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '2');
			// Type 3
			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[1]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '3');

		}else
		if ( ssParams->subsets_list[i].members[0].cost >= mid_cost && ssParams->subsets_list[i].members[1].cost >= mid_cost )	// type 3
		{
			if (rndreal(0, 1) < 0.5)
			{	// Type 1
				generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '1');
			}else
			{	// Type 3
				generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[1]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '3');
			}
			// Type 2
			generate_ind_candidate(ssParams, &(ssParams->subsets_list[i].members[0]), &(ssParams->candidates_set->members[ candidates_count++ ]), dists, '2');

		}


	}

	ssParams->candidates_set_size = candidates_count;

	free(dists);

}

/**
 * @brief      Generate a new candidate by applying some noise فخ an individual.
 *
 * @param      base       Individual to perturb
 * @param      candidate  New individual created by pertubing values of `base`
 * @param      dists      Size of perturbation
 * @param[in]  type       Type of perturbation
 */
void generate_ind_candidate(SSType *ssParams, individual *base, individual *candidate, double *dists, char type){

	int i;
	double new_value;
	double rnd;
	switch (type){
		case '0':
			rnd = rndreal(0, 1);
			for (i = 0; i < ssParams->nreal; ++i)
			{
				new_value = base->params[i] - (dists[i] * rnd);
				new_value = MIN(new_value, ssParams->max_real_var[i]);
				new_value = MAX(new_value, ssParams->min_real_var[i]);
				candidate->params[i] = new_value;
			}
			break;
		case '1':
			for (i = 0; i < ssParams->nreal; ++i)
			{
				new_value = base->params[i] - (rndreal(0,1) * dists[i]);
				new_value = MIN(new_value, ssParams->max_real_var[i]);
				new_value = MAX(new_value, ssParams->min_real_var[i]);
				candidate->params[i] = new_value;
			}
			break;
		case '2':
			for (i = 0; i < ssParams->nreal; ++i)
			{
				new_value = base->params[i] + (rndreal(0,1) * dists[i]);
				new_value = MIN(new_value, ssParams->max_real_var[i]);
				new_value = MAX(new_value, ssParams->min_real_var[i]);
				candidate->params[i] = new_value;
			}
			break;
		case '3':
			for (i = 0; i < ssParams->nreal; ++i)
			{
				new_value = base->params[i] + (rndreal(0,1) * dists[i]);
				new_value = MIN(new_value, ssParams->max_real_var[i]);
				new_value = MAX(new_value, ssParams->min_real_var[i]);
				candidate->params[i] = new_value;
			}
			break;
	}

}