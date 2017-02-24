/**
 * @file stats.c
 * @author Amir M. Abdol
 * @contact a.m.abdol@uva.nl
 * @date Created on May 2015
 *
 * @brief Implementation of simple frequency tracker.
 */

#include "ss.h"

/**
 * @brief      Keep track of number of values in each sub-regions of search space.
 * This can be used to plot a heatmap of a search space.
 *
 */
void update_frequency_matrix(SSType *ssParams, individual *ind){

	for (int i = 0; i < ssParams->nreal; ++i)
	{
		for (int j = 0; j < ssParams->p; ++j)
		{
			if (ind->params[i] > ssParams->min_boundary_matrix[i][j] && ind->params[i] < ssParams->max_boundary_matrix[i][j])
			{
				ssParams->freqs_matrix[i][j]++;
				break;
			}
		}
	}
}

/**
 * @brief	Calculate average cost of the individuals in the refset.
 */
double average_cost_refset(SSType *ssParams, Set *set, int set_size) {

	/* Naive implementation */
	double sum = 0.0;
	for (int i = 0; i < set_size; ++i)
	{
		sum += set->members[i].cost;
	}
	return sum / set_size;
}

/**
 * @brief	Calculate variation of individual's cost (in the refset).
 */
double var_cost_refset(SSType *ssParams, Set *set, int set_size) {

    /* https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance */
    double mean = 0.0, s = 0.0, delta = 0.0;
    for (int i = 0; i < set_size; ++i)
	{       
        delta = set->members[i].cost - mean;
        mean += delta / (i + 1);
        s += delta * (set->members[i].cost - mean);
	}
    return s / (set_size - 1);
}
