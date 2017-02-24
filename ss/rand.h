/**
 * @file rand.h
 * @author Amir M. Abdol
 * @contact a.m.abdol@uva.nl
 * @date Created on May 2015
 *
 * @brief Implementation of randomization function
 */

#include "ss.h"

/**
 * @brief      Return a real random number between `low` and `high`
 *
 * @param[in]  low   
 * @param[in]  high  
 */
double rndreal (double low, double high) 
{
    return (low + (high-low) * RandomReal() );
}


/**
 * @brief      Generate ::Set of random ::individual
 */
void generate_random_set(SSType *ssParams, Set *set, int set_size){

	for (int i = 0; i < set_size; ++i)
	{
		random_ind(ssParams, &(set->members[i]));
	}
}

/*
	Generate random values with respect to the min_real_var and max_real_var for every individual in a set
 */
/**
 * @brief      Generate an ::individual which each of the parameter are bounded by
 * its corresponding bounds (`min_real_var[i]`, `max_real_var[i]`)
 *
 */
void random_ind(SSType *ssParams, individual *ind){

	for (int i = 0; i < ssParams->nreal; ++i)
	{
		ind->params[i] = rndreal(ssParams->min_real_var[i], ssParams->max_real_var[i]);
	}
}