/**
 * @file evaluate.c
 * @author Amir M. Abdol
 * @contact a.m.abdol@uva.nl
 * @date Created on May 2015
 *
 * @brief Implementation of objective function and evaluation of 
 * ::individual and ::Set using Score() function
 */

#include "ss.h"

# include "score.h"
# include "zygotic.h"

/**
 * @brief      Updates values of parameters marked to be tweak in `inp` with
 * values corresponding values of passed as `s`; then calculate the score and
 * return `score+penalty` which will be set as a new `cost` for the passed ::individual
 *
 * @param      s         Array of parameters from an ::individual struct
 * @param      ssParams  
 * @param      inp       An ::Input variable that is being use to pass data to Score()
 * @param      out       An ::ScoreOutput variable uses to reterive data from Score()
 *
 * @return     `score + penalty` of the passed individual.
 */
double objective_function( double *s, SSType *ssParams, Input *inp, ScoreOutput *out ) {

	/* copy array of individual into another */
    for ( int i = 0; i < inp->tra.size; ++i ) {
        *( inp->tra.array[i].param  ) = s[i];
    }

    Score( inp, out, 0 );
    ssParams->n_function_evals++;
    return out->score + out->penalty;
}

/**
 * @brief      Evaluate the cost of an ::individual 
 */
void evaluate_ind(SSType *ssParams, individual *ind, Input *inp, 
	ScoreOutput *out) {

	ind->cost = objective_function(ind->params, ssParams, inp, out);
}

/**
 * @brief      Evaluate cost of each individual in a set
 * 
 * @param[in]  set_size  Set size
 */
void evaluate_set(SSType *ssParams, Set *set, int set_size, Input *inp, ScoreOutput *out) {

	for (int i = 0; i < set_size; ++i)
	{
		evaluate_ind(ssParams, &(set->members[i]), inp, out);
	}
}
