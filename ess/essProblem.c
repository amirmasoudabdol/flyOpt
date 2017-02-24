/**
 * @file essProblem.c
 * @author Amir M. Abdol
 * @contact a.m.abdol@uva.nl
 * @date Created on Sep 2015
 *
 * @brief Implementaion of Objective Function by linking eSS to fly simulator.
 */

 #include "ess.h"

#include "score.h"
#include "zygotic.h"

/**
 * This will call the problem simulator to the cost function for the individual.
 * The function should first copy ind->params to the inp->params and call the 
 * simulator like `simulator(inp, out)`, and finally write the `out` to `ind->cost`
 */
double objectiveFunction(eSSType *eSSParams, individual *ind, void *inp, void *out){

    for (int i = 0; i < ((Input *)inp)->tra.size; ++i){
        *( ((Input *)inp)->tra.array[i].param  ) = ind->params[i];
    }

	Score(inp, out, 0);
    return ((ScoreOutput*)out)->score + ((ScoreOutput*)out)->penalty;

}


double objfn(double x[]){
	return 0;
}


double nelder_objfn(const gsl_vector *x, void *data){
	
	return objfn(x->data);
}

int levermed_objfn(const gsl_vector *x, void *data, gsl_vector *f){
	return 0;
}