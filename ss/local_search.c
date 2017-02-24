/**
 * @file local_search.c
 * @author Amir M. Abdol, Anton Crombach, Damjan Cicin-Sain
 * @contact a.m.abdol@uva.nl
 * @date Created on May 2015
 *
 * @brief Implementation of Local Search routines.
 * 
 * Modified on August 2018
 */

#include "ss.h"
#include <gsl/gsl_vector_double.h>

// AC: probably added by Damjan
SSType *ssParamsLocal;
Input *inpLocal;
ScoreOutput *outLocal;
individual indLocal;
// \AC

/**
 * @brief      A link to GSL implementation of Nelder-Mead function for 
 * performing the local optimization on a individual.
 */
void nelder_mead(SSType *ssParams, individual *ind, Input *inp, 
	ScoreOutput *out) {

	/* AC: Damjan added 3 static vars (i.e. file wide accessible */
    outLocal = out;
    inpLocal = inp;
    ssParamsLocal = ssParams;

    int status;
    unsigned int iter = 0;
    const size_t p = ssParams->nreal;
    double size;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_multimin_function f;

    gsl_vector_view x = gsl_vector_view_array( ind->params, p );

/*#ifdef DEBUG
    print_ind( ssParams, ind, 1 );
    print_ind( ssParams, ssParams->best, 1 );

    gsl_vector_fprintf(stdout, &x.vector, "%f");
    printf("---\n");
#endif*/

    gsl_vector *ss = gsl_vector_alloc( p );
    gsl_vector_set_all( ss, 0.1 );

    gsl_rng_env_setup();
    const gsl_rng_type *type = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc ( type );

    f.f = &nelder_objfn;
    f.n = p;
    f.params = &x.vector;

    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc ( T, p );
    gsl_multimin_fminimizer_set( s, &f, &x.vector, ss );

    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate( s );

/*#ifdef DEBUG
        for ( int i = 0; i < s->x->size; ++i ) {
            printf( "%lf, ", gsl_vector_get( s->x, i ) );
        }

        printf( ": %lf\n", s->fval );
#endif*/

        if ( s->fval < ind->cost ) {
            /* Replace ind->params with newly generated params */
            ind->params = s->x->data;
            ind->cost = s->fval;
        }

        if ( status ) {
            break;
        }

        size = gsl_multimin_fminimizer_size( s );
        status = gsl_multimin_test_size( size, 1e-3 );
/*#ifdef DEBUG            
        printf("NM iter=%d\n", iter);
#endif*/
    } while ( status == GSL_CONTINUE && iter < ssParams->max_no_improve );

    gsl_vector_free( ss );
    gsl_rng_free ( r );
}

/* AC: probably added by Damjan */
/**
 * @brief      Evaluating the cost of an individual by updating `indLocal` values
 * by the new vector generated from `gsl_multimin_`
 *
 * @param[in]  x     gsl_vector containing values from corresponding individual
 * @param      data  The data
 *
 * @return     The cost of an individual
 */
double nelder_objfn( const gsl_vector *x, void *data ) {

    indLocal.params = x->data;
    evaluate_ind( ssParamsLocal, &( indLocal ), inpLocal, outLocal );
    return indLocal.cost;
}

/**
 * @brief      Simple Stochastic Hill Climbing routine. Climbs from `params` to `new_params`
 * and finally return `new_params`
 */
void take_step(SSType *ssParams, double *params, double *new_params) {

	double max_v, min_v;
	for (int i = 0; i < ssParams->nreal; ++i)
	{
		/* Determine lower and upper bound for drawing a random number. */
		min_v = MAX( ssParams->min_real_var[i], 
			params[i] - ssParams->step_size );
		max_v = MIN( ssParams->max_real_var[i], 
			params[i] + ssParams->step_size );

		new_params[i] = rndreal( min_v, max_v );
	}
}
