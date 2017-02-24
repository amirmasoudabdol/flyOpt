/**
 * @file essRand.c
 * @author Amir M. Abdol
 * @contact a.m.abdol@uva.nl
 * @date Created on Sep 2015
 *
 * @brief Implementation of simple randomization routines
 */

 #include "ess.h"

double rndreal(double low, double high){

	// printf("%lf\n", RandomReal());
    return (low + (high-low) * RandomReal() );	
}


void random_Set(eSSType *eSSParams, Set *set, double *low, double *high){
	for (int i = 0; i < set->size; ++i)
	{
		random_Ind(eSSParams, &(set->members[i]), low, high);
	}
}

void random_Ind(eSSType *eSSParams, individual *ind, double *low, double *high){

	for (int i = 0; i < eSSParams->n_Params; ++i){
		if (!eSSParams->logBound)
			ind->params[i] = rndreal(low[i], high[i]);
		else{
			double mx, mn, la;
			mn = low[i];
			mx = high[i];
			la = log10(mx) - log10(mn);

          // determine if linear or log scale
          if ((mn < 0.0) || (mx <= 0.0))
            ind->params[i] = mn + rndreal(0,1) * (mx - mn);
          else
            {
              if (la < 1.8)
                ind->params[i] = mn + rndreal(0, 1) * (mx - mn);
              else
                ind->params[i] = pow(10.0, log10(mn)) + la * rndreal(0, 1);
            }
		}

		ind->n_notRandomized = 0;
	}
}