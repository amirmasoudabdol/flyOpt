/**
 * @file ssTools.c
 * @author Amir M. Abdol
 * @contact a.m.abdol@uva.nl
 * @date Created on May 2015
 *
 * @brief Implementation of tools used in different functions.
 */

#include "ss.h"
#include "rand.h"
#include "string.h"
#include <math.h>

/*
	Compute the distance of a member from all the members of set
 */
// void distance_to_set(SSType *ssParams, Set *set, int set_size, individual *ind){

// 	double dist = 0;
// 	for (int i = 0; i < set_size; ++i)
// 	{
// 		dist += euclidean_distance(ssParams, &(set->members[i]), ind);
// 	}
// 	ind->distance = dist;
// }

/**
 *	@brief Compute the Eculidean Distance between two individuals
 */
double euclidean_distance(SSType *ssParams, individual *ind1, individual *ind2){

	double distance = 0;
	for (int i = 0; i < ssParams->nreal; ++i)
	{
		distance += (ind1->params[i] - ind2->params[i]) * (ind1->params[i] - ind2->params[i]);
	}
	return sqrt(distance);

}

/**
 * @brief      Compute the product of two matrices.
 *
 */
void matrix_product(SSType *ssParams, double **A, int a_row, int a_col, double **B, int b_row, int b_col, double **P, int p_row, int p_col){

	int sum  = 0;
    for (int i = 0 ; i < a_row ; i++ )
    {
      for (int j = 0 ; j < b_col ; j++ )
      {
        for (int k = 0 ; k < b_row ; k++ )
        {
          sum = sum + A[i][k]*B[k][j];
        }
 
        P[i][j] = sum;
        sum = 0;
      }
    }
    p_row = a_row;
    p_col = b_col;

}

/**
 * @brief      Return the index of closest member in a set to a selected individual.
 * The distance will be computed using Euclidean distance.
 *
 * @param      ind        Individual to find the closest member to.
 * @param[in]  ind_index  An index from Set indicating the closest member to `ind`
 *
 */
int closest_member(SSType *ssParams, Set *set, int set_size, individual *ind, int ind_index){
	double dist;
	double min;
	int min_index;

	if (ind_index == set_size - 1 ){
		min = euclidean_distance(ssParams, ind, &(set->members[set_size - 2]));
		min_index = set_size - 2;
	}
	else if (ind_index == 0){
		min = euclidean_distance(ssParams, ind, &(set->members[1]));
		min_index = 1;		
	}
	else{
		min = euclidean_distance(ssParams, ind, &(set->members[ind_index - 1]));
		min_index = ind_index - 1;		
	}


	for (int i = 0; i < set_size; ++i)
	{
		// printf("%d: ", i);
		// printf("%lf, %d\n", min, min_index);
		if ( i != ind_index ){
			// printf("%d++@#\n", ind_index);
			dist = euclidean_distance(ssParams, ind, &(set->members[i]));
			if (dist < min ){
				min = dist;
				min_index = i;
			}

		}
		

	}
	// printf("====%d\n", min_index);
	return min_index;
}


/**
 * @brief      Return an index to minimum value in first `length` member of an array.
 *
 * @param[in]  arr     The arr
 * @param[in]  length  Indicate the subset length of array to search.
 * @param      index   Index of the minimum value.
 *
 */
double min(const double *arr, int length, int *index) {

    double minimum = arr[0];
    *index = 0;
    for (int i = 1; i < length; ++i) {
        if (minimum > arr[i]) {
            minimum = arr[i];
            *index = i;
        }
    }
    return minimum;
}

/**
 * @brief      Return the maximum value of an array with its `index`
 *
 * @param[in]  arr     The arr
 * @param[in]  length  Indicate the subset length of array to search
 * @param      index   Index of the maximum value.
 *
 */
double max(const double *arr, int length, int *index) {

    double maximum = arr[0];
    *index = 0;
    for (int i = 1; i < length; ++i) {
        if (maximum < arr[i]) {
            maximum = arr[i];
            *index = i;
        }
    }
    return maximum;
}

/**
 * @brief      Delete an item at `index_to_delete` from a Set and shift all the members to left.
 *
 * @param      ssParams         The ss parameters
 * @param      set              The set
 * @param[in]  set_size         The set size
 * @param[in]  index_to_delete  The index to delete
 */
void delete_and_shift(SSType *ssParams, Set *set, int set_size, int index_to_delete){

	for (int i = index_to_delete; i < set_size - 1; ++i)
	{
		copy_ind(ssParams, &(set->members[i]), &(set->members[i + 1]));
	}

}

/**
 * @brief      Check if two individual are equal to each others. 
 * (being close enough to each other)
 *
 *
 * @return     True if equal, False otherwise.
 */
bool is_equal(SSType *ssParams, individual *ind1, individual *ind2){

	bool isEqual = false;
	if ( euclidean_distance(ssParams, ind1, ind2) < ssParams->dist_epsilon )
		isEqual |= 1;	

	return isEqual;
}

/**
 * @brief      Determines if an individual exist in a set or not
 *
 * @return     -1 if the member doesn't exist. `index` if it's exist in the set.
 */
int is_exist(SSType *ssParams, Set *set, int set_size, individual *ind){

	int index = -1;
	for (int i = set_size - 1; i >= 0; --i)
	{
		if ( is_equal(ssParams, &(set->members[i]), ind) ){
			index = i;
			break;
		}

	}

	return index;
}

/**
 * @brief      
 * @note Not implemented
 *
 * @param[in]  member_length      The member length
 *
 * @return     True if subset exist, False otherwise.
 */
bool is_subset_exist(SSType *ssParams, Set *subsets_list, int subsets_list_size, Set *subset, int subset_size, int member_length){

	return false;
}

/**
 * @brief      Determines if `ind1` are `ind2` together make a subset.
 *
 * @note       Probably not used
 *
 * @param      ssParams  The ss parameters
 * @param      ind1      The ind 1
 * @param      ind2      The ind 2
 *
 * @return     True if exist in subsets list, False otherwise.
 */
bool is_exist_in_subsets_list(SSType *ssParams, individual *ind1, individual* ind2){
	for (int i = 0; i < ssParams->subsets_list_size; ++i)
	{
		if ( (is_equal(ssParams, ind1, &(ssParams->subsets_list[i].members[0]))
			  	&& is_equal(ssParams, ind2, &(ssParams->subsets_list[i].members[1]))) 
			 ||
			 (is_equal(ssParams, ind2, &(ssParams->subsets_list[i].members[0]))
			  	&& is_equal(ssParams, ind1, &(ssParams->subsets_list[i].members[1])))
			 )
		{
				return true;
				break;
		}

	}

	return false;
}

/**
 * @brief      Copy one individual to another using memcpy. (`src` to `dest`)
 * 
 * @note       Assume that the memory of `dest` is already allocated!
 *
 * @param      ssParams  The ss parameters
 * @param      dest      The destination
 * @param      src       The source
 */
void copy_ind(SSType *ssParams, individual *dest, individual *src){
	memcpy(dest->params, src->params, ssParams->nreal*sizeof(double));
	dest->cost =  src->cost;
	// dest->distance = src->distance;
}

/**
 * @brief      Determines if an individual is in a flatzone.
 *
 *
 * @return     True if in flatzone, False otherwise.
 */
bool is_in_flatzone(SSType *ssParams, Set *set, int set_size, individual *ind){

	bool isInFlatzone = false;

	/* 
		The loop doesn't check the last item since it is the best sol and every good solution in
		comparison to that is in flatzone coverd by that!
	 */
	for (int i = 0; i < set_size; ++i)	// set: ref_set
	{
		/* code */
		// printf("%f, %f, %f\n", ind->cost , set->members[i].cost, set->members[i].cost * ( 1 - ssParams->fitness_epsilon) );
		// if ( ind->cost > set->members[i].cost * ( 1 - ssParams->fitness_epsilon) )
		if (ind->cost < set->members[i].cost + (set->members[i].cost * ssParams->fitness_epsilon) 
		 && ind->cost > set->members[i].cost - (set->members[i].cost * ssParams->fitness_epsilon))
		{
			ssParams->n_flatzone_detected++;
			isInFlatzone |= 1;
			break;
		}

	}
	return isInFlatzone;
}

/**
 * @brief      Parse a tabular row of double values.
 *
 * @param      ssParams  The ss parameters
 * @param      line      The line
 * @param      row       The row
 */
void parse_double_row(SSType *ssParams, char *line, double *row){

    int i = 0;
    const char* tok;
    for (tok = strtok(line, "\t"); tok && *tok; i++, tok = strtok(NULL, "\t\n"))
    {
        row[i] = atof(tok);
    }
}

/**
 * @brief      Parse a tabular row of integer values.
 *
 * @param      ssParams  The ss parameters
 * @param      line      The line
 * @param      row       The row
 */
void parse_int_row(SSType *ssParams, char *line, int *row){

    int i = 0;
    const char* tok;
    for (tok = strtok(line, "\t"); tok && *tok; i++, tok = strtok(NULL, "\t\n"))
    {
        row[i] = atoi(tok);
    }
}