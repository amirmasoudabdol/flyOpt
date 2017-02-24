/**
 * @file sort.c
 * @author Amir M. Abdol
 * @contact a.m.abdol@uva.nl
 * @date Created on May 2015
 *
 * @brief Simple implementation of QuickSort algorithm.
 */

#include "ss.h"

/**
 * @brief      Apply the Quick Sort algorithm on a set
 *
 * @param[in]  key       Meter to be used. 'c': cost, 'd': distance
 */
void quick_sort_set(SSType *ssParams, Set *set, int set_size/*, char key*/) {

	double *numbers = (double *)malloc(set_size * sizeof(double));

	for (int i = 0; i < set_size; ++i)
	{
		/*if (key == 'c')*/
		numbers[i] = set->members[i].cost;
		/*
		else if (key == 'd')
		 	numbers[i] = set->members[i].distance;
	 	*/
	}

	quick_sort(ssParams, set, set_size, numbers, 0, set_size - 1);
	free(numbers);
}


/**
 * @brief      Perform QuickSort in-place on array of numbers and rearrange Set based
 * on the result in-place.
 *
 */
void quick_sort(SSType *ssParams, Set *set, int set_size, double* numbers, int left, int right) {

	static individual pivot_ind;
	double pivot;
	int l_hold, r_hold;

	l_hold = left;
	r_hold = right;
	pivot = numbers[left];

	pivot_ind = set->members[left];

	while (left < right)
	{
		while ((numbers[right] >= pivot) && (left < right))
			right--;

		if (left != right)
		{
			numbers[left] = numbers[right];
			set->members[left] = set->members[right];
			left++;
		}

		while ((numbers[left] <= pivot) && (left < right))
			left++;

		if (left != right)
		{
			numbers[right] = numbers[left];
			set->members[right] = set->members[left];
			right--;
		}
	}

	numbers[left] = pivot;
	set->members[left] = pivot_ind;

	pivot = left;
	left = l_hold;
	right = r_hold;

	if (left < pivot)
		quick_sort(ssParams, set, set_size, numbers, left, pivot-1);

	if (right > pivot)
		quick_sort(ssParams, set, set_size, numbers, pivot+1, right);
}

/**
 * @brief      Remove the last item of the set and put the ind at it's correct place, assuming that 
	the list is already a sorted list. It's basically a reverse Insertion Sort!
 * @param[in]  key       Meter type, 'c': cost. 'd': distance. Currently not implemented, cost will be always used.
 */
void insertion_sort(SSType *ssParams, Set *set, int set_size /*, char key*/) {
        
	individual temp;
	/* Start at `set_size - 1` since the last item should be replaced anyway */
    int j =  set_size - 1;
 
    while (j > 0 && set->members[j].cost < set->members[j - 1].cost) 
    {
    	temp 				= set->members[j - 1];
    	set->members[j - 1] = set->members[j];
    	set->members[j] 	= temp;

     	j--;

    }  
}
