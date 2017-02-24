/**
 * @file ss.h
 * @author Amir M. Abdol, Anton Crombach
 * @contact a.m.abdol@uva.nl
 * @date Created on May 2015
 * 
 * Modified August 2018
 *
 * @brief Contains all necessary data type and function declaration for Scatter Search
 */

#ifndef SS_INCLUDED
#define SS_INCLUDED

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "maternal.h"
#include "../utils/random.h"

/* AC: Nelder-Mead local search */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit_nlin.h>

// #define STATS
#define DEBUG

// AC: Ansi control codes?
#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

#define eul  2.71828182845905
#define pi 3.14159265358979

#ifndef MAX
	#define MAX(x, y) (((x) > (y)) ? (x) : (y))
	#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif

typedef enum { false, true } bool;


typedef struct individual
{
	double* params;
	double cost;
	// double distance; 		// The distance of an individual to a set. Only uses in select_ref_set function.

} individual;

typedef struct Set
{
	individual *members;
	// int size;

} Set;

// typedef struct subSetItem
// {
// 	individual *pairs;
// 	struct subSetItem *next;
// } subSetItem;

/* Struct to keep track of files (originally from lsa.c) */
typedef struct Files {
    char *inputfile;            /* name of the input file */
    char *outputfile;           /* name of the output file */
} Files;

typedef struct SSType
{
	int seed;							//!< Random seed used to initialize the random generator
	int max_iter;						//!< Maximum number of iteration allowed
	int max_elite;						//!< Maximum number of elit members allowed in a set: `ref_set_size / 2`
	int n_iter;							//!< Counter for current iteration

	int nreal;							//!< Number of parameters to be optimized. ::individual size will be set to this.
	double sol;							//!< Possible guess or analytical value of the minima
	           							//! for instace, if it is minimization, 0 would be a guess.
	           							//! It will be used to avoid unneccessary local searches

	double* min_real_var;				//!< Array of values indicating lower bounds for each variables.
	double* max_real_var;				//!< Array of values indicating upper bounds for each variables.
	
	int p;								//!< The number of sub-regions for each boundary interval. This will be used to make a grid of p*n to 
	      								//! select paramters from to form the Scatter Set
	double **min_boundary_matrix;		//!< Matrix of lower boundaries for sub-regions
	double **max_boundary_matrix;		//!< Matrix of upper boundaries for sub-regions

	individual *best;					//!< Pointer to the best individual in a RefSet

	int ref_set_size;					//!< Size of Reference Set
	Set *ref_set;						//!< Reference Set

	int scatter_set_size;				//!< If it's not set manually it will be calculated as: `scatter_set_size = 10 * ref_set_size`
	Set *scatter_set;					//!< Scatter Set
	
	int pair_size;						//!< The size of the pairs in subSetItem; normally `2`
	int subsets_list_size;				//!< The size of the linked list containing the candidates
	Set *subsets_list;					//!< List of sets with size of two

	// Set *recombined_list;			// Basically the subsets_list, first will be recombined
	Set	*candidates_set;				//!< Candidate Set consists of individuals generated from linear combination of every pair of ::individual
	int candidates_set_size;			//!< Set of Candidate Set
	
	double dist_epsilon;				//!< Minimum Euclidean distance between two parameters vectors
	double fitness_epsilon;	 			//!< Minimum difference between the cost of two individuals

	int perform_ref_set_regen;			//!< Whether the Reference Set should be regenrated during the optimization or not.
	int ref_set_regen_freq;				//!< The frequency of performing regenration on Reference Set.


	/* Stats */
	int n_refinement;					//!< Number of local search performed
	int n_ref_set_update;				//!< Number of substitution in refSet
	int n_duplicates;					//!< Number of duplicates deteceted in candidateSet
	int n_flatzone_detected;			//!< Number of flatzone detected during the updating procedure
	int n_function_evals;				//!< Number of function objective function evaluations
	int n_regen;						//!< Number of refSet regenration performed
	int n_duplicate_replaced;			//!< Number of duplicate ::individual being detected and replaced during the process
	
	int **freqs_matrix;					//!< Frequencies of parameters being in sub-regions
	double **probs_matrix;				//!< Probabilities of parameters being in sub-regions

	int perform_stop_criteria;			//!< Whether to stop the optimization when the algorithm detects the difference between the
	                          			//! cost of the worst and best individual in Reference Set is less than certain value. 
	double stop_criteria;				//!< Minimum difference between best and worst individauls in a set for perfoming 'stop_criteria'

	int perform_warm_start;				//!< Whether to perform warm start or not


	int perform_flatzone_detection;		//!< Whether to detect if an individual is in a flat zone during the optimization or not


	/* Output */
	char *ref_set_final_filename;		//!< Filename to be used for exporting the final Reference Set
	char *freq_mat_final_filename;		//!< Filename to be used for exporting the final Frequency Matrix
	char *prob_mat_final_filename;		//!< Filename to be used for exporting the final Probablity Matrix


	/* Local Search Parameters */
	int perform_local_search;			//!< Whether to perform Local Search on individuals or not
	int local_search_freq;				//!< Frequency of applying Local Search. Since Local Search is usually very expensive
	                      				//! in comparison to the rest of procedure, one might one to perform it less frequently
	char local_search_method;			//!< Local Search Method:
	                         			//!	- 'n': Nelder-Mead
	                         			//!	- 't': Stochastic Hill Climbing

	int filter_good_enough;				//!< Flag to restrict local search to well-scoring individuals
	double good_enough_score_diff;  	//!< Score difference with the best solution that makes an individual good enough

	int filter_different_enough;		//!< Flag to restrict local search to individuals that are different enough from all others in the ref set
	double different_enough_param_dist;	//!< Minimal parameter distance to be different enough
	double different_cost_margin;       //!< Fraction of closest member cost that is added and subtracted to the closest member cost to create an interval outside of which the filtered individual should be

	int max_no_improve;					//!< Maximum number of attemps for Stochastic Hill Climbing algorithm
	double step_size;					//!< Step size for fluctuating paramters in Stochastic Hill Climbing

} SSType;


/*
				Output Files [Global]
 */
extern FILE *ref_set_history_file;
extern FILE *best_sols_history_file;
extern FILE *freqs_matrix_file;
extern FILE *freq_mat_final_file;
extern FILE *prob_mat_final_file;
extern FILE *ref_set_final_file;
extern FILE *can_set_history_file;
extern FILE *stats_file;

/*
				Functions Prototypes
 */

// input.c
void read_input(SSType *ssParams, int argc, char const *argv[]);

// ss.c
void InitSS(Input *inp, SSType *ssParams, Files *files);
void RunSS(Input *inp, SSType *ssParams, Files *files);

// init.c
void init_ssParams(SSType *ssParams);
void init_scatter_set(SSType *ssParams, Set *set);
void diversify(SSType *ssParams, Set *set, int set_size);

// recombine.c
void recombine_subsets_list(SSType *ssParams);
void recombine_subset(SSType *ssParams, Set *subset);
void recombine_ind(SSType *ssParams, individual *ind, double step);

// select.c
void select_subsets_list(SSType *ssParams, Set *set, int set_size/*, subSetItem *subsets_list*/);
void select_ref_set(SSType *ssParams, Set *set, int set_size);
bool select_best(SSType *ssParams);
void init_scatter_set(SSType *ssParams, Set *set);
void init_ref_set(SSType *ssParams);
void generate_candiates(SSType *ssParams);
void generate_ind_candidate(SSType *ssParams, individual *base, individual *candidate,  double *dists, char type);

// refine.c
bool is_in_flatzone(SSType *ssParams, Set *set, int set_size, individual *ind);
void update_ref_set(SSType *ssParams);
void replace(SSType *ssParams, individual *dest, individual *src /*, char sort_to_perform*/);
void compute_Mt(SSType *ssParams, Set *set, int set_size, double **M, int m_row, int m_col);
void re_gen_ref_set(SSType *ssParams, Set *set, int set_size, char type, Input *inp, ScoreOutput *out);

// allocate.c
void allocate_ind_memory(SSType *ssParams, individual *ind, int member_length);
void allocate_set_memory(SSType *ssParams, Set *set, int set_size, int member_length);
void deallocate_ind_memory(SSType *ssParams, individual *ind);
void deallocate_set_memory(SSType *ssParams, Set *set, int set_size);
void deallocate_subsets_list_memory(SSType *ssParams);
void deallocate_ssParam(SSType *ssParams);

// sort.c
void quick_sort_set(SSType *ssParams, Set *set, int set_size /*, char key*/);
void quick_sort(SSType *ssParams, Set *set, int set_size, double *numbers, int left, int right);
void insertion_sort(SSType *ssParam, Set *set, int set_size /*, char key*/);

// local_search.c
void nelder_mead(SSType *ssParams, individual *ind, Input *inp, ScoreOutput *out);
void refine_set(SSType *ssParams, Set *set, int set_size, char method, Input *inp, ScoreOutput *out);
void refine_individual(SSType *ssParams, Set *set, int set_size, individual *ind, char method, Input *inp, ScoreOutput *out);
void take_step(SSType *ssParams, double *params, double *new_params);
double nelder_objfn( const gsl_vector *x, void *data ); //Damjan

// ssTools.c
void distance_to_set(SSType *ssParams, Set *set, int set_size, individual *ind);
double euclidean_distance(SSType *ssParams, individual *ind1, individual *ind2);
double dist(double param1, double param2);
void random_ind(SSType *ssParams, individual *ind);
void generate_random_set(SSType *ssParams, Set *set, int set_size);
double rndreal(double low, double high);
void copy_params(SSType *ssParams, individual *ind1, individual *ind2);
void update_bestSet(SSType *ssParams, individual *best);
bool is_equal(SSType *ssParams, individual *ind1, individual *ind2);
int is_exist(SSType *ssParams, Set *set, int set_size, individual *ind);
bool is_exist_in_subsets_list(SSType *ssParams, individual *ind1, individual *ind2);
bool is_subset_exist(SSType *ssParams, Set *subsets_list, int subsets_list_size, Set *subset, int subset_size, int member_length);
double min(const double *arr, int length, int *index);
double max(const double *arr, int length, int *index);
void delete_and_shift(SSType *ssParams, Set *set, int set_size, int index_to_delete);
int closest_member(SSType *ssParams, Set *set, int set_size, individual *ind, int ind_index);
void copy_ind(SSType *ssParams, individual *src, individual *dest);
void parse_double_row(SSType *ssParams, char *line, double *row);
void parse_int_row(SSType *ssParams, char *line, int *row);
void init_warm_start(SSType *ssParams);
void matrix_product(SSType *ssParams, double **A, int a_row, int a_col, double **B, int b_row, int b_col, double **P, int p_row, int p_col);

// report.c
void init_report_files(SSType *ssParams, Files *files);
void write_set(SSType *ssParams, Set *set, int set_size, int member_length, FILE *fpt, int iter, char mode);
void write_ind(SSType *ssParams, individual *ind,  int member_length, FILE *fpt, int iter, char mode);
void print_set(SSType *ssParams, Set *set, int set_size, int member_length);
void print_ind(SSType *ssParams, individual *ind, int member_length);
void print_subsets_list(SSType *ssParams);
void print_double_matrix(SSType *ssParams, double **matrix, int row, int col);
void print_int_matrix(SSType *ssParams, int **matrix, int row, int col);
void write_int_matrix(SSType *ssParams, int **matrix, int row, int col, FILE *ftp, int iter, char mode);
void write_double_matrix(SSType *ssParams, double **matrix, int row, int col, FILE *ftp, int iter, char mode);
void loadBar(int x, int n, int r, int w);
void write_refset_eqparms(SSType *ssParams, Files *files, Input *inp);
void write_stats_header( FILE *fp );

// stats.c
void update_frequency_matrix(SSType *ssParams, individual *ind);
double average_cost_refset(SSType *ssParams, Set *set, int set_size);
double var_cost_refset(SSType *ssParams, Set *set, int set_size);

// evaluate.c
double objective_function(double *s, SSType *ssParams, Input *inp, ScoreOutput *out);
void evaluate_ind(SSType *ssParams, individual *ind, Input *inp, ScoreOutput *out);
void evaluate_set(SSType *ssParams, Set *set, int set_size, Input *inp, ScoreOutput *out);

#endif
