/**
 * @file report.c
 * @author Amir M. Abdol
 * @contact a.m.abdol@uva.nl
 * @date Created on May 2015
 *
 * @brief Implementaiton of pretty prints, debugging and file operation routines.
 */

#include "ss.h"
#include "fly_io.h"
#include <stdio.h>


/*
	Write the set to the file
 */
/**
 * @brief      Writing a set to the file
 * @param[in]  iter           Iteration number. The first column of a file refer indicate the current iteration.
 * @param[in]  mode           Either to overwrite or append to the file.
 */
void write_set(SSType *ssParams, Set *set, int set_size, int member_length, FILE *fpt, int iter, char mode){

	for (int i = 0; i < set_size; ++i)
	{
		write_ind(ssParams, &(set->members[i]), member_length, fpt, iter, mode);
	}
	// fflush(fpt);

}

/**
 * @brief      Write (append) the paramters of an individual to a file
 *
 */
void write_ind(SSType *ssParams, individual *ind, int member_length, FILE *fpt, int iter, char mode){

	// fprintf(fpt, "%d\t", iter);
	if (iter != -1)
		fprintf(fpt, "%d\t", iter);
	
	for (int i = 0; i < member_length; ++i)
	{
		fprintf(fpt, "%.5lf\t", ind->params[i]);
	}
	fprintf(fpt, "%lf\n", ind->cost);
	// fprintf(fpt, "\n");


}

/**
 * @brief      Print a Set to the terminal
 *
 */
void print_set(SSType *ssParams, Set *set, int set_size, int member_length){
	printf("-----------------------------------\n");
	for (int i = 0; i < set_size; ++i)
	{
		printf("%d: ", i); print_ind(ssParams, &(set->members[i]), member_length);
	}
	printf("\n");

}

/**
 * @brief      Print individual to the terminal
 *
 */
void print_ind(SSType *ssParams, individual *ind, int member_length){

	// for (int i = 0; i < member_length; ++i)
	// {
	// 	printf("%.5lf, ", ind->params[i]);
	// }	
	// printf("\t (cost: %lf)\t(distance: %lf)\n", ind->cost, ind->distance);
	printf("\t (cost: %lf)\n", ind->cost/*, ind->distance*/);

}

/**
 * @brief      Print the subset to the terminal
 *
 */
void print_subsets_list(SSType *ssParams){

	for (int i = 0; i < ssParams->subsets_list_size; ++i)
	{
		printf("[i: %d]\n", i);
		print_set(ssParams, &(ssParams->subsets_list[i]), ssParams->pair_size, ssParams->nreal);
	}
}

/**
 * @brief      Print a matrix of double
 *
 * @param      matrix    Pointer to a matrix
 * @param[in]  row       Number of row
 * @param[in]  col       Number of column
 */
void print_double_matrix(SSType *ssParams, double **matrix, int row, int col){
	for (int r=0; r<row; r++)
	{
		printf("%d: ", r);
	    for(int c=0; c<col; c++)
	         printf("%.4lf     ", matrix[r][c]);
	    printf("\n");
	 }
	 printf("===========================================\n");

}

/**
 * @brief      Print a matrix of integer
 *
 * @param      matrix    Pointer to a matrix
 * @param[in]  row       Number of row
 * @param[in]  col       Number of column
 */
void print_int_matrix(SSType *ssParams, int **matrix, int row, int col){
	for (int r=0; r<row; r++)
	{
		printf("%d: ", r);
	    for(int c=0; c<col; c++)
	         printf("%d\t", matrix[r][c]);
	    printf("\n");
	 }
	 printf("===========================================\n");

}


/**
 * @brief      Experimental implementation of progress bar in the terminal
 *
 */
void loadBar(int x, int n, int r, int w)
{
// Process has done i out of n rounds,
// and we want a bar of width w and resolution r.
    // Only update r times.
    if ( x % (n/r +1) != 0 ) return;
 
    // Calculuate the ratio of complete-to-incomplete.
    float ratio = x/(float)n;
    int   c     = ratio * w;
 
    // Show the percentage complete.
    printf("%3d%% [", (int)(ratio*100) );
 
    // Show the load bar.
    for (int x=0; x<c; x++)
       printf("=");
 
    for (int x=c; x<w; x++)
       printf(" ");
 
    // ANSI Control codes to go back to the
    // previous line and clear it.
    printf("]\n\033[F\033[J");
}

/**
 * @brief      Write a matrix of integer to a file
 *
 */
void write_int_matrix(SSType *ssParams, int **matrix, int row, int col, FILE *fpt, int iter, char mode){

	for (int r=0; r<row; r++)
	{
		if (iter != -1)
			fprintf(fpt, "%d", iter);
		
	    for(int c=0; c<col; c++)
	         fprintf(fpt, "\t%d", matrix[r][c]);
	    fprintf(fpt, "\n");
	 }
}

/**
 * @brief      Write a matrix of double to a file
 *
 */
void write_double_matrix(SSType *ssParams, double **matrix, int row, int col, FILE *fpt, int iter, char mode){

	for (int r=0; r<row; r++)
	{
		if (iter != -1)
			fprintf(fpt, "%d", iter);
		
	    for(int c=0; c<col; c++)
	         fprintf(fpt, "\t%lf", matrix[r][c]);
	    fprintf(fpt, "\n");
	 }
}

/**
 * @brief      Use WriteParameters() function to write paramters of Reference Set to 
 * individual output file. Output file names is genrated based on the ssParams->inname (input file name)
 * as follow: 'inname_parm_%6d.fout'
 *
 * @param      inp       Input struct to be modified and passed to WriteParamters()
 */
void write_refset_eqparms( SSType *ssParams, Files *files, Input *inp ){

	// Get the pointer to the 'inp' array.
    ParamList *ptab = inp->tra.array;

	/* File name */    
    char *out_fname= (char *) calloc(MAX_RECORD + 1, sizeof(char));
    char *shell_cmd= (char *) calloc(MAX_RECORD + 1, sizeof(char));

    for(int i = 0; i < ssParams->ref_set_size; i++ )
    {
    	// Updating parameters of `inp` with each individual
        for( int h = 0; h < ssParams->nreal; h++ ) 
        {
        	*(ptab[h].param) = ssParams->ref_set->members[i].params[h];
        }

        /* Copy input file to the one we want to change */
        sprintf( out_fname, "%s_ref_%02d", files->outputfile, i );
        sprintf( shell_cmd, "cp -f %s %s", files->outputfile, out_fname );
		if( -1 == system( shell_cmd ) )
	        error( "WriteParameters: error writing output file" );

        WriteParameters(out_fname, &(inp->zyg.parm), "eqparms", 8, inp->zyg.defs);

        // Replacing the seed in the output file with the random seed drawn.
        // Only if seed set to [random]
        sprintf( shell_cmd, "perl -0pi -e 's/$ss\nseed:\n\\[random\\]\n/$ss\nseed:\n%d\n/' %s", ssParams->seed, out_fname );
        if( -1 == system( shell_cmd ) )
	        error( "WriteParameters: could not replace random seed" );
    } //End of for loop
    
    free(out_fname);
    free(shell_cmd);
}

/**
 * @brief	Write header to statistics file.
 */
void write_stats_header( FILE *fp ) {

	fprintf( fp, "# %s %s %s %s %s %s %s %s %s %s\n", 
		"Iterations", 
		"Accumulated_function_evaluations",
		"Min_cost_refset",
		"Average_cost_refset",
		"Var_cost_refset",
		"Replacement_in_reference_set",
		"Local_searches",
		"Flatzones",
		"Duplicates",
		"Candidate_set_size" );
	fflush( fp );
}
