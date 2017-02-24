/**
 *                                                               
 *   @file global.h                                              
 *                                                               
 *****************************************************************
 *                                                               
 *   written by JR, modified by Yoginho                          
 *                                                               
 *****************************************************************
 *                                                               
 *   THIS HEADER CONTAINS *ONLY* STUFF THAT IS NEEDED IN COST
 *   FUNCTION-SPECIFIC FILES.                     
 *                                                               
 *****************************************************************
 *                                                               
 *   PLEASE THINK TWICE BEFORE PUTTING ANYTHING IN HERE!!!!      
 *                                                               
 *****************************************************************
 *                                                               
 *   global.h gets included automatically by other header files  
 *   if needed and is not included explicitly in .c files.       
 *                                                               
 */


#ifndef GLOBAL_INCLUDED
#define GLOBAL_INCLUDED

#include <float.h>


/*** A global in global for debugging **************************************/

int debug;                      /* debugging flag */
int proc_id;

/*** Constants *************************************************************/

/** IMPORTANT NOTE: the following used to be only 80 to be backward 
 *                compatible with punch cards that were of that maximum length.
*                 We have decided to abandon punch-card compatibility,    
*                 since few people are actually using them anymore...     
*                 It was a tough decision though.      JR & JJ, July 2001 
*/ 
#define MAX_RECORD 1024          /* max. length of lines read from file */

/** The following defines the maximum float precision that is supported by  
 * the code.
 */
extern const int MAX_PRECISION;

/** the following constant as a score tells the annealer to reject a move,  
* no matter what. It had better not be a number that could actually be a  
* score.
*/
extern const double FORBIDDEN_MOVE;     /* the biggest possible score, ever */

/** out of bound control in score */
extern const int OUT_OF_BOUND;  

/** The whole output */
typedef struct ScoreOutput {    
    int size_resid_arr;
    double score;
    double penalty;
    double *residuals;
    double *jacobian;
} ScoreOutput;


#endif
