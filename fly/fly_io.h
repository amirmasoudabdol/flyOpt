/**
 * @file fly_io.h
 * @author Damjan Cicin-Sain
 * @contact damjan.cicin@crg.es
 * @date Created on May 27, 2010, 12:07 PM
 *
 * @copyright Copyright (C) 1989-2003 John Reinitz, 2009-2013 Damjan Cicin-Sain, 
 * Anton Crombach and Yogi Jaeger
 * 
 * @brief I/O functions for the SimAnn code.
 */

#ifndef _FLY_IO_H
#define	_FLY_IO_H

/* this def needed for func. defs that refer to (* FILE) */
#include <stdio.h>
#include <ctype.h>

#include "global.h"
#include "ioTools.h"
#include "distributions.h"
#include "score.h"  // Added to resolve the implicit declarion of InitLimits and Penalty2Limits in fly_io.c

#ifdef SS
    #include "ss.h"
#elif ESS
    #include "ess.h"
#endif


#ifdef SS
    SSType ReadSSParameters( FILE * fp, Input *inp );
#elif ESS
    eSSType ReadeSSParameters( FILE * fp, Input *inp );
#endif


/** Optimize by executing scatter search or extended scatter search
 */
void Optimize( DistParms * distp, ScoreOutput * out );


/* A function that reads EqParms from the data file */

/** ReadParamters: reads the parameters for a simulation run from the 
 *                  eqparms or input section of the data file as indicated 
 *                  by the section_title argument and does the conversion  
 *                  of protein half lives into lambda parameters.          
 */
EqParms ReadParameters( FILE * fp, TheProblem defs, char *section_title );

/** ReadParametersX for the input passed by the mex file - needed for the scatter algorithm test */
EqParms ReadParametersX( double *x, TheProblem defs );

/** @brief ReadDivTimes: reead divison times from file */
/** 
     * MITOSIS SCHEDULE: hard-wired cell division tables ***********************                                                                         
     * The division schedule is taken from the $times section in the input file; 
     * If there is none, hard-coded times from maternal.c is taken.
     * 
     * We assume that ndivs is the number of divisions that we want to simulate:
     * 
     * i.e. 
     * 
     * - simulation starts from ccycle 14 means ndivs = 0
     * - simulation starts from ccycle 13 means ndivs = 1
     * - simulation starts from ccycle 12 means ndivs = 2
     * ...
     * 
     * the section entryes are defined as following:
     * 
     * <total_divisions> - this is the maximum number of divisions for that time 
     * schedule. total_divisions >= ndivs
     * 
     * <gastrulation_times> - each entry here corresponds to a gastrulation time
     * assuming that time 0 corresponds to the first division that we are 
     * simulating. 
     * 
     * i.e. 
     * gast_time = 50.0 if starting from ccycle 14 (ndivs = 0)
     * gast_time = 71.1 if starting from ccycle 13 (ndivs = 1)
     * 
     * <division_times> - in each line we have the times of all (total_divisions) 
     * divisions assuming that time 0 corresponds to the start of simulation, and
     * time entries are ordered in descending order of the cell cycles 
     * (cc14, cc13, cc12...).
     * 
     * i.e. 
     * for ndivs = 0 we will read the first line. We start the simulation 
     * from the cycle 14, so we assume that the first simulated cycle (cycle 14) 
     * started at time 0 and no other divisions will happen. So all the times 
     * corresponding to prevous (not simulated) divisions are negative.
     * for ndivs = 1 we will read the second line. We start the simulation 
     * from the cycle 13, so we assume that the first simulated cycle (cycle 13) 
     * started at time 0 and only one other division will happen (cc 13 -> cc 14)
     * at time 21.1. So (again) all the times corresponding to prevous 
     * (not simulated) divisions are negative.
     * 
     *  
     * division_durations - those are the division durations (in seconds) for
     * each division starting from the last one (cc 14).
     *                                                                        
     * Here's an example of the $times section
     *
     * $times
       total_divisions
       6
       gastrulation_times
       50.0 71.1 83.5 93 100.8
       division_times 
       0 -21.1 -33.5 -43.0 -51.8 -57.8
       21.1 0.0 -12.4 -21.9 -30.7 -36.7
       33.5 12.4 0.0 -9.5 -18.3 -24.3
       43.0 21.9  9.5 0.0 -8.8 -14.8
       50.8 29.7 17.3 7.8 -1.0 -7.0
       division_durations
       5.1  3.3  3.0 3.3 3.0 3.0
       $$   
     */   
Times ReadDivTimes( FILE * fp, TheProblem defs );

/* Functions that write or print EqParms */

/** WriteParameters: writes the out_parm struct into a new section in the 
 *                    file specified by filename; the new 'eqparms' sec-   
 *                    tion is inserted right after the 'input' section;    
 *                    to achieve this, we need to write to a temporary     
 *                    file which is then renamed to the output file name   
 *              NOTE: lambdas are converted back into protein half lives!! 
 */
void WriteParameters( char *filename, EqParms * p, char *title, int ndigits, TheProblem defs );

/** PrintParameters: prints an eqparms section with 'title' to the stream 
 *                    indicated by fp                                      
 */
void PrintParameters( FILE * fp, EqParms * p, char *title, int ndigits, TheProblem defs );

/** ReadGuts: reads the $gutsdefs section in a data file into an array 
 *             of strings which then need to get parsed                    
 */
char **ReadGuts( FILE * fp );

/** ReadData: reads in a data or bias section and puts it in a linked 
 *             list of arrays, one line per array; ndp is used to count
 *             the number of data points in a data file (ndp), which is
 *             used to calculate the root mean square (RMS) if required    
 *                                                                         
 *             ReadData allows comments that start with a letter or punc-  
 *             tuation mark, and treats lines starting with a number or    
 *             .num or -.num as data. It expects to read an int and        
 *             ngenes + 1 data points: lineage, time and ngenes protein    
 *             concentrations. It expects data to be in increasing spatial 
 *             and temporal order.                                         
 */
Dlist *ReadData( FILE * fp, char *section, int *ndp, TheProblem * defs );


/** ReadInterpData: reads the history section from a data file into an array 
  *  a dedicated structure
  */
Dlist *ReadInterpData( FILE * fp, char *section, int num_genes, int *ndp );

/** ReadTimes: reads a time table from a file and returns a DArrPtr 
 * FILE FORMAT: one time per line separated with newlines                  
 *        NOTE: max. times tab size is 10000                               
 */
DArrPtr ReadTimes( char *timefile, Zygote zyg );

/** ReadTheProblem: reads the problem section of a data file into the 
 *                   TheProblem struct.                                    
 */
TheProblem ReadTheProblem( FILE * fp );

/** ReadGenotypes: This function reads all the genotypes in a datafile & 
 *                  returns an SList with genotype number and pointers to  
 *                  the corresponding section titles for bias, bcd & facts 
 */
Slist *ReadGenotypes( FILE * fp, int ngenes );

/** ReadBicoid: reads the bcd section of a data file into a linked list; 
 *               also determines maxconc from the bicoid gradient          
 */
Blist *ReadBicoid( FILE * fp, char *section );

/* from score.h
   --------------- */

/**
 * ReadLimits: reads the limits section of a data file and returns the  
 *             approriate SearchSpace struct to the calling function     
 *
 Experimental feature: free limits
 @author Anton Crombach
 @date 2012, August
*/ 
SearchSpace *ReadLimits( FILE * fp, TheProblem defs );


/** ReadTweak: reads the tweak array passed from the calling 
 *              function through the mask array. This array has a value of  
 *              1 or 0 for each parameter in the model and is used by       
 *              Translate to create the array of pointers to the            
 *              parameters-to-be-tweaked. If mask array is NULL, it reads   
 *              that information from the input file                        
 */
Tweak ReadTweak( FILE * fp, int *mask, TheProblem defs );


/**
 * This routine reads the distribution parameters  from the input file
 * and stores them in DistP.xx from distributions.h                   
 * LG 03-02: need q for gen visiting distribution input file          
 * LG 05-02: set factors only dependent on q for general visiting     
 * distribution by calling qgt2_init or qlt2_init from distributions.c
 */
DistParms InitDistribution( FILE * fp );

/* from integrate.h
   ------------------ */

/** PrintBlastoderm: writes the output of the model to a stream specified 
 *                    by the fp file pointer. The Table is a solution of   
 *                    the model as returned by Blastoderm, the id speci-   
 *                    fies the title of the output and ndigits specifies   
 *                    the floating point precision to be printed.          
 *                    PrintBlastoderm adjusts its format automatically to  
 *                    the appropriate number of genes.                     
 */
void PrintBlastoderm( FILE * fp, NArrPtr table, char *id, int ndigits, Zygote * zyg );

/** WriteVersion: prints the version and the complete command line used 
 *                 to run fly_sa into the $version section of the data     
 *                 file                                                    
 */
void WriteVersion( char *filename, char *version, char *argvsave );


/** WriteTimes: writes the timing information to wherever it needs to be 
 *               written to at the end of a run                            
 */
void WriteTime( double time, char *input_fname );

/** PrintTimes: writes a line with the time in seconds */
void PrintTime( FILE * fp, double time );

#endif /* _FLY_IO_H */
