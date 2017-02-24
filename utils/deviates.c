/**
 * @file deviates.c  
 * Created 5-02 by: Lorraine Greenwald
 * Main wrapper for program that generates basic deviates for plotting and
 * sanity checks.  This program uses exact files and routines from the 
 * simulated annealing code.  Files distributions.c and distributions.h which 
 * are stored in the lam area on CVS. 
 *
 * takes distribution type, q value and theta_bar on command line
 * USAGE: 
 *  gen_deviates distribution qvalue theta_bar dist_deviate_output_file_name
 **************************************************************************
 *                                                               
 * Copyright (C) 1989-2003 John Reinitz                          
 * the full GPL copyright notice can be found in lsa.c          
 *                                                              
 **************************************************************************/

#ifdef ICC
#include <mathimf.h>
#else
#include <math.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include "distributions.h"      /* problem independent distributions */
/* need qgt2_init and qlt2_init prototypes */
#ifndef ERROR_INCLUDED
#include "error.h"              /* use fly_team standard error routine */
#endif
#include "random.h"             /* prototype for InitRand call */

const int BYTESIZE = 8;
const int LOWBITS = 0x330E;     /* for drand (32 bits) to erand (48 bits) compatibility */


// GLOBAL CONSTANTS

/** The following defines the maximum float precision that is supported by  
 * the code.
 */
const int MAX_PRECISION = 16;
                                                                  
/* the following constant as a score tells the annealer to reject a move,  */
/* no matter what. It had better not be a number that could actually be a  */
/* score.                                                                  */
const double FORBIDDEN_MOVE = DBL_MAX;  /* the biggest possible score, ever */

const int OUT_OF_BOUND = -1;

/**
 * This is main for the deviate generation code    
 * LG 05-02                                        */
int
main( int argc, char **argv ) { /* begin MAIN deviate */

    DistParms DistP;

    double theta_deviate;
    //unsigned short xsubj[3];    /* used for InitRand call */
    /* NOTE: xsubj gets changed with each call to RandomReal */

    /**************************************************
     * these are for the deviate distributions
     ***************************************************/
    double theta_bar = 1.;      /* mean for deviates accepted on command line input */

    FILE *outptr;               /*deviate output file */
    FILE *logptr;               /*log file */

    int i;
    int number_of_trials = 1000000;
    long seedval;
//    unsigned short left16, middle16;
//    int left;

    /*******************
     * initialization:  *
     *******************/

    if( argc != 5 ) {
        error( "gen_deviates: <dist:1=exp,2=uni,3=absnor,4=abslor,5=lor2,6=poi,7=gen,8=stdnorm,9=pareto,10=nor> <1<qvalue<3> <theta_bar> <deviates outputfile>\n" );
        exit( 1 );
    }

    /* distribution = argv[1]; */

    DistP.distribution = atoi( argv[1] );
    printf( "distribution type is %d \n", DistP.distribution );
    if( ( DistP.distribution > 10 ) || ( DistP.distribution < 1 ) ) {
        error( "gen_deviates: distribution must be int from 1 to 10 \n" );
    } else if( ( DistP.distribution == 5 ) || ( DistP.distribution == 8 ) || ( DistP.distribution == 10 ) ) {
        printf( "gen_deviates: distributions return negative values \n" );
    }

    /* q value = argv[2];  */

    DistP.q = atof( argv[2] );
    printf( "q value is %lf \n", DistP.q );

    /* theta_bar value = argv[3];  */

    theta_bar = atof( argv[3] );
    printf( "theta_bar value is %lf \n", theta_bar );


    logptr = stdout;
    outptr = fopen( argv[4], "w" );
    if( outptr == NULL ) {
        error( "gen_deviates: Error in  opening output file\n" );
        exit( 1 );
    }


    /* initialize the random number generator, now dSFMT() */

    seedval = 514804963;
/*
    xsubj[0] = LOWBITS;
    middle16 = ( unsigned short ) seedval;
    xsubj[1] = middle16;
    left = seedval >> ( BYTESIZE * sizeof( unsigned short ) );
    left16 = ( unsigned short ) left;
    xsubj[2] = left16;
*/
    //InitRand( xsubj );         /* makes the xsubj array static to random.c */
    InitRand( seedval );

    /* general visiting intialization */

    if( DistP.distribution == 7 ) {
        if( ( DistP.q >= 3.0 ) || ( DistP.q <= 1.0 ) ) {
            error( "gen_deviates: q must be between 1 and 3 \n" );
        } else if( DistP.q > 2.0 ) {
            qgt2_init( &DistP );
            /*        print_qgt2_visit(theta_bar);   * debug 1-13-03 */
        } else {                /* less than or equal to 2 */


            qlt2_init( &DistP );
        }
    }

    /* deviate generation loop */

    for( i = 0; i < number_of_trials; i++ ) {
        theta_deviate = generate_dev( theta_bar, &DistP );
        fprintf( outptr, " %g\n", theta_deviate );
        fflush( outptr );
    }

    /* printf("When q= %lf and trunc= %lf,  %16.15f moves are rejected.\n", DistP.q ,DistP.trunc,DistP.rejects);debug 1-13-03 */

    /* have done all generations here */


    fclose( logptr );
    fclose( outptr );
    return 0;

}                               /* end MAIN deviate */
