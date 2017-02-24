/**
 * @file fly.c                                                
 * @author JR, modified by Yoginho,                            
 *   -D landscape option by Lorraine Greenwald in Oct 2002,            
 *   -g option by Yousong Wang in Feb 2002,        
 *   -a option by Marcel Wolf in Apr 2002                          
 *
 * @copyright Copyright (C) 1989-2003 John Reinitz, 2009-2013 Damjan Cicin-Sain, 
 * Anton Crombach and Yogi Jaeger
 * 
 * It contains most of its problem specific code
 *
 */


#include <time.h>
#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>             /* for command line option stuff */

/*======
    Utils*/
#include "error.h"              /* error handling funcs */
#include "distributions.h"      /* DistP.variables and prototypes */
#include "integrate.h"
#include "random.h"             /* for InitRand() */

/*======
    Solver*/
#include "maternal.h"           /* for olddivstyle and such */
#include "score.h"              /* for init and Score funcs */
#include "solvers.h"            /* for name of solver funcs */
#include "zygotic.h"            /* for init, mutators and derivative funcs */
#include "fly_io.h"

/*=================
    Scatter Search
    */
#ifdef SS
    #include "ss.h"
#elif ESS
    #include "ess.h"
#endif

/*========
    defining the algorithm specific paramters...
    */
#ifdef SS
    SSType ssParams;
#endif

#ifdef ESS
    eSSType essParams;
#endif

// ScoreOutput which will be passed around for calculation
ScoreOutput out;    

// Variables moved from lsa.c
static Files files;
DistParms distp;                /*variable set by InitDistribution() */

/*** Constants *************************************************************/

/* command line option string */
const char *OPTS = ":a:b:Bc:C:De:Ef:g:hi:lLm:nNopQr:s:StTvw:W:y:";
/* D will be debug, like scramble, score */
/* must start with :, option with argument must have a : following */


/*** STATIC VARIABLES ******************************************************/

/* Help, usage and version messages */
static const char usage[] =
    "Usage: fly_X [-a <accuracy>] [-b <bkup_freq>] [-B] [-e <freeze_crit>] [-E]\n"
    "              [-f <param_prec>] [-g <g(u)>] [-h] [-i <stepsize>] [-l] [-L] \n"
    "              [-m <score_method>] [-n] [-N] [-p] [-Q] [-s <solver>] [-t] [-v]\n"
    "              [-w <out_file>] [-y <log_freq>]\n" "              <datafile>\n";

static const char help[] =
    "Usage: fly_X [options] <datafile>\n\n"
    "Argument:\n"
    "  <datafile>          input data file\n\n"
    "Options:\n"
    "  -a <accuracy>       solver accuracy for adaptive stepsize ODE solvers\n"
    "  -b <bkup_freq>      write state file every <bkup_freq> * tau moves\n" "  -B                  run in benchmark mode (only do fixed initial steps)\n"
    "  -D                  debugging mode, prints all kinds of debugging info\n"
    "  -e <freeze_crit>    set annealing freeze criterion to <freeze_crit>\n"
    "  -E                  run in equilibration mode\n"
    "  -f <param_prec>     float precision of parameters is <param_prec>\n"
    "  -g <g(u)>           chooses g(u): e = exp, h = hvs, s = sqrt, t = tanh\n"
    "  -h                  prints this help message\n"
    "  -i <stepsize>       sets ODE solver stepsize (in minutes)\n" "  -l                  echo log to the terminal\n"
    "  -m <score_method>   w = wls, o=ols score calculation method\n"
    "  -n                  nofile: don't print .log or .state files\n"
    "  -N                  generates landscape to .landscape file in equilibrate mode \n"
    "  -o                  use oldstyle cell division times (3 div only)\n" "  -p                  prints move acceptance stats to .prolix file\n"
    "  -s <solver>         choose ODE solver\n"
    "  -v                  print version and compilation date\n" "  -w <out_file>       write output to <out_file> instead of <datafile>\n"
    "  -y <log_freq>       write log every <log_freq> * tau moves\n\n" "Please report bugs to <yoginho@usa.net>. Thank you!\n";

static char version[MAX_RECORD];        /* version gets set below */
static char *argvsave;          /* static string for saving command line */

static double stepsize = 1.;    /* stepsize for solver */
static double accuracy = 0.001; /* accuracy for solver (not used yet) */
static int precision = 8;       /* precision for eqparms */
static int method = 0;          /* 0 for wls, 1 for ols */

// static int prolix_flag = 0;     /* to prolix or not to prolix */
// static int landscape_flag = 0;  /* generate energy landscape data */

/* the whole input - static not to read data from file at every loop */
static Input inp;


/*** FUNCTION POINTERS *****************************************************/

/* the following lines define a pointers to:                               */
/*            - pd:    dvdt function, Dvdt_sqrt or DvdtOrig in zygotic.c   */
/*            - pj:    Jacobian function, in zygotic.c                     */
/* This pointer needs to be static since both InitialMove and Restore-     */
/* State need it; it gets set in ParseCommandLine()                        */
/*                                                                         */
/* NOTE: ps (solver) is declared as global in integrate.h                  */
void ( *pd ) ( double *, double, double *, int, SolverInput *, Input * );
void ( *pj ) ( double, double *, double *, double **, int, SolverInput *, Input * );


/*** FUNCTIONS *************************************************************/

/*** COMMAND LINE OPTS ARE ALL PARSED HERE *********************************/

/**  ParseCommandLine: well, parses the command line and returns an index  
 *                     to the 1st argument after the command line options  
 */
int
ParseCommandLine( int argc, char **argv ) {

    int c, i;                   /* used to parse command line options */

    /* external declarations for command line option parsing (unistd.h) */
    extern char *optarg;        /* command line option argument */
    extern int optind;          /* pointer to current element of argv */
    extern int optopt;          /* contain option character upon error */
    
#ifdef SS
    sprintf( version, "fly_ss for journal Computation" );
#elif ESS
    sprintf( version, "fly_ess for journal Computation" );
#endif

    /* set default values for command line options */
    pd = DvdtOrig;              /* default derivative function */
    //pd = Dvdt_sqrt;
    pj = JacobnOrig;            /* default Jacobian function */
    ps = Rkck;                  /* default solver */
    dd = DvdtDelay;             /* delayed derivative fnuction */

    /* parse command line for options and their arguments */
    optarg = NULL;
    while( ( c = getopt( argc, argv, OPTS ) ) != -1 ) {
        switch ( c ) {
        case 'a':
            accuracy = atof( optarg );
            if( accuracy <= 0 )
                error( "fly_X: accuracy (%g) is too small", accuracy );
            break;
        case 'D':
            debug = 1;
            break;
        case 'f':
            precision = atoi( optarg ); /* -f determines float precision */
            if( precision < 0 )
                error( "fly_X: what exactly would a negative precision be???" );
            if( precision > MAX_PRECISION )
                error( "fly_X: max. float precision is %d!", MAX_PRECISION );
            break;
        case 'g':              /* -g choose g(u) function */
            pd = DvdtOrig;      //
            if( !( strcmp( optarg, "s" ) ) )
                gofu = Sqrt;
            else if( !( strcmp( optarg, "t" ) ) )
                gofu = Tanh;
            else if( !( strcmp( optarg, "e" ) ) )
                gofu = Exp;
            else if( !( strcmp( optarg, "h" ) ) )
                gofu = Hvs;
            else if( !( strcmp( optarg, "k" ) ) ) {
                gofu = Kolja;
            } else
                error( "fly_X: %s is an invalid g(u), should be e, h, s or t", optarg );
            break;
        case 'h':              /* -h help option */
            PrintMsg( help, 0 );
            break;
        case 'i':              /* -i sets the stepsize */
            stepsize = atof( optarg );
            if( stepsize < 0 )
                error( "fly_X: going backwards? (hint: check your -i)" );
            if( stepsize == 0 )
                error( "fly_X: going nowhere? (hint: check your -i)" );
            if( stepsize > MAX_STEPSIZE )
                error( "fly_X: stepsize %g too large (max. is %g)", stepsize, MAX_STEPSIZE );
            break;
        case 'm':              /* -m sets the score method: w for wls, o for ols */
            if( !( strcmp( optarg, "w" ) ) )
                method = 0;
            else if( !( strcmp( optarg, "o" ) ) )
                method = 1;
            break;
        case 'o':              /* -o sets old division style (ndivs = 3 only! ) */
            olddivstyle = 1;
            break;
        case 's':              /* -s sets solver to be used */
            if( !( strcmp( optarg, "a" ) ) )
                ps = Adams;
            else if( !( strcmp( optarg, "bd" ) ) )
                ps = BaDe;
            else if( !( strcmp( optarg, "bs" ) ) )
                ps = BuSt;
            else if( !( strcmp( optarg, "e" ) ) )
                ps = Euler;
            else if( !( strcmp( optarg, "h" ) ) )
                ps = Heun;
            else if( !( strcmp( optarg, "mi" ) ) || !( strcmp( optarg, "m" ) ) )
                ps = Milne;
            else if( !( strcmp( optarg, "me" ) ) )
                ps = Meuler;
            else if( !( strcmp( optarg, "r4" ) ) || !( strcmp( optarg, "r" ) ) )
                ps = Rk4;
            else if( !( strcmp( optarg, "r2" ) ) )
                ps = Rk2;
            else if( !( strcmp( optarg, "rck" ) ) )
                ps = Rkck;
            else if( !( strcmp( optarg, "rf" ) ) )
                ps = Rkf;
            else if( !( strcmp( optarg, "sd" ) ) )
                ps = SoDe;
            else if( !( strcmp( optarg, "kr" ) ) )
                ps = Krylov;
            /* else if (!(strcmp(optarg, "bnd")))
               ps = Band; */
            else
                error( "fly_X: invalid solver (%s), use: a,bs,e,h,kr,mi,me,r{2,4,ck,f}", optarg );
            break;
        case 'v':              /* -v prints version message */
            fprintf( stderr, "%s\n", version );
            exit( 0 );
        case 'w':              /* -w sets output file */
            files.outputfile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
            files.outputfile = strcpy( files.outputfile, optarg );
            // SetOutname( outname );      /* communicates outname to lsa.c */
            break;
        case ':':
            error( "fly_X: need an argument for option -%c", optopt );
            break;
        case '?':
        default:
            error( "fly_X: unrecognized option -%c", optopt );
        }
    }

    /* error checking here */
    if( ( ( argc - ( optind - 1 ) ) != 2 ) )
        PrintMsg( usage, 1 );

    argvsave = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    for( i = 0; i < argc; i++ ) {
        if( i > 0 )
            argvsave = strcat( argvsave, " " );
        argvsave = strcat( argvsave, argv[i] );
    }
    argvsave = strcat( argvsave, "\n" );

    return optind;
}

/*
    Initializing optimization specific variables and call the select optimization procedure.
*/
void
Optimize( DistParms *distp, ScoreOutput *out ) {

    /* all that info in various subsections of the structure 'inp' */
    FILE *infile;
    /* solver log file where Blastoderm function writes (debugging mode) */
    FILE *slogfile = NULL;

    // in case there's no -w, outname = inname
    if( !files.outputfile ) {
        files.outputfile = ( char * ) calloc( MAX_RECORD + 1, sizeof( char ) );
        files.outputfile = strcpy( files.outputfile, files.inputfile );
    }

    infile = fopen( files.inputfile, "r" );
    if( !infile )
        file_error( "fly_X error opening input file" );
    if( debug ) {
        slogfile = fopen( strcat( files.inputfile, ".slog" ), "w" );
        if( !slogfile )
            file_error( "fly_X error opening slog file" );
    }

    inp.zyg = InitZygote( infile, pd, pj, &inp, "input" );
    inp.sco = InitScoring( infile, method, &inp );
    inp.his = InitHistory( infile, &inp );  //It fills the polations vector
    inp.ext = InitExternalInputs( infile, &inp );
    inp.ste = InitStepsize( stepsize, accuracy, slogfile, files.inputfile );
    /* read the list of parameters to be tweaked */
    inp.twe = InitTweak( infile, NULL, inp.zyg.defs );
    inp.tra = Translate( &inp );

    /* reading optimization algorithm specific paramters */
    #ifdef SS
        ssParams = ReadSSParameters(infile, &inp);
    #elif defined(ESS)
        init_defaultSettings(&essParams);
        essParams = ReadeSSParameters(infile, &inp);
    #endif        

    /* input file read, copy parameters */
    fclose( infile );
    inp.lparm = CopyParm( inp.zyg.parm, &( inp.zyg.defs ) );
    /* write out command line to version string */
    WriteVersion( files.outputfile, version, argvsave );

    /* do actual optimization */
    #ifdef SS
        InitSS(&inp, &ssParams, &files);
        RunSS(&inp, &ssParams, &files);
    #elif defined(ESS)
        init_eSS( &essParams, &inp, out, files.inputfile);
        run_eSS( &essParams, &inp, out, files.inputfile);
    #endif

    /* Clean up */
    FreeMutant( inp.lparm );
}


/** 
    Starting point of the optimization algorithm. 

    Init output structure, read input file, and optimize.
*/
int
main( int argc, char **argv ) {

    /* timing the program */
    clock_t tic = clock();

#ifdef SS
    printf("Starting Scatter Search Algorithm\n");
#elif defined(ESS)
    printf("Starting Enhanced Scatter Search Algorithm\n");
#endif

    debug = 0;
    /* pointer to current argument of command line */
    int opt_index;

    /* 
        initialize output structure, score starts with large number that 
        we want to minimize
    */
    out.score = 1e38;
    out.penalty = 0;
    out.size_resid_arr = 0;
    out.jacobian = NULL;
    out.residuals = NULL;

    /* read input file */
    files.inputfile = ( char * ) calloc( MAX_RECORD, sizeof( char ) );
    opt_index = ParseCommandLine( argc, argv );
    strcpy( files.inputfile, argv[opt_index] );

    /* and now the magic happens! */
    Optimize(&distp, &out);

    /* how much time did it take? Re-using the path+input file */
    clock_t toc = clock();
    WriteTime( (double)(toc - tic) / CLOCKS_PER_SEC, files.inputfile );

    return 0;
}
