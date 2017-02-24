/**
 * @file distributions.h                                       
 * created 5-02 for function prototypes found in distributions.c   
 * used in tsp_sa.c (for qgt2_init and qlt2_init), and move.c      
 *           need pi for distributions                             
 *           used in input file to set factors, in distributions   */

#ifndef DISTRIBUTIONS_INCLUDED
#define DISTRIBUTIONS_INCLUDED

typedef struct my_distrib {
    int distribution;           /* move generation distribution type RO */
    /* 1 - exp; 2 - uni; 3 - absnor; 4 - abs lorentz */
    /*  LG: 07-05-00 formerly dist_type in lj code */
    double q;                   /* gen visiting distribution parameter RO */
    /* 1=guassian; 2=lor; but 1<q<3 */
    /*  LG: 03-02 need q and factors for GSA visit dist */
    /*  1   <  q < 2   uses qlt2_visit       */
    /*  2   <  q < 2.6 uses qgt2_visit       */
    /*  2.6 <= q < 3   uses binom_qgt2_visit */

    /*****variables that depend on q ********/
    double gam1;                /* gammaln of 1/(q-1)       */
    double gam2;                /* gammaln of 1/((q-1)-0.5) */
    double fact2;               /* exp(gam1-gam2)           */
    double alpha;               /* sqrt((q-1)/pi) * exp(gam1-gam2) NOT dependent on theta_bar */
    double alpha2;              /* (sqrt(2.)/ 2. )* sqrt((q-1)/pi) * exp(gam1-gam2) NOT dependent on theta_bar */
    double c;                   /*  1/(q-1)                  */
    double rejects;             /*  number of recursive calls; useful for efficiency analysis */
    double trunc;               /*  range for x is [-trunc, trunc]             */
    /*  2    <  q < 2.6   trunc is OK at 2 million */
    /*  2.6  <= q < 2.85  trunc is 99 and 16 zeros */
    /*  2.85 <= q < 3     trunc is 24 nines */

    /*****variables that depend on theta_bar ********/



} DistParms;


/* prototype functions for distributions */

/** gasdev returns a normally distributed 
 * deviate with zero mean and unit variance
 * use RandomReal() rather than ran1
 * taken from Numerical Recipes in C p.289
 */
double gasdev( void );   

/**
 * gammln returns the ln(gamma(xx)) for xx > 0 
 * uses ln because without it too many overflows 
 * taken from Numerical Recipes in C p.214
 */
double gammln( double xx );     

/** poidev returns as a float an integer value 
 * that is a random deviate from a Poisson distribution
 * with mean = xm.  Uses RandomReal() rather than ran1
 * uses ln because without it too many overflows 
 * taken from Numerical Recipes in C p.294
 */
double poidev( double xm );     

/** binom_qgt2_visit generates random deviates for 2.6 <= q < 3 
 * for Tsallis GSA q-distribution using the REJECTION METHOD. 
 * when q is 2.6 or bigger, the tails of the distribution are
 * significant, so it becomes necessary to use a larger 
 * range for x (the variable trunc).  When q >2.85, an even  
 * larger range for x is needed. This value is set in qgt2_init.  
 * These large values computationally challenge the computer, 
 * so the first two terms of the binomial approximation are used
 * in the calculation of kappa and krappa.  
 * This function may be called recursively because a 
 * deviate value must be returned.  Due to theoretic limit  
 * the comparison function is the limit of genvis at q=3. 
 * for values of q between 1 and 2, use qlt2_visit 
 * for values of q=2 use lorentz, q=1 use normal
 * for values of 2< q <2.6 use qgt2_visit 
 * for values of 2.6<= q < 3 use binom_qgt2_visit (this routine) 
 * trunc must be set to a larger number for q > 2.8 (see qgt2_init)
 */
double binom_qgt2_visit( double theta_bar, DistParms * DistP ); /* 2.6 <= q < 3 */

/** qgt2_visit generates random deviates for q between 2 and 2.6 
 * for Tsallis GSA q-distribution using the REJECTION METHOD. 
 * This function may be called recursively because a 
 * deviate value must be returned.  Due to theoretic limit  
 * the comparison function is the limit of genvis at q=3. 
 * for values of q between 1 and 2, use qlt2_visit 
 * for values of q=2 use lorentz, q=1 use normal
 * for values of 2< q <2.6 use qgt2_visit (this routine) 
 * for values of 2.6<= q < 3 use binom_qgt2_visit 
 */
double qgt2_visit( double theta_bar, DistParms * DistP );       /* 2 < q < 2.6 */

/** qlt2_visit generates random deviates for values of 
 * q between 1 and 2,  distributed with Tsallis
 * GSA q-distribution using the REJECTION METHOD. 
 * This function may be called recursively because a deviate 
 * value must be returned.  
 * This uses a lorentzian bounding function.  
 * for values of q between 1 and 2, use qlt2_visit (this routine)
 * for q between 2 and 3 use qgt2_visit
 * for values of q=2 use lorentz, q=1 use normal
 * for values of 2< q <2.6 use qgt2_visit  
 * for values of 2.6<= q < 3 use binom_qgt2_visit 
 */
double qlt2_visit( double theta_bar, DistParms * DistP );       /* q between 1 and 2 */

/** print_qgt2_visit created for calculation debugging */
void print_qgt2_visit( double theta_bar, DistParms * DistP );   /* debug 1-13-03 */

/** generate_dev generates basic deviates for any distribution  
 * The distribution type and theta_bar are input as parameters 
 * deviate is the value returned.                              
 */
double generate_dev( double theta_bar, DistParms * DistP );

/**qgt2_init initializes the factors used for the  general 
 * visiting distribution for q between 2 and 3, but are    
 * not dependent on theta_bar. These get calculated only 
 * once during initialization of the SA run and remain     
 * during the entire run.                                  
 */
void qgt2_init( DistParms * DistP );    /* this routine in distributions.c */

/** qlt2_init initializes the factors used for the  general 
 * visiting distribution for q between 1 and 2, but are    
 * not dependent on theta_bar. These get calculated only 
 * once during initialization of the SA run and remain     
 * during the entire run.                                  
 */
void qlt2_init( DistParms * DistP );    /* this routine in distributions.c */
#endif
