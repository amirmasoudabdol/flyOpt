package Lsa::Fly;

=head1 NAME

Lsa::Fly - a Perl toolbox for fly LSA datafiles

=head1 SYNOPSIS

    use Lsa::Fly;
    ($genes,$geneID,$ndivs,$nnucs,$diff_schedule)= get_problem($datafile);

    use Lsa::Fly;
    $lin = get_lin($datafile);

    use Lsa::Fly;
    @genotypes = get_genotypes($datafile);

    use Lsa::Fly;
    @facts = get_facts($genotype, $datafile);

    use Lsa::Fly;
    @bias = get_bias($genotype, $datafile);

    use Lsa::Fly;
    @bcd  = get_bcd($genotype, $datafile);

    use Lsa::Fly;
    @penalty = get_penalty($genotype, $datafile);

    use Lsa::Fly;
    @mat_penalty = get_maternal_penalty($genotype, $datafile);

    use Lsa::Fly;
    @parm = get_eqparms($datafile);

    use Lsa::Fly;
    @limits = get_limits($datafile);

    use Lsa::Fly;
    @limits = get_penalty_limits($datafile);

    use Lsa::Fly;
    ($mmax, $vmax) = get_max($datafile);

    use Lsa::Fly;
    @divtimes = get_divtimes($ndivs, $oldstyle);

    use Lsa::Fly;
    @duration = get_duration($ndivs, $oldstyle);

    use Lsa::Fly;
    $gasttime = get_gasttime($ndivs, $oldstyle);

    use Lsa::Fly;
    @datatimes = get_datatimes(\@facts);

    use Lsa::Fly;
    @timedata = get_time_data($time,\@facts);

    use Lsa::Fly;
    @genedata = get_gene_data(\@garray,\@facts);

    use Lsa::Fly;
    @cyclebcd = get_cyclebcd(\$ccycle, \@bcd);

    use Lsa::Fly;
    %opts = get_flysa_opts($datafile);

    use Lsa::Fly;
    %out_opts = make_opts($datafile, %in_opts);

    use Lsa::Fly;
    @genes = parse_geneID($geneID);

    use Lsa::Fly;
    $geneID = make_geneID(@genes);

    use Lsa::Fly;
    @key = parse_keyline($keyline);

    use Lsa::Fly; 
    $ccycle = get_ccycle($time, $ndivs, $oldstyle);

    use Lsa::Fly;
    $time = t2time($temp_class, $ndivs);

    use Lsa::Fly;
    $temp_class = time2t($time, $ndivs);

=head1 DESCRIPTION

This module contains functions for retrieving specific information from a 
fly lsa datafile and functions that do tasks related to fly problems (e.g.
converting temporal classes into explicit times and vice versa or parsing 
gene ID strings or lineage numbers). It also contains some useful fly-spe-
cific data, like a list of all the genes for which we have data, all the 
diffusion schedules and the canonical times for temporal classes. We also 
have old- and newstyle division schedules here. If you have to edit any of 
this, edit it here and use the following variables in your Perl scripts:

$Lsa::Fly::genestring     contains a geneID string of all genes in the model
$Lsa::Fly::divisions      contains all valid numbers of cell divisions
$Lsa::Fly::dschedules     contains all valid diffusion schedules as a string

To retrieve values in the other structures mentioned above, use the appro-
priate subroutines described below.

get_problem() reads the problem section in a $datafile and returns an array 
with the problem information: number of genes ($genes), gene ID string 
($geneID), number of nuclear divisions ($ndivs), number of nuclei at clea-
vage cycle 14 ($nnucs) and diffusion schedule ($diff_schedule).

get_lin() returns the lineage number of the first nucleus at cycle 14 in a 
$datafile; this is needed, for instance, for plotting the correct range of 
nuclei for all cell cycles in the various viewer scripts.

#get_genotypes() returns an array with all the genotypes in $datafile.
get_genotypes() returns total number of genotypes in $datafile.

get_facts() returns the facts section for a specific $genotype in $datafile.

get_bias() returns the bias section for a specific $genotype in $datafile.

get_bcd() returns the bcd gradients for a specific $genotype in $datafile.

get_penalty() returns the penalty section for $genotype in $datafile.
CAUTION: get_penalty() simply ignores missing penalty sections and returns
         an empty array!

get_maternal_penalty() returns the maternal penalty section for $genotype
in $datafile.
CAUTION: get_maternal_penalty() simply ignores missing penalty sections and
         returns an empty array!

get_eqparms() returns an array of references to arrays, which contain the 
equation parameters in $datafile (after the annealer has been run). Note 
that some diffusion schedules (such as A and C) only have one diffusion 
parameter, no matter how many genes there are in $datafile. 

get_limits() returns the following array of limits in $datafile:

               Penalty                 Explicit Limits

     [0]       $Lambda                 "N/A"
     [1]       "N/A"                   ref to @Tlim
     [2]       "N/A"                   ref to @mlim
     [3]       "N/A"                   ref to @hlim
     [4]       ref to @dlim            ref to @dlim
     [5]       ref to @lambdalim       ref to @lambdalim

The refereces point to arrays of references which point to ranges of the 
following format: ( $lower_limit, $upper_limit ).

Note that certain diffusion schedules (such as A or C) only have one set of 
limits for the diffusion parameter d.

get_penalty_limits() returns limits just as get_limits() does, but converts
penalty limits into explicit limits. Make sure you do not use this on files
with explicit limits.
CAUTION: currently, this only works properly for sqrt-g(u)!

get_max() returns maximum values for bcd ($mmax) and a reference to an 
array of maximum values for zygotic genes ($vmax). If penalty sections are 
present, they are also used to evaluate maximum values.

get_divtimes() returns the times at which mitosis is finished for the spe-
cified number of cell divsions ($ndivs). If $oldstyle is true, the oldstyle
division schedule is returned (for 3 cell divisions only).

get_duration() returns the durations of the cell divisions for the specified
number of cell divsions ($ndivs). If $oldstyle is true, the oldstyle divi-
sion schedule is returned (for 3 cell divisions only).

get_gasttime() returns the gastrulation time for the specified number of 
cell divsions ($ndivs). If $oldstyle is true, the oldstyle division sche-
dule is returned (for 3 cell divisions only).

get_datatimes() takes an array which contains a @facts section (use 
extract() to extract it from a datafile) and returns all the times for 
which there is data.

get_time_data() takes a $time and an array containing datafile @facts or 
unfold output and returns the data for the indicated time.

get_gene_data() takes a list of gene indices (@garray) to be extracted 
from a data array (@facts) and then produces output with data for only 
these genes. The gene indices correspond to the columns of the genes in 
the @facts array.

get_cyclebcd() returns the bcd gradient for a specific cleavage cycle.

get_flysa_opts() looks for the fly_sa option line in the $version section 
of the $datafile and returns a hash with elements for each option. For 
example: if -o is set, %opts{"o"} will be 1, or if -s e is chosen, 
%opts{"s"} will be set to "e".

make_opts() checks the version section of the datafile for command line
options and adds specific command line options that are passed by the
caller in the %in_opts hash. %in_opts always override options read from
the version section.

parse_geneID() takes a $geneID string and returns an array with all the
corresponding 3-letter gene abbreviations. Valid abbreviations are:

       B       Bicoid*	       Bcd
       C       Caudal	       Cad
       E       Even-skipped    Eve
       F       Fushi tarazu    Ftz
       G       Giant	       Gt
       H       Hunchback       Hb
       I       Hairy           H
       K       Kruppel         Kr
       N       Knirps          Kni
       O       Odd-paired      Odd
       P       Paired          Prd
       Q       Huckebein       Hkb
       R       Runt            Run
       S       Sloppy-paired   Slp
       T       Tailless        Tll
       V       Bicoid*         Bcd

*  V is time-variable Bcd implemented as an external input as opposed to
   B, which is Bcd fixed across each cleavage cycle

The same letters in small represent data for the corresponding genes.

make_geneID() does exactly the reverse of parse_geneID(). It takes an array
of 3-letter gene abbreviations (@genes) and returns the corresponding gene 
ID string.

parse_keyline() parses a guts keyline which must be a section title returned
by unfold -G. Legal abbreviations are all genes listed above, plus:

       A       promoter threshold term in u (parameter h)
       U       u in g(u)
       Z       g(u)
       J       the whole synthesis term, R * g(u)
       L       the decay term, -lambda * v(a)
       X       left diffusion D * (v(a,i-1) - v(a,i))
       Y       right diffusion D * (v(a,i+1) - v(a,i))
       D       the whole derivative, i.e. the whole right-hand side (rhs)

get_ccycle() returns the cleavage cycle number (10-14) for a given $time
and a division schedule (determined by $ndivs and $oldstyle).

t2time() converts a temporal class ($temp_class) to an explicit time (in 
minutes), valid temporal classes are:

       c10     cell cycle 10
       c11     cell cycle 11
       c12     cell cycle 12
       c13     cell cycle 13
       t1      cell cycle 14, temporal class 1
       t2      cell cycle 14, temporal class 2
       t3      cell cycle 14, temporal class 3
       t4      cell cycle 14, temporal class 4
       t5      cell cycle 14, temporal class 5
       t6      cell cycle 14, temporal class 6
       t7      cell cycle 14, temporal class 7
       t8      cell cycle 14, temporal class 8

The second argument of t2time() must be the number of divisions ($ndivs) in
the current division schedule.

time2t() is the reverse function of t2time(). It returns the temporal class 
for a given $ndivs and $time. If it cannot find any time class, it simply 
returns the time given in its argument without changing it in any way.
 
=cut 

use strict;

use Exporter;
use vars qw(@ISA @EXPORT @EXPORT_OK $VERSION);

use Carp;

use Lsa::Data 9.50;

$VERSION = 9.50;
@ISA = qw(Exporter);

@EXPORT = qw(&get_problem &get_genotypes
	     &get_lin
             &get_facts &get_bias &get_bcd
	     &get_penalty &get_maternal_penalty
	     &get_eqparms &get_limits &get_penalty_limits &get_max
	     &get_divtimes &get_duration &get_gasttime
	     &get_datatimes &get_time_data &get_cyclebcd 
             &get_gene_data
	     &get_flysa_opts &make_opts
	     &parse_geneID &make_geneID &parse_keyline
             &get_ccycle &t2time &time2t);

@EXPORT_OK = qw($genestring $divisions $dschedules);

# here be limits to what's allowed in the problem section
# numbers of genes MUST be <100, numbers of nucs <1000

use vars qw($genestring $divisions $dschedules);
$genestring = "ABCEFGHIKNOPQRSTVg";     # valid gene abbrevs (A is param_h)
$divisions  = "01234";                 # valid numbers of cell divisions
$dschedules = "ABCDE";                 # valid diff schedules

# this is the gene ID hash:
# note that the 'A' is used for guts and stands for the parameter h
# promoter threshold and not the gene hairy!
my(%id)=("A" => "Param_h",
         "B" => "Bcd",
         "C" => "Cad",
         "E" => "Eve",
         "F" => "Ftz",
         "G" => "Gt" ,
         "H" => "Hb" ,
         "I" => "H"  ,
         "K" => "Kr" ,
         "N" => "Kni",
         "O" => "Odd",
         "P" => "Prd",
	 "Q" => "Hkb",
         "R" => "Run",
         "S" => "Slp",
         "T" => "Tll",
	 "V" => "Bcd",
	 "g" => "Gt_mRNA");

# here are the cell division schedules and gastrulation times

my(@old_divtimes) = ( 52.0,28.0,12.0 );
my(@divtimes1)    = ( 21.1 );
my(@divtimes2)    = ( 33.5, 12.4 );
my(@divtimes3)    = ( 43.0,21.9, 9.5 );
my(@divtimes4)    = ( 50.8,29.7,17.3, 7.8 );

my(@old_duration) = (  4.0, 4.0, 4.0 );
my(@duration1)    = (  5.1 );
my(@duration2)    = (  5.1, 3.3 );
my(@duration3)    = (  5.1, 3.3, 3.0 );
my(@duration4)    = (  5.1, 3.3, 3.0, 3.3 );

my($old_gasttime) =  88.0;
my($gasttime0)    =  50.0;
my($gasttime1)    =  71.1;
my($gasttime2)    =  83.5;
my($gasttime3)    =  93.0;
my($gasttime4)    = 100.8;

# the following hashes define temporal classes for various div schedules

my(%tc0)=("start" =>  0.000,
	  "t1"    =>  3.125,
	  "t2"    =>  9.275,
	  "t3"    => 15.625,
	  "t4"    => 21.875,
	  "t5"    => 28.125,
	  "t6"    => 34.475,
	  "t7"    => 40.625,
	  "t8"    => 46.875 );

my(%tc1)=("start" =>  0.000,
	  "c13"   => 10.550,
	  "t1"    => 24.225,
	  "t2"    => 30.475,
	  "t3"    => 36.725,
	  "t4"    => 42.975,
	  "t5"    => 49.225,
	  "t6"    => 55.475,
	  "t7"    => 61.725,
	  "t8"    => 67.975 );

my(%tc2)=("start" =>  0.000,
	  "c12"   =>  6.200,
	  "c13"   => 22.950,
	  "t1"    => 36.625,
	  "t2"    => 42.875,
	  "t3"    => 49.125,
	  "t4"    => 55.375,
	  "t5"    => 61.625,
	  "t6"    => 67.875,
	  "t7"    => 74.125,
	  "t8"    => 80.375 );

my(%tc3)=("c11"   =>  0.000,   
          "c12"   => 15.700,   
          "c13"   => 32.450,   
          "t1"    => 46.125,   
          "t2"    => 52.375,   
          "t3"    => 58.625,   
          "t4"    => 64.875,   
          "t5"    => 71.125,   
          "t6"    => 77.375,   
          "t7"    => 83.625,   
          "t8"    => 89.875 );

my(%tc4)=("c10"   =>  0.000,   
          "c11"   => 12.550,   
          "c12"   => 23.500,   
          "c13"   => 40.250,   
          "t1"    => 53.925,   
          "t2"    => 60.175,   
          "t3"    => 66.425,   
          "t4"    => 72.675,   
          "t5"    => 78.925,   
          "t6"    => 85.175,   
          "t7"    => 91.425,   
          "t8"    => 97.675 );


### Fly-specific functions: ################################################
#   All the funcs below will extract some fly specific information from    #
#   a datafile or perform fly-specific tasks, like parsing gene ID strings #
#   or converting temporal classes into numerical time and vice versa. All #
#   of these function do moderate error checking.                          #

### @problem = get_problem($datafile) ######################################
#   returns the elements of the problem structure in $datafile; these ele- #
#   ments are: number of genes ($genes), gene ID string ($geneID), number  #
#   of nuclear divisions ($ndivs), number of nuclei at cleavage cycle 14   #
#   ($nnucs) and the diffusion schedule ($diff_schedule)                   #
############################################################################

sub get_problem {

    my($i, $genes, $egenes, $geneID, $egeneID, $ndivs, $nnucs, $diff_schedule); 
    my($datafile) = $_[0];
    my(@data)     = load_data($datafile);

    if ( ($i=find_section("problem", $datafile)) == -1 ) {
	croak "no problem section in $datafile"; }

    ($genes = $data[$i+2])=~s/\s+//;    
    if ($genes !~ /^\d{1,2}$/) {   
	croak "error reading ngenes ($genes)"; }

    # KLUDGE ALARM: this is ugly and probably not robust!
    # here we have to check whether the datafile has 
    # external inputs or not 
    if ( $data[$i+3] =~ /^geneIDs/ ) {
	($geneID = $data[$i+4])=~s/\s+//;
	if ($geneID !~ /^[$genestring]{1,256}$/) {
	    croak "error reading geneID ($geneID)"; }
	($ndivs = $data[$i+6])=~s/\s+//;
	if ($ndivs !~ /^[$divisions]$/) {
	    croak "error reading ndivs ($ndivs)"; }
	($nnucs = $data[$i+8])=~s/\s+//;
	if ($nnucs !~ /^\d{1,3}$/) {
	    croak "error reading nnucs ($nnucs)"; } 
	($diff_schedule = $data[$i+10])=~s/\s+//;
	if ($diff_schedule !~ /^[$dschedules]$/) {
	    croak "error reading diff_schedule ($diff_schedule)"; }
    } elsif ( $data[$i+3] =~ /^number_of_external_inputs/ ) {
	($egenes = $data[$i+4])=~s/\s+//;    
	if ($egenes !~ /^\d{1,2}$/) {   
	    croak "error reading egenes ($egenes)"; }
	($geneID = $data[$i+6])=~s/\s+//;
	if ($geneID !~ /^[$genestring]{1,256}$/) {
	    croak "error reading geneID ($geneID)"; }
	($egeneID = $data[$i+8])=~s/\s+//;
	if ($egeneID !~ /^[$genestring]{1,256}$/) {
	    croak "error reading egeneID ($egeneID)"; }
	($ndivs = $data[$i+10])=~s/\s+//;
	if ($ndivs !~ /^[$divisions]$/) {
	    croak "error reading ndivs ($ndivs)"; }
	($nnucs = $data[$i+12])=~s/\s+//;
	if ($nnucs !~ /^\d{1,3}$/) {
	    croak "error reading nnucs ($nnucs)"; } 
	($diff_schedule = $data[$i+14])=~s/\s+//;
	if ($diff_schedule !~ /^[$dschedules]$/) {
	    croak "error reading diff_schedule ($diff_schedule)"; }
    } else {
	croak "error reading 3rd line of the problem section"; 
    }

    return($genes, $geneID, $ndivs, $nnucs, $diff_schedule, $egenes, $egeneID);
}

### $lin = get_lin($datafile) ##############################################
#   returns the lineage number of the most anterior nucleus at cleavage    #
#   cycle 14                                                               #
############################################################################

# FIXME: this is based on the BCD section, need to change this!

sub get_lin {

    my($datafile) = $_[0];
    my(@bcd, @line);

    #my $gt = get_genotypes($datafile);
    @bcd = get_bcd(0, $datafile);
    foreach ( @bcd ) {
	s/^\s+//;
	unless ( /^\d/ ) { next; }
	@line = split /\s+/;
	if ( $line[0] >= 8192 ) { last; }
    }

    if ( !(($line[0] >= 8192) && ($line[0] < 8292)) ) {
	croak "error reading lineage numbers (1st num: $line[0]) in get_lin()"; }

    return $line[0];
}

### $lin = get_lin($datafile) ##############################################
#   returns the lineage number of the most anterior nucleus at cleavage    #
#   cycle 14                                                               #
############################################################################

sub get_genotypes {

    my($datafile)=$_[0];
    my(@genotypes, @gtypes, @line, $i);

    @genotypes = extract("genotypes", $datafile);
    
#  The following lines are added/disabled by Damjan Cicin-Sain to make it work with our new genotype identification by order number

    my $num_gen = @genotypes;
    return $num_gen;
#    for ($i=1; $i<(@genotypes-1); $i++) {   # find genotype string
#	if ( $genotypes[$i] =~ /\w/ ) {
#	if ($genotypes)
#	    @line = split /\s+/, $genotypes[$i];
	
	#
	#    if ( @line == 4 ) {             # no external inputs/history
	#	$gtypes[$i-1]=$line[4];
	#   } elsif ( @line == 6 ) {        # external inputs/history
	#	$gtypes[$i-1]=$line[6];
	#   } else {
	#	croak "not a valid format of the genotypes section";
	#   }
#	}
#    }

   # return @gtypes;
}

### @facts = get_facts($genotype, $datafile) ###############################
#   returns the facts section for $genotype in $datafile as an array       #
############################################################################

sub get_facts {

    my($genotype)=$_[0];
    my($datafile)=$_[1];
    my(@genotypes, @line, $factsection, @facts);

    @genotypes = extract("genotypes", $datafile);

    my $size = @genotypes;
    my $i;
    foreach $i (1.. ($size-2)) {   # find bcd section titles
	if ($genotypes[$i] =~ /^\w/ ) { 
	    @line = split(/\s+/, $genotypes[$i]);
	     if ( $genotype == ($i-1) ) {$factsection = $line[1];}

#    foreach ( @genotypes ) {   # find facts section titles
#	if ( /^\w/ ) { 
#	    @line = split /\s+/;
#	    if ( @line == 4 ) {       # no external inputs/history
#		if ( $genotype eq $line[4] ) { $factsection = $line[1]; }
#	    } elsif ( @line == 6 ) {  # external inputs/history
#		if ( $genotype eq $line[6] ) { $factsection = $line[1]; }
#	    }

	}
    }
	
    # go get the facts
    @facts = extract($factsection, $datafile);
    return @facts;
}

### @bias = get_bias($genotype, $datafile) #################################
#   returns the bias section for $genotype in $datafile as an array        #
############################################################################

sub get_bias {

    my($genotype)=$_[0];
    my($datafile)=$_[1];
    my(@genotypes, @line, $biassection, @bias);

    @genotypes = extract("genotypes", $datafile);

    my $size = @genotypes;
    my $i;
    foreach $i (1.. ($size-2)) {   # find bcd section titles
	if ($genotypes[$i] =~ /^\w/ ) { 
	    @line = split(/\s+/, $genotypes[$i]);
	     if ( $genotype == ($i-1) ) {$biassection = $line[0];}

#    foreach ( @genotypes ) {   # find bias section titles
#	if ( /^\w/ ) { 
#	    @line = split /\s+/;
#	    if ( @line == 4 ) {       # no external inputs/history
#		if ( $genotype eq $line[4] ) { $biassection = $line[0]; }
#	    } elsif ( @line == 6 ) {  # external inputs/history
#		if ( $genotype eq $line[6] ) { $biassection = $line[0]; }
#	    }
	}
    }

    # go get the bias
    @bias = extract($biassection, $datafile);
    return @bias;
}

### @bcd = get_bcd($genotype, $datafile) ###################################
#   returns the bcd section for $genotype in $datafile as an array         #
############################################################################

sub get_bcd {

    my($genotype)=$_[0];
    my($datafile)=$_[1];
    my(@genotypes, @line, $bcdsection, @bcd);

    @genotypes = extract("genotypes", $datafile);
	
	my $size = @genotypes;
	my $i;
    	foreach $i (1.. ($size-2)) {   # find bcd section titles
		if ($genotypes[$i] =~ /^\w/ ) { 
		    @line = split(/\s+/, $genotypes[$i]);
	            if ( $genotype == ($i-1) ) {
			$bcdsection = $line[2];
		    }
#	    if ( @line == 4 ) {       # no external inputs/history
#		if ( $genotype eq $line[4] ) { $bcdsection = $line[2]; }
#	    } elsif ( @line == 6 ) {  # external inputs/history
#		if ( $genotype eq $line[6] ) { $bcdsection = $line[2]; }
#	    }
		}
    }
                               # go get the bcd gradients

    @bcd = extract($bcdsection, $datafile);
    return @bcd;
}

### @penalty = get_penalty($genotype, $datafile) ###########################
#   returns the penalty data section for $genotype in $datafile            #
#   CAUTION: get_penalty() will not complain if the penalty section is not #
#            present, it will simply return an empty array                 #
############################################################################

sub get_penalty {

    my($genotype) = $_[0];
    my($datafile) = $_[1];
    my(@penalty,);

    unless ( (find_section("penalty_data.$genotype", $datafile)) == -1 ) {
	@penalty = extract("penalty_data.$genotype", $datafile); }

    return @penalty;
}

### @penalty = get_maternal_penalty($genotype, $datafile) ##################
#   returns the maternal pentaly data section for $genotype in $datafile   #
#   CAUTION: get_penalty() will not complain if the penalty section is not #
#            present, it will simply return an empty array                 #
############################################################################

sub get_maternal_penalty {

    my($genotype) = $_[0];
    my($datafile) = $_[1];
    my(@mat_penalty);

    unless ( (find_section("maternal_penalty_data.$genotype", 
			    $datafile )) == -1 ) {
	@mat_penalty = extract("maternal_penalty_data.$genotype", 
			    $datafile); }

    return @mat_penalty;
}

### @parm = get_eqparms($datafile) #########################################
#   returns an array of references to arrays, which contain the equation   #
#   parameters in $datafile (after the annealer has been run). Note that   #
#   some diffusion schedules (such as A and C) only have one diffusion     #
#   parameter, no matter how many genes there are in $datafile.            #
############################################################################

sub get_eqparms {

    my($datafile) = $_[0];
    my(@eqparms, $ngenes, @parm, @T, $i);

    $ngenes = (get_problem($datafile))[0];

    @eqparms = extract("eqparms", $datafile);
    foreach ( @eqparms ) { s/^\s+//; chomp; }

    $parm[0] = [ split /\s+/, $eqparms[2] ];          # ref to @R
    for ( $i=4; $i<4+$ngenes; $i++ ) {
	$T[$i-4] = [ split /\s+/, $eqparms[$i] ]; }   # refs to @Ts
    $parm[1] = [ @T ];                                # T is 2D
    $parm[2] = [ split /\s+/, $eqparms[ 5+$ngenes] ]; # ref to @m
    $parm[3] = [ split /\s+/, $eqparms[ 7+$ngenes] ]; # ref to @h
    $parm[4] = [ split /\s+/, $eqparms[ 9+$ngenes] ]; # ref to @d
    $parm[5] = [ split /\s+/, $eqparms[11+$ngenes] ]; # ref to @lambda

    return @parm;
}

### @limits = get_limits($datafile) ########################################
#   returns a limits structure (see beginning of this file for details)    #
#   containing the limits in $datafile                                     #
############################################################################

sub get_limits {

    my($datafile) = $_[0];
    my(@limits, $ngenes, @lim, @Tlim, $i);

    $ngenes = (get_problem($datafile))[0];

    @limits = extract("limits", $datafile);
    foreach ( @limits ) { chomp; s/^\s+//; s/[(),]//g; }

    $lim[0] = $limits[2];
    $lim[1] = make_limits($limits[4]);
    if ($lim[0] eq "N/A") {   # only read Tlim, mlim, hlim if no penalty
	for ( $i=6; $i<6+$ngenes; $i++ ) {
	    $Tlim[$i-6] = make_limits($limits[$i]); }
	$lim[2] = [ @Tlim ];
	$lim[3] = make_limits($limits[ 7+$ngenes]);
	$lim[4] = make_limits($limits[ 9+$ngenes]);
	$lim[5] = make_limits($limits[11+$ngenes]);
	$lim[6] = make_limits($limits[13+$ngenes]);
    } else {                  # otherwise these limits are N/A
	$lim[2] = "N/A";
	$lim[3] = "N/A";
	$lim[4] = "N/A";
	$lim[5] = make_limits($limits[12]);
	$lim[6] = make_limits($limits[14]);
    }

    return @lim;
}

### make_limits allocates a limits array for a specific kind of parameter ##
#   and returns a reference to the allocated structure; this is used by    #
#   get_limits above, but should not be visible to the outside world       #
############################################################################

sub make_limits {

    my($line)=$_[0];
    my(@line, @limline, $i);

    @line = split /\s+/, $line;
    for ($i=0;$i<(@line/2);$i++) {
	$limline[$i] = [ $line[2*$i], $line[2*$i+1] ]; }

    return \@limline;
}


### @limits = get_penalty_limits($datafile) ################################
#   returns a limits structure (see beginning of this file for details)    #
#   containing the limits in $datafile; the difference to get_limits() is  #
#   that penalty limits are converted to explicit limits                   #
#   CAUTION: as of yet, this only works for sqrt-g(u)!                     #
############################################################################

sub get_penalty_limits {

    my($datafile) = $_[0];
    my(@limits, $ngenes, $i, $j);
    my($gu_lower, $gu_upper, $x_lower, $x_upper);
    my( $y_lower,  $y_upper, $u_lower, $u_upper);

    $ngenes = (get_problem($datafile))[0]; 

    @limits = get_limits($datafile);

    if ( $limits[0] eq 'N/A' ) { croak "$datafile has explicit limits!"; }

    $gu_lower = $limits[0];     # limits[0] is Lambda       
    $gu_upper = 1 - $limits[0];
    $x_lower  = ( 2 * $gu_lower - 1 );
    $x_upper  = ( 2 * $gu_upper - 1 );
    $y_lower  = sqrt( 1 - $x_lower * $x_lower );
    $y_upper  = sqrt( 1 - $x_upper * $x_upper );
    $u_lower  = $x_lower / $y_lower;
    $u_upper  = $x_upper / $y_upper;
    $u_lower  = $u_lower / sqrt($ngenes);
    $u_upper  = $u_upper / sqrt($ngenes);
    $u_lower  = sprintf("%.3f", $u_lower);  # round floats
    $u_upper  = sprintf("%.3f", $u_upper);

    $limits[2] = [];            # reinitialize pointers: Tlim
    $limits[3] = [];            #                        mlim
    $limits[4] = [];            #                        hlim

    for ($i=0;$i<$ngenes;$i++) {
	for ($j=0;$j<$ngenes;$j++) {
	    $limits[2]->[$i][$j] = [ $u_lower, $u_upper ];
 	}
	$limits[3]->[$i] = [ $u_lower, $u_upper ];
	$limits[4]->[$i] = [ $u_lower, $u_upper ];
    }
    $limits[0] = 'N/A';
    return @limits;
}


### ($mmax, \$vmax) = get_max($datafile) ###################################
#   returns the maximum value of bcd ($mmax) and a reference to an array   #
#   of maximum values for all other genes (\$vmax); includes data and pe-  #
#   nalty sections in its search                                           #
############################################################################

sub get_max {

    my($datafile) = $_[0];
    my($ngenes, $mmax, @vmax, $i, @line);
    my($genotypes, $genotype);
    my(@bcd, @mat_penalty, @facts, @penalty);

    $ngenes = (get_problem($datafile))[0];  # initialize max variables to
    $mmax = -999999;                        # ridiculously small values
    for ($i=0;$i<$ngenes;$i++) { $vmax[$i] = -999999; }

    my $size = get_genotypes($datafile);  # look for maximum values in
    my $j;

    foreach $j (1.. ($size-2)) {   # find bcd section titles
	$genotype = $j-1;
#    foreach $genotype ( @genotypes ) {      # data and penalty sections
	@bcd = get_bcd($genotype, $datafile);
	foreach ( @bcd ) {
	    s/\s+//;
	    if ( /^[0-9.-]/ ) {
		@line = split /\s+/;
		if ( $line[1] > $mmax ) { $mmax = $line[1] }; } }

	@mat_penalty = get_maternal_penalty($genotype, $datafile);
	foreach ( @mat_penalty ) {
	    s/\s+//;
	    if ( /^[0-9.-]/ ) {
		@line = split /\s+/;
		if ( $line[1] > $mmax ) { $mmax = $line[1] }; } }

	@facts = get_facts($genotype, $datafile);
	foreach ( @facts ) {
	    s/\s+//;
	    if ( /^[0-9.-]/ ) {
		@line = split /\s+/;
		for ($i=0;$i<$ngenes;$i++) {    
		    if ( $line[$i+2] > $vmax[$i] ) {
			$vmax[$i] = $line[$i+2]; } } } }

	@penalty = get_penalty($genotype, $datafile);
	foreach ( @penalty ) {
 	    s/\s+//;
	    if ( /^[0-9.-]/ ) {
		@line = split /\s+/;
		for ($i=0;$i<$ngenes;$i++) {    
		    if ( $line[$i+2] > $vmax[$i] ) {
			$vmax[$i] = $line[$i+2]; } } } }
    }

    return ($mmax, \@vmax);
}

### @divtimes = get_divtimes($ndivs, $oldstyle) ############################
#   returns an array of division times for the division schedule defined   #
#   by $ndivs and $oldstyle (e.g. 4-div newstyle)                          #
############################################################################

sub get_divtimes {

    my($ndivs) = $_[0];
    my($old)   = $_[1];

    if ($old) {
	if    ($ndivs==3) { return @old_divtimes; }
	else { croak "only 3 cell divisions allowed for oldstyle"; }
    } else {
	if    ($ndivs==0) { croak "no cell divs for 0 div schedule!"; }
	elsif ($ndivs==1) { return @divtimes1; }
	elsif ($ndivs==2) { return @divtimes2; }
	elsif ($ndivs==3) { return @divtimes3; }
	elsif ($ndivs==4) { return @divtimes4; }
	else { croak "can't handle $ndivs cell divisions"; }
    }
}

### @duration = get_duration($ndivs, $oldstyle) ############################
#   returns an array of division durations for the division schedule defi- #
#   ned by $ndivs and $oldstyle (e.g. 4-div newstyle)                      #
############################################################################

sub get_duration {

    my($ndivs) = $_[0];
    my($old)   = $_[1];

    if ($old) {
	if    ($ndivs==3) { return @old_duration; }
	else { croak "only 3 cell divisions allowed for oldstyle"; }
    } else {
	if    ($ndivs==0) { croak "no cell divs for 0 div schedule!"; }
	elsif ($ndivs==1) { return @duration1; }
	elsif ($ndivs==2) { return @duration2; }
	elsif ($ndivs==3) { return @duration3; }
	elsif ($ndivs==4) { return @duration4; }
	else { croak "can't handle $ndivs cell divisions"; }
    }
}

### $gasttime = get_gasttime($ndivs, $oldstyle) ############################
#   returns the time of gastrulation for the division schedule defined by  #
#   $ndivs and $oldstyle (e.g. 4-div newstyle)                             #
############################################################################

sub get_gasttime {
    my($ndivs) = $_[0];
    my($old)   = $_[1];

    if ($old) {
	if    ($ndivs==3) { return $old_gasttime; }
	else { croak "only 3 cell divisions allowed for oldstyle"; }
    } else {
        if    ($ndivs==0) { return $gasttime0; }
	elsif ($ndivs==1) { return $gasttime1; }
	elsif ($ndivs==2) { return $gasttime2; }
	elsif ($ndivs==3) { return $gasttime3; }
	elsif ($ndivs==4) { return $gasttime4; }
	else { croak "can't handle $ndivs cell divisions"; }
    }
}

### @datatimes = get_datatimes(\@facts) ####################################
#   returns an array of times for which there is data in @facts            #
############################################################################

sub get_datatimes {

    my(@data)    = @{$_[0]};
    my($oldtime) = -1;         # initialized to any impossible time < 0
    my(@times, @line, $count);

    foreach (@data) {
	s/^\s+//;
	if ( /^\d/ ) {         # check for times for which there's data
	    @line = split /\s+/;
	    if ($line[1] != $oldtime) { 
		$times[$count++]=$line[1];
		$oldtime=$line[1]; } } }

    return @times;
}

### @timedata = get_time_data($time, \@facts) ##############################
#   returns data for $time in @facts                                       #
############################################################################

sub get_time_data {

    my($time)    = $_[0];
    my(@in_data) = @{$_[1]};
    my(@out_data, @line);
    my($i)=2;

    $time = sprintf("%.3f", $time);

    foreach (@in_data) {
	s/^\s+//;
	if ( /^\d/ ) {
	    @line = split /\s+/;
	    if ($line[1] == $time) { $out_data[$i++]=$_; }
	}
    }
    if ( scalar(@out_data) > 0 ) {  # if data was found: wrap it up
	$out_data[0]=$in_data[0];
	$out_data[1]="\n";
	$out_data[$#out_data+1]="\n";
	$out_data[$#out_data+1]="\$\$\n";
    }
    return @out_data;
}

### @genedata = get_gene_data(@genes, \@facts) #############################
#   returns data for genes with indices stored in @genes                   #
############################################################################

sub get_gene_data {

    my(@genes)   = @{$_[0]};
    my(@in_data) = @{$_[1]};
    my(@out_data,@line,$newline);
    my($i)=1;
    
    foreach (@in_data) {
	s/^\s+//;
	if ( $_ eq "" ) {
	    $out_data[$i++]="\n";   # keep empty lines
	} elsif ( /^lin/ ) {     
	    @line = split /\s+/;    # extract title
	    $newline = " lin     t ";
	    foreach (@genes) {
		unless ($line[$_+2] eq "") {
		    $newline .= sprintf("     %s", $line[$_+2]);
		} else {
		    carp "no title for gene no. $_";
		}
	    }
	    $out_data[$i++]=$newline . "\n";
	} elsif ( /^\d/ ) {
	    @line = split /\s+/;    # extract data for relevant genes
	    $newline = ' ' . sprintf("%4.d %7.3f", $line[0], $line[1]);
	    foreach (@genes) {
		unless ($line[$_+2] eq "") { 
		    $newline .= sprintf("  %6.2f", $line[$_+2]);
		} else {
		    croak "gene no. $_: data missing";
		}
	    }
	    $out_data[$i++]=$newline . "\n";
	}
    }
    if ( @out_data > 0 ) {  # if data was found: wrap it up
	$out_data[0]="\$output\n";
	$out_data[$#out_data+1]="\$\$\n";
    }
    return @out_data;    
}

### @cyclebcd = get_cyclebcd($ccycle, \@bcd) ###############################
#   returns bcd for a specific cell cycle                                  #
############################################################################

sub get_cyclebcd {

    my($ccycle) = $_[0];
    my(@in_bcd) = @{$_[1]};
    my(@out_bcd, @line, $exp);
    my($i,$j)=0;

    for ($i=0; $i<@in_bcd; $i++) {
	$in_bcd[$i] =~ s/^\s+//;
	if ( $in_bcd[$i] =~ /^\d/ ) {
	    @line = split /\s+/, $in_bcd[$i];
	    if ($line[0] >= (2 ** ($ccycle-1))) { last; }
	}
    }
    
    while ( $in_bcd[$i] =~ /^\d/ ) {
	@line = split /\s+/, $in_bcd[$i++];
	$out_bcd[$j++] = $line[1];
	$in_bcd[$i] =~ s/^\s+//
    }

    return @out_bcd;
}

### %opts = get_flysa_opts($datafile) ######################################
#   returns fly_sa options stored in the $version section of a datafile    #
############################################################################

sub get_flysa_opts {

    my($datafile) = $_[0];
    my(@data) = load_data($datafile);;
    my(@cmd, $cmd, $i, %opts);

    for ( $i=0; $i<=$#data; $i++ ) {
	if ( $data[$i] =~ /^\$version/ ) {
	    $cmd = $data[$i+2];
	    last; } }

    if ( $i>=$#data ) { croak "version could not be found in $datafile"; }

    @cmd = split /\s+/, $cmd;
    unless ($cmd[0] =~ /fly_sa/) {
	croak "fly_sa command line not found in $datafile!"; }
    for ( $i=0; $i<=$#cmd; $i++ ) {
	if ( $cmd[$i] =~ /^\-a/ ) {
	    $opts{"a"}=$cmd[$i+1];
	    $i++;
	} elsif ( $cmd[$i] =~ /^\-b/ ) {
	    $opts{"b"}=$cmd[$i+1];
	    $i++;
	} elsif ( $cmd[$i] =~ /^\-C/ ) {
	    $opts{"C"}=$cmd[$i+1];
	    $i++;
	} elsif ( $cmd[$i] =~ /^\-d/ ) {
	    $opts{"d"}=$cmd[$i+1];
	    $i++;
	} elsif ( $cmd[$i] =~ /^\-e/ ) {
	    $opts{"e"}=$cmd[$i+1];
	    $i++;
	} elsif ( $cmd[$i] =~ /^\-E/ ) {
	    $opts{"E"}=1;
	} elsif ( $cmd[$i] =~ /^\-f/ ) {
	    $opts{"f"}=$cmd[$i+1];
	    $i++;
        } elsif ( $cmd[$i] =~ /^\-g/ ) {
	    $opts{"g"}=$cmd[$i+1];
	    $i++;
	} elsif ( $cmd[$i] =~ /^\-i/ ) {
	    $opts{"i"}=$cmd[$i+1];
	    $i++;
	} elsif ( $cmd[$i] =~ /^\-l/ ) {
	    $opts{"l"}=1;
	} elsif ( $cmd[$i] =~ /^\-L/ ) {
	    $opts{"L"}=1;
	} elsif ( $cmd[$i] =~ /^\-o/ ) {
	    $opts{"o"}=1;
	} elsif ( $cmd[$i] =~ /^\-p/ ) {
	    $opts{"p"}=1;
	} elsif ( $cmd[$i] =~ /^\-Q/ ) {
	    $opts{"Q"}=1;
 	} elsif ( $cmd[$i] =~ /^\-r/ ) {
	    $opts{"r"}=$cmd[$i+1];
	    $i++;
	} elsif ( $cmd[$i] =~ /^\-s/ ) {
	    $opts{"s"}=$cmd[$i+1];
	    $i++; 
	} elsif ( $cmd[$i] =~ /^\-S/ ) {
	    $opts{"S"}=1;
	} elsif ( $cmd[$i] =~ /^\-t/ ) {
	    $opts{"t"}=1;
	} elsif ( $cmd[$i] =~ /^\-T/ ) {
	    $opts{"T"}=1;
	} elsif ( $cmd[$i] =~ /^\-w/ ) {
	    $opts{"w"}=$cmd[$i+1];
	    $i++; 
	} elsif ( $cmd[$i] =~ /^\-W/ ) {
	    $opts{"W"}=$cmd[$i+1];
	    $i++; 
	} elsif ( $cmd[$i] =~ /^\-y/ ) {
	    $opts{"y"}=$cmd[$i+1];
	    $i++; 
	} else {
	    if ( $cmd[$i] =~ /^\-/ ) { carp 
	        "unknown fly_sa opt ($cmd[$i]) in \$version of $datafile"; }
	}
    }

    return %opts;
}

### %out_opts = make_opts($datafile, %in_opts) #############################
#   combines version section options with overriding opts given by caller  #
############################################################################

sub make_opts {

    my($datafile)=$_[0];
    my(%in_opts)=%{$_[1]};
    my(%out_opts);

    # if $version section present: get defaults from there
    if ( -1 != find_section("version", $datafile) ) {
	%out_opts = get_flysa_opts($datafile); }

    # if explicit v options set -> override defaults
    if ( $in_opts{"a"} ne "N/A" ) { $out_opts{"a"}=$in_opts{"a"}; } 
    if ( $in_opts{"g"} ne "N/A" ) { $out_opts{"g"}=$in_opts{"g"}; } 
    if ( $in_opts{"i"} ne "N/A" ) { $out_opts{"i"}=$in_opts{"i"}; } 
    if ( $in_opts{"q"} ne "N/A" ) { $out_opts{"q"}=$in_opts{"q"}; }
    if ( $in_opts{"r"} ne "N/A" ) { $out_opts{"r"}=$in_opts{"r"}; }
    if ( $in_opts{"s"} ne "N/A" ) { $out_opts{"s"}=$in_opts{"s"}; }
    if ( $in_opts{"z"} ne "N/A" ) { $out_opts{"z"}=$in_opts{"z"}; }
 
    return %out_opts;
}

### @genes = parse_geneID($geneID) #########################################
#   returns an array of 3-letter gene names after parsing $geneID          #
############################################################################

sub parse_geneID {

    my($id)=$_[0];                      
    my($pos, $char, @chars, @genes); 

    @chars = split //, $id;         # split geneID string into characters
    $pos = 0;
    foreach $char ( @chars ) {
	unless ( $char =~ /^[$genestring]$/io ) {
	    croak "$char is not a valid 1-letter gene abbreviation"; }
	$genes[$pos++] = $id{$char};
    }

    return @genes;
}

### $geneID = make_geneID(@genes) ##########################################
#   returns a geneID according to the @genes array of 3-letter abbrevs     #
#   CAUTION: the whole geneID string will be in CAPS!                      #
############################################################################

sub make_geneID {

    my(@genes)=@{$_[0]};
    my($geneID, %myid);

    %myid = reverse(%id);
    foreach (@genes) {
	if ( defined $myid{$_} ) {
	    $geneID .= $myid{$_};
	} else {
	    croak "$_ is not a valid 3-letter gene abbreviation";
	}
    }

    return $geneID;
}

### @key = parse_keyline($keyline) #########################################
#   parses a guts keyline (which must be a section title line returned by  #
#   unfold -G); returns an array with the guts keys                        #
############################################################################

sub parse_keyline {

    my($line)=$_[0];
    my(@key, @line, @word, $i, $j);

    @line = split /\s+/, $line;
    for ( $i=1; $i<@line; $i++ ) {	
	@word = split //, $line[$i];
	if ( "ABCEFGHIKNOPQRSTV" =~ /$word[0]/ ) {
	    $key[$i-1] = (parse_geneID($word[0]))[0];
	    for ( $j=1; $j<@word; $j++ ) {
		$key[$i-1] .= ' + ';
		$key[$i-1] .= (parse_geneID($word[$j]))[0];
	    }
	} elsif ( "UZJLXYD" =~ /$word[0]/ ) {
	    if ( @word > 2 ) { die "gutv: error parsing guts title\n"; }
	    if ( $line[$i] eq "U" ) {
		$key[$i-1] = "u";
	    } elsif ( $line[$i] eq "Z" ) {
		$key[$i-1] = "g(u)";
	    } elsif ( $line[$i] eq "J" ) {
		$key[$i-1] = "tot. sythesis";
	    } elsif ( $line[$i] eq "L" ) {
		$key[$i-1] = "tot. decay";
	    } elsif ( $line[$i] eq "X" ) {
		$key[$i-1] = "left transport";
	    } elsif ( $line[$i] eq "Y" ) {
		$key[$i-1] = "right transport";
	    } elsif ( $line[$i] eq "XY" ) {
		$key[$i-1] = "tot. transport";
	    } elsif ( $line[$i] eq "D" ) {
		$key[$i-1] = "dv/dt";
	    }
	} else {
	    croak "unknown abbreviation ($line[$i]).\n";
	}
    }
    return @key;
}

### $ccycle = get_ccycle($time, $ndivs, $oldstyle) #########################
#   returns the cleavage cycle for $time and a given division schedule     #
############################################################################

sub get_ccycle {
    my($time)  = $_[0];
    my($ndivs) = $_[1];
    my($old)   = $_[2];

    if ($old) {
	if ($ndivs==3) {
	    if    ( ($time>$old_divtimes[0]) ) { return 14; }
	    elsif ( ($time>$old_divtimes[1]) ) { return 13; }
	    elsif ( ($time>$old_divtimes[2]) ) { return 12; }
	    else                               { return 11; }
	} else { croak "only 3 cell divisions allowed for oldstyle"; }
    } else {
	if ($ndivs==0) { 
	    return 14;
	} elsif ($ndivs==1) { 
	    if    ( ($time>$divtimes3[0]) ) { return 14; }
	    else                            { return 13; }
	} elsif ($ndivs==3) { 
	    if    ( ($time>$divtimes3[0]) ) { return 14; }
	    elsif ( ($time>$divtimes3[1]) ) { return 13; }
	    else                            { return 12; }
	} elsif ($ndivs==3) { 
	    if    ( ($time>$divtimes3[0]) ) { return 14; }
	    elsif ( ($time>$divtimes3[1]) ) { return 13; }
	    elsif ( ($time>$divtimes3[2]) ) { return 12; }
	    else                            { return 11; }
	} elsif ($ndivs==4) { 
	    if    ( ($time>$divtimes4[0]) ) { return 14; }
	    elsif ( ($time>$divtimes4[1]) ) { return 13; }
	    elsif ( ($time>$divtimes4[2]) ) { return 12; }
	    elsif ( ($time>$divtimes4[3]) ) { return 11; }
	    else                            { return 10; }
	} else { croak "can't handle $ndivs cell divisions"; }
    }
}

### $time = t2time($temp_class) ############################################
#   returns the time which corresponds to a temporal class (c10-13, t1-8)  #
############################################################################

sub t2time {

    my %tc;
    my($tin)  =$_[0];
    my($ndivs)=$_[1];

    if ( $ndivs == 0 ) {
	%tc=%tc0;
    } elsif ( $ndivs == 1 ) {
	%tc=%tc1;
    } elsif ( $ndivs == 2 ) {
	%tc=%tc2;
    } elsif ( $ndivs == 3 ) {
	%tc=%tc3;
    } elsif ( $ndivs == 4 ) {
	%tc=%tc4;
    } else { 
	croak "$ndivs not a valid number of cell divisions"; 
    }

    unless ( $tc{$tin} eq "" ) {
	return $tc{$tin};    
    } else {
	croak "no time class $tin for $ndivs divisions";
    }
} 

### $temp_class = time2t($time) ############################################
#   returns the temporal class (c10-13, t1-8) which corresponds to a time  #
############################################################################

sub time2t {

    my %tc;
    my($timein)=$_[0];
    my($ndivs) =$_[1];

    if ( $ndivs == 0 ) {
	%tc=reverse %tc0;
    } elsif ( $ndivs == 1 ) {
	%tc=reverse %tc1;
    } elsif ( $ndivs == 2 ) {
	%tc=reverse %tc2;
    } elsif ( $ndivs == 3 ) {
	%tc=reverse %tc3;
    } elsif ( $ndivs == 4 ) {
	%tc=reverse %tc4;
    } else { 
	croak "$ndivs not a valid number of cell divisions"; 
    }

    unless ( $tc{$timein} eq "" ) {
	return $tc{$timein};
    } else {
	return $timein;
    }
}

1;
