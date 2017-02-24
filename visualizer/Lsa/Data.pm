package Lsa::Data;

=head1 NAME

Lsa::Data - a Perl toolbox for LSA datafile handling

=head1 SYNOPSIS

    use Lsa::Data;
    @data = load_data($infile);

    use Lsa::Data;
    write_data($outfile, \@data);

    use Lsa::Data;
    get_filelist($orig_datafile);

    use Lsa::Data;
    @section = extract($section_title, $datafile);

    use Lsa::Data;
    replace($section_title, \@replacement_section, $datafile);

    use Lsa::Data;
    $i = find_section($section, $datafile);

    use Lsa::Data;
    $i = kill_section($section, $datafile);

=head1 DESCRIPTION

This module provides basic functions for reading, writing and manipulating
datafiles which are in the dataformat used in the Reinitz lab for all things
related to Lam simulated annealing (LSA). The functions in this module are 
problem-independent and can be used on any file that sticks to the basic data 
format rules defined in the doc/dataformatX.X documents that come with the 
code distributions..

load_data() loads the contents of the $infile into the @data array.

write_data() writes an array containing Lsa @data to the $outfile.

get_filelist() returns a sorted list of datafiles derived from $orig_datafile
according to the scheme $orig_filename_### where ### is some number

extract() extracts a section with $section_title from a $datafile.

replace() opens the $datafile and replaces the section with $section_title 
with a section stored in the @section array.

find_section() finds a section with $section_titlel in a $datafile and 
returns the line number $i of the $section_title (first line in $datafile 
is 0); if the section cannot be found, find_section() returns -1.

kill_section() removes $section from $datafile.

=cut 

use strict;
use vars qw(@ISA @EXPORT $VERSION);

use Carp;
use Exporter;

use File::Basename;                     # funcs to parse pathnames

$VERSION = 9.50;
@ISA = qw(Exporter);

@EXPORT = qw(&load_data &write_data 
             &get_filelist
             &extract &replace
             &find_section
             &kill_section
             &strip_time);

### @data = load_data($infile): ############################################
#   reads an lsa datafile and returns an array with its contents           #
############################################################################

sub load_data {

    my(@data);                    # data array to be returned 
    my($infile) = $_[0];          # name of input datafile

    open(IN, $infile) or croak "could not open $infile: $!";
    @data = <IN>;
    close(IN);

    return @data;
}

### write_data($outfile, \@data): ############****##########################
#   writes an array which contains lsa data to a data file                 #
############################################################################

sub write_data {

    my($outfile)=  $_[0];         # name of output data file
    my(@data)   =@{$_[1]};        # data array to be written to $outfile

    open(OUT, ">$outfile") or croak "could not write $outfile: $!";
    foreach (@data) { print OUT; }
    close(OUT);
}

### @stripped_data = strip_time(\@data) #################################### 
#   strips data of the time (second) column                                #
############################################################################

sub strip_time {

    my(@indata)=@{$_[0]};
    my(@inline,$outline);
    my @outdata;
    
    foreach ( @indata ) {
	s/^\s+//; chomp;
	if ( /^\d/ ) {
	    @inline = split /\s+/;
	    $outline = join ' ', @inline[0,2..$#inline];
	} else {
	    $outline = $_;
	}
	push @outdata, "$outline\n";
    }
	    
    return @outdata;
}

### @filelist = get_filelist($orig_datafile): ##############################
#   finds all datafiles derived from an original datafile in a directory   #
#   according to the naming scheme $orig_datafile_###, where ### is some   #
#   arbitrary number                                                       #
#   returns a sorted list of full paths to these files                     #
############################################################################

sub get_filelist {

    my($infile) = $_[0];                 # original datafile
    my(@filelist,$file,$base,$dir);      

    $base = basename($infile);           # parse path name
    $dir  = dirname($infile); 

    opendir(DIR, $dir) || die("Couldn't open $dir!\n");;
    @filelist = grep { /${base}_\d+\z/ } readdir(DIR);
    closedir(DIR);

    @filelist = sort @filelist;

    foreach $file ( @filelist ) {        # we want full file paths
	$file = $dir . "/" . $file;
    }    
    return @filelist;
}    

### @section = extract($section_title, $datafile): #########################
#   reads a specific section of a datafile and returns it as an array      #
#   returns 0 if replace failed, 1 otherwise                               # 
############################################################################

sub extract {

    my($section)  = $_[0];              # name of section to be extracted 
    my($datafile) = $_[1];
    my(@in_dat) = load_data($datafile); # array with contents of $datafile
    my($i, $j, @out_dat);               # @out_dat: array to be returned

    if ( ($i=find_section($section, $datafile)) == -1 ) {
	croak "section \$$section not found in $datafile";
    }

    while ($in_dat[$i] !~ /\$\$/) {
	$out_dat[$j++]=$in_dat[$i++];
    }
    $out_dat[$j]=$in_dat[$i];           # don't forget the terminating $$

    return @out_dat;
}
    
### $i = $replace($section_title, \@replacement_section, $datafile): #######
#   replaces a section of a datafile with the contents of an array;        #
#   returns 0 if replace failed, 1 otherwise                               # 
############################################################################

sub replace {

    my($section)  =   $_[0];            # section to be replaced
    my(@section)  = @{$_[1]};           # replacment section array
    my($datafile) =   $_[2];
    my(@data) = load_data($datafile);   # content of $datafile
    my($pos, $length);                  # pos and length of $section

    if ( ($pos=find_section($section, $datafile)) == -1 ) {
	return 0;
    }

    $length = (extract($section, $datafile));
    splice @data, $pos, $length, @section;
    write_data($datafile, \@data);

    return 1;
}

### $i = find_section($section, $datafile): ################################
#   finds a section in a datafile; returns the line number of the section  #
#   in the datafile (first line is 0) or -1 if it could not be found       #
############################################################################

sub find_section {

    my($section) = $_[0];               # section to be found 
    my(@data)    = load_data($_[1]);    # contents of $datafile
    my $i;

    for ($i=0;$i<=$#data; $i++) {
	if ($data[$i] =~ /^\$$section\s/) { return $i; }
    }

    return -1;                          # if not found -> return -1
}

### $i = kill_section($section, $datafile): ################################
#   deletes $section from $datafile; returns 0 if error, 1 otherwise       #
############################################################################

sub kill_section {

    my($section)  =   $_[0];            # section to be replaced
    my($datafile) =   $_[1];
    my(@data) = load_data($datafile);   # content of $datafile
    my($pos, $length);                  # pos and length of $section

    if ( ($pos=find_section($section, $datafile)) == -1 ) {
	return 0;
    }

    $length = (extract($section, $datafile));
    splice @data, $pos, $length, ();
    if ( $data[$pos] !~ /^\$/ ) {     # remove superfluous whitespace
	splice @data, $pos, 1, (); }
    write_data($datafile, \@data);

    return 1;
}

1;




