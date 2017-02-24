package Lsa::Embryo;

use strict;

use Exporter;
use vars qw(@ISA @EXPORT @EXPORT_OK $VERSION);

use Carp;

$VERSION = 9.50;
@ISA = qw(Exporter);

@EXPORT = qw(&convert_lin &convert_linnum);


### embryo functions #######################################################
#   all functions below deal with data format and model independent    
#   aspects of the Drosophila embryo

### @converted_data = convert_lin(@data) ###################################
#   takes an array of data and converts the lineage numbers in the first   #
#   column into percent A-P position                                       #
############################################################################

sub convert_lin {

    my $i;
    my @line;
    my(@data_array)=@{$_[0]};
    
    for ( $i=0; $i<@data_array; $i++ ) {
	$data_array[$i] =~ s/^\s+//;
	unless ( $data_array[$i] =~ /^\d/ ) { next; }
	@line = split /\s+/, $data_array[$i];

	$line[0] = convert_linnum($line[0]);
	$data_array[$i] = join ' ', @line, "\n";    
    }
    return @data_array;
}


### $converted_lin = convert_linnum($lin) ##################################
#   takes a lineage number and converts it into percent A-P position       #
############################################################################

sub convert_linnum {

    my($lin) = $_[0];

# the data points are set such that they are in the same relative position
# for all cell cycles; to achieve this we need to do the following steps:
# - subtract the number of the first lineage number from the current one
# - multiply the result with the relative distance between nucs at that ccycle
# - add an offset to compensate for the fact that we can assume that the
#   daughter nuclei are shifted the same distance anterior and posterior 
#   from their mother; this isn't stated explicitly in the model but makes 
#   the graphs look nicer
# this is the best I could do; it still leads to missing patterns at the an-
# terior border of the graphs whenever we 'loose' a daughter nucleus there 
# after a nuclear division; dunno what I could do about this...

    if ( ($lin >= 8192) && ($lin < 8292) ) {
	$lin = ($lin-8192)+0.5;
    } elsif ( ($lin >= 4096) && ($lin < 4146) ) {
	$lin = ($lin-4096)*2+1;
    } elsif ( ($lin >= 2048) && ($lin < 2073) ) {
	$lin = ($lin-2048)*4+2;
    } elsif ( ($lin >= 1024) && ($lin < 1037) ) {
	$lin = ($lin-1024.)*8+4;
    } elsif ( ($lin >=  512) && ($lin <  519) ) {
	$lin = ($lin-512.)*16+8;
    } else {
	croak "lineage number ($lin) out of range!\n";
    }

    return $lin;
}

1;



