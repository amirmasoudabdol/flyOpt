package GPlot;

=head1 NAME

GPlot - a very simple Perl package for plotting with gnuplot

=head1 SYNOPSIS

    use Lsa::GPlot;
    gplot(%arghash);

=head1 DESCRIPTION

The gplot() function provides an very simple interface between Perl 
and gnuplot.

gplots() plots any array of data of the following format:

<X-coordinate> <Column1> <Column2> ... <ColumnN>

The data can contain comments which must start with anything but
a number. Only lines beginning with a number are plotted (white 
space at the beginning of the line is ignored).

gplot() takes a hash of arguments which contains the following
elements:

    $arghash{"data"}         reference to the data array

    $arghash{"xlim"}         ref to x-axis limits [$xmin, $xmax]
    $arghash{"ylim"}         ref to y-axis limits [$ymin, $ymax]

    $arghash{"title"}        title of the plot (may be empty)
    $arghash{"xlabel"}       label for x-axis (optional)
    $arghash{"ylabel"}       label for y-axis (optional)
    $arghash{"keypos"}       position of the key
    $arghash{"keytitle"}     title line for the key (optional)
    $arghash{"key"}          reference to [key strings] (optional)

    $arghash{"mode"}         output mode, valid modes are: 
                             - disp    display on screen
                             - fig     write to XFig file (.fig)
                             - png     write to PNG format (.png)
    $arghash{"size"}         valid sizes for plots are: small, 
                             medium and large
    $arghash{"plotmode"}     valid pmodes are lines, points or
                             linespoints
    $arghash{"nogrid"}       flag for not plotting the grid

    $arghash{"savefile"}     filename of the file to save plot in,
                             plot_out() appends the appropriate
                             file extension automatically

    $arghash{"dbugflag"}     flag for debugging mode


=head1 SOFTWARE REQUIREMENTS

Obviously, gnuplot must be installed for this to work.

GPlot.pm was developed using gnuplot 4.0. It may or may not work
with other gnuplot versions.

=head1 KNOWN BUGS

The key title is not placed correctly in X11 and other output ter-
minals. I think this is due to a bug introduced in gnuplot 4.0.
I have no idea how to fix this so don't complain.

A long time ago, when GPlot was still Flygraph, I fixed a bug that 
displayed all cycle 14 nuclei displaced by 0.5 percent A-P position.

Please send reports of other bugs to the email address below. I may
or may not fix them.

=head1 AUTHOR

Johannes Jaeger (yoginho@usa.net)

=cut

use strict;
use Exporter;                           
use vars qw(@ISA @EXPORT $VERSION);

use Carp;                               

use File::Basename;                     # for file basename

$VERSION = 0.90;
@ISA = qw(Exporter);

@EXPORT = qw(&gplot);


# plot_out: plots data formatted in unfold output format
# ------------------------------------------------------

sub gplot {

### parse arguments and initalize local variables ##############################

#-- parsing the argument hash --------------------------------------------------

    my(%args)        = %{$_[0]};           # arguments are passed in a hash


    my(@data)        = @{$args{"data"}};     # data array

    my($xmin, $xmax) = @{$args{"xlim"}};   # x-range
    my($ymin, $ymax) = @{$args{"ylim"}};   # y-range

    my($title)       = $args{"title"};     # plot title
    my($xlabel)      = $args{"xlabel"};    # label for x-axis
    my($ylabel)      = $args{"ylabel"};    # label for y-axis
    my($keytitle)    = $args{"keytitle"};  # key title
    my(@key)         = @{$args{"key"}};    # key entries
    my($keypos)      = $args{"keypos"};    # key position

    my($mode)        = $args{"mode"};      # output mode: X11, fig etc.
    my($size)        = $args{"size"};      # plot size: tine to huge
    my($pmode)       = $args{"plotmode"};  # plot mode: lines etc.
    my($nogrid)      = $args{"nogrid"};    # grid or no grid?

    my($sfile)       = $args{"savefile"};  # rel path of file to save

    my($debug)       = $args{"dbugflag"};  # debugging flag
      
#-- other local variables ------------------------------------------------------

    my $i;                                 # local loop counter
    my @line;                              # array used to split lines

    my $gscript;                           # the gnuplot script
    my $dirname;                           # directory in which output is saved
    my($ncolumns, $plotcolumn);            # columns to be plotted

    my($geom, $fn, $fsize);                # X window geometry (size) and font
    my($thick) = $args{"thick"};                             # line thickness for .fig output
    my $term;
    my($xsize, $ysize);                    # plot sizes for .fig & .png output

### error checking #############################################################

    if ( @data == 0 ) {
	croak "data array is empty!"; }

    if ( $mode !~ /disp|eps|epslatex|fig|latex|png/ ) {
	croak "invalid output mode: $mode (use disp, eps, epslatex, fig, latex or png)"; }

    if ( $pmode !~ /lines|points|linespoints/ ) {
	croak "invalid plot mode: $pmode (use lines, points or linespoints)"; }

    if ( $sfile eq "" && $mode =~ /eps|fig|latex|png/ ) {
	croak "save-file name empty (required for eps, epslatex, fig, latex or png mode"; }

    if ( $size !~ /small|medium|large/ ) {
	croak "invalid size: $size (use small, medium or large)"; }

#-- remove leading white space and comment lines --------------------------------

    @data = grep { /^\s*\d/ } @data;
    foreach ( @data ) { s/^\s+//; }

#-- determine number of columns to be plotted -----------------------------------

    $ncolumns = (split /\s+/, $data[0])-1;    # CAUTION: clobbers subroutine args
    if ( $ncolumns < 1 ) { croak "no data columns found for plotting"; }

#-- if writing output file: make directories if necessary -----------------------

    unless ( $mode eq "disp" ) {
	$dirname = dirname($sfile);
	unless ( -e "$dirname" ) { mkdir "$dirname"; }
    }

### construct the gnuplot script ################################################

#-- construct the output-related part of the gnuplot script ---------------------

# case 1: we want output on the screen (terminal = X11)

    if ( $mode eq "disp" ) {
	if    ( $size eq "small" ) { $geom = '192x135'; $fsize = 8; }
	elsif ( $size eq "medium") { $geom = '320x225'; $fsize = 10; }
	elsif ( $size eq "large")  { $geom = '540x430'; $fsize = 12; }
        $fn = "-*-helvetica-medium-r-normal--$fsize-*-*-*-*-*-*-*";

	($gscript = <<END_WRITE) =~ s/^\s+//gm;
	        set output
                set terminal X11
END_WRITE
        
# case 2: we want to save the output as .fig

    } elsif ( $mode eq "fig" ) {
	if ( $size eq "small" ) { 
	    $xsize =  4; $ysize = 2; $thick = 2; $fn = 10; }
	elsif ( $size eq "medium" ) { 
	    $xsize =  6; $ysize = 3; $thick = 2; $fn = 11; }
	elsif ( $size eq "large" ) { 
	    $xsize =  8; $ysize = 4; $thick = 2; $fn = 12; }
  
        ($gscript = <<END_WRITE) =~ s/^\s+//gm;
	  set output "$sfile\.fig"
	  set terminal fig color fontsize $fn thickness $thick \\
	                   inches size $xsize $ysize
END_WRITE

# case 3: we want to save the output as .png

    } elsif ( $mode eq "png" ) {
        if    ( $size eq "small" ) { 
	    $xsize = 240; $ysize = 180; $fsize = "tiny"; }
	elsif ( $size eq "medium") { 
	    $xsize = 350; $ysize = 250; $fsize = "small"; }
	elsif ( $size eq "large") { 
	    $xsize = 560; $ysize = 440; $fsize = "medium"; }
 
	($gscript = <<END_WRITE) =~ s/^\s+//gm;
            set output "$sfile\.png"
        # set terminal png $fsize size $xsize,$ysize \\
        #                      xffffff x000000 x000000 \\
        #                      x6b8e23 x191970 xff0000 xffa500 \\
        #                      xffd700 x66cdaa xffdab9 x9acd32 \\
 	   	  #            x00ced1 xdda0dd xd02090 x000000

          # set term png background "x000000"

# settings for inverted graphs (black background) 
          # set terminal png $fsize size $xsize,$ysize \\
	         #             x000000 xffffff xffffff \\
          #                  xffd700 x9acd32 x00ced1 xff0000 xffdab9 
END_WRITE

# case 4: we want to save the output as .eps
# NOTE: size options will be ignored

    } elsif ( $mode eq "eps" ) {

	$thick = 2;             # line thickness
	$fsize = 18;            # fontsize * 2 (18 will give you fontsize 9)
        $fn = '"Helvetica"';    # font, don't forget to include double quotes

        ($gscript = <<END_WRITE) =~ s/^\s+//gm;
            set output "$sfile\.eps"
            set terminal postscript eps \\
                                    color solid linewidth $thick \\
                                    $fn $fsize
END_WRITE

# case 5: we want to save the output as latex (picture environment)
# NOTE: size options will be ignored
# CAUTION: this mode is UGLY!!! DO NOT USE FOR SERIOUS WORK!

    } elsif ( $mode eq "latex" ) {
 
        ($gscript = <<END_WRITE) =~ s/^\s+//gm;
	  set output "$sfile\.tex"
	  set terminal latex roman 10
END_WRITE

# case 6: we want to save plot as .eps and .tex 
# NOTE: size options will be ignored

    } elsif ( $mode eq "epslatex" ) {

	$fsize = 10;            # fontsize 
        $fn = '"default"';      # font, don't forget to include double quotes

        ($gscript = <<END_WRITE) =~ s/^\s+//gm;
            set output "$sfile\.eps"
            set terminal epslatex \\
                         color solid \\
			 $fn $fsize
END_WRITE

    }

#-- here we create the rest of the  gnuplot script -----------------------------
    
# title or not title?

    unless ( $title eq "" ) { $gscript .= "set title \"$title\"\n"; }

# key or no key?

    if ( @key == 0 && $keytitle eq "" ) {
	$gscript .= "set nokey\n";
    } else {
	$gscript .= "set key $keypos title \"$keytitle\"\n";
    }

# grid or no grid?

    if ( $nogrid ) { $gscript .= "set nogrid\n"; } 
    else           { $gscript .= "set grid\n"; }

# axes labels or no axes labels?

    unless ( $xlabel eq "" ) { $gscript .= "set xlabel \"$xlabel\"\n"; }
    unless ( $ylabel eq "" ) { $gscript .= "set ylabel \"$ylabel\"\n"; }

    if ( $thick eq ".5"){
        $term = "set term png lw .5"
    }else{
        $term = "set term png lw 2"
    }

# more gnuplot options; these are all set to default, change as you see fit

    ($gscript .= <<END_SCRIPT) =~ s/^\s+//gm;
    $term
        set border
	set noclip points
	set clip one
        set noclip two
        set dummy x,y
        set noarrow
        set nologscale
        set nooffsets
	set nopolar
	set noparametric
        set format x "%g"
        set format y "%g"
        set xrange [$xmin:$xmax]
	set yrange [$ymin:$ymax]
END_SCRIPT

    if ( $ncolumns == 1 ) {
	$gscript .= "plot '-' title '$key[0]' with $pmode\n";
    } else {
	$gscript .= "plot '-' title '$key[0]' with $pmode,\\\n";	
	for ( $i=1; $i<($ncolumns-1); $i++ ) {
	    $gscript .= "     '-' title '$key[$i]' with $pmode, \\\n"; }
	$gscript .= "     '-' title '$key[$i]' with $pmode\n";
    }


# attach inline data to the gnuplot script
# lines not beginning with a number are automatically filtered out
# this is all a bit complicated because gnuplot is too stupid to be
# able to use columns in inline data

    for ($i=1; $i<=$ncolumns; $i++) {
    foreach ( @data ) { 
        s/^\s+//; chomp;
        if ( /^\d/ ) { 
        @line = split /\s+/;
        $gscript .= "$line[0] $line[$i]\n";
        }
    }
    $gscript .= "end\n";
    }



    # debug mode
    if ( $debug ) {
	print "gnuplot script:\n\n";
	print "$gscript\n";
    }

# pipe script & data into gnuplot (no more temporary files required)

    if ( $mode eq "disp" ) {
	open ( GNUPLOT, "| gnuplot -geometry $geom -fn $fn -persist") 
	    or croak("could not open connection to gnuplot");  
    } else {
	open ( GNUPLOT, "| gnuplot") 
	    or croak("could not open connection to gnuplot");
    }
    print GNUPLOT $gscript;
    close GNUPLOT;

}

1;
