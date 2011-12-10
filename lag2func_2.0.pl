#!/usr/bin/env perl -w
use strict;
use Cwd qw(cwd abs_path realpath);
use File::Temp qw(tempfile);
use Getopt::Long;
use List::Util qw(min max);

################################################################################
# Name and Version
################################################################################
#
# lag2func.pl version 2.0
#
################################################################################
# Licence: CC-BY-NC-SA
################################################################################
#
# This script was written by Timo Ruppell (timo[dot]ruppell[at]gmail[dot]com)
# and can be found at:
#
# http://www.helsinki.fi/~ruppel/code/
#
# It is distributed under the Creative Commons Attribution, Noncommercial, and
# ShareAlike licence. This means You can copy, edit, modify, etc. to Your hearts 
# content as long as You give me credit, don't make money and share Your work in
# kind.
#
# http://creativecommons.org/
# http://en.wikipedia.org/wiki/Creative_Commons_licenses
#
################################################################################
# Description
################################################################################
#
# This is a small perl script for transforming lanhep/calchep files to increase
# performance by substituting excessively long "Lorentz" parts in lgrng.mdl
# files with placeholders and evaluating these parts in func.mdl instead.
#
# The way this script works is fairly simple:
#  - Lines are read in from the lgrng.mdl file.
#  - The lines are passed to Mathematica for collation and extraction of terms
#    proportional to (real,imaginary)x(scalar,G5,G(3),p_mu,m_i).
#  - A placeholder for each of these collections of terms in the form of A01,
#    B01 etc is made, where the numeral corresponds to the original line in the
#    lgrng.mdl file.
#  - The long expression in the lgrng.mdl file is replaced by, e.g., A01+B01*i
#  - The func.mdl file is expanded to include the placeholders and the long
#    expressions they stand for.
#  - In lines proportional to only a single Lorentz factor the placeholder A##  
#    is moved into the "Factor" part of lgrng.mdl 
#
# The reason for doing this is that objects in func.mdl are evaluated only once
# and the ensuing simplification in the calculation of the amplitudes from
# the vertices in lgrng.mdl speeds up the overall calculation of more complex
# models significantly. 
#
# The formatting required by the func.mdl file places some restrictions on this
# operation:
#   1) The "Expression" field is truncated at around 2000 characters. This may
#      be simply a convention which could be changed in the source code of
#      calchep. In any case we work around it by simply splitting up longer
#      substitutions.
#   2) The "Name" field is exactly 7 characters long. Using the file line number
#      and upper case letters to uniquely identify a substitution means that
#      we can only accomodate
#      - 999.999 lines at lengths < ~50000 characters (26*1950) or
#      -  99.999 lines at lengths < ~1.3M characters (26^2*1950)
# Even if it were the case that these are not simply the result of a coding
# convention in calchep, it does not seem likely for either of these limits to
# apply.
#
# In addition to the above "limitations" there are a couple of caveats:
#   1) This script uses Mathematica, so problems with incompatible versions
#      may appear. This script was developed and tested using Mathematica 7
#      but it should run on older versions since we are not doing anything
#      fancy.
#   2) If the lanhep/calchep files contain negative exponents e.g. X**-2 or X^-2
#      this script will malfunction. I may get around to fixing this but as it
#      is a very unlikely occurrence (X*Y**-2 is written as X/Y**2) the added
#      complexity is not worth it right now.
#   3) For some reason Mathematica in some instances leaves parentheses in place
#      after ComplexExpand. If these are in of form (X*Y*Z) they are removed,
#      but if they are of form (X+Y) this script will abort.
#   4) You can overwrite Your files with this script. Default behavior is to
#      avoid and always check, but for convenience if e.g. both input and output
#      files are given it is assumed that overwriting is intended.
#
# Substitution of ** --> ^ is performed on the input lgrng.mdl file's Lorentz
# part to avoid confusion with Mathematica's NonCommutativeMultiply. In fact,
# by default, this script also parses all other parts of the lgrng and func
# files and makes that substitution.
#
# With respect to the coding: It is very easy to make perl code completely
# unreadable, but faced with either learning pyhton or just commenting a lot
# and trying to keep things simple, I chose the latter. If you want this script
# without the metric ton of comments use perltidy --dsc or run something like
# cat lag2func.pl | sed -E 's/[ ]+(\# |\(\*).+//g.
#
# If You do release modified / new code based on this work I'd love know!
# Needless to say, if you find bugs or think of specific improvements,
# please don't hesitate to contact me.
#
# ++Timo Ruppell
#
################################################################################
# Change log:
################################################################################
# 
# v2.0: -Multiple coefficients of the same Lorentz factor get combined to reduce
#        lgrng.mdl "Lorentz" part's complexity even more. 
#       -Coefficients of single Lorentz factors get put in the "Factor" part of
#        lgrng.mdl.
#       -All Lorentz factors are now taken into consideration (v1.0 only sorted
#        scalar and G5)
#       -The test of the Mathematica kernel is commented out and need to be 
#        redone to test all the different Lorentz factors.
#
################################################################################
# To Do:
################################################################################
#
# -variable naming conventions are all over the place due to hare brained coding
#  process, needs streamlining :)
# -better help file
#
################################################################################
# Acknowledgements:
################################################################################
#
# I got the idea for this script and lifted a couple of code snippets from
# Daniel Reeves' MASH (Mathematica Scripting Hack) script.
# Daniel's Site: http://ai.eecs.umich.edu/people/dreeves/
#
################################################################################

################################################################################
# Set up short usage notification and some more detailed instructions
#
# (This bit needs some more work but suffices for now)
#
################################################################################

sub erm {                                                                       # Error notification to be used instead of die, which is too verbose
    my $name = ( $0 =~ m/([^\/]+)$/ ) ? $1 : $0;                                # Match the name of this script in a (possibly) full path
    print "$_[1]\n" if defined $_[1];                                           # Print error message if supplied
    print <<"USE"   if ( $_[0] > 0 );                                           # Provide small blurb on error code larger than 0

$1 replaces the Lorenz parts in calchep/lanhep lgrng.mdl files with
placeholders and moves them into the func.mdl files.

USE
    print <<"USE" if ( $_[0] >= 0 );                                            # Print short usage guide on non-negative error codes
Usage: $name in_file [out_file] [OPTIONS]
       $name -m model_number [OPTIONS]
USE
    print <<"USE" if ( $_[0] == 1 );                                            # Provide detailed explanations on error code 1
 
--------------------------------------------------------------------------------
 Option                Explanation
--------------------------------------------------------------------------------
    -c  --convert      Convert ** -> ^ in all parts of the func/lgrng file.
                       This is on by default, use --noconvert to disable.
    -d  --digits=n     Force n number of digits (2-5) for line ID's
    -h  --help         Print this help text.
    -i  --input=dir    Input directory, defaults to pwd.
    -l  --length=n     Maximum allowed length of "Lorentz" part. This script 
                       will skip all lines shorter than this. Defult is n=100.
    -m  --model=n      Model file number, e.g., 13 for lgrng13.mdl. Superceded 
                       by input / output files if they are given. If both input
                       file and -m are given this will allow overwriting.
    -o  --output=dir   Output directory, defaults to pwd.
    -v  --verbose      Be more verbose durng operation
    -C  --copy         Make copies of the other .mdl files, not just lgrng and 
                       func, if the input/output model numbers are different
    -I  --inplace      Force overwriting of input .mdl files, -o and output 
                       files will be ignored.           
    -L  --long         Allow very long lines (~1.3 million characters) in the
                       lgrng.mdl file. Due to formatting constraints this works
                       only if there are less than 100.000 lines total in the 
                       lgrng.mdl file. 
    -M  --math=d/b     Accepts either the directory d where the Mathematica 
                       kernel binary is, or the kernel binary b itself. For 
                       OSX the path to Mathematica.app is sufficient.
    -X  --overwrite    Allow the overwriting of existing files. If an output
                       file is given, this option is also switched on.
                       
--------------------------------------------------------------------------------
 Examples
--------------------------------------------------------------------------------

        $name lgrng10.mdl
        
    Will look for a file lgrng10.mdl and func10.mdl in the current directory and
    attempt to make new ones named lgrng11.mdl and func11.mdl but abort if they 
    already exist.

        $name -Im 11
        
    Will look for lgrng11.mdl and func11.mdl and overwrite them inplace.
    
        $name lgrng12.mdl lgrng99.mdl --copy
        
    Will look for lgrng12.mdl and func12.mdl and create / overwrite lgrng99.mdl
    and func99.mdl if they already exist. Also, will attempt to copy 
    extlib12.mdl, vars12.mdl and prtcls12.mdl to extlib99.mdl etc.

        $name -m 13 -i NMSSM/ -o NMSSM-NEW/ -CX -l 50 --math /usr/opt/math
        
    Model no. 13 from NMSSM/ to NMSSM-NEW/ with all associated files, allowing 
    overwriting of existing files, chopping the Lorentz part at 50 characters
    and using a custom location for the Math kernel.
USE
    exit;
}

################################################################################
# Parse command line arguments and take immediate actions (-s -h)
################################################################################

Getopt::Long::Configure("bundling");                                            # Configure Getopt

my $pathIn  = cwd();                                                            # Input dir, defaults to pwd
my $pathOut = cwd();                                                            # Output dir, defaults to pwd
my $model;                                                                      # Model number
my $inplace;                                                                    # Flag for overwriting model files in place
my $chopL   = 100;                                                              # Max allowed length for Lorentx parts
my $digits  = 2;                                                                # Number of digits in line ID's, defaults to 2 (this is also the minium)
my $usrMath = "";                                                               # User specified math kernel
my $overwrite;                                                                  # Flag to allow overwriting
my $help;                                                                       # Flag to disply usage tutotrial
my $long;                                                                       # Make line ID's start with two letters insted of one
my $copy;                                                                       # Flag to copy over associated model files as well
my $verbose;                                                                    # Flag to be more verbose
my $l2c = 1;                                                                    # Flag to make ** -> ^ changes in the model files, default is yes

GetOptions(
    "i|input=s"   => \$pathIn,                                                  # Get the options from the command line
    "o|output=s"  => \$pathOut,
    "m|model=i"   => \$model,
    "l|length=i"  => \$chopL,
    "v|verbose"   => \$verbose,
    "C|copy"      => \$copy,
    "I|inplace"   => \$inplace,
    "M|math=s"    => \$usrMath,
    "L|long"      => \$long,
    "X|overwrite" => \$overwrite,
    "h|help"      => \$help,
    "d|digits=i"  => \$digits,
    "c|convert!"  => \$l2c
) or erm(1);                                                                    # Print help on errors
if ($help) { erm(1) }                                                           # Print help on -h/--help

################################################################################
# Preliminary checks on the command line arguments
################################################################################

if ( not $model and not $ARGV[0] ) { erm( 0, "Input?" ) }                       # User needs to provide at least a model no. or input file
if ( $digits > 6 ) { erm( -1, "Can only force up to 6 digit line ID's." ) }     # This is due to the func.mdl file formatting
if ( $digits == 6 and $long ) { erm( -1, "-d 6 and -L are incompatible." ) }    # Check that --long isn't used when it shouldn't
if ( -d $usrMath ) {                                                            # If the user specifies a directory
    $usrMath =~ s/(\w)\s*$/$1\//;                                               #    add slash if missing and
    if ( $^O =~ /darwin/ ) {                                                    #    depending on the OS, put in a best
        $usrMath .= "Contents/MacOS/MathKernel";                                #    guess for the kernel name
    }
    else {
        $usrMath .= "math";
    }
    if ( not -e $usrMath ) {
        erm( -1, "User defined Math kenrel missing: $usrMath $!" );
    }
}
if ( not -d "$pathIn" ) { erm( -1, "Can't find model directory $pathIn: $!" ) } # Check input dir
if ( not -d "$pathOut" ) { mkdir $pathOut }                                     # Check output dir

################################################################################
# Set up the input and output files
################################################################################

my @outputF;                                                                    # The func.mdl output aggregator
my @outputL;                                                                    # The lgrng.mdl output aggregator
my $modelIn;                                                                    # The number of the input model file
my $modelOut;                                                                   # The number of the output model file

my @IOF = (@ARGV) ? map( abs_path($_), @ARGV ) : undef;                         # Any command line arguments left will be assumed to be input/output files

if ( $IOF[0] and $IOF[0] =~ m/^(.*)lgrng([1-9][0-9]*)\.mdl/ ) {                 # If an input file is given pick up the model no. and path from it.
    $pathIn  = $1;
    $modelIn = $2;
}
elsif ( $IOF[0] ) {
    erm( -1, "$ARGV[0] is not a valid input file name." );
}
if ( $IOF[1] and $IOF[1] =~ m/^(.*)(lgrng|func)([1-9][0-9]*)\.mdl/ ) {          # If an output file is given pick up the model no. and path from it and
    $overwrite = 1;                                                             #   also assume that overwriting is OK
    $pathOut   = $1;
    $modelOut  = $3;
}
elsif ( $IOF[1] ) {
    erm( -1, "$ARGV[1] is not a valid output file name." );
}

$pathIn  = realpath($pathIn)."/";                                               # Clean up the in/out paths (e.g. collapse x/../y -> y) and append a
$pathOut = realpath($pathOut)."/";                                              # trailing slash for convenient concatenation later

if ( not $modelIn ) {                                                           # If no input file is supplied by the user and
    if ( $pathIn =~ $pathOut ) {                                                #   the input and output paths are the same
        $modelIn  = $model;                                                     #       set input model no. to user supplied model no.
        $modelOut = $model + 1;
    }                                                                           #       and the output model no. to one higher
    else {                                                                      #   the input and output paths are different
        $modelOut = $modelIn = $model;                                          #       set output model no. and the input model no. to user supplied one
    }
}
elsif ( not $modelOut and not $model ) {                                        # If only an input file is supplied and
    if ( $pathIn =~ $pathOut ) {                                                #   the input and output paths are the same
        $modelOut = $modelIn + 1;
    }                                                                           #       set the output model no. to the input model no. + 1
    else {                                                                      #   the input and output paths are different
        $modelOut = $modelIn;                                                   #       set output model no. to the input model no.
    }
}
elsif ( not $modelOut and $model ) {                                            # If only an input _file_ but also a model no. is supplied and
    $modelOut = $model;                                                         #   set output model no. to the user supplied model no.
}

my $inL  = $pathIn  . "lgrng" . $modelIn  . ".mdl";                             # Lagrangian input file name
my $outL = $pathOut . "lgrng" . $modelOut . ".mdl";                             # Lagrangian output file name
my $inF  = $pathIn  . "func"  . $modelIn  . ".mdl";                             # Func input file name
my $outF = $pathOut . "func"  . $modelOut . ".mdl";                             # Func output file name

if ($inplace) {                                                                 # Scratch all the previous if the user specifies in-place substitution
    $outL      = $inL;
    $outF      = $inF;
    $overwrite = 1;
}

print <<"FILES" if ($verbose);
Input Lagrangian file:  $inL
Output Lagrangian file: $outL
Input Functions file:   $inF
Output Functions file:  $outF
FILES

if ( not -r "$inF" ) { erm( -1, "Can't read/locate model file $inF: $!" ) }
if ( not -r "$inL" ) { erm( -1, "Can't read/locate model file $inL: $!" ) }
if ( not $overwrite ) {                                                         # Make sure not to overwrite files. Works only if files are read/write!!
    if ( "$outL" eq "$inL" ) { erm( -1, "$inL would be overwritten!" ) }
    if ( "$outF" eq "$inF" ) { erm( -1, "$inF would be overwritten!" ) }
    if ( -e "$outL" ) { erm( -1, "Output file $outL already exists!" ) }
    if ( -e "$outF" ) { erm( -1, "Output file $outF already exists!" ) }
}

################################################################################
# Set up the Mathematica kernel, script and callMath() subroutine
################################################################################

my ( $scrfh, $script ) = tempfile( "math-XXXX", SUFFIX => '.m', UNLINK => 1 );  # Temp file to store the Mathematica
print $scrfh <<'MATH';                                                          #    bit that does the grunt work
SetOptions[$Output,PageWidth->Infinity];                                        (* PageWidth set to infinite     *)
i=I;                                                                            (* to avoid Mathematica's weird  *)
line=ToExpression[Drop[$CommandLine, 4][[1]]];                                  (* Fortran line break fromatting *)
cLine = ComplexExpand@{Re@line, Im@line};                                       (* Split into re and im parts    *)
factors = {G5, G3, G3 G5,                                                       (* List of factorts to expect    *)
           m1p2, m1p3, m1p4, m2p1, m2p3, m2p4,                                  (* must be same list as later    *)
           m3p1, m3p2, m3p4, m4p1, m4p2, m4p3, 
           m1m2, m1m3, m1m4, m2m3, m2m4, m3m4};       
cList = Join[{Fold[Coefficient[#1, #2, 0] &, cLine, factors]},                  (* Serious Fold[] mojo to pick   *)
              Fold[Coefficient[#1, #2, 0] &,                                    (* out only coeffs of factors    *)
                   Coefficient[cLine, #, 1] & /@ factors, factors]];            (* and not any cross terms       *)
cSingle = (First@Dimensions@# - Count[#, {0, 0}] - {1}) &@cList;                (* Check if only one coefficient *)
cRemainder = {line - Plus @@ ((#1 + I #2) #3) & @@                              (* Check for cross terms         *)
              Join[Transpose@cList, {Join[{1}, factors]}] // Expand};
Print /@ FortranForm /@ Flatten@Join[cList, cSingle, cRemainder];               (* Printout, cRemainder should=0 *)
Exit[0];                                                                        (* Exit cleanly                  *)
MATH

my @mathpath = (                                                                # Path to mathematica kernel
    "$usrMath",                                                                 # Possible user defined entry
    "/usr/bin/math",                                                            # *nix
    "/usr/local/bin/math",                                                      # *nix
    "/Applications/Mathematica.app/Contents/MacOS/MathKernel"                   # Mac OS X
);

my $math;
for (@mathpath) { if ( -e $_ ) { $math = $_; last; } }                          # Use the first kernel to be found
if ( not $math ) { erm( -1, "Mathematica kernel not found.\n" ) }

sub callMath {                                                                  # Subroutine for calling Mathematica
    my $code = $_[0];
    $code =~ s/\*{2}/\^/g;                                                      # Fix ** -> ^ in the input (in case of calchep files)
    my $cmd = qq{$math -noprompt -run '<<$script' '$code'};                     # Set up the command line call for Mathematica
    my $pid = open( F, "$cmd |" ) or erm( -1, "Can't open pipe from @_: $!" );  # Open pipe to Mathematica and note the process number
    $SIG{INT}  = sub { kill( 'INT',  $pid ); print "\n"; exit; };               # Pass along some interrupt signals to make
    $SIG{TERM} = sub { kill( 'TERM', $pid ); print "\n"; exit; };               #   sure the kernel quits as well if
    $SIG{QUIT} = sub { kill( 'QUIT', $pid ); print "\n"; exit; };               #   the script is stopped via interrupt
    $SIG{ABRT} = sub { kill( 'ABRT', $pid ); print "\n"; exit; };               #   For some reason, very rarely this doesn't work and
    $SIG{HUP}  = sub { kill( 'HUP',  $pid ); print "\n"; exit; };               #   the Kernel is left hanging using 100%CPU :-(.
    my @result = <F>;                                                           # Run Mathematica
    close(F);
    map s/\*{2}/\^/g, @result;                                                  # Fix ** -> ^ in the output (Math FortranForm uses **)
    return \@result;                                                            # Return array ref, since sub() only returns scalars
}

################################################################################
# Read in the lgrng.mdl and func.mdl file
################################################################################

open( INF, "<$inF" ) or erm( -1, "Can't open model file $inF: $!" );            # Open the Functions file
@outputF = <INF>;                                                               # Put the original func file into the output
close(INF);
if ($l2c) { map s/\*{2}/\^/g, @outputF }                                        # Fix ** -> ^

my @lines;
open( INL, "<$inL" ) or erm( -1, "Can't open model file $inL: $!" );            # Open the Lagrangian model file
@lines = <INL>;                                                                 # Read the file
close(INL);

################################################################################
# Main routine preamble
################################################################################

$" = "";                                                                        # $LIST_SEPARATOR need to be empty string for formatting reasons
my $line;                                                                       # A line in the lgrng file
my @parts;                                                                      #   the split up parts of that line (split at |'s)
my $part;                                                                       #   the "Lorentz part" of that line
my @forms;                                                                      #   that part's coefficient list
my @coeff = ("","*G5","*G(m3)","*G(m3)*G5",
             "*m1.p2","*m1.p3","*m1.p4","*m2.p1","*m2.p3","*m2.p4",
             "*m3.p1","*m3.p2","*m3.p4","*m4.p1","*m4.p2","*m4.p3",
             "*m1.m2","*m1.m3","*m1.m4","*m2.m3","*m2.m4","*m3.m4");
my $lineNum      = 0;                                                           # Line no. from start of file, not start of substitutions
my $maxLenF      = 0;                                                           # For keeping track of the longest substitution going into func.mdl
my $maxLenL      = 0;                                                           # For keeping track of the longest line left in lgrng.mdl
my @defLine      = split /\|/, $lines[2];                                       # The line that defines the various widths for entries in lgrng.mdl
my $factorLength = length $defLine[4];                                          #   the width of the "Factor" part therein
my $maxLineNum   = $#lines + 1;                                                 # Largest line number
my $cutoff       = 1950;                                                        # The hard cutoff is probably(?) 1992 (2000 - 8). Change at your own risk!!

$digits = max( length $maxLineNum, $digits, 2 );                                # How many digits does the line ID need (min 2 to avoid mixup with e.g. G5)
if ( $digits > 6 ) { erm( -1, "Too many lines!" ); }                            # More than 999999 lines cannot be accommodated
if ( $digits == 6 and $long ) { erm( 0, " Too many lines for -L option" ); }    # Check that --long isn't used when it shouldn't
my @letters = split //, "abcdefghijklnopqrtsuvwxyzABCDEFGHIJKLNOPQRSTUVWXYZ";   # Make an array of letters to use in the line ID (M is illegal)
if ( $digits < 6 and $long ) {                                                  # If there are only <100k lines and the user wants to
    my @temp = @letters;                                                        #   we can make two letter ID's so we can
    my @ids;                                                                    #   accomodate much, much longer lines (~1.2 million characters)
    while (@letters) {                                                          # For each letter
        my $t = shift @letters;                                                 #   select it
        push( @ids, map( $t.$_, @temp ) );                                      #   concatenate with all letters and store in $ids
    }
    @letters = @ids;
}

#my $test = callMath("a**2+b**3*i+(c1-c2)*G5+d*i*G5");                           # Test that the Math kernel works correctly
#if ( "@$test" !~ m/a\^2\nb\^3\nc1\s-\sc2\nd\n/ ) {
#    erm( -1, "Mathematica not working properly." );
#}

################################################################################
# Main routine proper
################################################################################

foreach $line (@lines) {
    $lineNum++;
    if ($l2c) { $line =~ s/\*{2}/\^/g }                                         # Fix ** -> ^
    my $status = sprintf "%".length($maxLineNum)."d/%d (%3d%%)", 
        $lineNum, $maxLineNum, $lineNum / ($maxLineNum) * 100;
    print "\rProcessing line $status";                                          # Small bit to print the progress
    if ( $lineNum < 4 ) {                                                       # Skip the first three lines of the model file
        push @outputL, $line;                                                   #   but move them into the output
        next;
    }
    @parts = split /(\|)/, $line;                                               # Split the line
    if ( $#parts != 10 ) {                                                      # Simple check for Lagrangian file formatting
        erm( -1, "\n$inL: Error at line $lineNum\n$#parts\n" );
    }
    $parts[8] .= " " x ( $factorLength - length $parts[8] );                    # The ** -> ^ transform shortens the lines, so we fill them up
    $part = pop @parts;                                                         # Take the "Lagrangian" part
    if ( length $part < "$chopL" ) {                                            # Check for too short form to bother with replacement
        push @outputL, "@parts".$part;                                          #   move the line into the output if yes
        $maxLenL = max( $maxLenL, length $part );                               #   keep track of the longest line
        next;
    }
    if ( $part =~ m/[pmM][1234]\.[pmM][1234]\*[pmM][1234]\.[pmM][1234]/ ) {     # Check for quadruple indices e.g. m1.m2*m3.m4
        push @outputL, "@parts".$part;                                          #   move the line into the output if yes
        $maxLenL = max( $maxLenL, length $part );                               #   keep track of the longest line
        next;
    }
    $part   =~ s/G\(m3\)/G3/g;                                                  # Do some substitutions so that factoring in
#    $part   =~ s/\(1-G5\)/PL/g;                                                 #   Mathematica is much easier, e.g., remove dot 
#    $part   =~ s/\(1\+G5\)/PR/g;                                                #   product and parens
    $part   =~ s/\.//g;
    my $ref = callMath($part);                                                  # Make Mathematica do the chopping up
    @forms  = @$ref;                                                            # Pass of the array via reference
    chomp @forms;                                                               # Pesky newlines
    if ( pop @forms ne 0 ) {                                                    # Check that the routine was able to get all the factors
       erm( -1, "\nCould not isolate all factors on line $lineNum.");           # i.e. that no weird cross term remainder was left over
    }
    my $manyCoeff =             pop @forms;                                     # If there is only one coefficient this is = 0
    map s/\s//g,                    @forms;                                     # Mathematica uses unnecessary whitespace
    map s/\(([^)+-]+)\)/$1/g,       @forms;                                     # Sometimes there are parentheses leftover from ComplexExpand ??
    my @check = map m/(\([^)]+\))/, @forms;                                     #   this is very disturbing !! So we double check that all such parentheses
    if (@check) {                                                               #   are removed correctly and exit if this is not possible
        $" = "\n";                                                              #   i.e. if they are of the form (X+Y) or (X-Y) and not (X*Y)
        erm( -1, "\nSuspicious parentheses found on line $lineNum:\n@check" );
    }
    my @lineID    = map sprintf( '%s%0*d', $_, $digits, $lineNum ), @letters;   # The line ID numbers in AA000## format
    my $wP        = ( ( $long and $digits == 5 ) or $digits == 6 ) ? "" : " ";  # Really hacky way of prepending correct whitespace to lineID later
    my $wA       .= " " x ( 7 - length $wP.$lineID[0] );                        # Ditto for appending
    my $plus      = "";
    my $reimForm  = "";
    my $finalForm = "";
    my $superID   = "";
    my $IDcount   = 0;
	my $nCoeff    = 2 * ( $#coeff + 1 );                                        # Number of coefficients (incl. re and im parts)
    for ( my $i = 0 ; $i < $nCoeff ; ++$i ) {                                   # Start putting the pieces together
        my $k     = int $i / 2;
        my $oddQ  = $i / 2 - $k;                                                # Used to test for even (re part) or odd (im part) iterations of the loop
        my $oCB   = ( $oddQ ) ? "" : "(";                                       # Open parens for the whole coefficient
        $reimForm = ( $oddQ ) ? $reimForm : "";                                 # Initialize $reimForm only on even (re) loop cycles
        if ( $forms[$i] ) {                                                     # Don't make substitution ID's for non-existent parts
            my $protoForm = "";                                                 # The first variable we collect lump IDs into
            my $manyParts = 0;                                                  # To check if we have only one lump
            my @bits      = split /([+-]?[A-Za-z0-9\/*^]+)/, $forms[$i];        # Chop the substitution into monomial parts [+-]?[^+-] works but is dangerous
            while (@bits) {
                my $newLine = "";                                               # The variable we collect the new func file line into
                if ( length $bits[0] >= $cutoff ) {                             # This should never, ever happen, but let's check just in case
                    erm( -1, "\nAssembly error at line: $lineNum" );            #   code that could handle the situation would go here
                }
                while ( @bits and length $newLine.$bits[0] < $cutoff ) {        # Reassemble the substitutions in ~< 1950 character long lumps
                    $newLine .= shift @bits;
                }
                $maxLenF = max( $maxLenF, length $newLine );                    # Keep track of the longest line
                push @outputF, $wP.$lineID[$IDcount].$wA."|".$newLine."\n";     # Fill up the .func output with the lump
                $protoForm .= $plus.$oCB.$lineID[$IDcount];                     # Fill up the protoForm of the .lgrng line with the lump ID
                if ( not @bits ) {                                              # If we run out of stuff for the lumps
                    if ( $manyParts ) {                                         #   and filled more than one lump
                        $protoForm =~ s/^(\(|\+)//;                             #   we remove the opening parens / "+" (for re / im respectively)
                        my $sID = ( $oddQ ) ? $#lineID : $#lineID - 1;          #   and fetch a lumpID number from the end of the ID list
                        my $ppp = ( $oddQ ) ? "+" : "(";                        # For prepending the new lumpID with correct parens / "+"
                        push @outputF,$wP.$lineID[$sID].$wA."|".$protoForm."\n";# Put old lumpIDs into the func file and
                        $protoForm = $ppp.$lineID[$sID];                        #    append the new collective lumpID to the .lgrng line
                    }
                    $reimForm .= $protoForm ;                                   # Move the single lumpID to a new variable that collects the re + im part
                    $reimForm .= ( $oddQ ) ? "*i" : "" ;                        #   and append a factor *i if neccessary
                }
                $plus = "+";                                                    # A "+" between later parts
                $oCB = "";                                                      #   but no parens
                $manyParts++;                                                   #   since we now fill up many lumps
                $IDcount++;                                                     # Increment the ID number
                if ( $IDcount == $#letters - 1) {                               # We can only make as many "lumps" as there are separate indices per line
                    if ( $digits < 6 ) {                                        #   note that the last two (Y&Z or YZ&ZZ) Ids are reserved, so no ( $IDcount == $#letters +1 )
                        erm( -1, "\nToo long line: $lineNum \nTry -L option" );
                    }
                    else {                                                      # This should practically never happen (knock on wood)
                        my $lineLength = length $line;
                        erm( -1, <<"ERROR");                                    # Very ugly formatting, maybe just print a "look at source" msg here
\nWow!
It is unfortunately impossible to accomodate Your model due to the combination
of having very many (>100k) and very long lines. Assembly of the placeholder
ID's ran out of ID's on line $lineNum ($lineLength characters long).\n
Please see the source code for a more detailed explanation.
ERROR
                    }
                }
            }    
        }
        if ( $oddQ ) {                                                          # After every odd (im) part we need to make some checks and cleanup
            if ( $forms[$i] and $forms[$i-1] and $i ne 1) {                     # If both re AND im part were present AND the coefficient is not "1"
               $reimForm .= ")";                                                #   close the parens, so you get (A01+B01*i)
            }
            else {
               $reimForm =~ s/\(([^\(]*)/$1/;                                   #   otherwise strip the leading parens
            }
            if ( ( $forms[$i] xor $forms[$i-1] ) and not $manyCoeff ) {         #  If there is only a single coefficient and only a real or only an imaginary part
                $reimForm =~ s/^\+//;                                           #   strip the sign
                my $imQ = "";                                                   #
                if ( $forms[$i] ) {                                             #   if the part was imaginary
                    $reimForm =~ s/\*i$//;                                      #   strip the *i
                    $imQ = "*i"                                                 #   for adding *i later
                }
                $parts[8] =~ s/\s*$/\*$reimForm/;                               #   add the ID into the Factor part of lgrng.mdl
                $finalForm  = ( $coeff[$k] ) ? substr $coeff[$k],1 : "1";       #   and put just the coeffcient in the Lorentz part  
                $finalForm  .= $imQ;                                            #   with *i if it was imaginary
            } 
            elsif ( $forms[$i] or $forms[$i-1] ) {                              # In all other cases
                $finalForm .= $reimForm.$coeff[$k];                             #   add the coefficient and append the re+im form to the finalForm
            }
        }
    }
    $parts[8] .= " " x ( $factorLength - length $parts[8] );                    # Our process shortens the Factor part, so we fill it up with whitespace
    push @outputL, "@parts".$finalForm."\n";                                    # Fill up the .lgrng output
    $maxLenL = max( $maxLenL, length $finalForm );                              # Keep track of the longest line
}

################################################################################
# Main routine cleanup
################################################################################

my $blankF .= " " x $maxLenF;                                                   # Make an empty pad as long as the longest substitution
$outputF[2] =~ s/(\|>\sExpression)\s+(<&>)/$1$blankF$2/g;                       # Pad out the defining line in the func file

# Lisää Factor osaan lisää tilaa

my $blankL .= " " x $maxLenL;                                                   # Make an empty pad as long as the longest line left in lgrng.mdl
$outputL[2] =~ s/\s+(<\|$)/$blankL$1/g;                                         # Pad out the defining line in the lgrng file

################################################################################
# Write out the new lgrng and func files
################################################################################

open( OUTF, ">$outF" );
print OUTF "@outputF";
close(OUTF);
open( OUTL, ">$outL" );
print OUTL "@outputL";
close(OUTL);
print " Done!\n";

################################################################################
# Copy associated model files
################################################################################

if ($copy) {
    my @fnames = ( "extlib", "prtcls", "vars" );
    for ( 0 .. 2 ) {
        system "cp " . $pathIn  . $fnames[$_] . $modelIn  . ".mdl "
                     . $pathOut . $fnames[$_] . $modelOut . ".mdl";
        if ($verbose) { 
            print "Copied " . $pathIn  . $fnames[$_] . $modelIn  . ".mdl"
                  . " --> " . $pathOut . $fnames[$_] . $modelOut . ".mdl\n";
        }
    }
}
