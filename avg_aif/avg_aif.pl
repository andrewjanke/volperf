#! /usr/bin/env perl
#
# Andrew Janke - rotor@cmr.uq.edu.au
# Center for Magnetic Resonance
# The University of Queensland
# http://www.cmr.uq.edu.au/~rotor
#
# Copyright Andrew Janke, The University of Queensland.
# Permission to use, copy, modify, and distribute this software and its
# documentation for any purpose and without fee is hereby granted,
# provided that the above copyright notice appear in all copies.  The
# author and the University of Queensland make no representations about the
# suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
#
# Mon Aug  4 14:17:05 EST 2003 - initial version

$| = 1;

use strict;
use warnings "all";
use Getopt::Tabular;
use File::Basename;

my($Help, $Usage, $me, @opt_table, $history, %opt);
my(@args, @infiles, $outfile);

chomp($me = `basename $0`);
%opt = ('clobber' => 0,
        'verbose' => 0,
        'nsamples' => 50,
        'mask' => '',
        );

$Help = <<HELP;
| $me creates an average AIF from an input time series of
|        files and a mask.
|
| Problems or comments should be sent to: rotor\@cmr.uq.edu.au
HELP

$Usage = "Usage: $me [options] -mask <mask.mnc> <in1.mnc> [<in2.mnc> [...]] <avg.aif>\n".
         "       $me -help to list options\n\n";

@opt_table = (
     ["-verbose",    "boolean", 0, \$opt{'verbose'},
        "be verbose"                                             ],
     ["-clobber",    "boolean", 0, \$opt{'clobber'},
        "clobber existing files"                                 ],
     ["-nsamples", "integer", 1, \$opt{'nsamples'},
        "The number of samples to choose from the mask."         ],
     ["-mask", "string",  1, \$opt{'mask'},
        "mask image to use"                                      ],
     );


# get history string
chomp($history = `date`);
$history .= '>>>> ' . join(' ', $me, @ARGV);

# Check arguments
&Getopt::Tabular::SetHelp($Help, $Usage);
&GetOptions(\@opt_table, \@ARGV) || exit 1;
die $Usage if ($#ARGV < 0);

# get infiles and outfile
$outfile = pop(@ARGV);
@infiles = @ARGV;

# check for output file
if(-e $outfile && !$opt{'clobber'}){
   die "$me: $outfile exists, use -clobber to overwrite\n\n";
   }
   
# slurp in the data and build the averages
my($c, @totals);

print STDOUT "Getting Samples:";
@args = ('mincsample', '-clobber',
         '-mask', $opt{'mask'},
         '-random_samples', $opt{'nsamples'},
         '-ascii',
         @infiles);
open(FH, join(' ', @args) . " | ");

$c = 0;
foreach (<FH>){
   chomp;
   
   $totals[$c] += $_;
   
   $c++;
   
   if($c == ($#infiles + 1)){
      $c = 0;
      
      print STDOUT '.';
      }
   }
print STDOUT "Done\n";
close(FH);

# create the average time-series
foreach (@totals){
   $_ /= $opt{'nsamples'};
   }

# add the closing ';'
$totals[$#totals] .= ';';

# write out the average AIF
open(AIF, ">$outfile");
print AIF "MINC Vector File\n" .
          "%\n" .
          "% Created by $me\n" .
          "% $history \n" .
          "\n" .
          "MINC_Vector_Type = Normal_MINC_Vector;\n" .
          "MINC_Vector =\n" .
          join("\n", @totals) .
          "\n";
close(AIF);
