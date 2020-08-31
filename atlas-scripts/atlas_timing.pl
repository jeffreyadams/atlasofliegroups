#!/usr/bin/perl

#perl script for timing atlas calculations
#see timing_A.pl, timing_B.pl, ..., timing_SO.pl, timing_exceptional.pl

#atlas_timing.pl "Sp(4,R)" writes a temporary file
#timing_file.at which contains this:

#<all.at
#is_unitary(Sp(4,R).trivial)
#quit

#it then launches atlas and executes this file, which quits with
#an error at the quit statemt. Actually it runs this with time
#like this:
#`time  -p ../atlas timing_file.at`

#NEXT

#a command like this is found in timing_A.pl, or something similar:

#./atlas_timing.pl "SL(2,R)" >  timing_A.orig.txt 2>&1`;

#will run this computation, and direct both the name of the group, an
#error message, and the output of time to the desired file.

($G)=@ARGV[0];
print("\n$G\n",);
open(OUT,">timing_file.at");
print OUT "<all.at\n";
print OUT "is_unitary($G.trivial)\n";
print OUT "quit";
close(OUT);
`time  -p ../atlas timing_file.at`;





