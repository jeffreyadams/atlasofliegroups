#!/usr/bin/perl
`./atlas_timing.pl "SO(2,2)" >  timing_D.orig.txt 2>&1`;
`./atlas_timing.pl "SO(3,3)" >> timing_D.orig.txt 2>&1`;
`./atlas_timing.pl "SO(4,4)" >> timing_D.orig.txt 2>&1`;
`./atlas_timing.pl "SO(5,5)" >> timing_D.orig.txt 2>&1`;
#`./atlas_timing.pl "SO(6,6)" >> timing_D.orig.txt 2>&1`;
#`./atlas_timing.pl "SO(7,7)" >> timing_D.orig.txt 2>&1`;

