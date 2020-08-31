#!/usr/bin/perl
$label="B";
`./atlas_timing.pl "SO(3,2)" >  timing_B.orig.txt 2>&1`;
`./atlas_timing.pl "SO(4,3)" >> timing_B.orig.txt 2>&1`;
`./atlas_timing.pl "SO(5,4)" >> timing_B.orig.txt 2>&1`;
`./atlas_timing.pl "SO(6,5)" >> timing_B.orig.txt 2>&1`;
`./atlas_timing.pl "SO(7,6)" >> timing_B.orig.txt 2>&1`;
`./atlas_timing.pl "SO(8,7)" >> timing_B.orig.txt 2>&1`;
`./atlas_timing.pl "SO(9,8)" >> timing_B.orig.txt 2>&1`;

`grep -v "Cannot\\|double\\|Abandon\\|Command" $label.orig.txt > $label.txt`

