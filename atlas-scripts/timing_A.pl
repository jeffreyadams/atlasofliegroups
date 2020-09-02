#!/usr/bin/perl
#perl script for timing atlas calculations
#calls atlas_timing.pl
$label="timing_A";

`./atlas_timing.pl "SL(2,R)" >  $label.orig.txt 2>&1`;
`./atlas_timing.pl "SL(3,R)" >> $label.orig.txt 2>&1`;
`./atlas_timing.pl "SL(4,R)" >> $label.orig.txt 2>&1`;
`./atlas_timing.pl "SL(5,R)" >> $label.orig.txt 2>&1`;
`./atlas_timing.pl "SL(6,R)" >> $label.orig.txt 2>&1`;
`./atlas_timing.pl "SL(7,R)" >> $label.orig.txt 2>&1`;
`./atlas_timing.pl "SL(8,R)" >> $label.orig.txt 2>&1`;
`./atlas_timing.pl "SL(9,R)" >> $label.orig.txt 2>&1`;

`grep -v "Cannot\\|double\\|Abandon\\|Command" $label.orig.txt > $label.txt`





