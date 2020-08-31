#!/usr/bin/perl
$label="C";

`./atlas_timing.pl "Sp(4,R)" >  timing_C.orig.txt 2>&1`;
`./atlas_timing.pl "Sp(6,R)" >> timing_C.orig.txt 2>&1`;
`./atlas_timing.pl "Sp(8,R)" >> timing_C.orig.txt 2>&1`;
`./atlas_timing.pl "Sp(10,R)" >> timing_C.orig.txt 2>&1`;
#`./atlas_timing.pl "Sp(12,R)" >> timing_C.orig.txt 2>&1`;
#`./atlas_timing.pl "Sp(14,R)" >> timing_C.orig.txt 2>&1`;

`grep -v "Cannot\\|double\\|Abandon\\|Command" $label.orig.txt > $label.txt`
