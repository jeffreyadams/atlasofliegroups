#!/usr/bin/perl
$label="timing_exceptional";

`./atlas_timing.pl "G2_s" > $label.orig.txt 2>&1`;
#`./atlas_timing.pl "F4_s" >> $label.orig.txt 2>&1`;
#`./atlas_timing.pl "E6_s" >> $label.orig.txt 2>&1`;
#`./atlas_timing.pl "E6_q" >> $label.orig.txt 2>&1`;

`grep -v "Cannot\\|double\\|Abandon\\|Command" $label.orig.txt > $label.txt`









    
