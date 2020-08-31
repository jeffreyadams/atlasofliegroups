#!/usr/bin/perl
$label="A";
`./atlas_timing.pl "SO(4,0)" >  timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(3,1)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(2,2)" >> timing_SO.txt 2>&1`;

`./atlas_timing.pl "SO(5,0)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(4,1)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(3,2)" >> timing_SO.txt 2>&1`;

`./atlas_timing.pl "SO(6,0)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(5,1)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(4,2)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(3,3)" >> timing_SO.txt 2>&1`;

`./atlas_timing.pl "SO(7,0)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(6,1)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(5,2)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(4,3)" >> timing_SO.txt 2>&1`;

`./atlas_timing.pl "SO(9,0)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(8,1)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(7,2)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(6,3)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(5,4)" >> timing_SO.txt 2>&1`;

`./atlas_timing.pl "SO(10,0)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(9,1)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(8,2)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(7,3)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(6,4)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(5,5)" >> timing_SO.txt 2>&1`;

`./atlas_timing.pl "SO(11,0)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(10,1)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(9,2)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(8,3)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(7,4)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(6,5)" >> timing_SO.txt 2>&1`;

`./atlas_timing.pl "SO(12,0)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(11,1)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(10,2)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(9,3)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(8,4)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(7,5)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(6,6)" >> timing_SO.txt 2>&1`;

`./atlas_timing.pl "SO(14,0)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(12,1)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(11,2)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(10,3)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(9,4)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(8,5)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(7,6)" >> timing_SO.txt 2>&1`;
`./atlas_timing.pl "SO(6,6)" >> timing_SO.txt 2>&1`;


`grep -v "Cannot\\|double\\|Abandon\\|Command" $label.orig.txt > $label.txt`
