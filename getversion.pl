#!/usr/bin/perl
use strict;
use warnings;

my $version_file = "sources/version.h";
open(my $in, "<", $version_file) or die ("no file $version_file: $!");
while(<$in>) {  print $1 if /VERSION += +"(\S+)"/; }
