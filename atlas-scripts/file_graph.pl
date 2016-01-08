#!/usr/bin/perl
#generate graph of inclusions of at files
#usage: file_graph.pl files
# generates files and files.ps

use strict;
use List::MoreUtils qw/uniq/;

my $out=$ARGV[0];

open(OUT,">$out")||die("Can't open $out for output");

#my $files=`grep -h "^<.*at" *`;
my $files=`ls *at`;
my @files=split "\n", $files;

select(OUT);
print "
digraph G {
ratio=\"1.5\"
size=\"7.5,10!\"";
foreach my $file (@files){
    my $string=$file;
    $string=~ s/\.at//;
    print("\"$string\"\n");
}

foreach my $file (@files){
    my $targets=`grep "^<.*at" $file`;
    $file =~ s/\.at//g;
    my @targets=split "\n", $targets;
    foreach my $target (@targets){
	$target=~ s/<//;
	$target=~ s/{.*//g;
	$target=~ s/\s//g;
	$target=~ s/\.at//;
	print("\n\"$file\"->\"$target\"");
    }
}
print("\n}");
`dot -Tps -o$out.ps $out`;
close(OUT);








