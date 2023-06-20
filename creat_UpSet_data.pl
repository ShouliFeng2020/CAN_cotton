#!/usr/bin/perl
use warnings;
use strict;
unless(@ARGV == 2){
	print "Usage: perl $0 <synteny> <blast>\n";
	exit;
}

open my $f1,"$ARGV[0]" or die;
open my $f2,"$ARGV[1]" or die;

my %f1;
my %f2;
my %hash;

while(<$f1>){
chomp;
next if /^#/;
my @mem = split;
my $lst = "$mem[0]_$mem[1]";
$f1{$lst}++;
$hash{$lst}++;
}

while(<$f2>){
chomp;
next if /^#/;
my @mem = split;
my $lst = "$mem[0]_$mem[1]";
$f2{$lst}++;
$hash{$lst}++;
}

print "Homologous\tsynteny\tRBblast\n";

foreach(sort keys %hash){

my $o1 = exists $f1{$_}?1:0;
my $o2 = exists $f2{$_}?1:0;

print "$_\t$o1\t$o2\n";
}


