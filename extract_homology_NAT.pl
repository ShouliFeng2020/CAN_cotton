#!/usr/bin/perl
use warnings;
use strict;
unless(@ARGV == 3){
	print "Usage:perl $0 <TM1> <H7> <homology.lst>\n";
	exit;
}

open my $f1,"$ARGV[0]" or die;
open my $f2,"$ARGV[1]" or die;
open my $ho,"$ARGV[2]" or die;

my %gh;
my %gb;

while(<$f1>){
chomp;
my @mem = split;

$mem[11] =~ s/_[ekq]+//;
$gh{$mem[11]} = $mem[0];
}

while(<$f2>){
chomp;
my @mem = split;

$mem[11] =~ s/_[ekq]+//;
$gb{$mem[11]} = $mem[0];
}

while(<$ho>){
chomp;
my @mem = split;
if(exists $gh{$mem[0]} && exists $gb{$mem[1]}){
	print "$gh{$mem[0]}\t$gb{$mem[1]}\n";
}
}

