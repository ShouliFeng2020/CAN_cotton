#!/usr/bin/perl
use warnings;
use strict;
unless(@ARGV == 2){
	print "Usage: perl $0 <blast1> <blast2>\n";
	exit;
}

open my $in1,"$ARGV[0]" or die;
open my $in2,"$ARGV[1]" or die;

my %blast1;
my %blast2;

while(<$in1>){
chomp;
my @mem = split;
if(exists $blast1{$mem[0]}){
	my @ar=split/\t/,$blast1{$mem[0]};
	$blast1{$mem[0]} = $_ if $mem[11] > $ar[11];
}else{
	$blast1{$mem[0]} =$_;
}

}


while(<$in2>){
chomp;
my @mem = split;
if(exists $blast2{$mem[0]}){
        my @ar=split/\t/,$blast2{$mem[0]};
        $blast2{$mem[0]} = $_ if $mem[11] > $ar[11];
}else{
        $blast2{$mem[0]} =$_;
}

}


foreach(sort keys %blast1){
 my @mem = split/\t/,$blast1{$_};
 if(exists $blast2{$mem[1]}){
	my @ar = split/\t/,$blast2{$mem[1]};
	print "$mem[0]\t$mem[1]\n" if $mem[0] eq $ar[1];
}

}


