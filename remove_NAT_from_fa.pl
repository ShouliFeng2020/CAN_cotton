#!/usr/bin/perl
use warnings;
use strict;
unless(@ARGV == 1){
	print "\nUsage:perl $0 <cds.fa>\n";
	exit;
}


my %hash;
my $key;
while(<>){
chomp;
if(/>(\w+)/){
 $key = $1;

}else{

if(exists $hash{$key}){
	$hash{$key} .= $_;
}else{
	$hash{$key} = $_;
}

}

}


foreach(sort keys %hash){
	print ">$_\n$hash{$_}\n" unless /NAT/;


}

