#!/usr/bin/perl
use warnings;
use strict;
unless(@ARGV == 2){
	print "Usage:perl $0 <lncRNA.lst>  <novel.longRNA.gtf>\n";
	exit;
}
open my $tmap,"$ARGV[0]" or die;
open my $gtf,"$ARGV[1]" or die;
my %hash;
my $pos = 0;
while(<$tmap>){
chomp;
my @mem = split;
my $id = $mem[$pos];
$hash{$id}++;
}


while(<$gtf>){
chomp;
	if (/^#/){
		print "$_\n";
	}else{
		if(/gene_id "([\w\.\-]+)"/){
		print "$_\n" if exists $hash{$1};
		}
	}

}



