#!/usr/bin/perl
use warnings;
use strict;
unless(@ARGV == 3){
	print "Usage:perl $0 <tmap> <pos> <merged.gtf>\n";
	exit;
}
open my $tmap,"$ARGV[0]" or die;
open my $gtf,"$ARGV[2]" or die;
my %hash;
my $pos = $ARGV[1];
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
		if(/transcript_id "([\w\.\-]+)"/){
		print "$_\n" if exists $hash{$1};
		}
	}

}



