#!/usr/bin/perl
use warnings;
use strict;

unless(@ARGV == 2){
	print "Usage:perl $0 <lncRNA.gtf> <taxon_code>\n";
	exit;
	}
open my $in,"$ARGV[0]" or die;
my %hash;
my %chr;
my $n = 0;
my $taxon = $ARGV[1];
while(<$in>){
chomp;
next if /^#/;
my @mem = split;

unless(exists $chr{$mem[0]}){
$n = 0;
$chr{$mem[0]}++;
}

/gene_id "(\w+\.\w+)"/;
my $gene = $1;
if(exists $hash{$gene}){
	my $lst = $hash{$gene};
	s/$gene/$lst/g;
	print "$_\n";
}else{
	$n += 1;
my $new_name;	
	if($n < 10){
	 $new_name = "${taxon}_$mem[0]LNC000$n";
	}elsif($n < 100){
	 $new_name = "${taxon}_$mem[0]LNC00$n";
	}elsif($n < 1000){
	 $new_name = "${taxon}_$mem[0]LNC0$n";
	}else{
	 $new_name = "${taxon}_$mem[0]LNC$n";
	}
	$hash{$gene} = $new_name;
	s/$gene/$new_name/g;
	print "$_\n";

}


}



