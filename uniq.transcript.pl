#!/usr/bin/perl
use warnings;
use strict;

unless(@ARGV ==2 ){
	print "Usage: perl $0 <lst> <gtf>\n";
	exit;
}

open my $lst,"$ARGV[0]" or die;
open my $gtf,"$ARGV[1]" or die;

my %hash;
while(<$lst>){
/transcript_id "([\w\.:]+)"/;
$hash{$1}++;
}

while(<$gtf>){
/transcript_id "([\w\.:]+)"/;
print if exists $hash{$1};

}


