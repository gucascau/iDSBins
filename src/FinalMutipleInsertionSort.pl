#!/usr/bin/perl
use strict;
use warnings;
my $version="1.0 version";
use Getopt::Long;

# Copyright (c) 2021 Dr. Kaifu lab    
# Author:Xin Wang 
# email: xin.wang@childrens.harvard.edu
# PI: Kaifu Chen

# Function: This script is to update the sample ID for the large insertion for multiple donors

my %opts;
GetOptions(\%opts,"i:s","o:s","h:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{o}  ||!defined $opts{i}||defined $opts{h}) {
	die "************************************************************************
	Usage: $0.pl  -i FinalMultipleDonors -o FinalSortedMultipleDonor -o Blast_output 
	
	### 

	Request Parameters:
	-i Final insertion with multiple donors or unknown donors
	-o Final insertion with multiple donors or unknown donors that are clearly sorted
	
	Optional Restriction Parameters:
	-h Help
	
************************************************************************\n";
}



### Read the raw insertion events

my $request=$opts{i};

open REQ,"$request" or die $!;

#my %finaltype; 
my %finalstring; my %finalnum;
my $output=$opts{o};
open OUT,">$output" or die $!;
while (<REQ>){
	chomp;
	my ($indexid,$id,$type,$wholeseq,$insertseq,$length,$num,$iden,$cov,$qual,$chr,$start,$end,$strand,$feature,$distance,$partial,$alter,$upstream,$downstream) =(split/\t/,$_);
	if ($indexid eq "CaseID"){
		print OUT "$_\n";
		next;
	}
	$finalstring{$id}.=$_."\n";
	next if($num =~/^T/);
	$finalnum{$id}= $num;	

}

foreach my $t (sort {$finalnum{$a} cmp $finalnum{$b}} keys %finalnum){	
	print OUT "$finalstring{$t}";

}

close OUT;
close REQ;