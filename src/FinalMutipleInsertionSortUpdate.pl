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
my %hash; my $length0=0; my $id0; my $counts0; 
#my %finaltype; 
my %finalstring; my %finalnum;
my $output=$opts{o}; my %delete;
open OUT,">$output" or die $!;

while (<REQ>){
	chomp;
	my @array=(split/\t/,$_);
	my ($indexid,$id,$type,$wholeseq,$insertseq,$length,$num,$iden,$cov,$qual,$chr,$start,$end,$strand,$feature,$distance,$partial,$alter,$upstream,$downstream,$status,$upshift,$downshift) =@array;

	if ($indexid eq "CaseID"){
		print OUT "$_\n";
		next;
	}
	

	$finalstring{$id}.=$_."\n";
	if($upshift eq "Unknown"){
		$finalnum{$id}=$num;
		next;
	}
	
	
	if(abs($length-$length0)<=5 && $upshift <=10 && $downshift<=10 && $feature=~/^gene/){
		my $fcov=$cov+$counts0;
		if ($cov >$counts0){
			#$finalcount{$id}=$fcov;
			$delete{$id0}++;
			print "$finalstring{$id0}\n$finalstring{$id}\n";
	
			
		}else{
			#$finalcount{$id0}=$fcov;
			$delete{$id}++;
			print "$finalstring{$id0}\n$finalstring{$id}\n";
		}		
		$length0=$length;
		$id0=$id;
		$counts0=$cov;
		
	}else{
		$length0=$length;
		$id0=$id;
		$counts0=$cov;
	}
	
	next if($num =~/^T/);
	$finalnum{$id}= $num;	
	
}

foreach my $t (sort {$finalnum{$a} cmp $finalnum{$b}} keys %finalnum){
	next if (exists $delete{$t});
	my @array=split/\n/,$finalstring{$t};
	if($#array==1){
		print OUT "$finalstring{$t}";
	}else{
		my @sortedArray= sort (@array);
		
		my $string=join "\n",@sortedArray;
		print OUT "$string\n";
		
	}


}

close OUT;
close REQ;
