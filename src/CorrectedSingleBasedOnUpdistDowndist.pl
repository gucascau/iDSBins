#!/usr/bin/perl
use strict;
use warnings;
my $version="1.0 version";
use Getopt::Long;

# Copyright (c) 2021 Dr. Kaifu lab    
# Author:Xin Wang 
# email: xin.wang@childrens.harvard.edu
# PI: Kaifu Chen

# Function: This script is to get rid of potential duplicates in following way. Check the distance between NEIGHBOR insertions. As you group them at the end by chromosome and location within chromosome it is very easy. All 10 of the fake insertions have endpoints located  +/- 0 to 10 from each other and in all cases one of two insertions are very low read count – 1 typically and up to 5. 



my %opts;
GetOptions(\%opts,"i:s","o:s","g:s","h:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i}|| !defined $opts{o} ||defined $opts{h}) {
	die "************************************************************************
	Usage: $0.pl -i index -g Cut off shiftsize -o Blast_output 
	
	### 

	Request Parameters:
	-i Final raw insertion events
	-o Final high quality events, we sorted the single donor insertions, then multiple donor insertions and final unclear insertions.
	
	Optional Restriction Parameters:
	-h Help
	
************************************************************************\n";
}



### Read the insertion events

my $request=$opts{i};

open REQ,"$request" or die $!;

my %hash; my $length0=0; my $id0; my $counts0=0; 

my $output=$opts{o}; 

while (<REQ>){
	chomp;
	my @array=(split/\t/,$_);
	my ($indexid,$id,$type,$wholeseq,$insertseq,$length,$num,$iden,$cov,$qual,$chr,$start,$end,$strand,$feature,$distance,$partial,$alter,$upstream,$downstream,$status,$upshift,$downshift) =@array;
	next if ($indexid eq "CaseID");
	if(abs($length-$length0)<=5 && $upshift <=10 && $downshift<=10 && $feature=~/^gene/){
		my $fcov=$cov+$counts0;
		if ($cov >$counts0){
			my @array1=@array;
			splice(@array1,8,1,$fcov);
			#$array1[8]=$fcov;
			$hash{$id}=join "\t",@array1;
			print "$hash{$id0}\n$hash{$id}\n";
			delete $hash{$id0};
			
		}else{
			my @array0=split/\t/,$hash{$id0};
			#$array0[8]=$fcov;
			splice(@array0,8,1,$fcov);
			$hash{$id0}=join "\t",@array0;	
			print "$hash{$id0}\n$_\n";	
		}
		$length0=$length;
		$id0=$id;
		$counts0=$cov;
		
	}else{
		$hash{$id}=join "\t",@array;
		$length0=$length;
		$id0=$id;
		$counts0=$cov;
	}

	
}
close REQ;

open ONE,">$output.OneCorrected.txt" or die $!; 
print ONE "CaseID\tRepresentativeRead\tSampleID\tWhole sequence\tInsertion Sequence\tInsertion Size (bp)\tDonor Number\tRepresentive Read Identity\tTotalReadCount\tRead Quality\tDonor  Chromosome\tDonor Start\tDonor End\tDonor Strand\tDonor Feature\tDistance To Feature\tDistance Description\tAlternative Feature\tJunction_5\tJunction_3\tGapOfReads(3kb)\tUpdist\tDowndist\n";

foreach my $i (keys %hash){
	print ONE "$hash{$i}\n";
}

close ONE;


