#!/usr/bin/perl
use strict;
use warnings;
my $version="1.0 version";
use Getopt::Long;

# Copyright (c) 2021 Dr. Kaifu lab    
# Author:Xin Wang 
# email: xin.wang@childrens.harvard.edu
# PI: Kaifu Chen

# Function: This script is to trim the low quality events and correct the read coverage, then updated the "1orMore", "2orMore" or "3orMore" if they have unclear alignments
#	 requirement: 1
#		trim criterion: the low coverage with low quality (<=2 && <25)


my %opts;
GetOptions(\%opts,"i:s","o:s","g:s","h:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{g}  ||!defined $opts{i}|| !defined $opts{o} ||defined $opts{h}) {
	die "************************************************************************
	Usage: $0.pl -g Assembled Reads -i index -c Clustering Files With Coverage -o Blast_output 
	
	### 

	Request Parameters:
	-g Final clustering with read and read coverage
	-i Final raw insertion events
	-o Final high quality events, we sorted the single donor insertions, then multiple donor insertions and final unclear insertions.
	
	Optional Restriction Parameters:
	-h Help
	
************************************************************************\n";
}


## Read the clustering files 
my $cluter=$opts{g}; my %Rcount;  

open ASS,"$cluter" or die $!;
while (<ASS>){
	chomp;
	my ($count,$Rid)=(split/\t/,$_)[2,3];
	$Rcount{$Rid}=$count;
}
close ASS;



### Read the raw insertion events

my $request=$opts{i};

open REQ,"$request" or die $!;

my %finaltype; 

my $output=$opts{o};
open ONE,">$output.One.txt" or die $!;
open MU,">$output.Multiple.txt" or die $!;
print ONE "CaseID\tRepresentativeRead\tSampleID\tWhole sequence\tInsertion Sequence\tInsertion Size (bp)\tDonor Number\tRepresentive Read Identity\tTotalReadCount\tRead Quality\tDonor  Chromosome\tDonor Start\tDonor End\tDonor Strand\tDonor Feature\tDistance To Feature\tDistance Description\tAlternative Feature\tJunction_5\tJunction_3\tGapOfReads(3kb)\n";
print MU "CaseID\tRepresentativeRead\tSampleID\tWhole sequence\tInsertion Sequence\tInsertion Size (bp)\tDonor Number\tRepresentive Read Identity\tTotalReadCount\tRead Quality\tDonor  Chromosome\tDonor Start\tDonor End\tDonor Strand\tDonor Feature\tDistance To Feature\tDistance Description\tAlternative Feature\tJunction_5\tJunction_3\tGapOfReads(3kb)\n";
while (<REQ>){
	chomp;
	my ($indexid,$id,$type,$wholeseq,$insertseq,$length,$num,$iden,$cov,$qual,$chr,$start,$end,$strand,$feature,$distance,$partial,$alter,$upstream,$downstream) =(split/\t/,$_);
	next if ($id eq "ID");
	my $fcov=(exists $Rcount{$id})?$Rcount{$id}:$cov;
	next if ($fcov <=2 && $qual<25);
	
	my $sub=substr $wholeseq, 3,-3;
	
	next if (($sub!~/^GCA/ || $sub!~/TTG$/) && $fcov<=5);
	
	if ($num eq "1orMore"){
		$num= ($feature=~/^gene\|/)?"2orMore":$num;
	}
	
	my $describe= ($indexid =~/^E/)?"YES":"NO";
	my $fstring=join"\t", ($indexid,$id,$type,$wholeseq,$insertseq,$length,$num,$iden,$fcov,$qual,$chr,$start,$end,$strand,$feature,$distance,$partial,$alter,$upstream,$downstream,$describe);

	$finaltype{$num}.=$fstring."\n";
	
	if($num eq "1" && $iden ne "Unknown"){
		print ONE "$fstring\n" ;
	}else{
		print MU "$fstring\n";
	}
	
}

close ONE;
close REQ;
close MU;
