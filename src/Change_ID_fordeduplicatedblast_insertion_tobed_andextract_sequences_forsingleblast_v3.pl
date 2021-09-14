#!/usr/bin/perl
#author:wangxin
### #### caluculate the side effect on the final associated genes results
### Here is to calculate the sample size effect


#!/usr/bin/perl
use strict;
use warnings;
my $version="1.0 version";
use Getopt::Long;

#Author:Xin Wang 
#email: xin.wang@childrens.harvard.edu
#PI: Kaifu Chen

### This script is to generate the large insertion table and bedfile for further annotation


my %opts;
GetOptions(\%opts,"i:s","o:s","g:s","q:s","h:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{g} || !defined $opts{q} ||!defined $opts{i}||!defined $opts{o} ||defined $opts{h}) {
	die "************************************************************************
	Usage: $0.pl -g Assembled fasta -q Reads quality file -i final single insertion -o Output string name
	
	Request Parameters:
	 
	-g: Assembled fasta file
	-q: Quality of each read file
	-i: Final insertion with single donor information
	-o: Output string files, including insertion assembled fasta, assembled bed file for further annotation, final primilary table with statistic results
	
	Optional Parameters:
	-h Help
	
************************************************************************\n";
}


#### The index of the reads sequence information

my $fasta=$opts{g};

open FASTA, "$fasta" or die "cannot open file $fasta";

my %str; my $id3;
while (<FASTA>){
	chomp;	
	if ($_=~/>(\S+)/){
		$id3=$1;
	}else{
		$str{$id3}.=$_;
	}  
}
close FASTA;


#### The index of the read quality information

my $readquality=$opts{q};
my %qual;
open RQ, "$readquality" or die $!;
while (<RQ>){
	chomp;
	my ($id,$quality)=(split/\t/,$_)[0,5];
	$qual{$id}=$quality;
}
close RQ;

my $out=$opts{o};
my $insert=$opts{i};

open INS,"$insert" or die $!;
open BED,">$out.Fsingle.bed" or die $!;
open FA,">$out.Fsingle.fasta" or die $!;
open FT, ">$out.Fsingle.table.txt" or die $!;

my $n=0;
while (<INS>) {
	chomp;
	s/\r//g;
	my @array=split/\t/,$_;
	my $identity=$array[3];	
	my $string=$str{$array[1]};
	$n++;
	my $type=($array[11]>$array[10])?"+":"-";
	my $start=($array[11]>$array[10])?$array[10]:$array[11];
	my $end=($array[11]>$array[10])?$array[11]:$array[10];
	my $chro=$array[2];
	my $length1=$array[7];
	my $length2=$array[9]-$array[8]+1;
	my $coverage=$array[15];
	my $Fquality=$qual{$array[1]};
	
	my $Donorstart=$array[7];
	my $Donorend=$array[8];
	my $Donorlength=$Donorend-$Donorstart+1;
	
	my $Donorseq=substr $string, $Donorstart, $Donorlength;
	
	print BED "$chro\t$start\t$end\t$array[1]\t$coverage\t$type\n";
	
	print FA ">$array[1]\n$string\n";
	

	if ($length1 <=55 && $length2<=55){
		print FT "$array[1]\t$string\t1\t$coverage\t$identity\t$Fquality\t$Donorseq\t$Donorlength\t$chro\t$start\t$end\t$type\n";
	}else{
		print FT "$array[1]\t$string\t1ormore\t$coverage\t$identity\t$Fquality\t$Donorseq\t$Donorlength\t$chro\t$start\t$end\t$type\n";
	}
}
close INS;
close FT;
close FA;
close BED;


