#!/usr/bin/perl
use strict;
use warnings;
my $version="1.0 version";
use Getopt::Long;

#Author:Xin Wang 
#email: xin.wang@childrens.harvard.edu
#PI: Kaifu Chen

### This script is to simply to cluster the identical sequences to improve the speed.
### requirement: 1
#####


my %opts;
GetOptions(\%opts,"i:s","o:s","g:s","b:s","f:s","r:s","h:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{g} ||!defined $opts{f} || !defined $opts{r} ||!defined $opts{i}|| !defined $opts{o} ||defined $opts{h}) {
	die "************************************************************************
	Usage: $0.pl -g Assembled Reads -f Unassembled R1 -r Unassembled R2 -i index -c Clustering Files With Coverage -o Blast_output 
	
	### This script is to simply to cluster the identical sequences to improve the speed.

	
	Request Parameters:
	-g Assembled fastq
	-f Unassembled forward reads
	-r Unassembled reverese reads (Attention: this read will be different and normally reverse the original reverse read after PEAR merge)
	-i Detected reads with quality, identify rmatches and Relength type (JKM139-1_S2_detected_combined.highqual.txt)
	-o Simple cluster file that output with quality, identify rmatches Relength type and Read counts supports
	
	Optional Restriction Parameters:
	-h Help
	
************************************************************************\n";
}


### read from assembled reads
my $assembly=$opts{g}; my %hash0;  my %string0;
my $Aid;
open ASS,"$assembly" or die $!;
while (<ASS>){
	chomp;
    if ($. % 4 == 1)  {
		$Aid=(split/\s+/,$_)[0];
		$Aid=~s/^@//;
		#print "$id1\n";
		$hash0{$Aid}=$_."\n";
	}elsif($. % 4 ==2){
		$hash0{$Aid}.=$_."\n";
		$string0{$Aid}=$_;
	}else{
		$hash0{$Aid}.=$_."\n";
	}
	
}
close ASS;


### read forward reads and reverse reads
my $Fread=$opts{f};
open FASTQ1,"$Fread" or die $!;

my $Fid;  my %hash1;  my %string1;
while (<FASTQ1>) {
	chomp;
    # print "$_" ;

    if ($. % 4 == 1)  {
		$Fid=(split/\s+/,$_)[0];
		$Fid=~s/^@//;
		#print "$id1\n";
		$hash1{$Fid}=$_."\n";
	}elsif($. % 4 ==2){
		$hash1{$Fid}.=$_."\n";
		$string1{$Fid}=$_;
	}else{
		$hash1{$Fid}.=$_."\n";
	}
	
}
close FASTQ1;

my $Rread=$opts{r};

open FASTQ2,"$Rread" or die $!;

my $Rid;  my %hash2; my %string2;
while (<FASTQ2>) {
	chomp;
    # print "$_" ;

    if ($. % 4 == 1)  {
		$Rid=(split/\s+/,$_)[0];
		$Rid=~s/^@//;
		#print "$id1\n";
		$hash2{$Rid}=$_."\n";
	}elsif($. % 4 ==2){
		$hash2{$Rid}.=$_."\n";
		$string2{$Rid}=$_;
		
	}else{
		$hash2{$Rid}.=$_."\n";
	}
	
}


## Read the clustering files and seperated these reads into two categories: assembled and unassembled:

my $request=$opts{i};


my $output=$opts{o};


open REQ,"$request" or die $!;

open OUT,">$output" or die $!;

my %num; my $string; my %final; my %inf;
while (<REQ>){
	chomp;
	my $id =(split/\t/,$_)[0];
	next if ($id eq "ID");
	$inf{$id}=$_;
	
	if (exists $hash0{$id}){
		$string=$string0{$id};
		$final{$string}.="$id;";
		$num{$string}++;
		# print ASEM "$hash0{$id}";
# 		print ASEMA ">$id\n$string0{$id}\n";
# 		print OUT "Assembly\t$_\n";
		
	}else{
	
		$string=join '', ($string1{$id},$string2{$id});
		$final{$string}.="$id;";
		$num{$string}++;
	}
}

print OUT "ID\tReadQuality\tRIdentity\tRMatches\tRLength\tType\tRcounts\n";
foreach my $m (keys %final){
	my @array=split/;/,$final{$m};
	my $information=$inf{$array[0]};
	my $finalnum=$num{$m};
	print OUT "$information\t$num{$m}\n";
	
}
close OUT;
