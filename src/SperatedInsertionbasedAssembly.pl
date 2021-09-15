#!/usr/bin/perl
use strict;
use warnings;
my $version="1.0 version";
use Getopt::Long;

#Author:Xin Wang 
#email: xin.wang@childrens.harvard.edu
#PI: Kaifu Chen

### This script is to divided all the insertion into two type based on the assembled information.
### requirement: 1
#####1.	We didived them into one donor, two donor, three donor and without donors


my %opts;
GetOptions(\%opts,"i:s","o:s","g:s","b:s","f:s","r:s","h:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{g} ||!defined $opts{f} || !defined $opts{r} ||!defined $opts{i}|| !defined $opts{o} ||defined $opts{h}) {
	die "************************************************************************
	Usage: extract_fastq.pl -g Assembled Reads -f Unassembled R1 -r Unassembled R2 -i index -c Clustering Files With Coverage -o Blast_output 
	
	### Here divide all the representive reads into assembled reads, unassembled reads. We setup assembled as the first prioxity.

	
	Request Parameters:
	-g Assembled fastq
	-f Unassembled forward reads
	-r Unassembled reverese reads (Attention: this read will be different and normally reverse the original reverse read after PEAR merge)
	-i Clustering File with coverage, read quality and identity. This provide the coverage, read quality, identity of each represent insertions generated before.
	-o Output substring of mutiple insertion types
	
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

### print the output of fastq files for all the insertions.
open ASEM,">$output.Insassembled.fastq";
open UASEM,">$output.InsUnassembled.R1.fastq";
open UBSEM,">$output.InsUnassembled.R2.fastq";

### print the output of fastq files for all the insertions.
open ASEMA,">$output.Insassembled.fasta";
open UASEMA,">$output.InsUnassembled.R1.fasta";
open UBSEMA,">$output.InsUnassembled.R2.fasta";
open OUT,">$output.CoverageWithAssembly.txt";
print OUT "AsemStatus\tInsertionID\tReadNumber\tRepRead\tRepQual\tRepIden\tClsReads\tClsIden\tClsQual\tRankNum\n";

while (<REQ>){
	chomp;
	my $id =(split/\t/,$_)[2];
	next if ($id eq "RepRead");
	
	if (exists $hash0{$id}){
		print ASEM "$hash0{$id}";
		print ASEMA ">$id\n$string0{$id}\n";
		print OUT "Assembly\t$_\n";
		
	}else{
		print UASEM "$hash1{$id}";
		print UBSEM "$hash2{$id}";
		
		print UASEMA ">$id\n$string1{$id}\n";
		print UBSEMA ">$id\n$string2{$id}\n";
		print OUT "UnAssembly\t$_\n";
		
	}
}

close UASEM;
close UBSEM;
close OUT;
close REQ;
close UASEMA;
close UBSEMA;



