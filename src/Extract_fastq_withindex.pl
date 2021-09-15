#!/usr/bin/perl
use strict;
use warnings;
my $version="1.0 version";
use Getopt::Long;

# Copyright (c) 2021 Dr. Kaifu lab    
# Author:Xin Wang 
# email: xin.wang@childrens.harvard.edu
# PI: Kaifu Chen

# Function: This script is to eliminate the reads that have error index


my %opts;
GetOptions(\%opts,"f:s","r:s","o:s","i:s","g:s","h:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{f} || !defined $opts{r} || !defined $opts{i} ||!defined $opts{o}||!defined $opts{g} || defined $opts{h}) {
	die "************************************************************************
	Usage: $0.pl -f ffastq -r  Rfastq -i forward_index -g reverse_index -o Fastq_withindex
	
	### 

	Request Parameters:
	-f Forward reads (fastq format)
	-r Reverse reads (fastq format)
	-i Forward index
	-g Reverse index
	-o Output files with corrected index, including corrected index fastq files and error index fastq files.
	
	Optional Parameters:
	-h Help
************************************************************************\n";
}

my $output=$opts{o};
my $fastq1=$opts{f};
my $fastq2=$opts{r};

my $str1=$opts{i};
my $str2=$opts{g};

open FASTQ1,"$fastq1" or die "cannot open file $fastq1";

my $id1; my %string1; my %hash1;
while (<FASTQ1>) {
	chomp;
    # print "$_" ;

    if ($. % 4 == 1)  {
		$id1=(split/\s+/,$_)[0];
		$id1=~s/^@//;
		print "$id1\n";
		$hash1{$id1}=$_."\n";
	}elsif($.%4 == 2){
		$string1{$id1}=$_;
		$hash1{$id1}.=$_."\n";
		
	}else{
		$hash1{$id1}.=$_."\n";
	}

}
close FASTQ1;

my $id2; my %string2; my %hash2; 
open FASTQ2, "$fastq2" or die "cannot open file $fastq2";

while (<FASTQ2>) {
	chomp;
    # print "$_" ;

    if ($. % 4 == 1)  {
		$id2=(split/\s+/,$_)[0];
		$id2=~s/^@//;
		$hash2{$id2}=$_."\n";
	}elsif($.%4 == 2){
		$string2{$id2}=$_;
		$hash2{$id2}.=$_."\n";
		#$str2{$id2}=$_;
		
	}else{
		$hash2{$id2}.=$_."\n";
	}
}
close FASTQ2;

open PAIR1,">$output.R1.fastq" or die $!;
open PAIR2 ,">$output.R2.fastq" or die $!;

open UPAIR1,">$output.unindex.R1.fastq" or die $!;
open UPAIR2 ,">$output.unindex.R2.fastq" or die $!;
foreach my $id (keys %hash1){
	#print "$string1{$id}\n";
	if ($string1{$id}=~/^$str1/m && $string2{$id}=~/^$str2/m){
		print PAIR1 "$hash1{$id}";
		print PAIR2 "$hash2{$id}";		
	}else{
		print UPAIR1 "$hash1{$id}";
		print UPAIR2 "$hash2{$id}";
	}
	
}	
	


close PAIR1;
close PAIR2;
close UPAIR1;
close UPAIR2;


