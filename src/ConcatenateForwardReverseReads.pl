#!/usr/bin/perl
use strict;
use warnings;
#use Data::Dump qw(dump);
#Author:Xin Wang 
#email: xin.wang@childrens.harvard.edu
#PI: Kaifu Chen

#use Bio::Seq;
#use Bio::SeqIO;


#### the script is to concacanate the forward reads and reverse reads

my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"f:s","o:s","r:s","h:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{f} || !defined $opts{o}  || !defined $opts{r} || defined $opts{h}) {
	die "************************************************************************
	Usage: extract_fastq.pl -f forward fasta -r reverse fastq  -o Output of concatenarate sequences
	
	Request Parameters:
	-f Forward fastq file
	-r Reverse fastq file
	-o Concatenated forward and reverse reads
	
	-h Help
************************************************************************\n";
}

my $output=$opts{o};
my $ffastq=$opts{f};
my $rfastq=$opts{r};


open FFASTQ,"$ffastq" or die $!;

my $id1;  my %Fhash;  
while (<FFASTQ>) {
	chomp;
    # print "$_" ;

    if ($. % 4 == 1)  {
		$id1=(split/\s+/,$_)[0];
		$id1=~s/^@//;
		#print "$id1\n";
		
	}elsif($. % 4 == 2){
		$Fhash{$id1}=$_;
	}
	
}
close FFASTQ;


open RFASTQ,"$rfastq" or die $!;

my $id2;  my %Rhash;  
while (<RFASTQ>) {
	chomp;
    # print "$_" ;

    if ($. % 4 == 1)  {
		$id2=(split/\s+/,$_)[0];
		$id2=~s/^@//;
		#print "$id1\n";
		
	}elsif($. % 4 == 2){
		$Rhash{$id2}=$_;
	}
	
}
close RFASTQ;




open OUT,">$output.Concatenated.fasta" or die $!;
#open REV,">$output.rev.fasta" or die $!;
foreach my $i (keys %Fhash){
	
	if(length $Rhash{$i} >=290 && length $Fhash{$i} >=290){
		
		my $string=join "",($Fhash{$i},$Rhash{$i});
		print OUT ">$i\n$string\n";	
	}
	
}
close OUT;


