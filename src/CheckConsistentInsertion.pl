#!/usr/bin/perl
use strict;
use warnings;
my $version="1.0 version";
use Getopt::Long;

#Author:Xin Wang 
#email: xin.wang@childrens.harvard.edu
#PI: Kaifu Chen

### This script is to confirm whether the read are identified as inserted reads and double confirm
### requirement: Final results are those inserted reads identified consistently.
#####1.	The reads are identified consistently from both method, we will export the read as confident inserted reads.  
#####2. The reads are identified from either of the method, otherwide we will have a further look at the inserted quality and identity. Here we eliminated those reads.

#### The inconsistent result came from low quality, we could easily used the information from both read and eliminate the inconsistent result due to low quality.
### The result indicated that specifically inserted reads identified from  first method showed lower quality. 



my %opts;
GetOptions(\%opts,"i:s","o:s","g:s","c:s","h:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{g}  ||!defined $opts{o} ||defined $opts{h}) {
	die "************************************************************************
	Usage: extract_fastq.pl -i InsertedReadsFromFirst -g InsertedReadsFromSecond -o FinalString
	
	Request Parameters:
	-i Inserted reads detected from the first method that did not assemble read 
	-g Inserted assembled reads detected from the second method that assemble read.
	-o Output substring of final inserted reads.
	
	Optional Parameters:
	-c The min quality requirements for inconsistent inserted reads
	-h Help
	
************************************************************************\n";
}



my $output=$opts{o};

my $FIns=$opts{i};
my $SInsA=$opts{g};
my $SInsUA=$opts{f};

my $minQ=(defined $opts{c})?$opts{c}:25;


my %fInsIndex; my %num;
open FIN,"$FIns" or die $!;
while (<FIN>){
	chomp;
	my ($id,$quality,$identity,$matches,$length,$type)=split/\t/,$_;
	$fInsIndex{$id}->{name}=$id;
	$fInsIndex{$id}->{quality}=$quality;
	$fInsIndex{$id}->{identity}=$identity;
	$fInsIndex{$id}->{matches}=$matches;
	$fInsIndex{$id}->{length}=$length;
	$fInsIndex{$id}->{type}=$type;
	$num{$id}++;
}
close FIN;

my %SInsIndex;
open SFINA,"$SInsA" or die $!;
while (<SFINA>){
	chomp;
	my ($id,$quality,$identity,$matches,$length,$type)=split/\t/,$_;
	$SInsIndex{$id}->{name}=$id;
	$SInsIndex{$id}->{quality}=$quality;
	$SInsIndex{$id}->{identity}=$identity;
	$SInsIndex{$id}->{matches}=$matches;
	$SInsIndex{$id}->{length}=$length;
	$SInsIndex{$id}->{type}=$type;
	$num{$id}++;
}
close SFINA;


open OUT, ">$output.overall.txt" or die $!;
open HIGH, ">$output.highqual.txt" or die $!;

foreach my $i (keys %num){
	if(exists $fInsIndex{$i}->{name} && exists $SInsIndex{$i}->{name}){		
		print OUT "$i\t$fInsIndex{$i}->{quality}\t$fInsIndex{$i}->{identity}\t$fInsIndex{$i}->{matches}\t$fInsIndex{$i}->{length}\tConsistent\n";
		print HIGH "$i\t$fInsIndex{$i}->{quality}\t$fInsIndex{$i}->{identity}\t$fInsIndex{$i}->{matches}\t$fInsIndex{$i}->{length}\tConsistent\n";
	}elsif(exists $fInsIndex{$i}->{name} && !exists $SInsIndex{$i}->{name}){
	
		print OUT "$i\t$fInsIndex{$i}->{quality}\t$fInsIndex{$i}->{identity}\t$fInsIndex{$i}->{matches}\t$fInsIndex{$i}->{length}\tFirst\n";
		#print HIGH "$i\t$fInsIndex{$i}->{quality}\t$fInsIndex{$i}->{identity}\t$fInsIndex{$i}->{matches}\t$fInsIndex{$i}->{length}\tFirst\n" if ($fInsIndex{$i}->{quality} >$minQ);
		
	}elsif(!exists $fInsIndex{$i}->{name} && exists $SInsIndex{$i}->{name}){
		print OUT "$i\t$SInsIndex{$i}->{quality}\t$SInsIndex{$i}->{identity}\t$SInsIndex{$i}->{matches}\t$SInsIndex{$i}->{length}\tSecond\n";
		#print HIGH "$i\t$SInsIndex{$i}->{quality}\t$SInsIndex{$i}->{identity}\t$SInsIndex{$i}->{matches}\t$SInsIndex{$i}->{length}\tSecond\n" if ($SInsIndex{$i}->{quality} >$minQ);
	
	}	
	
}
close OUT;
close HIGH;
