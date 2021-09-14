#!/usr/bin/perl
use strict;
use warnings;


#### This scripts is to add the extral Genomic feature with ENTIRELY, PARTIALLY, and PROXIMITY




my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","b:s","o:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{b} || !defined $opts{o}  ) {
	die "************************************************************************
	Usage: $0.pl -i Requested insertion file
				-b Annotation Bed file
				-o Output of insertion file with Proximity/Partially/ Entirely
************************************************************************\n";
}



my $insert=$opts{i};


my $annotaton=$opts{b};
my $output=$opts{o};


my $n=0;
my %hash;

open AN,"$annotaton" or die "cannot open file $annotaton";
while (<AN>){
	chomp;
	
	my @arrayA=split/\t/,$_;
	$hash{$arrayA[6]}->{Start}=$arrayA[1];
	$hash{$arrayA[6]}->{End}=$arrayA[2];
	
}


open IN,"$insert" or die "cannot open file $insert";
open OUT,">$output" or die $!;
while (<IN>){
	chomp;
	$n++;
	my @array=split/\t/,$_;
	if ($n==1){
		print OUT "Case\tNID\tGenotype\tWholeinsertion\tInseredSeq\tInsertionSize\tDonorNum\tIdentity\tReadCounts\tDonorChromosome\tDonorStart\tDonorEnd\tStrand\tAnnotationFeature\tDistanceToFeature\tFeatureDiscription\tAlterFeature\tQuality\n";
		next;
	}
	my $stringA=join "\t",@array[0..8]; 
	my $stringB=join "\t",@array[11..15];
	# my $stringB=join "\t",@array[11..17];

	my $quality=($array[22]+$array[24])/2;
	my $feature;

	my $StartE=$hash{$array[15]}->{Start};
	my $EndE=$hash{$array[15]}->{End};
	
	my $StartI=$array[12];
	my $EndI=$array[13];
	
	next if ($array[8]<=2 && $quality<=20);
	
	
	if ($array[12] eq "NO"){
		
		### here is for two donor intergrated without any information	
		$feature="NO";
	}elsif ($array[16] != 0){
		$feature="PROXIMITY";
	}elsif($StartI>=$StartE && $EndI<=$EndE){
		$feature="ENTIRELY";
	}else{
		$feature="PARTIALLY";
	}
	my $stringC= join "\t",($feature,$array[16],$array[17],$quality);
	
	print OUT "$stringA\t$stringB\t$stringC\n";
	
	
}

