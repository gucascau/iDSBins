#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#use Bio::Seq;
#use Bio::SeqIO;


my $version="1.0 version";

#### This is to create the index of final blast result for further annotation

### Transfer blast result to bedfile
#### create bedfile for blast

my %opts;
GetOptions(\%opts,"i:s","o:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{o} ) {
	die "************************************************************************
	Usage: $0.pl -i Blast result (with mutiple insertions)
				-o Further annotation for mutiple blast
************************************************************************\n";
}

my $output=$opts{o};
my $input=$opts{i};



open I,"$input" or die $!;

my $f=0;

open OUT,">$output" or die $!;

while (<I>) {
	chomp;
    # print "$_" ;
	my @array=split/\s+/,$_;
	my $information=$_;
	#$basic1{$qid}=$information;
	#next if ($array[6] <5 || ($array[8]-$array[7]) <5);
	#next if ($identity <90);
	
	my $chr2=$array[1];
	my ($identity,$bitscore,$pvalue)=@array[2,12,13];
	my $type=($array[10]>$array[9])?"+":"-";
	my $start=($array[10]>$array[9])?$array[9]:$array[10];
	my $end=($array[10]>$array[9])?$array[10]:$array[9];
	$f++;
			
	print OUT "$chr2\t$start\t$end\t$array[0]\t0\t$type\n";
	
}

close OUT;
close I;
	
