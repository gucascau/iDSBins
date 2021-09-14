#!/usr/bin/perl
#author:wangxin
### #### caluculate the side effect on the final associated genes results
### Here is to calculate the sample size effect
use strict;
use warnings;
#use Statistics::Test::WilcoxonRankSum;
#use List::Util qw(sum);
#use List::Util 'shuffle';

my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","s:s","f:s","g:s");
print "*************\n*$version*\n*************\n";
if ( !defined $opts{i}|| !defined $opts{s} ||!defined $opts{f}||!defined $opts{g}  ) {
       	die "************************************************************************
       	Usage: $0.pl
				-i: Yeast id index
				-g: insert sequences
			-f: Blast results
			-s: Blast bed file with sequence and coverage information
************************************************************************\n";
}


my $input=$opts{i};

my $insert=$opts{f};
### the sequence index
my $fasta=$opts{g};

my %alter;

##### Here is the yeast genomic information Chromosme information
open IN,"$input" or die $!;
while (<IN>){
    chomp;
	s/\r//g;
    my ($id,$chr) =(split /\t/,$_);
    #next if (exists $exist{$id});
	
	$alter{$id}=$chr;
    #$exist{$id}++;
    #$print PH1 ">$id\n$string1{$id}\n";
    #print PH2 ">$id\n$string2{$id}\n";
    #$phix{$_}++;
}
close IN;


#### Here is the sequence of each insert ID


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


my $out=$opts{s};


open INS,"$insert" or die "cannot open file $insert";
open OUT,">$out" or die $!;
my $n=0;

while (<INS>) {
	chomp;
	s/\r//g;
    # print "$_" ;

	my @array=split/\s+/,$_;
	my $identity=$array[2];	
	my $string=$str{$array[0]};
	$n++;
	my $type=($array[10]>$array[9])?"+":"-";
	my $start=($array[10]>$array[9])?$array[9]:$array[10];
	my $end=($array[10]>$array[9])?$array[10]:$array[9];
	my $chro=$alter{$array[1]};

	my $length1=$array[6];
	my $length2=$array[8]-$array[7]+1;
	my $coverage=$array[14];
	if ($length1 <=55 && $length2<=55){
		print OUT "$chro\t$start\t$end\tIM$n\t0\t$type\t$coverage\t$identity\t1\t$string\n";
	}elsif ($length1 >55 && $length2>55){
		print OUT "$chro\t$start\t$end\tIM$n\t0\t$type\t$coverage\t$identity\t3\t$string\n";
	}else{
		print OUT "$chro\t$start\t$end\tIM$n\t0\t$type\t$coverage\t$identity\t2\t$string\n";
	}
}
close INS;
close OUT;


