#!/usr/bin/perl
#author:wangxin
### #### Transfer estimated donor into bed file with sequence information
#
use strict;
use warnings;
#use Statistics::Test::WilcoxonRankSum;
#use List::Util qw(sum);
#use List::Util 'shuffle';

my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","s:s","f:s","g:s","m:s","n:s");
print "*************\n*$version*\n*************\n";
if ( !defined $opts{i}|| !defined $opts{s}|| !defined $opts{m} || !defined $opts{n}  ||!defined $opts{f}||!defined $opts{g}  ) {
       	die "************************************************************************
       	Usage: $0.pl
				-i: Yeast id index
			-g: Yeast genome sequence
			-m: forward reads
			-n: reverse reads
			-f: Estimated two donor Blast results
			-s: Blast bed file with sequence and coverage information
************************************************************************\n";
}


my $input=$opts{i};

my $estimate=$opts{f};
### the sequence index
my $fasta=$opts{g};

my $fastq1=$opts{m};
my $fastq2=$opts{n};

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



open FASTQ1, "$fastq1" or die "cannot open file $fastq1";

my %strf; my $id4;
while (<FASTQ1>){
	chomp;
	
    if ($. % 4 == 1)  {
  		$id4=(split/\s+/,$_)[0];
  		$id4=~s/^@//;
  	}elsif($.%4 == 2){
  		$strf{$id4}=$_;
  		#$hash1{$id1}.=$_."\n";
	}
    
}
close FASTQ1;

open FASTQ2, "$fastq2" or die "cannot open file $fastq2";

my %strr; my $id5;
while (<FASTQ2>){
	chomp;
	
    if ($. % 4 == 1)  {
  		$id5=(split/\s+/,$_)[0];
  		$id5=~s/^@//;
  	}elsif($.%4 == 2){
  		$strr{$id5}=$_;
  		#$hash1{$id1}.=$_."\n";
	}
    
}
close FASTQ2;


my $out=$opts{s};




open EST,"$estimate" or die "cannot open file $estimate";
open OUT,">$out" or die $!;
my $n=0;

while (<EST>) {
	chomp;
	s/\r//g;
    # print "$_" ;
	my @array=split/\s+/,$_;
	#my $string=$str{$array[0]};
	$n++;
	my $type=($array[7]>$array[4])?"+":"-";
	my $start=($array[7]>$array[4])?$array[4]:$array[7];
	my $end=($array[7]>$array[4])?$array[7]:$array[4];
	my $chro=$alter{$array[2]};
	my $coverage=$array[0];
	my $length=$end-$start+1;
	
	my $length1=$array[8];
	my $string1=substr $strf{$array[1]},0,$length1;
	
	my $string2=substr $str{$chro},$start,$length;
	
	my $length3=$array[12]-$array[11]+1;
	my $string3=substr $strr{$array[1]},-$length3;
	
	my $string=join '', ($string1,$string2,$string3);
	
	if ($length1<=55 && $length3 <=55){
		print OUT "$chro\t$start\t$end\tES$n\t0\t$type\t$coverage\t100\t1\t$string\n";
	}elsif ($length1>55 && $length3 >55 ){
		print OUT "$chro\t$start\t$end\tEM$n\t0\t$type\t$coverage\t100\t3\t$string\n";		
	}else{
		print OUT "$chro\t$start\t$end\tEM$n\t0\t$type\t$coverage\t100\t2\t$string\n";
	}
}
close EST;
close OUT;


