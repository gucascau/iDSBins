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
GetOptions(\%opts,"i:s","o:s","g:s","m:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{g} || !defined $opts{m} || !defined $opts{o} ) {
	die "************************************************************************
	Usage: $0.pl -i Blast result (with mutiple insertions)
		 		-g chromosome id changes
				-m sequence with id infomation
				-o Further annotation for mutiple blast
************************************************************************\n";
}

my $output=$opts{o};
my $input=$opts{i};
my $alteration=$opts{g};
my $sequence=$opts{m};


#### index for chromosome id ####

my %alter;

open A,"$alteration" or die "cannot open file $alteration";
while (<A>){
	chomp;
	s/\r//g;
	my ($i,$ch)=split /\t/,$_;
	
	$alter{$i}=$ch;
	
}
close A;




### Index for sequence id #####

my %seq;
my $name;
open S,"$sequence" or die "cannot open file $sequence";
while (<S>){
	chomp;
	s/\r//g;
	if (/>(\S+)/){
		$name=$1;
	}else{
		$seq{$name}.=$_;
	}
	
}

close S;



open I,"$input" or die $!;

my $f=0;

open OUT,">$output" or die $!;

while (<I>) {
	chomp;
    # print "$_" ;
	my @array=split/\s+/,$_;
	my $information=$_;
	#$basic1{$qid}=$information;
	next if ($array[6] <5 || ($array[8]-$array[7]) <5);
	#next if ($pid eq "ref-NC_001135-" && $pstart >= 13650 && $pend<= 13850 );
	#next if ($pid eq "ref-NC_001135-" && $pstart >= 200750 && $pend <= 201000);
	#next if ($identity <90);
	my $chr=$array[1];
	
	my $chr2=$alter{$chr};
	my ($identity,$bitscore,$pvalue)=@array[2,12,13];
	my $type=($array[10]>$array[9])?"+":"-";
	my $start=($array[10]>$array[9])?$array[9]:$array[10];
	my $end=($array[10]>$array[9])?$array[10]:$array[9];
	$f++;
	
	my $sequence=$seq{$array[0]};
		
	print OUT "$chr2\t$start\t$end\tID$f\t0\t$type\t$identity\t$bitscore\t$pvalue\t$sequence\n";
	
	
}

close OUT;
close I;
	
