#!/usr/bin/perl
use strict;
use warnings;



#use Bio::Seq;
#use Bio::SeqIO;

###  This script is to extract the clustering sequence (fasta) and blast results for each cluster
#### 

my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"f:s","b:s","o:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{f} ||!defined $opts{b}|| !defined $opts{o}) {
	die "************************************************************************
	Usage: $0.pl 
	-f original fasta(MAT and unmappable sequence with coverage information) 
	-b blast information for each sequence (Here we only used all BLAST results of each ID) 
	-o Output string to seperate MAT and unmappable sequences
************************************************************************\n";
}


#### Told the input data
#my $output=$opts{o};
my $fasta=$opts{f};
#my $dir=$opts{d};
#my $FCL=$opts{c};
#my $start=$opts{t};
my $blast=$opts{b};
#my $input=$opts{d};
my $out=$opts{o};



## read orginal fasta file:
open FASTA,"$fasta" or die $!;
my $id1;  my %hash;  my %cov; my $qual;
while (<FASTA>) {
	chomp;
    # print "$_" ;
    if (/>(\S+)\s+(\S+)\s+(\S+)/)  {
	
		$id1=$1;
		$cov{$id1}=$2;
		$qual{$id1}=$3;
	}else{
		$hash{$id1}.=$_;
	
	}
}
close FASTA;

#### index of blast results 
open BLAST,"$blast" or die $!;

my %str;
while (<BLAST>) {
	chomp;
    # print "$_" ;
	my ($id2,$rid,$start,$end,$length,$rstart,$rend)=(split/\t/,$_)[0,1,6,7,8,9,10];
	next unless (exists $hash{$id2});
	if($start >30 && ($length-$end)>30){
		$str{$id2}=$_;
	}
}
close BLAST;

open MAT,">$out.mat.txt" or die $!;
open NMAP,">$out.unmappable.fasta" or die $!;
foreach my $i (keys %hash){
	if (exists $str{$i}){
		print MAT "$str{$i}\t$cov{$i}\t$qual{$i}\n";
	}else{
		print NMAP ">$i\t$cov{$i}\t$qual{$i}\n$hash{$i}\n";
	}
}

close MAT;
close NMAP;







