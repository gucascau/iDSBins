#!/usr/bin/perl
use strict;
use warnings;

#### 
### This script is to extract Uniq blast results (Single donors) and statistic their information.
### requirement: 1
#####1

my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s","a:s","b:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{o}|| !defined $opts{a} || !defined $opts{b}   ) {
	die "************************************************************************
	Usage: $0.pl -i Single_insertion.txt -a gapsize of reads -b gapsize of mapping regions -o Blast_statatistic 
************************************************************************\n";
}


### write a report for first round of insertion files
### number of reads in different chromsomes including duplciates reads
### number of uniq events in different chromosome using the uniq reads

### inserted length distribution
### uniq inserted length distribution
my $input=$opts{i};
my $output=$opts{o};

my $gapA=$opts{a};
my $gapB=$opts{b};

my %num; my %match; my %hash;
open I,"$input" or die $!;
my $n=1;
while (<I>){
	chomp;
	s/\r//g;
	my ($cov,$qid,$pid,$identity,$length,$mismatches,$gapopen,$qstart,$qend,$qlength,$pstart,$pend,$plength,$score,$pvalue)=split/\s+/,$_;
	
	#$Fqstart,$Rqstart,$Fqend,$Rqend,$Fpstart,$Rpstart,$Fpend,$Rpend,$Flength,$Rlength)=split/\s+/,$_;
	my $information=$_;
	### split into two categories: one had no other insertions
	###  							another had 
	next if ($qid eq "qID");
	#next if (exists $hash{$qid});
	$pstart=($pstart<$pend)?$pstart:$pend; 
	$pend=($pstart<$pend)?$pend:$pstart;
	# my $nstart=$Fqstart;
# 	my $nend=$Rqend;
	my $m=0;
	my $q=$n-1;
	#my $pid_max=$pid;
	foreach my $i (1..$q){	
		if ($pid eq $hash{$i}->{pid} && (abs($qstart-$hash{$i}->{qstart}) <=$gapA) && (abs($qend-$hash{$i}->{qend}) <=$gapA)  && (abs($pstart-$hash{$i}->{pstart}) <=$gapB) && (abs($pend-$hash{$i}->{pend}) <=$gapB)){
			$hash{$i}->{num} +=$cov;
			$hash{$i}->{string}.="\t".$qid;
			$m++;
			print "$qid\t$m\n";
			
			$hash{$i}->{inf}=$information if ($cov > $hash{$i}->{cov} || ($cov == $hash{$i}->{cov} && $qual{$pid} >$hash{$i}->{qual}));
			$hash{$i}->{cov}=$cov if ($cov > $hash{$i}->{cov} || ($cov == $hash{$i}->{cov} && $qual{$pid} >$hash{$i}->{qual}));
		}
	
	}
	
	if ($m==0){

		#print "$n\n";
		$hash{$n}->{cov}=$cov;
		$hash{$n}->{pid}=$pid;
		$hash{$n}->{qstart}=$qstart;
		$hash{$n}->{qend}=$qend;
		$hash{$n}->{pend}=$pend;
		$hash{$n}->{pstart}=$pstart;
		$hash{$n}->{inf}=$information;
		$hash{$n}->{num}=$cov;
		$hash{$n}->{string}=$qid;
		$n++;
		
	}
}

close I;

#### 
open OUT,">$output.uniq.txt" or die $!;
open CLT,">$output.cls.txt" or die $!;
foreach my $num(keys %hash){

	print  OUT "$hash{$num}->{inf}\t$hash{$num}->{num}\n";
	print CLT "$hash{$num}->{num}\t$hash{$num}->{string}\n";
}

close OUT;
close CLT;


