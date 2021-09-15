#!/usr/bin/perl
use strict;
use warnings;

#### 
### This script is to extract Uniq blast results (single donors) and statistic their information.
### requirement: 1
#####1

my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s","a:s","b:s","q:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{o} || !defined $opts{q} ) {
	die "************************************************************************
	Usage: $0.pl -i Single_insertion.txt -q quality file of each reads -o Unique insertions output string
	
	Request Parameters for single donor deduplication:
	-i Insertion with detected one donor
	-q The quality information of each read
	-o Output of unique insertion that only detected with one donor, this including the clusting file and unique insertion event file
	
	Optional Parameters:
	-a The difference allowed between two reads, start site, end start site on the inserted reads (default: 6)
	-b The difference allowed between two reads, start site, end start site on the chromosome locus of donor (default: 6)
	-h Help
************************************************************************\n";
}


### write a report for first round of insertion files
### number of reads in different chromsomes including duplciates reads
### number of uniq events in different chromosome using the uniq reads

### inserted length distribution
### uniq inserted length distribution
my $input=$opts{i};
my $output=$opts{o};

## the gap size between 
my $gapA=(defined $opts{a})?$opts{a}:6;
my $gapB=(defined $opts{b})?$opts{b}:6;


### read the quality file
my $Quality=$opts{q};

my %qual;
open QUAL,"$Quality" or die $!;
while (<QUAL>){
	chomp;
	my ($id,$qu)=(split/\t/,$_)[0,5];
	$qual{$id}=$qu;
}
close QUAL;

#########################################################
#### perform the clustering based on the 
#########################################################
my %num; my %match; my %hash;my %Rcount; my %Rident; my %inf;
open I,"$input" or die $!;
my $n=1;
while (<I>){
	chomp;
	s/\r//g;
	my ($cov,$qid,$pid,$identity,$length,$mismatches,$gapopen,$qstart,$qend,$qlength,$pstart,$pend,$plength,$score,$pvalue)=split/\s+/,$_;
	$Rcount{$qid}=$cov;
	$Rident{$qid}=$identity;
	#$Fqstart,$Rqstart,$Fqend,$Rqend,$Fpstart,$Rpstart,$Fpend,$Rpend,$Flength,$Rlength)=split/\s+/,$_;
	my $information=$_;
	$inf{$qid}=$information;
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
			#print "$qid\t$m\n";
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
		$hash{$n}->{qual}=$qual{$qid};
		$n++;
		
	}
}

close I;

#### 
open OUT,">$output.uniq.txt" or die $!;
open CLT,">$output.cls.txt" or die $!;

my %Fnum;
foreach my $num(keys %hash){
	
	my @array=split/\t/,$hash{$num}->{string};
	
	### Here we also considered the quality. If the coverage the same, we only considered the high quality one
	### identify the representive reads
	my $fqual; my $fiden; my $fcov; my $fid; my $FRcount;
	
	my %Lqual=(); my %Lcov=(); my %Liden=();
	foreach my $i (@array){
		$Lqual{$i}=$qual{$i};
		$Lcov{$i}=$Rcount{$i};
		$Liden{$i}=$Rident{$i};
		
		$FRcount += $Rcount{$i};
		$fid.=$i.";";
		$fcov.=$Rcount{$i}.";";
		$fiden.=$Rident{$i}.";";
		$fqual.= $qual{$i}.";";
	}
	
	my @keyF = sort {$Lcov{$b} <=> $Lcov{$a} || $Liden{$b} <=> $Liden{$a} || $Lqual{$b} <=> $Lqual{$a} } keys %Lcov;
	

	
	my $RepRead=$keyF[0];
	next if (exists $Fnum{$RepRead});
	$Fnum{$RepRead}++;
	
	print  OUT "$inf{$RepRead}\t$FRcount\n";

	print CLT "SecondFinal\tID$num\t$FRcount\t$RepRead\t$qual{$RepRead}\t$Rident{$RepRead}\t$fid\t$fqual\t$fiden\t$fcov\n";
}

close OUT;
close CLT;


