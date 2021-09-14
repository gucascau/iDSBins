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
GetOptions(\%opts,"i:s","o:s","a:s","b:s","q:s","s:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{o}|| !defined $opts{a} || !defined $opts{q} || !defined $opts{b}  ) {
	die "************************************************************************
	Usage: $0.pl -i Single_estimated_insertion.txt -a gapsize of reads -b gapsize of mapping regions -s aligned shift allowed at end of R1 or the beging of R2 (default 5bp) -q quality file of each reads -o Blast_statatistic 
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
my $shift= (defined $opts{s})?$opts{s}:5;

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


my %num; my %match; my %hash;
open I,"$input" or die $!;
my $n=1;
while (<I>){
	chomp;
	s/\r//g;
	
	
	#### here four locus take into consideration : refercence locus: pstart, pend, and reads locus: qstart, qend
	my ($cov,$qid,$pid,$identity,$length,$pstart,$pmiddle1,$pmiddle2,$pend,$qstart,$qmiddle1,$qmiddle2,$qend,$qlengthf,$qlengthr,$quality,$strand)=split/\s+/,$_;
	
	
	
	#my ($cov,$qid,$pid,$identity,$length,$pstart,$pmiddle1,$pmiddle2,$pend,$qstart,$qmiddle1,$qmiddle2,$qend,$quality)=split/\s+/,$_;
	
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


### we set up stringent parameter for these reads that if they have 5bp inner unknown sequence, we will not consider these event
	next unless ($qmiddle1<=$shift && (($qlengthf-$qmiddle2) <=$shift));
	my $m=0;
	my $q=$n-1;
	#my $pid_max=$pid;
	foreach my $i (1..$q){	
		if (($pid eq $hash{$i}->{pid}) && (abs($qstart-$hash{$i}->{qstart}) <=$gapA) && (abs($qend-$hash{$i}->{qend}) <=$gapA)  && (abs($pstart-$hash{$i}->{pstart}) <=$gapB) && (abs($pend-$hash{$i}->{pend}) <=$gapB)){
			$hash{$i}->{num} +=$cov;
			$hash{$i}->{string}.="\t".$qid;
			$m++;
			#print "$qid\t$m\n";
				
			### Here we also considered the quality. If the coverage the same, we only considered the high quality one
			$hash{$i}->{inf}=$information if ($cov > $hash{$i}->{cov} || ($cov == $hash{$i}->{cov} && $qual{$qid} > $hash{$i}->{qual}));
			$hash{$i}->{cov}=$cov if ($cov > $hash{$i}->{cov} || ($cov == $hash{$i}->{cov} && $qual{$qid} >$hash{$i}->{qual}));
			$hash{$i}->{qual}=$qual{$qid} if ($cov > $hash{$i}->{cov} || ($cov == $hash{$i}->{cov} && $qual{$qid} >$hash{$i}->{qual}));
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
foreach my $num(keys %hash){

	print  OUT "$hash{$num}->{inf}\t$hash{$num}->{num}\n";
	
	my @array=split/\t/,$hash{$num}->{string};
	my $fid=join ";",@array;
	my ($Reid,$Rident)=(split/\t/,$hash{$num}->{inf})[1,3];
	my $Requal=$qual{$Reid};
	my $fqual='';
	foreach my $i (@array){
		$fqual.= $qual{$i}.";";
		
	}
	print CLT "SecondFinal\tID$num\t$hash{$num}->{num}\t$Reid\t$Requal\t$Rident\t$fid\t$fqual\n";
}

close OUT;
close CLT;


