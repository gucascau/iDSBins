#!/usr/bin/perl
use strict;
use warnings;


### This script is to extract the reads that consistently showed non insertion events:
### requirement: 1
#####1.	Completely blasted against reference genome from start to end of the reads (allow for 2 gaps and 2 mismatches at most, Identify is more than 95%);
#####2.	Forward and Reward regards consistently mapped against the identical reference location.


#use Bio::Seq;
#use Bio::SeqIO;


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s","g:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{o} || !defined $opts{g}) {
	die "************************************************************************
	Usage: extract_fastq.pl -i BLAST -g allow_gap_size_for_PE -o Blast_categories
************************************************************************\n";
}

my $index=$opts{i};
my $output=$opts{o};
my $gap=$opts{g};

open G,"$index" or die $!;
open O,">$output.oneinsert.txt" or die $!;
open CO,">$output.oneinsert_withgap.txt" or die $!;
open TI, ">$output.twoinsert.txt" or die $!;

print O "qID\tPID\tInsertsize\tFstart\tFend\tRstart\tRend\tFqstart\tRqstart\tFqend\tRqend\n";
my %inf2; my %confInsert2;

while (<G>) {
	chomp;
    # print "$_" ;
	my ($qid,$pid,$identity,$match,$qgap,$mismatch,$qstart,$qend,$qlength,$pstart,$pend,$plength,$score,$evlaue)=split/\s+/,$_;
	my $information=$_;
	my $start=($pstart<$pend)?$pstart:$pend;
	my $end=($pstart<$pend)?$pend:$pstart;

	if (exists $confInsert2{$qid}){
		
		# if ($pid eq $inf2{$qid}->{chr}){
		# 	$size=$inf2{$qid}->{end}-$start+1;
		# }
		if ($pid eq $inf2{$qid}->{chr} && (($start>=$inf2{$qid}->{start} && $start<=$inf2{$qid}->{end}) || ($inf2{$qid}->{start}>=$start && $inf2{$qid}->{start} <=$end)) ){
			my @array= sort { $a <=> $b }($start,$inf2{$qid}->{start},$end,$inf2{$qid}->{end});
			my $length=$array[3]-$array[0];
			my $str=join "\t",@array;
			print O "$qid\t$pid\t$length\t$str\t$inf2{$qid}->{qstart}\t$qstart\t$inf2{$qid}->{qend}\t$qend\t$qlength\n";
		}elsif ($pid eq $inf2{$qid}->{chr} && $start > $inf2{$qid}->{end} && ($start - $inf2{$qid}->{end}) < $gap){
			my @array= sort { $a <=> $b }($start,$inf2{$qid}->{start},$end,$inf2{$qid}->{end});
			my $length=$array[3]-$array[0];
			my $str=join "\t",@array;
			print CO "$qid\t$pid\t$length\t$str\t$inf2{$qid}->{qstart}\t$qstart\t$inf2{$qid}->{qend}\t$qend\t$qlength\n";
        }elsif ($pid eq $inf2{$qid}->{chr} && $inf2{$qid}->{start} > $end && ($inf2{$qid}->{start}-$end) <$gap){
        	
			my @array= sort { $a <=> $b }($start,$inf2{$qid}->{start},$end,$inf2{$qid}->{end});
			my $length=$array[3]-$array[0];
			my $str=join "\t",@array;
			print CO "$qid\t$pid\t$length\t$str\t$inf2{$qid}->{qstart}\t$qstart\t$inf2{$qid}->{qend}\t$qend\t$qlength\n";
		
		
		}else{
			print TI "$confInsert2{$qid}\n$information\n";
		}
		
		$inf2{$qid}->{chr}=$pid;
		$inf2{$qid}->{start}=$start;
		$inf2{$qid}->{end}=$end;
		$inf2{$qid}->{qstart}=$qstart;
		$inf2{$qid}->{qend}=$qend;
		$confInsert2{$qid}=$information;
	}else{
		$inf2{$qid}->{chr}=$pid;
		$inf2{$qid}->{start}=$start;
		$inf2{$qid}->{end}=$end;
		$inf2{$qid}->{qstart}=$qstart;
		$inf2{$qid}->{qend}=$qend;
		$confInsert2{$qid}=$information;
		#print "$information\n";
	}

}
close G;
close CO;
close O;
close TI;

