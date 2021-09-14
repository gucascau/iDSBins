#!/usr/bin/perl
#author:wangxin
### function: This script is to extract the potential elimanated insertions and added the removed insertion coverage into previous insertions
### The script is based on the previous self treatment 


use strict;
use warnings;

my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","g:s","o:s");
print "*************\n*$version*\n*************\n";
if ( !defined $opts{i}|| !defined $opts{g} || !defined $opts{o}) {
       	die "************************************************************************
       	Usage: $0.pl
	       		-i: Potential removed duplicated list (yWH475-1_S10_all.potential.txt)
			-g: Orignal insertion list (yWH475-1_S10_all.fasta)
			-o: Insertion list after removing duplicated insertion(Output.txt). We also added the removed coverage into the orignal insertions
************************************************************************\n";
}

my $duplicates=$opts{i};
my $origin=$opts{g};
my $output=$opts{o};


my %remove; my %con; my %contain; my %hash; my %name; my $n=0;
open D, "$duplicates" or die $!;
while (<D>){
	chomp;
	s/\r//;
	my ($id1,$id2)=(split/\t/,$_)[0,1];
	if (!exists $name{$id1} && !exists $name{$id2}){
		$n++;
		
		
		#push @{$name{$id1}},($id1,$id2);
		my $string=join "\t",($id1,$id2);
		$name{$id1}=$n;
		$name{$id2}=$n;
		$hash{$n}=$string;
	}elsif (exists $name{$id1} && !exists $name{$id2} ){
		my $num=$name{$id1};
		#push @{$name{$id1}},$id2;
		$name{$id2}=$num;
		$hash{$num}.="\t$id2";
	}elsif (exists $name{$id2} && !exists $name{$id1} ){
		my $num=$name{$id2};
		#push @{$name{$id2}},$id1;
		$hash{$num}.="\t$id1";
		$name{$id1}=$num;
	}
	
	

	
}
close D;

open CLS,">$output.cls" or die $!;
foreach my $i (sort keys %hash){
	my $max; my $n1Tp=0; my $ST; my $SA; my $SB; my $SC; my $SD;
	print CLS "$hash{$i}\n";
	my @array=split /\t/,$hash{$i};
	
	foreach my $p (@array){
		my ($n1,$n1T)=(split/\|/,$p)[0,3];
	
		if ($n1T>$n1Tp){
			$max=$p;
			$n1Tp=$n1T;
		}
		
		$ST += $n1T;

		
	}
	

	$contain{$max}=$ST;
		
	foreach my $q (@array){
		## my ($n1,$n1T)=(split/\|/,$q)[0,3];
		
		if ($q eq $max){
			$con{$q}++;
		}else{
			$remove{$q}++;
		}
		
	}

	
}
	
close CLS;
	
	


close D;

### index of fasta 
my %str; my $idN;
open G,"$origin" or die $!;
open OUT,">$output" or die $!;
while (<G>){
	chomp;
	s/\r//;
	if (/>(\S+)/){
		$idN=$1;
		
	}else{
		$str{$idN}=$_;
	}
	
}
close G;

foreach my $f (keys %str){
	next if (exists $remove{$f});
	if (exists $con{$f}){
		my @inf=split/\|/,$f;
		$inf[3]=$contain{$f};
		my $string=join "|",@inf;
		print OUT "$string\t$str{$f}\n";
		
	}else{
		print OUT "$f\t$str{$f}\n";
	}
	
}

close G;
close OUT;



