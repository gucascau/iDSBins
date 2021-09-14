#!/usr/bin/perl
#author:wangxin
### function:  The script is to add both quality for each annotated insertions
use strict;
use warnings;


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","a:s","b:s","o:s","r:s");
print "*************\n*$version*\n*************\n";
if ( !defined $opts{i}|| !defined $opts{a} || !defined $opts{b}  || !defined $opts{o} || !defined $opts{r}) {
       	die "************************************************************************
       	Usage: $0.pl
	       		-i: the final insertion events (yWH475-1_S10_combined_dedupe_div.single_annotation.txt)
			-r: reads with the quality score (yWH475-1_S10_QC.highquality.stat)
			-a: aligned raw insertions with reads inf (All_mappable_eval.txt)
			-b: unaligned raw insertions with reads inf (All_unmappable_eval.txt)
			-o: the final insertion events adding reads quality
************************************************************************\n";
}


my $input=$opts{i};

my $aligned=$opts{a};
my $unaligned=$opts{b};

my $Readq=$opts{r};

my $output=$opts{o};

my %inf;

open QU, "$Readq" or die "cannot open file $Readq";
my %qual;

### add the read quality inforamtion;
while (<QU>){
	chomp;
	my ($id,$Fupstream,$Fdownstream,$Rupstream,$Rdownstream)=(split/\t/,$_);
	
	my $string=join"\t",($Fupstream,$Fdownstream,$Rupstream,$Rdownstream);
	
	$qual{$id}=$string;
    
}
close QU;




open AN,"$aligned" or die $!;

while (<AN>){
	chomp;
	my @array=split/\t/,$_;
	
	### find the read information
	my $read=$array[10];
	
	### calculate the read qulaity
	my $quality=$qual{$read};
	my $id=$array[3];
	
	### concatenate read and quality together
	my $str=join "\t",($read,$quality);
	
	
	$inf{$id}=$str;
}

close AN;

open UA,"$unaligned" or die $!;
while (<UA>){
	chomp;
	my ($id1,$reads,$seq)=split/\t/,$_;
	my $id2=(split/\|/,$id1)[0];
	
	my $quality_un=$qual{$reads};
	
	my $str2=join "\t",($reads,$quality_un);
	
	$inf{$id2}=$str2;
	
	
}

close UA;

open IN, "$input" or die "cannot open file $input";
open OUT,">$output" or die $!;
while (<IN>){
	chomp;
	my @array=(split/\t/,$_);
	my $insertion=$_;
	
	my ($id,$detail,$sample,$cov)=split/\|/,$array[1];
	
	if($detail eq "Unmappable"){
		
		my $information=$inf{$id};
		print OUT "$insertion\t$information\n";
	}else{
		my @array2=split/\:/,$detail;
		
		my $information=$inf{$array2[3]};
		my $string=join "\t",@array;
		
		print OUT "$insertion\t$information\n"
	
	}
    
}
close IN;
close OUT;
