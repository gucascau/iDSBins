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
GetOptions(\%opts,"o:s","f:s","g:s","m:s","n:s","t:s","h:s");
print "*************\n*$version*\n*************\n";
if ( !defined $opts{o}|| !defined $opts{m} || !defined $opts{n}  ||!defined $opts{f}||!defined $opts{g} ||!defined $opts{t} ||defined $opts{h}) {
       	die "************************************************************************
       	Usage: $0.pl
			-g: Yeast genome sequence
			-m: forward reads
			-n: reverse reads
			-f: Estimated two donor Blast results
			-t: Sample name
			-o: Output of the estimated files
			
			-h: help
************************************************************************\n";
}



my $estimate=$opts{f};
### the sequence index
my $fasta=$opts{g};

my $fastq1=$opts{m};
my $fastq2=$opts{n};

my $sample=$opts{t};

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


my $out=$opts{o};


open EST,"$estimate" or die "cannot open file $estimate";


open OUT,">$out.txt" or die $!;
open ESTW,">$out.wholeseq.fasta" or die $!;
open ESTINS,">$out.inserted.fasta" or die $!;
open BED,">$out.bed" or die $!;
my $n=0;

while (<EST>) {
	chomp;
	s/\r//g;
    # print "$_" ;
	my ($ocov,$id,$chr,$iden,$Inslength,$Rstart,$RM1,$RM2,$Rend,$fstart,$rstart,$fend,$rend,$rlength,$rlengthr,$rQual,$strand,$cov)=split/\s+/,$_;
	#my $string=$str{$array[0]};
	$n++;
	
	
	### For estimated one we set up strigent criterion that there should be no more than 10bp gaps for 
	
	my $startF= ($fstart<=50)?$fstart:45;
	my $endF=(($rlengthr-$rend)<=60)?$rend:($rlengthr-50);
	
	my $length=$Inslength;
	
	print "$length\n";
	### Calculate the donor number of insertions

	my $ndonor=1;
	if ($fstart >60 || ($rlengthr-$rend)>60){
		$ndonor="1orMore";

	}
	

	### Upstream sequence:
	
	my $string1=$strf{$id};
	
	my $InsUplength=$rlength-$fstart+1;
	my $InsUpstring=substr $strf{$id},($fstart-1),$InsUplength;
	
	
	### Middle estimated sequence
	my $lengthM=($RM2-$RM1+1);
	my $string2=substr $str{$chr},$RM1,$lengthM;
	
	
	## we need to consider the strand
	
	my $fMidString;
	if ($strand eq "+"){
		$fMidString=$string2;
	}else{
		$string2=~tr/atcgATCG/tagcTAGC/;
		$fMidString=reverse $string2;
	}
	
	print "$fMidString\n";
	### Downstream sequence:
	
	my $string3=$strr{$id};
	my $InsDownstring=substr $strf{$id}, 0,$rend;
	
	#print "$string3\n";
	### Add the junction sequence
	
	# ### The inserted junction information (30bp upstream and downstream)
# 	my $start_3=$startF-15;
# 	my $start_5=$endF-15;
# 	my $upjunction=substr $string{$id}, $start_3,30;
# 	my $downjunction=substr $string{$id}, $start_5,30;
#
#
	
	
	### The final insertion sequence and fragment with inserted DNA
	my $ffragment=join "",($string1,$fMidString,$string3);
	my $finsertion=join "",($InsUpstring,$fMidString,$InsDownstring);
	
	
	### add the junction sequence
	my $start_3=$startF-15;
	my $start_5=$endF-15;
	my $upjunction=substr $strf{$id}, $start_3,30;
	my $downjunction=substr $strr{$id}, $start_5,30;
	
	
	print ESTW ">E$n\t$id\n$ffragment\n";
	print ESTINS ">E$n\t$id\n$finsertion\n";
	print BED "$chr\t$Rstart\t$Rend\t$id\t0\t$strand\n";
	
	#print OUT "$n\t$id\t$sample\t$string{$id}\t$insertion\t$length\t$ndonor\t$iden\t$fcov\t$fqual\t$chr\t$dstart\t$dend\t$strand\t$annotation\t$dist\t$feature\t$FMstring\t$upjunction\t$downjunction\n";
	print OUT "E$n\t$id\t$sample\t$ffragment\t$finsertion\t$length\t$ndonor\t$iden\t$cov\t$rQual\t$chr\t$Rstart\t$Rend\t$strand\tNO\tNO\tNO\tNO\t$upjunction\t$downjunction\n";

	
}
close EST;
close OUT;
close ESTINS;
close ESTW;
close BED;


