#!/usr/bin/perl
use strict;
use warnings;
#use Data::Dump qw(dump);
#Author:Xin Wang 
#email: xin.wang@childrens.harvard.edu
#PI: Kaifu Chen

#use Bio::Seq;
#use Bio::SeqIO;


#### extract two edge (30+30bp) of insertions from two PE reads for further deduplication

my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"f:s","o:s","u:s","t:s","e:s","r:s","h:s","i:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{f} || !defined $opts{o}  || !defined $opts{r} || !defined $opts{i} || defined $opts{h}) {
	die "************************************************************************
	Usage: extract_fastq.pl -f forward fastq -r reverse fastq -u cut-off of upstream (33bp) -t start of upstream -e the start site of reverse reads -o Output of concatenarate sequences
	
	Request Parameters:
	-f Forward fastq file
	-r Reverse fastq file
	-i Inserted reads with quality file
	-o Output substring of quality control Reads
	
	Optional Parameters:
	-u the cut-off of upstream and downstream (default 30bp, 11+19bp)
	-t The start site of Forward reads (default 33)
	-e The start site of Reverse reads (default 39)
	-h Help
************************************************************************\n";
}

my $output=$opts{o};
my $ffastq=$opts{f};
my $rfastq=$opts{r};
my $insert=$opts{i};

my $cutoff=(defined $opts{u})?$opts{u}:30;
my $fstart=(defined $opts{t})?$opts{t}:33;
my $rstart=(defined $opts{e})?$opts{e}:39;


open FFASTQ,"$ffastq" or die $!;

my $id1;  my %Fhash;  
while (<FFASTQ>) {
	chomp;
    # print "$_" ;

    if ($. % 4 == 1)  {
		$id1=(split/\s+/,$_)[0];
		$id1=~s/^@//;
		#print "$id1\n";
		
	}elsif($. % 4 == 2){
		$Fhash{$id1}=$_;
	}
	
}
close FFASTQ;


open RFASTQ,"$rfastq" or die $!;

my $id2;  my %Rhash;  
while (<RFASTQ>) {
	chomp;
    # print "$_" ;

    if ($. % 4 == 1)  {
		$id2=(split/\s+/,$_)[0];
		$id2=~s/^@//;
		#print "$id1\n";
		
	}elsif($. % 4 == 2){
		$Rhash{$id2}=$_;
	}
	
}
close RFASTQ;



open IN,"$insert" or die $!;
open OUT,">$output" or die $!;
my %num;
while (<IN>) {
	chomp;
    # print "
	next if ($_=~/^$/);
	my $i =(split/\t/,$_)[0];
	next unless (exists $Fhash{$i});
	#print "$i\n";
	#print "$Fhash{$i}\n$Rhash{$i}\n\n";
	next if (exists $num{$i});
	$num{$i}++;


	my $Fqual_str=$Fhash{$i};
	
	my $Ftotal_length=length $Fqual_str;
	
	#my $Rtcutoff=$cutoff;
	my $str1= substr ($Fqual_str,$fstart,$cutoff);
	
	#my $str1= ($Ftotal_length >$cutoff)?(substr $Fqual_str,33,$cutoff):$Fqual_str;
	
	
	# my $Rqual_str=($start==1)?$Rhash{$i}:(substr $Rhash{$i},($start-1),(-$start+1));
	#my $qual_str=$hash{$i};
	
	my $Rqual_str=$Rhash{$i};
	my $Rtotal_length=length $Rqual_str;

	my $str2= substr ($Rqual_str,$rstart,$cutoff);
	
	### Reverse Complement Of Dna for R2
	$str2 =~ tr/ATGCatgc/TACGtacg/;
	my $revcomp = reverse $str2;
	
	my $str_long=join "",($str1,$revcomp);

	#my $str= ($total_length >$tcutoff)?$str_long:$qual_str;

	print OUT ">$i\n$str_long\n";

}

close OUT;
close IN;





