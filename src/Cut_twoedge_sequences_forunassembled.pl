#!/usr/bin/perl
use strict;
use warnings;



#use Bio::Seq;
#use Bio::SeqIO;


#### extract two edge (100bp) of insertions from two PE reads

my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"f:s","o:s","u:s","t:s","i:s","r:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{f} || !defined $opts{o} ||!defined $opts{u} ||!defined $opts{t} || !defined $opts{r}||!defined $opts{i}) {
	die "************************************************************************
	Usage: extract_fastq.pl -f forward fasta -r reverse fasta -i Insertion id -u cut-off of upstream (30bp) -t start of upstreamdownstream (1bp) -o Output of concatenarate sequences
************************************************************************\n";
}

my $output=$opts{o};
my $ffasta=$opts{f};
my $rfasta=$opts{r};
my $insert=$opts{i};

my $cutoff=$opts{u};
my $start=$opts{t};

open FFASTA,"$ffasta" or die $!;

my $id;  my %Fhash; 
while (<FFASTA>) {
	chomp;
    # print "$_" ;

    if (/>(\S+)/)  {
	
		$id=$1;
	}else{
		$Fhash{$id}.=$_;
	
	}
}
close FFASTA;


open RFASTA,"$rfasta" or die $!;

my $id2;  my %Rhash; 
while (<RFASTA>) {
	chomp;
    # print "$_" ;

    if (/>(\S+)/)  {
	
		$id2=$1;
	}else{
		$Rhash{$id2}.=$_;
	
	}
}
close RFASTA;

open IN,"$insert" or die $!;
open OUT,">$output" or die $!;
my %num;
while (<IN>) {
	chomp;
    # print "
	next if ($_=~/^$/);
	my $i =(split/\t/,$_)[0];
	#print "$i\n";
	#print "$Fhash{$i}\n$Rhash{$i}\n\n";
	next if (exists $num{$i});
	$num{$i}++;

	my $Fqual_str=($start==1)?$Fhash{$i}:(substr $Fhash{$i},($start-1),(-$start+1));
	#my $qual_str=$hash{$i};
	my $Ftotal_length=length $Fqual_str;
	#my $Rtcutoff=$cutoff;
	my $str1= ($Ftotal_length >$cutoff)?(substr $Fqual_str,33,$cutoff):$Fqual_str;
	
	
	my $Rqual_str=($start==1)?$Rhash{$i}:(substr $Rhash{$i},($start-1),(-$start+1));
	#my $qual_str=$hash{$i};
	my $Rtotal_length=length $Rqual_str;
	
	#my $Rtcutoff=2*$cutoff;
	my $str2= ($Ftotal_length >$cutoff)?(substr $Rqual_str,-69,$cutoff):$Rqual_str;
	
	my $str_long=join "",($str1,$str2);

	#my $str= ($total_length >$tcutoff)?$str_long:$qual_str;

	print OUT ">$i\n$str_long\n";

}

close OUT;
close IN;





