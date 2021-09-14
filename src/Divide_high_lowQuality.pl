#!/usr/bin/perl
use strict;
use warnings;



#use Bio::Seq;
#use Bio::SeqIO;


#### perform the filteration of upstream of 100bp using quality control (>30) and the left downstream using quality control (>15)


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"f:s","r:s","o:s","u:s","g:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{f} || !defined $opts{r} ||!defined $opts{o} ||!defined $opts{u} ||!defined $opts{g}) {
	die "************************************************************************
	Usage: extract_fastq.pl -f Ffastq -r  Rfastq -u cut-off of upstream (100bp) -g cut-off of downstream quality(leftbp) -o Output of Quality Control Reads
************************************************************************\n";
}

my $output=$opts{o};
my $fastq1=$opts{f};
my $fastq2=$opts{r};
my $cutoff=$opts{u};
my $dcutoff=$opts{g};

open FASTQ1,"$fastq1" or die $!;

my $id1;  my %hash1;  my %qual1;
while (<FASTQ1>) {
	chomp;
    # print "$_" ;

    if ($. % 4 == 1)  {
		$id1=(split/\s+/,$_)[0];
		$id1=~s/^@//;
		#print "$id1\n";
		$hash1{$id1}=$_."\n";
	}elsif( $.%4==0) {
		my $score;
		my $qual_str=$_;
		my $total_length=length $qual_str;
		
		#### extract the upstream 100bp of sequences (reads longer than 100bp)
		my $str= ($total_length >100)?substr($qual_str,0,100):$qual_str;
		#print "$id1\n$str\n";
		my $length= ($total_length >100)?100:($total_length);
		map{ $score += ord($_)-33} split("",$str);
		#print "$id1\t$score\n";
		my $av_qual=$score/$length;
		$qual1{$id1}->{upstream}=$av_qual;
		
		#### extract the left downstream sequences (reads shorter than 100bp retain the quality of upstream)
		my $dscore;
		
		my $dlength= ($total_length >100)? ($total_length -100):($total_length);
		my $dstr= ($total_length >100)?substr($qual_str,100,$dlength):$qual_str;
		
		 map{ $dscore += ord($_)-33} split("",$dstr);
		 my $dav_qual=$dscore/$dlength;
		$qual1{$id1}->{downstream}=$dav_qual;
		
		
		
		$hash1{$id1}.=$_."\n";
	}elsif($.%4 == 2){
		#$string1{$id1}=$_;
		$hash1{$id1}.=$_."\n";
	}else{
		$hash1{$id1}.=$_."\n";
	}
}
close FASTQ1;

my $id2; my %hash2; my %qual2;
open FASTQ2, "$fastq2" or die $!;

while (<FASTQ2>) {
	chomp;
    # print "$_" ;

    if ($. % 4 == 1)  {
		$id2=(split/\s+/,$_)[0];
		$id2=~s/^@//;
		$hash2{$id2}.=$_."\n";
	}elsif( $.%4==0) {
		my $score;
		my $qual_str=$_;
		my $total_length=length $qual_str;
		
		#### extract the upstream 100bp of sequences (reads longer than 100bp)
		my $str= ($total_length >100)?substr($qual_str,0,100):$qual_str;
		#print "$id1\n$str\n";
		my $length= ($total_length >100)?100:($total_length);
		 map{ $score += ord($_)-33} split("",$str);
		#print "$id1\t$score\n";
		my $av_qual=$score/$length;
		$qual2{$id2}->{upstream}=$av_qual;
		
		#### extract the left downstream sequences (reads shorter than 100bp retain the quality of upstream)
		my $dscore;
		
		my $dlength= ($total_length >100)? ($total_length -100):($total_length);
		my $dstr= ($total_length >100)?substr($qual_str,100,$dlength):$qual_str;
		
		 map{ $dscore += ord($_)-33} split("",$dstr);
		 my $dav_qual=$dscore/$dlength;
		$qual2{$id2}->{downstream}=$dav_qual;
		
		
		
		$hash2{$id2}.=$_."\n";
	}elsif($.%4 == 2){
		#$string2{$id2}=$_;
		$hash2{$id2}.=$_."\n";
	}else{
		$hash2{$id2}.=$_."\n";
	}
}
close FASTQ2;




open OUT1,">$output.highquality.R1.fastq" or die $!;
open OUT2,">$output.highquality.R2.fastq" or die $!;
open OUTH,">$output.highquality.stat" or die $!;

open OUT3,">$output.lowquality.R1.fastq" or die $!;
open OUT4,">$output.lowquality.R2.fastq" or die $!;
open OUTL,">$output.lowquality.stat" or die $!;

open OUTA, ">$output.Allquality.stat" or die $!;

my $hnum=0;
my $lnum=0;

print OUTH "ID\tFQUpstream\tFQDownstream\tRQUpstream\tRQDownstream\n";
print OUTL "ID\tFQUpstream\tFQDownstream\tRQUpstream\tRQDownstream\n";

foreach my $i (keys %qual1){
	if ($qual1{$i}->{upstream} >=$cutoff && $qual2{$i}->{upstream} >=$cutoff && $qual1{$i}->{downstream} >$dcutoff && $qual2{$i}->{downstream} >$dcutoff ){
		print OUT1 "$hash1{$i}";
		print OUT2 "$hash2{$i}";
		$hnum++;
		print OUTH "$i\t$qual1{$i}->{upstream}\t$qual1{$i}->{downstream}\t$qual2{$i}->{upstream}\t$qual2{$i}->{downstream}\n";
		print OUTA "$i\t$qual1{$i}->{upstream}\t$qual1{$i}->{downstream}\t$qual2{$i}->{upstream}\t$qual2{$i}->{downstream}\n";
	}else{
		print OUT3 "$hash1{$i}";
		print OUT4 "$hash2{$i}";
		print OUTL "$i\t$qual1{$i}->{upstream}\t$qual1{$i}->{downstream}\t$qual2{$i}->{upstream}\t$qual2{$i}->{downstream}\n";
		print OUTA "$i\t$qual1{$i}->{upstream}\t$qual1{$i}->{downstream}\t$qual2{$i}->{upstream}\t$qual2{$i}->{downstream}\n";
		$lnum++;
	}
	
}

print OUTH "High quality number: $hnum\n";

print OUTL "Low quality number: $lnum\n";

close OUT1;
close OUT2;
close OUT3;
close OUT4;

close OUTH;
close OUTL;

