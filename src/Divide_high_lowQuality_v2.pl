#!/usr/bin/perl
use strict;
use warnings;

#use Bio::Seq;
#use Bio::SeqIO;
#Author:Xin Wang 
#email: xin.wang@childrens.harvard.edu
#PI: Kaifu Chen


### This script is used to filter the low quality reads and measure the read quality of inserted regions using both read quality.


#### Criterion:
### 1. The quality of two edge of MATA region is higher than 25, you could set up according if possible.
### 2. The quality of inserted region\/or left region is higher than 15, due to the low aboundance or quality of Miseq with increased read length.
### 3. Output for the seperation of read into high quality and low quality. Meanwhile, we also measure the number and put the read quality informaiton into different files.


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"f:s","r:s","o:s","u:s","g:s","h:s","t:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{f} || !defined $opts{r} ||!defined $opts{o} || defined $opts{h}) {
	die "************************************************************************
	Usage: extract_fastq.pl -f Ffastq -r Rfastq -o Output of Quality Control Reads
	
	Request Parameters:
	-f Ffastq (Forward reads)
	-r Rfastq (Reverse reads)
	-o Output substring of quality control Reads
	
	Optional Parameters:
	-t Set up the while size of MAT region(this region was set up for the rough size of MATA region, default (90bp) )
	-u Cut-off of upstream quality(this regions mainly MAT regions, therefore we set up with high quality, default (25)) 
	-g Cut-off of downstream quality(this regions mainly inserted regions, therefore we set up with a low quality, default(15)) 
	-h Help
************************************************************************\n";
}

my $output=$opts{o};
my $fastq1=$opts{f};
my $fastq2=$opts{r};

### set up the parameters:

# set up the upstream quality (two sides of MAT region, half of whole MAT)
my $cutoff=(defined $opts{u})?$opts{u}:25;
# set up the downstream quality(the left part (Inserted region or MAT region), this might contain other side of MAT region if the reads are very short )
my $dcutoff=(defined $opts{g})?$opts{g}:15;
# set up the whole size of MAT region.
my $SMAT=(defined $opts{t})?$opts{t}:90;

open FASTQ1,"$fastq1" or die $!;

my $id1;  my %hash1;  my %qual1; my %Flength; my %seq1;
while (<FASTQ1>) {
	chomp;
    # print "$_" ;
    if ($. % 4 == 1)  {
		
		### put the ID into hash
		$id1=(split/\s+/,$_)[0];
		$id1=~s/^@//;
		#print "$id1\n";
		$hash1{$id1}=$_."\n";
		
	}elsif( $.%4==0) {
		my $score;
		my $qual_str=$_;
		my $total_length=length $qual_str;
		
		$Flength{$id1}=$total_length;
		### set up the upstream size of whole MAT region
		my $halfMAT=$SMAT/2;
		#print "$halfMAT\n";
		
		### Here some reads will be eliminated if they are shorter than half of MAT (45bp);
		next if ($total_length<=$halfMAT);
	
		
		my $str=substr($qual_str,0,$halfMAT);
		####
		# my $str= ($total_length >100)?substr($qual_str,0,100):$qual_str;
	# 	#print "$id1\n$str\n";
	# 	my $length= ($total_length >100)?100:($total_length);
		map{ $score += ord($_)-33} split("",$str);
		#print "$id1\t$score\n";
		my $av_qual=$score/$halfMAT;
		$qual1{$id1}->{upstream}=$av_qual;
		
		#### extract the left downstream sequences (This region might contain the other side of MAT region or inserted region based on the legnth of whole read)
		my $dscore;
		
		my $dlength=$total_length - $halfMAT;
		my $dstr= substr($qual_str, -$dlength);
		
		# #my $dlength= ($total_length >100)? ($total_length -100):($total_length);
# 		my $dstr= ($total_length >100)?substr($qual_str,100,$dlength):$qual_str;
#	
		 map{ $dscore += ord($_)-33} split("",$dstr);
		 my $dav_qual=$dscore/$dlength;
		$qual1{$id1}->{downstream}=$dav_qual;

		### Put all the information into hash
		$hash1{$id1}.=$_."\n";

	}elsif($.%4 == 2){
		#$string1{$id1}=$_;
		$hash1{$id1}.=$_."\n";
		$seq1{$id1}=$_;
	}else{
		$hash1{$id1}.=$_."\n";
	}
}
close FASTQ1;

my $id2; my %hash2; my %qual2; my %Rlength; my %seq2;
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

		$Rlength{$id2}=$total_length;
		### set up the upstream size of whole MAT region
		my $halfMAT=$SMAT/2;
		
	
		
		### Here some reads will be eliminated if they are shorter than half of MAT (45bp);
		next if ($total_length<=$halfMAT);
		
		
		my $str=substr($qual_str,0,$halfMAT);
		####
		# my $str= ($total_length >100)?substr($qual_str,0,100):$qual_str;
	# 	#print "$id1\n$str\n";
	# 	my $length= ($total_length >100)?100:($total_length);
		map{ $score += ord($_)-33} split("",$str);
		#print "$id1\t$score\n";
		my $av_qual=$score/$halfMAT;
		$qual2{$id2}->{upstream}=$av_qual;
		
		#### extract the left downstream sequences (This region might contain the other side of MAT region or inserted region based on the legnth of whole read)
		my $dscore;
		
		my $dlength=$total_length - $halfMAT;
		my $dstr= substr($qual_str, -$dlength);
		
		# #my $dlength= ($total_length >100)? ($total_length -100):($total_length);
# 		my $dstr= ($total_length >100)?substr($qual_str,100,$dlength):$qual_str;
#	
		 map{ $dscore += ord($_)-33} split("",$dstr);
		 my $dav_qual=$dscore/$dlength;
		$qual2{$id2}->{downstream}=$dav_qual;

		
		$hash2{$id2}.=$_."\n";
	}elsif($.%4 == 2){
		#$string2{$id2}=$_;
		$seq2{$id2}=$_;
		$hash2{$id2}.=$_."\n";
	}else{
		$hash2{$id2}.=$_."\n";
	}
}
close FASTQ2;




open OUT1,">$output.highquality.R1.fastq" or die $!;
open OUT2,">$output.highquality.R2.fastq" or die $!;
open HQFFA,">$output.highquality.R1.fasta" or die $!;
open HQRFA,">$output.highquality.R2.fasta" or die $!;


open OUTH,">$output.highquality.stat" or die $!;

open OUT3,">$output.lowquality.R1.fastq" or die $!;
open OUT4,">$output.lowquality.R2.fastq" or die $!;

open LQFFA,">$output.lowquality.R1.fasta" or die $!;
open LQRFA,">$output.lowquality.R2.fasta" or die $!;

open OUTL,">$output.lowquality.stat" or die $!;

open OUTA, ">$output.Allquality.stat" or die $!;

my $hnum=0;
my $lnum=0;

print OUTH "ID\tFQUpstream\tFQDownstream\tRQUpstream\tRQDownstream\tAverageInserted\n";
print OUTL "ID\tFQUpstream\tFQDownstream\tRQUpstream\tRQDownstream\tAverageInserted\n";
print OUTA "ID\tFQUpstream\tFQDownstream\tRQUpstream\tRQDownstream\tAverageInserted\n";

foreach my $i (keys %qual1){
	
	## mainly eliminated the shorter reads (shorter than half MATA)
	next if (!exists $qual2{$i}->{upstream});
	
	
	### Here we require the high-quality reads have a upstream higher than 25 and a downstream higher than 15.
	if ($qual1{$i}->{upstream} >=$cutoff && $qual2{$i}->{upstream} >=$cutoff && $qual1{$i}->{downstream} >$dcutoff && $qual2{$i}->{downstream} >$dcutoff ){
		print OUT1 "$hash1{$i}";
		print OUT2 "$hash2{$i}";
		
		print HQFFA ">$i\n$seq1{$i}\n";
		print HQRFA ">$i\n$seq2{$i}\n";
		
		$hnum++;
		my $averageDown=($qual1{$i}->{downstream}+$qual2{$i}->{downstream})/2;
		print OUTH "$i\t$qual1{$i}->{upstream}\t$qual1{$i}->{downstream}\t$qual2{$i}->{upstream}\t$qual2{$i}->{downstream}\t$averageDown\n";
		print OUTA "$i\t$qual1{$i}->{upstream}\t$qual1{$i}->{downstream}\t$qual2{$i}->{upstream}\t$qual2{$i}->{downstream}\t$averageDown\n";
	}else{
		print OUT3 "$hash1{$i}";
		print OUT4 "$hash2{$i}";
		
		print LQFFA ">$i\n$seq1{$i}\n";
		print LQRFA ">$i\n$seq2{$i}\n";
		
		my $averageDown2=($qual1{$i}->{downstream}+$qual2{$i}->{downstream})/2;
		print OUTL "$i\t$qual1{$i}->{upstream}\t$qual1{$i}->{downstream}\t$qual2{$i}->{upstream}\t$qual2{$i}->{downstream}\t$averageDown2\n";
		print OUTA "$i\t$qual1{$i}->{upstream}\t$qual1{$i}->{downstream}\t$qual2{$i}->{upstream}\t$qual2{$i}->{downstream}\t$averageDown2\n";
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
close OUTA;

