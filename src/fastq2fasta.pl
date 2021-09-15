#!/usr/bin/perl
use strict;
use warnings;

### 
### Fastq to fasta file

#use Bio::Seq;
#use Bio::SeqIO;


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i}  ||!defined $opts{o}) {
	die "************************************************************************
	Usage: extract_fastq.pl -i Rfastq  -o Rfasta
************************************************************************\n";
}

my $output=$opts{o};
my $fastq=$opts{i};

# my $index1=$opts{r};
# my $index2=$opts{f};

open FASTQ1,"$fastq" or die $!;
open OUT,">$output" or die $!;
my $id1; 
while (<FASTQ1>) {
	chomp;
    # print "$_" ;

    if ($. % 4 == 1)  {
		$id1=(split/\s+/,$_)[0];
		$id1=~s/^@//;
		print OUT ">$id1\n";
		#$hash1{$id1}=$_."\n";
	}elsif($.%4 == 2){
		print OUT "$_\n";
		#$hash1{$id1}.=$_."\n";

	}

}
close OUT;
close FASTQ1;
