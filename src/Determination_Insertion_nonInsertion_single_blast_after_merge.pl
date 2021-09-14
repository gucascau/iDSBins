#!/usr/bin/perl
use strict;
use warnings;


### This script is to extract the reads that consistently showed non insertion events:
### requirement: 1
#####1.	Completely blasted against reference genome from start to end of the reads (allow for 2 gaps and 2 mismatches at most, Identify is more than 95%);
#####2.	Forward and Reward regards consistently mapped against the identical reference location.
#### If only they have one read mapping information

#use Bio::Seq;
#use Bio::SeqIO;


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{o} ) {
	die "************************************************************************
	Usage: extract_fastq.pl -i R1_Blast -o Blast_categories
************************************************************************\n";
}

my $output=$opts{o};
my $input=$opts{i};
my $index=$opts{g};
my $yeast=$opts{m};

my %conf1; my %inf1; my %int1; my %confInsert1; my %potential1_two; my %potential1_chr; my %check_all1; my %basic1;

open I,"$input" or die "cannot open file $input";



open OUT,">$output.noninsert.confident.txt" or die $!;

open CI,">$output.insert.confident.txt" or die $!;

open PI,">$output.potential.insertion.txt" or die $!;
open GID,">$output.Pdeletion.txt" or die $!;
#open GII, ">$output.Pinsertion.txt" or die $!;
open GIU, ">$output.Pinsertion_mat.txt" or die $!;
open GIUT, ">$output.unclassified.txt" or die $!;




while (<I>) {
	chomp;
    # print "$_" ;
	my ($qid,$pid,$identity,$match,$gap,$mismatch,$qstart,$qend,$qlength,$pstart,$pend,$plength,$score,$evlaue)=split/\s+/,$_;
	my $information=$_;
	$basic1{$qid}=$information;
	next if ($pid eq "ref-NC_001135-" && $pstart >= 13650 && $pend<= 13850 );
	next if ($pid eq "ref-NC_001135-" && $pstart >= 200750 && $pend <= 201000);
	next if (exists $conf1{$qid});
	next if (exists $confInsert1{$qid});
	next if (exists $potential1_chr{$qid});
	
	
	
	my $cov=$match/$qlength;
	print "$cov\n";
	

	#### confident from other chromosome insertion
	if ($pid ne "ref-NC_001135-"){
		$confInsert1{$qid}=$information;
		
		print CI "$information\n";
		
		
		###### confident non-insertion:
	}elsif($pstart>=294300  && $pstart <=294500 && $pend>=294300 && $pend<=294500 && $match>=70){
		$conf1{$qid}=$information;
		
		print OUT "$information\n";
		##### confident insertion from other location of chr3
	}elsif($pstart>=294500 || $pstart<=13600 || ($pstart>=13800 && $pstart<=294300)){
			$potential1_chr{$qid}=$information;
			#### potential insertion or deletion but need to check original data  
			print CI "$information\n";
			
		
		#### confident deletions due to short fragements (<=80 && matches <=60)
			
	}elsif($pstart>=294300  && $pstart <=294500 && $pend>=294300 && $pend<=294500 && $match<=60 && $qlength<80){
		print GID "$information\n";
		
		
		
		#### potential split the fragments with insertion (telemere repetitive elements) 
		#### undertermine the insertion or deletion due to the last one 
		#### potential split the fragments with insertion from MAT elements
	}elsif($pstart>=294300  && $pstart <=294500 && $pend>=294300 && $pend<=294500 && $match<=60 && $qlength >90){
		
		####################################	
		if(exists ($potential1_two{$qid})){
			#my $length=abs($end -$inf1{$qid}->{start})
			my @array=sort { $a <=> $b }($qstart,$qend,$inf1{$qid}->{start},$inf1{$qid}->{end});
			if ($array[0]<=5 && $qlength-$array[3] <5){
				print PI "$potential1_two{$qid}\n$information\n";
			}else{
			
				#### Potential insertion that came from MAT elements

				print GIU "$potential1_two{$qid}\n$information\n"
			}

		}else{
			$potential1_two{$qid}=$information;
				$inf1{$qid}->{start}=$qstart;
				$inf1{$qid}->{end}=$qend;
				$inf1{$qid}->{match}=$match;
				$inf1{$qid}->{length}=$qlength;
				############
				#potential insertion
		}

		
	######## determine the insert or deletion events:	
	}else{
		print GIUT "$information\n";
	}
	

	
}

close GID;
close GIU;
close GIUT;
close CI;
close OUT;
close CI;





