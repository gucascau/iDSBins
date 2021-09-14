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
if (!defined $opts{i} || !defined $opts{g} ||!defined $opts{o} ) {
	die "************************************************************************
	Usage: extract_fastq.pl -i R1_Blast -g R2_Blast  -o Blast_categories
************************************************************************\n";
}

my $output=$opts{o};
my $input=$opts{i};
my $index=$opts{g};
#my $yeast=$opts{m};

my %conf1; my %inf1; my %int1; my %confInsert1; my %potential1_two; my %potential1_chr; my %check_all1; my %basic1;

open I,"$input" or die "cannot open file $input";
while (<I>) {
	chomp;
    # print "$_" ;
	my ($qid,$pid,$identity,$match,$gap,$mismatch,$qstart,$qend,$qlength,$pstart,$pend,$plength,$score,$evlaue)=split/\s+/,$_;
	my $information=$_;
	$basic1{$qid}=$information;
	next if (exists $conf1{$qid});
	next if (exists $confInsert1{$qid});
	
	my $cov=$match/$qlength;
	print "$cov\n";
	

	#### potential other chromosome insertion
	if ($pid ne "ref-NC_001135-"){
		$confInsert1{$qid}=$information;
		###### confident non-insertion:
	}elsif($pstart>=294300  && $pstart <=294500 && $pend>=294300 && $pend<=294500 && $match>=70){
		$conf1{$qid}=$information;
		
		##### potential insertion from other location of chr3
	}elsif($pstart>=294500 || $pstart<=13600 || ($pstart>=13800 && $pstart<=294300)){
		$potential1_chr{$qid}=$information;
		
		
		#### potential split the fragment with insertion (telemere repetitive elements)
	}elsif($pstart>=294300  && $pstart <=294500 && $pend>=294300 && $pend<=294500 && $match<=50){
		$potential1_two{$qid}=$information;
		####################################	
			$inf1{$qid}->{start}=$qstart;
			$inf1{$qid}->{match}=$match;
			$inf1{$qid}->{length}=$qlength;
			#############
		

		#### potential insertion or deletion but need to check original data  
		
	}# elsif($pstart<13800 && $pstart>13600){
# 		$check_all1{$qid}=$information;
#
# 		$int1{$qid}->{start}=$qstart;
# 		$int1{$qid}->{match}=$match;
# 		$int1{$qid}->{length}=$qlength;
#
	# }
}
close I;


my %conf2; my %inf2;my %int2; my %confInsert2;my %potential2_two; my %potential2_chr; my %check_all2; my %basic2;
open G,"$index" or die "cannot open file $index";
while (<G>) {
	chomp;
    # print "$_" ;
	my ($qid,$pid,$identity,$match,$gap,$mismatch,$qstart,$qend,$qlength,$pstart,$pend,$plength,$score,$evlaue)=split/\s+/,$_;
	my $information=$_;
	$basic2{$qid}=$information;
	next if (exists $conf2{$qid});
	next if (exists $confInsert2{$qid});
	my $cov=$match/$qlength;
	

		#### potential other chromosome insertion
	if ($pid ne "ref-NC_001135-"){
		$confInsert2{$qid}=$information;
		###### confident non-insertion:
	}elsif($pstart>=294300  && $pstart <=294500 && $pend>=294300 && $pend<=294500 && $match>=70){
		$conf2{$qid}=$information;
		##### confident insertion from other location of chr3
	}elsif($pstart>=294500 || $pstart<=13600 || ($pstart>=13800 && $pstart<=294300)){
		$potential2_chr{$qid}=$information;
		#### potential insertion or deletion but need to check original data  
		#### potential split the fragment with insertion (telemere repetitive elements)
	}elsif($pstart>=294300  && $pstart <=294500 && $pend>=294300 && $pend<=294500 && $match<=50){
		$potential2_two{$qid}=$information;
	####################################	
		$inf2{$qid}->{start}=$qstart;
		$inf2{$qid}->{match}=$match;
		$inf2{$qid}->{length}=$qlength;
		#############
	# }elsif($pstart<13800 && $pstart>13600){
# 		$check_all2{$qid}=$information;
#
# 		$int2{$qid}->{start}=$qstart;
# 		$int2{$qid}->{match}=$match;
# 		$int2{$qid}->{length}=$qlength;
	}

}
close G;







open OUT,">$output.noninsert.confident.txt" or die $!;

open R1,">$output.noninsert.R1.txt" or die $!;
open R2,">$output.noninsert.R2.txt" or die $!;

### including the insertion from other chromosome and chromosome:
open CI,">$output.insert.two_confident.txt" or die $!;
#### one with unknown insertion another with insertion from other location.
open OI,">$output.insert.potential_otherchr.txt" or die $!;
open OIC,">$output.insert.potential_chr3.txt" or die $!;

open GID,">$output.Pdeletion.txt" or die $!;

#### Potential insertion only showed two edge of insertion events 4-44 for two PE, and reads length is longer than 90bp
open GII, ">$output.Pinsertion.txt" or die $!;
#### Potential insertion only showed two edge of insertion events 4-44 for two PE, and reads length is longer than 90bp
open GIU, ">$output.Pdeletion_uncertain.txt" or die $!;
open GIUT, ">$output.Pinsertion_telemore.txt" or die $!;


#open PI,">$output.Pinsertordeletiondoublecheck.txt" or die $!;



###### if we removed the bias location 
####next if ($pid eq "ref-NC_001135-" && $pstart<13800 && $pstart>13600);####next if ( $pid eq "ref-NC_001135-" && $pstart<200950 && $pstart>200800);

#open PII,">$output.Pinsertiondoublecheck.txt" or die $!;
#open PID,">$output.Pdeletiondoublecheck.txt" or die $!;
#open PIU,">$output.Pdeletion_uncertaindoublecheck.txt" or die $!;

open UN,">$output.unclassified.txt" or die $!;


#open LIST,"$yeast" or die $!;

foreach my $i (keys %basic1){

	
	##### print confident no-insertion
	if (exists $conf2{$i} && exists $conf1{$i}){
		print OUT "$conf1{$i}\n$conf2{$i}\n";
	}elsif(exists $conf1{$i} && !exists $conf2{$i}){
		print R1 "$conf1{$i}\n$basic2{$i}\n";
	}elsif(exists $conf2{$i} && !exists $conf1{$i}){
		print R2 "$basic1{$i}\n$conf2{$i}\n";
		
		#### print confident insertion (insertion from other chromomsome)
	}elsif(exists $confInsert1{$i} && exists $confInsert2{$i}){
		print CI "$confInsert1{$i}\n$confInsert2{$i}\n";
		#### print insertion from other chromosome
	}elsif(exists $confInsert1{$i} && !exists $confInsert2{$i}){
		print OI "$confInsert1{$i}\n$basic2{$i}\n";
	}elsif(exists $confInsert2{$i} && !exists $confInsert1{$i}){
		print OI "$basic1{$i}\n$confInsert2{$i}\n";
		#### print insertion from Chromosome 3 other location
	}elsif(exists $potential1_chr{$i} && exists $potential2_chr{$i}){
		print CI "$potential1_chr{$i}\n$potential2_chr{$i}\n";
		
		
		##### print insertion from Chromosome 3 other location but other elements are also inside
	}elsif(exists $potential1_chr{$i} && !exists $potential2_chr{$i}){
		print OIC "$basic1{$i}\n$basic2{$i}\n";	
		
	}elsif(exists $potential2_chr{$i} && !exists $potential1_chr{$i}){
		print OIC "$basic1{$i}\n$basic2{$i}\n";
			
		#### print unconfirmed insertion or deletions of exactly 294345-294384 
		### normally this happened with small insertion or small deletion based on their correpsonding size
		### Here we perform and function in order to reconsider the best last that showed aligned to 13600- 13700
	}elsif(exists $potential1_two{$i} && exists $potential2_two{$i} ){
		#print GI "$basic1{$i}\n$basic2{$i}\n";
		###
		### this categories are deletions:
		print GID "$basic1{$i}\n$basic2{$i}\n" if ($inf1{$i}->{length} <= 90 || $inf2{$i}->{length}<=90);
		
		### this categories are insertions that the overlapping regions are unknown and cannot defined the length
		
		### unless the both reads length are less than 300bp
		
		print GII "$basic1{$i}\n$basic2{$i}\n" if ($inf1{$i}->{length} >90 && $inf2{$i}->{length} >90 && $inf1{$i}->{start} < 10 && $inf2{$i}->{start} < 10);
		
		### This categories should be small deletions
		print GIU "$basic1{$i}\n$basic2{$i}\n" if ($inf1{$i}->{length} >90 && $inf2{$i}->{length} >90 && ($inf1{$i}->{start} < 10 && $inf2{$i}->{start} > 10 && $inf2{$i}->{start} < 50));
		print GIU "$basic1{$i}\n$basic2{$i}\n" if ($inf1{$i}->{length} >90 && $inf2{$i}->{length} >90 && ($inf2{$i}->{start} < 10 && $inf1{$i}->{start} > 10 && $inf1{$i}->{start} < 50));
		print GIUT "$basic1{$i}\n$basic2{$i}\n" if ($inf1{$i}->{length} >90 && $inf2{$i}->{length} >90 && ($inf1{$i}->{start} > 50 || $inf2{$i}->{start} > 50));
		
		#### print unconfirmed insertion or deletions of unconfirmed location
		### normally this happened with small insertion or small deletion based on their correpsonding si
	# }elsif(exists $check_all1{$i} && exists $check_all2{$i} ){
# 		###print PI "$basic1{$i}\n$basic2{$i}\n";
#
# 		print PID "$basic1{$i}\n$basic2{$i}\n" if ($int1{$i}->{length} <= 90 || $int2{$i}->{length}<=90);
# 		print PII "$basic1{$i}\n$basic2{$i}\n" if ($int1{$i}->{length} >90 && $int2{$i}->{length} >90 && $int1{$i}->{start} < 10 && $int2{$i}->{start} < 10);
# 		print PIU "$basic1{$i}\n$basic2{$i}\n" if ($int1{$i}->{length} >90 && $int2{$i}->{length} >90 && ($int1{$i}->{start} > 10 || $int2{$i}->{start} > 10));
#
		### print unknown information
	}else{
		print UN "$basic1{$i}\n$basic2{$i}\n";
	}

}

close R1;
close R2;
close OI;
close PI;
close CI;
close UN;
close GID;
close GII;
close OUT;
close GIU;
# close PID;
# close PII;
# close PIU;
