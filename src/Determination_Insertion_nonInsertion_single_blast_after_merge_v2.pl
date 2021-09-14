#!/usr/bin/perl
use strict;
use warnings;


### This script is to divide all the reads into different categories: noninserted, confident-inserted, potential inserted, potential deleted, or unclassified.
### requirement: 1
#####1.	Completely blasted against reference genome from start to end of the reads (allow for 2 gaps and 2 mismatches at most);
#####2.	Forward and Reward regards consistently mapped against the identical reference location.
#### If only they have one read mapping information

#use Bio::Seq;
#use Bio::SeqIO;


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s","q:s","f:s","h:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{o} || !defined $opts{q} || !defined $opts{f} || defined $opts{h}) {
	die "************************************************************************
	Usage: $0.pl -i ReadBLast -q Read_Quality -f Assembled_fasta -o Blast_categories
	
	Request Parameters:
	-i The blast result of each reads
	-q The quality of each reads, here we specific focused on the average read quality that eliminated the MAT.
	-f The assembled fasta file
	-o Output of different categories of the reads
	
	Optional Parameters:
	
 	-c The chromosme of MATA reference position (chrIII:294345-294428, default chrIII)
 	-t The start site of MATA reference position (chrIII:294345-294428, default 294300)
 	-e The end site of MATA reference position (chrIII:294345-294428, default 294500)
 	-m The mininus length of large insertions (default 10bp)
 	-n The size of MATA region (default 90bp)
	-h Help
************************************************************************\n";
}

my $output=$opts{o};
my $input=$opts{i};
my $RQual=$opts{q};

my $FA=$opts{f};



my $chr=(defined $opts{c})?$opts{c}:"chrIII";
my $Rstart= (defined $opts{t})?$opts{t}:294300;
my $Rend= (defined $opts{e})?$opts{e}:294500;
my $Isize= (defined $opts{m})?$opts{m}:10;
my $Msize= (defined $opts{n})?$opts{n}:90;

my $MinRegionSize=$Isize+$Msize;


### Index of fasta sequence

open FA,"$FA" or die $!;

my $id;  my %hash;  my %sequence; my %len;
while (<FA>) {
	chomp;
    # print "$_" ;
	if (/^>(\S+)/){
		$id=$1;
		$hash{$id}++;
	}else{
		$sequence{$id}.=$_;
		$len{$id}=length $_;
	}

}
close FA;



my %conf1; my %inf1; my %int1; my %confInsert1; my %potential1_two; my %potential1_chr; my %check_all1; my %basic1;

### Make the index for the quality of reads, here we did not take into consideration the R1 upstream and R2 30bp due to they mainly came from MAT region with high quality.

my %qual;
open QU, "$RQual" or die $!;
while (<QU>){
	chomp;
	my ($id,$readqual)=(split/\t/,$_)[0,1];
	$qual{$id}=$readqual;
}

close QU;

unless (-e $output){
 	system ("mkdir $output");
 }else{
	 print "$output have already existed, no need to create! Continue ...\n";
 }


## Store all the insertion events and assigned them with the insertion or non-insertion events

open I,"$input" or die $!;

while (<I>) {
	chomp;
    # print "$_" ;
	my ($qid,$pid,$identity,$match,$gap,$mismatch,$qstart,$qend,$qlength,$pstart,$pend,$plength,$score,$evlaue)=split/\s+/,$_;
	my $information=$_;
	$basic1{$qid}=$information;
	
	#### Here due to mutiple alignments of MAT right side HMR, HML.
	next if ($pid eq "chrIII" && $pstart >= 13650 && $pend<= 13850 );
	next if ($pid eq "chrIII" && $pstart >= 200750 && $pend <= 201000);
	
	### Usually the best matches were identified from the first one, so if we could identifed them as non-inserted, potential insertion, or confident insertion, we will ignore the following blast results.
	
	next if (exists $conf1{$qid});
	next if (exists $confInsert1{$qid});
	next if (exists $potential1_chr{$qid});
	
	### We tried to measure the mapping coverage of the whole read, this would like to give brief information whether they have the large insertion events.
	my $cov=$match/$qlength;
	#print "$cov\n";
	
	### We also put the quality of reads 
	my $quality=$qual{$qid};

	#### Firstly, confident insertion from other chromosome insertion
	### if the reads have blast results from other chrosome, we confidently considered them as read with large insertion events
	if ($pid ne "chrIII"){
		$confInsert1{$qid}=$information;	
	###### confident non-insertion, here might contain the short indels. Due to the concentration of large insertion, we will take these result into consideration.
 	}elsif($pstart>=$Rstart  && $pstart <=$Rend && $pend>=$Rstart && $pend<=$Rend && $match>=70 && $qlength <=$MinRegionSize ){
		$conf1{$qid}=$information;	

		
	##### confident insertion from other location of chr3
	}elsif($pstart>=$Rend  || $pstart<=$Rstart){
		$confInsert1{$qid}=$information;
			
	 #### potential insertion or deletion but need to check original data  
	 #### potential split the fragment with insertion (telemere repetitive elements or MAT, therefore these regions located with two sides of MATA)
	}elsif($pstart>=$Rstart  && $pstart <=$Rend && $pend>=$Rstart && $pend<=$Rend && $match<70 && $qlength >$MinRegionSize){
		
		$potential1_chr{$qid}=$information;
	######## This might contain lots of short deletions. Since the short indel we applied other strategies, we will not focus on this part.
	}else{
		my $non++;
		#print "$information\n"  ;
	}
		
}

close I;


open OUT,">$output/$output.noninsert.confident.txt" or die $!;


### including the insertion from other chromosome and chromosome:
open CI,">$output/$output.insert.blast.txt" or die $!;
#### one with unknown insertion another with insertion from other location.

#### Potential insertion only showed two edge of insertion events 4-44 for two PE, and reads length is longer than 90bp

open UN,">$output/$output.unclassified.txt" or die $!;

open OUTA,">$output/$output.Ainsertion.fasta" or die $!;

open OUTS,">$output/$output.Ainsertion.quality.txt" or die $!;

open OUTU,">$output/$output.Anoinsertion.quality.txt" or die $!;

print OUTS "ID\tReadQuality\tRIdentity\tRMatches\tRLength\n";


foreach my $m (keys %hash){

	if (exists $confInsert1{$m}){
		print CI "$confInsert1{$m}\n";
		
		my ($identity,$match,$qlength)=(split/\t/,$confInsert1{$m})[2,3,8];
		print OUTS "$m\t$qual{$m}\t$identity\t$match\t$qlength\n";
		print OUTA ">$m\n$sequence{$m}\n";

	}elsif (exists $potential1_chr{$m}){
		print CI "$potential1_chr{$m}\n";
		
		my ($identity,$match,$qlength)=(split/\t/,$potential1_chr{$m})[2,3,8];
		print OUTS "$m\t$qual{$m}\tNA\tNA\t$qlength\n";
		print OUTA ">$m\n$sequence{$m}\n";

	}elsif(exists $conf1{$m}){
		print OUT "$conf1{$m}\n";
		my ($identity,$match,$qlength)=(split/\t/,$conf1{$m})[2,3,8];
		print OUTU "$m\t$qual{$m}\t$identity\t$match\t$qlength\n";
		
	# the left if the length longer than 100bp, we considered them as with insertion	
	}elsif($len{$m} > $MinRegionSize){
		
		print CI "$m\tNA\n";
		print OUTS "$m\t$qual{$m}\tNA\tNA\t$len{$m}\n";
		print OUTA ">$m\n$sequence{$m}\n";
	
	}elsif(exists $basic1{$m}){
		print UN "$basic1{$m}\n";
		my ($identity,$match,$qlength)=(split/\t/,$basic1{$m})[2,3,8];
		print OUTU "$m\t$qual{$m}\t$identity\t$match\t$qlength\n";
		
	}else{
		print UN "$m\tUnaligned\n";
		print OUTU "$m\t$qual{$m}\tNA\tNA\t$len{$m}\n";
	}
	
}




close OUTS;
close CI;
close OUT;
close OUTA;
close I;
close UN;






