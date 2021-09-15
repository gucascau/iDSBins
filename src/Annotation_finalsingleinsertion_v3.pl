#!/usr/bin/perl
use strict;
use warnings;
my $version="1.0 version";
use Getopt::Long;

#Author:Xin Wang 
#email: xin.wang@childrens.harvard.edu
#PI: Kaifu Chen
### #### caluculate the side effect on the final associated genes results
### Here is to calculate the sample size effect

### This script is to generate the large insertion table and bedfile for further annotation
my %opts;
GetOptions(\%opts,"i:s","o:s","g:s","q:s","t:s","h:s","e:s","d:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{g} || !defined $opts{q} ||!defined $opts{i}||!defined $opts{o}||!defined $opts{t} ||!defined $opts{d}||defined $opts{h}) {
	die "************************************************************************
	Usage: $0.pl -g Assembled fasta -q Reads quality file -i final single insertion -d All assembled blast results -e Annotation of each locus -o Output string name
	
	Request Parameters:
	 
	-g: Assembled fasta file
	-q: Quality of each read file
	-d: Blast with all possible hits
	-e: Annotation of different elements based on the chr,start
	-i: Final insertion with single donor information
	-t: sample name
	-o: Output string files, including insertion assembled fasta, assembled bed file for further annotation, final primilary table with statistic results
	
	Optional Parameters:
	-h Help
	
************************************************************************\n";
}


#### The index of the reads sequence information

my $fasta=$opts{g};

open FASTA, "$fasta" or die "cannot open file $fasta";

my %seq; my $id3;
while (<FASTA>){
	chomp;	
	if ($_=~/>(\S+)/){
		$id3=$1;
	}else{
		$seq{$id3}.=$_;
	}  
}
close FASTA;


#### The index of the read quality information

my $readquality=$opts{q};
my %qual;
open RQ, "$readquality" or die $!;
while (<RQ>){
	chomp;
	my ($id,$quality)=(split/\t/,$_)[0,5];
	$qual{$id}=$quality;
}
close RQ;




####### this file is corresponding annnotation of each chr,start,end


my %element;

my $annot=$opts{e};

open E, "$annot" or die $!;

while (<E>){
	chomp;
	
	my ($chra,$starta,$enda,$elementa,$dista)=(split/\t/,$_)[0,1,2,12,13];
	# chrI	64	320	ID23025	0	+	chrI	63	336	ID2	0	-	X_element_combinatorial_repeat|TEL01L_X_element_combinatorial_repeat	0
	# chrI	427	675	ID33022	0	+	chrI	335	649	ID3	0	+	gene|YAL069W	0
	# chrI	530	675	ID4835	0	+	chrI	335	649	ID3	0	+	gene|YAL069W	0
	# chrI	1048	1197	ID3396	0	+	chrI	337	801	ID4	0	-	X_element|TEL01L_X_element	-248
	my $stringa=join ":",($elementa,$dista);
	
	my $indexa=join ":",($chra,$starta,$enda);
	
	$element{$indexa}=$stringa;
	
	
}

close E;


#### index for the mutiple blast for single donor ####

my $all_blast=$opts{d};

my %hash; my %str;
open I,"$all_blast" or die "cannot open file $all_blast";
while (<I>) {
	chomp;
    # print "$_" ;
	my ($qid,$pid,$identity,$match,$gap,$mismatch,$qstart,$qend,$qlength,$pstart,$pend,$plength,$score,$evlaue)=split/\s+/,$_;
	my $information=$_;
	#$basic1{$qid}=$information;
	next if ($qstart <5 || ($qlength-$qend) <5);
	#next if ($pid eq "ref-NC_001135-" && $pstart >= 13650 && $pend<= 13850 );
	#next if ($pid eq "ref-NC_001135-" && $pstart >= 200750 && $pend <= 201000);
	#next if ($identity <90);
	#next unless (exists $single{$qid});

	if (!exists $hash{$qid}){
		$str{$qid}->{chr}=$pid;
		$str{$qid}->{iden}=$identity;
		$str{$qid}->{matches}=$match;
		$str{$qid}->{score}=$score;
		#$str{$qid}->{evalue}=$evlaue;

		$str{$qid}->{qstart}=$qstart;
		$str{$qid}->{qend}=$qend;

		$hash{$qid}++;
		$str{$qid}->{string}=$information;
	}else{

		if ($identity/$str{$qid}->{iden} >= 0.99 && $score/$str{$qid}->{score} >= 0.95 ){
			$str{$qid}->{string}.="\n".$information;
			$hash{$qid}++;
		}
	}
}
close I;


my $out=$opts{o};
my $insert=$opts{i};

open INS,"$insert" or die $!;
open BED,">$out.Fsingle.bed" or die $!;
open FA,">$out.Fsingle.fasta" or die $!;
open FT, ">$out.Fsingle.table.txt" or die $!;

my $n=0;
my $sample=$opts{t};
while (<INS>) {
	chomp;
	s/\r//g;
	my ($cov0,$id,$chr,$iden,$match,$mismatch,$gapsize,$start,$end,$total,$Rstart,$Rend,$Rlength,$bitscore,$pvalue,$coverage)=split/\t/,$_;
	
#
# 1	M06255:18:000000000-JLRH5:1:1101:20211:4373	chrmt	93.939	99	4	2	44	140	171	47255	47353	85779	126	1.01e-29	1
# 1	M06255:18:000000000-JLRH5:1:1101:19370:2977	chrmt	84.615	130	6	13	44	165	197	77895	78018	85779	117	9.32e-27	1
# 1	M06255:18:000000000-JLRH5:1:1101:23060:18468	chrmt	96.296	54	2	0	44	97	136	37838	37785	85779	79.1	1.89e-15	1

	## annotation with different components:
	
	### the chromosome information
	### here is $chr
	
	### the start, end , and strand information;
	my $strand=($Rstart<$Rend)?"+":"-";
	my $dstart=($Rstart<$Rend)?$Rstart:$Rend;
	my $dend=($Rstart<$Rend)?$Rend:$Rstart;
			
	### the annotation of this insertion for elemement and distance
	my $inse=join ":",($chr,$dstart,$dend);
	
	my ($annotation,$dist)=split ":",$element{$inse};
	
	
	### The length 
	
	my $length=$end-$start+1;
	
	## the whole sequence information
	

	
	## The inserted sequence inforamtion:
	my $insertion=substr $seq{$id}, $start,$length;
	
	
	### The inserted junction information (30bp upstream and downstream)
	my $start_3=$start-15;
	my $start_5=$end-15;
	my $upjunction=substr $seq{$id}, $start_3,30;
	my $downjunction=substr $seq{$id}, $start_5,30;
	
	
	### calculate the donor:
	my $ndonor=1;
	if ($start>60 || ($total-$end)>60){
		$ndonor="1orMore";
	}
	

	### the coverage inforamtion of one sample
	## $coverage
		
	
	$n++;
	
	##### the alterative annotation from the blast;
	
	if (!exists $str{$id}->{string}){
		
		print FT "$n\t$id\t$sample\t$seq{$id}\t$insertion\t$length\t$ndonor\t$iden\t$coverage\t$chr\t$dstart\t$dend\t$strand\t$annotation\t$dist\tNA\n";
		print BED "$chr\t$dstart\t$dend\t$id\t$coverage\t$strand\n";
		print FA ">$id\n$seq{$id}\n";
		
		next;
	}
	
	my @array=split/\n/,$str{$id}->{string};
	
	my $FMstring;
	if ($#array>=1){
		foreach my $i (1..$#array){
			my ($Mchr,$Midentity,$Mmatch,$Mstart,$Mend)=(split/\t/,$array[$i])[1,2,3,9,10];
		
			my $Mchr2=$chr;
			my $Mtype2=($Mstart<$Mend)?"+":"-";
			my $Mstart2=($Mstart<$Mend)?$Mstart:$Mend;
			my $Mend2=($Mstart<$Mend)?$Mend:$Mstart;
			
			my $com=join ":",($Mchr2,$Mstart2,$Mend2);
			##### if need to add the Melement ### we need to generate the file
			my $Melement=$element{$com};		

			my $Mstring=join":",($Mchr2,$Midentity,$Mmatch,$Mtype2,$Mstart2,$Mend2,$Melement);
		
			$FMstring.=$Mstring.";";
		
		
		}
	}else{
		$FMstring="NO";
	}
	
	print FT "$n\t$id\t$sample\t$seq{$id}\t$insertion\t$length\t$ndonor\t$iden\t$coverage\t$chr\t$dstart\t$dend\t$strand\t$annotation\t$dist\t$FMstring\n";
	print BED "$chr\t$dstart\t$dend\t$id\t$coverage\t$strand\n";
	print FA ">$id\n$seq{$id}\n";
	#print FT "$array[1]\t$string\t$ndonor\t$coverage\t$identity\t$Fquality\t$Donorseq\t$Donorlength\t$chro\t$start\t$end\t$type\n";
	
}
close INS;
close FT;
close FA;
close BED;


