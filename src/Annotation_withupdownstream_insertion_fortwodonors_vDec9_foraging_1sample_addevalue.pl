#!/usr/bin/perl
use strict;
use warnings;


#### This scripts is to divide into single insertion


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"a:s","b:s","o:s","t:s","c:s","d:s","e:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{a} || !defined $opts{b} || !defined $opts{o} || !defined $opts{t}|| !defined $opts{c} || !defined $opts{d} || !defined $opts{e}  ) {
	die "************************************************************************
	Usage: $0.pl	 -a fasta file
				-b Two blast result
				-c chromosome changes
				-d Blast with all possible hits
				-e Annotation of different elements based on the chr,start
				-t sample name
				-o Output with blast result
************************************************************************\n";
}



my $fasta=$opts{a};
my $blast=$opts{b};
my $sample=$opts{t};

my $alteration=$opts{c};
my $all_blast=$opts{d};
my $anno=$opts{e};


my $output=$opts{o};

my %string;
open A,"$fasta" or die $!;
my $id;

while (<A>){
	chomp;
	if (/^>(\S+)/){
		$id=$1;
	}else{
		$string{$id}.=$_;
	}
}

close A;

#### index for chromosome id ####

my %alter;

open A,"$alteration" or die "cannot open file $alteration";
while (<A>){
	chomp;
	s/\r//g;
	my ($i,$ch)=split /\t/,$_;

	$alter{$i}=$ch;

}
close A;



#### index for the mutiple blast for mutiple donor ####

#### Here is specific for the alterantive information

open I,"$all_blast" or die "cannot open file $all_blast";

my %str; my %hash;


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
		$str{$qid}->{evalue}=$evlaue;
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


####### this file is corresponding annnotation of each chr,start,end
open E,"$anno" or die $!;

my %element; my %inf;

while (<E>){
	chomp;
	my ($chra,$starta,$enda,$strand,$iden,$bitscore,$evalue,$elementa,$dista)=(split/\t/,$_)[0,1,2,5,6,7,8,16,17];
	my $stringa=join ":",($elementa,$dista);
	
	my $indexa=join ":",($chra,$starta,$enda);
	
#	$element{$indexa}=$stringa;
	
	$inf{$indexa}->{strand}=$strand;
	$inf{$indexa}->{identity}=$iden;
	$inf{$indexa}->{element}=$stringa;
	$inf{$indexa}->{evalue}=$evalue;
	$inf{$indexa}->{bitscore}=$bitscore;
}

close E;




##### ##### here for two donor ###

open B,"$blast" or die $!;
open OUT,">$output" or die $!;

my $n=0;
print OUT "Case\tID\tGenotype\tWholeinsertion\tInsertionSequence\tInsertionSize\tDonorNumber\tIdentity\tTotalCoverage\tBitscore\tEvalue\tDonorChromosome\tDonorStart\tDonorEnd\tStrand\tAnnotationFeature\tDistanceToFeature\tAlterFeature\tJunction_5\tJunction_3\n";
while(<B>){
	chomp;
	
	my ($lengthall,$id,$OchrA,$qstartA,$qendA,$rstartA1,$rendA1,$length1,$OchrB,$qstartB,$qendB,$rstartB2,$rendB2,$length2)=split/\t/,$_;
	
	my $length=$qendB-$qstartA+1;
	my $insertion=substr $string{$id},$qstartA, $length;
	
	#### modified the start and end for two insertions
	my $rstartA=($rendA1>$rstartA1)?$rstartA1:$rendA1;
	my $rendA=($rendA1>$rstartA1)?$rendA1:$rstartA1;
	
	my $rstartB=($rendB2>$rstartB2)?$rstartB2:$rendB2;
	my $rendB=($rendB2>$rstartB2)?$rendB2:$rstartB2;
	
	
	my $start_3=$qstartA-15;
	my $start_5=$qendB-15;
	my $upjunction=substr $string{$id}, $start_3,30;
	my $downjunction=substr $string{$id}, $start_5,30;
	
	$n++;
	
	my $chrA=$alter{$OchrA};
	my $chrB=$alter{$OchrB};
	
	##### output for the id
	#my @array=split/\n/,$str{$id}->{string};
	#my ($chr2,$dstart,$dend,$im,$nons,$strand,$tcov,$identity,$ndonor,$cov,$covA,$covB,$annotation,$dist)=split/\:/,$id;
	### DNA2 cause had four samples
	#my ($chr2,$dstart,$dend,$im,$nons,$strand,$tcov,$identity,$ndonor,$cov,$covA,$covB,$covC,$covD,$annotation,$dist)=split/\:/,$id;
	


	my $cov=(split/\|/,$id)[3];

	
	##### find the insertion of the first Part
	
	
	my $lengthA=$qendA-$qstartA+1;
	my $insertionA= substr $string{$id},$qstartA,$lengthA;

	my $partA=join ":",($chrA,$rstartA,$rendA);

	my $strandA=$inf{$partA}->{strand};
	#my $annoA=$inf{$partA}->{element};
	my $idenA=$inf{$partA}->{identity};
	my $evalueA=$inf{$partA}->{evalue};
	my $bitscoreA=$inf{$partA}->{bitscore};
	
	
	my ($annoA,$distA)=split ":",$inf{$partA}->{element};

	
	##### find the insertion of the second Part

	
	
	
	my $lengthB=$qendB-$qstartB+1;
	my $insertionB= substr $string{$id},$qstartB,$lengthB;

	my $partB=join ":",($chrB,$rstartB,$rendB);

	my $strandB=$inf{$partB}->{strand};
	my ($annoB,$distB)=split ":",$inf{$partB}->{element};
	
	my $idenB=$inf{$partB}->{identity};
	my $evalueB=$inf{$partB}->{evalue};
	my $bitscoreB=$inf{$partB}->{bitscore};
	
	### Calculate the donor number of insertions

	my $ndonor=2;
	if ($qstartA >60 || abs($qstartB-$qendA) >20 || ($lengthall-$qendB)>60){
		$ndonor="2orMore";

	}

	
	#### print the overall insertion ##
	my $chrall=join";",($chrA,$chrB);
	my $idenall=join ";",($idenA,$idenB);
	
	my $aveevalue= ($evalueA+$evalueB)/2;
	my $avebitscore= ($bitscoreA+$bitscoreB)/2;
	
	print OUT "T$n\t$id\t$sample\t$string{$id}\t$insertion\t$length\t$ndonor\t$idenall\t$cov\t$avebitscore\t$aveevalue\t$chrall\tNO\tNO\tNO\tNO\tNO\tNO\t$upjunction\t$downjunction\n";
	
	### print the first element of insertion #
	
	
	##### how to add the alterative feature 
	print OUT "T$n\.A\t$id\t$sample\t$string{$id}\t$insertionA\t$lengthA\tT1\t$idenA\t$cov\t$bitscoreA\t$evalueA\t$chrA\t$rstartA\t$rendA\t$strandA\t$annoA\t$distA\tNO\tNO\tNO\n";
	
	##### how to add the alterative feature 
	### print the second element of insertion #
	print OUT "T$n\.B\t$id\t$sample\t$string{$id}\t$insertionB\t$lengthB\tT2\t$idenB\t$cov\t$bitscoreA\t$evalueA\t$chrB\t$rstartB\t$rendB\t$strandB\t$annoB\t$distB\tNO\tNO\tNO\n";
	
	
	
	#print OUT "$n\t$sample\t$string{$id}\t$insertion\t$length\t$ndonor\t$identity\t$cov\t$covA\t$covB\t$chr2\t$dstart\t$dend\t$strand\t$annotation\t$dist\t$FMstring\t$upjunction\t$downjunction\n";
	
		#
	#
	# #### output for the alteration sequnence;
	# my @array=split/\n/,$str{$id}->{string};
	#
	#
	#
	#
	#
	#
	#
	# my ($id,$chr1,$iden,$match,$start,$end)=(split/\t/,$_)[0,1,2,3,6,7];
	# my $length=$end-$start+1;
	# my $insertion=substr $string{$id}, $start,$length;
	#
	#
	# my $start_3=$start-15;
	# my $start_5=$end-15;
	# my $upjunction=substr $string{$id}, $start_3,30;
	# my $downjunction=substr $string{$id}, $start_5,30;
	#
	#
	#
	# my ($chr2,$dstart,$dend,$im,$nons,$strand,$tcov,$identity,$ndonor,$cov,$covA,$covB,$annotation,$dist)=split/\:/,$id;
	# $n++;
	#
	#
	# ##### output for the alteration sequence;
	#
	# my @array=split/\n/,$str{$id}->{string};
	#
	# my $FMstring;
	# if ($#array>=1){
	# 	foreach my $i (1..$#array){
	# 		my ($Mchr,$Midentity,$Mmatch,$Mstart,$Mend)=(split/\t/,$array[$i])[1,2,3,9,10];
	#
	# 		my $Mchr2=$alter{$Mchr};
	# 		my $Mtype2=($Mstart<$Mend)?"+":"-";
	# 		my $Mstart2=($Mstart<$Mend)?$Mstart:$Mend;
	# 		my $Mend2=($Mstart<$Mend)?$Mend:$Mstart;
	#
	# 		my $com=join ":",($Mchr2,$Mstart2,$Mend2);
	# 		##### if need to add the Melement ### we need to generate the file
	# 		my $Melement=$element{$com};
	#
	#
	#
	# 		my $Mstring=join":",($Mchr2,$Midentity,$Mmatch,$Mtype2,$Mstart2,$Mend2,$Melement);
	#
	# 		$FMstring.=$Mstring.";";
	#
	#
	# 	}
	# }else{
	# 	$FMstring="NO";
	# }
	#
	#
	#
	#
	# print OUT "$n\t$sample\t$string{$id}\t$insertion\t$length\t$ndonor\t$identity\t$cov\t$covA\t$covB\t$chr2\t$dstart\t$dend\t$strand\t$annotation\t$dist\t$FMstring\t$upjunction\t$downjunction\n";
	#
}

close B;
close OUT;


