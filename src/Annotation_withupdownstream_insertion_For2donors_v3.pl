#!/usr/bin/perl
use strict;
use warnings;


#### This script is generate the final table for the single donor insertions


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"a:s","b:s","o:s","t:s","c:s","d:s","e:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{a} || !defined $opts{b} || !defined $opts{o} || !defined $opts{t}|| !defined $opts{c} || !defined $opts{d} || !defined $opts{e}  ) {
	die "************************************************************************
	Usage: $0.pl -a AssembledFasta -b QualityFile -c FinalTwodonorInsertion	-d BLASThit -e BLASTAnnotation -t Sample name -o Output string
				-a Assembled fasta file (Test_detected_assemblysep.Insassembled.fasta)
				-b Reads with Quality results (Test_QC.Allquality.stat)
				-c Final unique two donor insertions (Test_detected_second_twoassembly.uniq.txt, this contain the reads coverage) 
				-d Blast with all possible hits 
				-e Annotation of different elements based on the chr,start
				-t Sample name
				-o Output with blast result
************************************************************************\n";
}



my $fasta=$opts{a};
my $Quality=$opts{b};
my $sample=$opts{t};

my $blast=$opts{c};
my $all_blast=$opts{d};
my $anno=$opts{e};


my $output=$opts{o};
my %string;

### Read the fasta file
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


### Read the read's quality file
my %qual;
open QUAL,"$Quality" or die $!;
while (<QUAL>){
	chomp;
	my ($id,$qu)=(split/\t/,$_)[0,5];
	$qual{$id}=$qu;
}
close QUAL;

#### index for the mutiple blast for single donor ####

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
my %element;

while (<E>){
	chomp;
	my ($chra,$starta,$enda,$chrb,$startb,$endb, $elementa,$dista)=(split/\t/,$_)[0,1,2,6,7,8,12,13];
	
	my ($start,$end)=($enda>$starta)?($starta,$enda):($enda,$starta);
	
	my $stringa=join ":",($chrb,$startb,$endb, $elementa,$dista);
	
	my $indexa=join ":",($chra,$start,$end);
	
	$element{$indexa}=$stringa;
	
	
}

close E;


##### 
open B,"$blast" or die $!;
open OUT,">$output.Ftable" or die $!;
open INS,">$output.2donor.fasta" or die $!;
open INSSEP,">$output.2donorsep.fasta" or die $!;
my $n=0;
print OUT "Case\tNID\tGenotype\tWholeinsertion\tInsertionSequence\tInsertionSize\tDonorNumber\tIdentity\tTotalCoverage\tQuality\tDonorChr\tDonorStart\tDonorEnd\tStrand\tDonorFeature\tDistanceToFeature\tDistanceDescript\tAlterFeature\tJunction_5\tJunction_3\n";
while(<B>){
	chomp;
	
	my ($fcov,$id,$cov,$lengthall,$id2,$chrA,$startA,$endA,$RstartA,$RendA,$lengthA0,$identA,$chrB,$startB,$endB,$RstartB,$RendB,$lengthB0,$identB)=split/\t/,$_;
	
	$n++;
	
	
	### The whole inserted region: Here we allow 10bp overhang region
	
	my $startF= ($startA<=50)?$startA:45;
	
	my $endF=(($lengthall-$endB)<=60)?$endB:($lengthall-50);
	my $length=$endF-$startF+1;
	my $insertion=substr $string{$id},($startF-1), $length;
	
	### The first inserted junction information (30bp upstream and downstream of inserted elements)
	my $start_3=$startF-15;
	
	my $start_5=$endF-15;
	my $upjunction=substr $string{$id}, $start_3,30;
	my $downjunction=substr $string{$id}, $start_5,30;
	
	

	
	##### find the insertion of the first Part
	
	

	my $strandA=($RstartA<$RendA)?"+":"-";
	my $dstartA=($RstartA<$RendA)?$RstartA:$RendA;
	my $dendA=($RstartA<$RendA)?$RendA:$RstartA;
	
	my $lengthA=$endA-$startA+1;
	
	# the first inserted sequence
	my $insertionA=substr $string{$id},($startA-1),$lengthA;
	
	
	### output with the PROXIMITY, ENTIRELY, PARTIALLY
	### the annotation of distance, approximately,overlapping this insertion,
	my $partA=join ":",($chrA,$dstartA,$dendA);
	my ($EchrA,$EstartA,$EendA,$EannotationA,$distA)=split ":",$element{$partA};
	
	my $featureA;
	if ($distA != 0){
		$featureA="PROXIMITY";
	}elsif($dstartA>=$EstartA && $dendA<=$EendA){
		$featureA="ENTIRE";
	}else{
		$featureA="PARTIALLY";
	}
	
	
	
	##### find the insertion of the second Part


	my $strandB=($RstartB<$RendB)?"+":"-";
	my $dstartB=($RstartB<$RendB)?$RstartB:$RendB;
	my $dendB=($RstartB<$RendB)?$RendB:$RstartB;
	
	
	my $lengthB=$endB-$startB+1;
	my $insertionB= substr $string{$id},($startB-1),$lengthB;

	my $partB=join ":",($chrB,$dstartB,$dendB);
	my ($EchrB,$EstartB,$EendB,$EannotationB,$distB)=split ":",$element{$partB};

	my $featureB;
	if ($distB != 0){
		$featureB="PROXIMITY";
	}elsif($dstartB>=$EstartB && $dendB<=$EendB){
		$featureB="ENTIRE";
	}else{
		$featureB="PARTIALLY";
	}
	
	
	### Calculate the donor number of insertions

	my $ndonor=2;
	if ($startA >60 || ($startB-$endA) >20 || ($lengthall-$endB)>60){
		$ndonor="2orMore";

	}

	
	#### print the overall insertion ##
	my $chrall=join";",($chrA,$chrB);
	my $idenall=join ";",($identA,$identB);
	
	my $quality=$qual{$id};
	print INSSEP ">Tw$n.A\t$id\n$insertionA\n>Tw$n.B\t$id\n$insertionB\n";
	print INS ">$id\n$string{$id}\n";
	
	print OUT "Tw$n\t$id\t$sample\t$string{$id}\t$insertion\t$length\t$ndonor\t$idenall\t$fcov\t$quality\t$chrall\tNO\tNO\tNO\tNO\tNO\tNO\tNO\t$upjunction\t$downjunction\n";
	
	### print the first element of insertion #
	
	
	##### how to add the alterative feature 
	print OUT "Tw$n\.A\t$id\t$sample\t$string{$id}\t$insertionA\t$lengthA\tT1\t$identA\t$fcov\t$quality\t$chrA\t$dstartA\t$dendA\t$strandA\t$EannotationA\t$distA\t$featureA\tNO\tNO\tNO\n";
	
	##### how to add the alterative feature 
	### print the second element of insertion #
	print OUT "Tw$n\.B\t$id\t$sample\t$string{$id}\t$insertionB\t$lengthB\tT2\t$identB\t$fcov\t$quality\t$chrB\t$dstartB\t$dendB\t$strandB\t$EannotationB\t$distB\t$featureB\tNO\tNO\tNO\n";
	
	
	
}

close B;
close OUT;
close INS;
close INSSEP;



