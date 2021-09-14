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
	Usage: $0.pl 
				-a Assembled fasta file (Test_detected_assemblysep.Insassembled.fasta)
				-b Reads with Quality results (Test_QC.Allquality.stat)
				-c Final unique single donor insertions (this contain the reads coverage) 
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

		### the identity and score value are important value to define whether they do have same elements
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
open INS,">$output.fasta" or die $!;
my $n=0;
print OUT "Case\tNID\tGenotype\tWholeinsertion\tInsertionSequence\tInsertionSize\tDonorNumber\tIdentity\tTotalCoverage\tQuality\tDonorChr\tDonorStart\tDonorEnd\tStrand\tDonorFeature\tDistanceToFeature\tDistanceDescript\tAlterFeature\tJunction_5\tJunction_3\n";
while(<B>){
	chomp;
	my ($cov,$id,$chr,$iden,$match,$start,$end,$lengthall,$Rstart,$Rend,$fcov)=(split/\t/,$_)[0,1,2,3,4,7,8,9,10,11,15];
	
	$n++;
	### Annotate with the representive reads quality information:
	
	my $fqual=$qual{$id};
	
	### the coverage, identity inforamtion of one sample -- $fcov,$iden,
	
	
	### the start, end , and strand information;
	my $strand=($Rstart<$Rend)?"+":"-";
	my $dstart=($Rstart<$Rend)?$Rstart:$Rend;
	my $dend=($Rstart<$Rend)?$Rend:$Rstart;
			

	
	### The length :
	
	my $startF= ($start<=50)?$start:44;
	
	my $endF=(($lengthall-$end)<=60)?$end:($lengthall-50);
	my $length=$endF-$startF+1;
	my $insertion=substr $string{$id},$startF, $length;
	
	
	### The inserted junction information (30bp upstream and downstream)
	my $start_3=$startF-15;
	my $start_5=$endF-15;
	my $upjunction=substr $string{$id}, $start_3,30;
	my $downjunction=substr $string{$id}, $start_5,30;
	
	
	### calculate the donor:
	my $ndonor=1;
	if ($start>60 || ($lengthall-$end)>60){
		$ndonor="1orMore";
	}
	
	
	### output with the PROXIMITY, ENTIRELY, PARTIALLY
	### the annotation of distance, approximately,overlapping this insertion,
	my $inse=join ":",($chr,$dstart,$dend);
	my ($chrA,$startA,$endA,$annotation,$dist)=split ":",$element{$inse};
	
	my $feature;
	if ($dist != 0){
		$feature="PROXIMITY";
	}elsif($dstart>=$startA && $dend<=$endA){
		$feature="ENTIRE";
	}else{
		$feature="PARTIALLY";
	}


	##### the alterative annotation from the blast;
	
	my @array=split/\n/,$str{$id}->{string};
	
	my $FMstring;
	if ($#array>1){
		foreach my $i (1..$#array){
			my ($Mchr,$Midentity,$Mmatch,$Mstart,$Mend)=(split/\t/,$array[$i])[1,2,3,9,10];
		
			
			my $Mtype2=($Mstart<$Mend)?"+":"-";
			my $Mstart2=($Mstart<$Mend)?$Mstart:$Mend;
			my $Mend2=($Mstart<$Mend)?$Mend:$Mstart;
			
			my $com=join ":",($Mchr,$Mstart2,$Mend2);
			##### if need to add the Melement ### we need to generate the file
			my $Melement=$element{$com};
			my $Mstring=join":",($Mchr,$Midentity,$Mmatch,$Mtype2,$Mstart2,$Mend2,$Melement);
		
			$FMstring.=$Mstring.";";
		}
	}else{
		$FMstring="NO";
	}
	
	print OUT "$n\t$id\t$sample\t$string{$id}\t$insertion\t$length\t$ndonor\t$iden\t$fcov\t$fqual\t$chr\t$dstart\t$dend\t$strand\t$annotation\t$dist\t$feature\t$FMstring\t$upjunction\t$downjunction\n";
	print INS ">$id\n$string{$id}\n";
	
	#print OUT "$n\t$sample\t$string{$id}\t$insertion\t$length\t$ndonor\t$identity\t$cov\t$covA\t$covB\t$chr2\t$dstart\t$dend\t$strand\t$annotation\t$dist\t$FMstring\t$upjunction\t$downjunction\n";
	
}

close B;
close OUT;
close INS;



