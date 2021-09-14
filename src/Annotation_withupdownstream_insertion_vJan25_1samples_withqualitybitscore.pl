#!/usr/bin/perl
use strict;
use warnings;


#### This scripts is to divide into single insertion for further annotation with mutiple alignments


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"a:s","b:s","o:s","t:s","c:s","d:s","e:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{a} || !defined $opts{b} || !defined $opts{o} || !defined $opts{t}|| !defined $opts{c} || !defined $opts{d} || !defined $opts{e}  ) {
	die "************************************************************************
	Usage: $0.pl -a fasta file
				-b single blast result
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
	my ($chra,$starta,$enda,$elementa,$dista)=(split/\t/,$_)[0,1,2,16,17];
	my $stringa=join ":",($elementa,$dista);
	
	my $indexa=join ":",($chra,$starta,$enda);
	
	$element{$indexa}=$stringa;
	
	
}

close E;







##### 
open B,"$blast" or die $!;
open OUT,">$output" or die $!;

my $n=0;
print OUT "Case\tNID\tGenotype\tWholeinsertion\tInsertionSequence\tInsertionSize\tDonorNumber\tIdentity\tTotalCoverage\tBitscore\tPvalue\tRead\tFupstreamq\tFdownstreamq\tRupstreamq\tRdownstreamq\tDonorChromosome\tDonorStart\tDonorEnd\tStrand\tAnnotationFeature\tDistanceToFeature\tAlterFeature\tJunction_5\tJunction_3\n";
while(<B>){
	chomp;
	my ($id,$chr1,$iden,$match,$start,$end,$total,$Rstart,$Rend,$bitscore,$pvalue)=(split/\t/,$_)[0,1,2,3,6,7,8,9,10,12,13];
	
	### Annotate with the representive reads quality information:
	
	
	
	
	
	### the chromosome information
	my $chr2=$alter{$chr1};
	
	### the start, end , and strand information;
	my $strand=($Rstart<$Rend)?"+":"-";
	my $dstart=($Rstart<$Rend)?$Rstart:$Rend;
	my $dend=($Rstart<$Rend)?$Rend:$Rstart;
			
	### the annotation of this insertion
	my $inse=join ":",($chr2,$dstart,$dend);
	
	my ($annotation,$dist)=split ":",$element{$inse};
	
	
	### The length :
	
	my $length=$end-$start+1;
	
	## The inserted sequence inforamtion:
	my $insertion=substr $string{$id}, $start,$length;
	
	
	### The inserted junction information (30bp upstream and downstream)
	my $start_3=$start-15;
	my $start_5=$end-15;
	my $upjunction=substr $string{$id}, $start_3,30;
	my $downjunction=substr $string{$id}, $start_5,30;
	
	
	### calculate the donor:
	my $ndonor=1;
	if ($start>60 || ($total-$end)>60){
		$ndonor="1orMore";
	}
	
	
	### the coverage inforamtion of one sample
	my ($Nid,$annot,$st,$cov)=split/\|/,$id;
	
	
	# my ($chr2,$dstart,$dend,$im,$nons,$strand,$tcov,$identity,$ndonor,$cov,$covA,$covB,$annotation,$dist)=split/\:/,$id;
	$n++;
	
	
	##### the alterative annotation from the blast;
	
	my @array=split/\n/,$str{$id}->{string};
	
	my $FMstring;
	if ($#array>=1){
		foreach my $i (1..$#array){
			my ($Mchr,$Midentity,$Mmatch,$Mstart,$Mend)=(split/\t/,$array[$i])[1,2,3,9,10];
		
			my $Mchr2=$alter{$Mchr};
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
	
	
	print OUT "$n\t$id\t$sample\t$string{$id}\t$insertion\t$length\t$ndonor\t$iden\t$cov\t$bitscore\t$pvalue\t$chr2\t$dstart\t$dend\t$strand\t$annotation\t$dist\t$FMstring\t$upjunction\t$downjunction\n";
	
	
	#print OUT "$n\t$sample\t$string{$id}\t$insertion\t$length\t$ndonor\t$identity\t$cov\t$covA\t$covB\t$chr2\t$dstart\t$dend\t$strand\t$annotation\t$dist\t$FMstring\t$upjunction\t$downjunction\n";
	
}

close B;
close OUT;



