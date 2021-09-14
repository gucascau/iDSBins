#!/usr/bin/perl
use strict;
use warnings;


#### We set up stringent criterion for the estimated insertions: 
#### 1. There are no other extra sequence at the end of R1 and the begining of R2. This means that the completely aligened to the similar locus of the genome
#### 2. There are no alternative feature of these two elements. For example, they consistently came from Ty region, this would result in the inaccuracy of the large insertion events.
#### This scripts is to add addtional features for estimated insertion. Here we igonred the insertions that mapped against mutiple regions or mitochondrion genome.


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","b:s","o:s","c:s","h:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{b} || !defined $opts{o}  || !defined $opts{c} || defined $opts{h} ) {
	die "************************************************************************
	Usage: $0.pl 
			-i Estimated annotation but not contain alternative feature
			-b Blast results for estimated donors
			-c Annotation of different elements based on the chr,start,end
			-o Output of Estimated annotation with alternative features
			-h help
************************************************************************\n";
}



my $multiple=$opts{i};

my $blast=$opts{b};

my $anno=$opts{c};
my $output=$opts{o};


open B,"$blast" or die "cannot open file $blast";

my %str; my %hash;

while (<B>) {
	chomp;
    # print "$_" ;
	my ($qid,$pid,$identity,$match,$gap,$mismatch,$qstart,$qend,$qlength,$pstart,$pend,$plength,$score,$evlaue)=split/\s+/,$_;
	my $information=$_;
	#$basic1{$qid}=$information;
	#next if ($qstart <5 || ($qlength-$qend) <5);
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
close B;

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


####### this file is corresponding annnotation of each chr,start,end
open I,"$multiple" or die $!;
open OUT,">$output" or die $!;
#my %element; my %inf;

while (<I>){
	chomp;
	my @array=split /\t/,$_;
	my ($id,$chr,$dstart,$dend)=@array[0,10,11,12];

	### add the final feature of this line:
	my $inse=join ":",($chr,$dstart,$dend);
	### add the annotation, distance feature for this line:
	my ($chrA,$startA,$endA,$annotation,$dist)=split ":",$element{$inse};
	### add the appproximate feature
	my $feature;
	if ($dist != 0){
		$feature="PROXIMITY";
	}elsif($dstart>=$startA && $dend<=$endA){
		$feature="ENTIRE";
	}else{
		$feature="PARTIALLY";
	}
	
	### check and determine whether they do have alterative feature:
	my @inf=split/\n/,$str{$id}->{string};
	
			
	### Add the alterative feature ($FMstring)
	my $FMstring;

		
	### No mutiple feature we just export the orignal feature
	if ($#inf<1){
		$FMstring="NO";
	}else{
		foreach my $i (1..$#inf){
			my ($Mchr2,$Midentity,$Mmatch,$Mstart,$Mend,$score,$evlaue)=(split/\t/,$inf[$i])[1,2,3,9,10,12,13];

			my $Mtype2=($Mstart<$Mend)?"+":"-";
			my $Mstart2=($Mstart<$Mend)?$Mstart:$Mend;
			my $Mend2=($Mstart<$Mend)?$Mend:$Mstart;
			#
			my $com=join ":",($Mchr2,$Mstart2,$Mend2);
			# 		if need to add the Melement ### we need to generate the file
			my $Melement=$element{$com};

			my $Mstring=join":",($Mchr2,$Midentity,$score,$evlaue,$Mmatch,$Mtype2,$Mstart2,$Mend2,$Melement);
			#
			 $FMstring.=$Mstring.";";
			
		}
		
	}
	
	## We ignore the insertion that have the alterative features 	
	my $infstring1=join "\t",@array[0..13];
	my $infstring2=join "\t",($annotation,$dist,$feature,$FMstring);
	my $infstring3=join "\t",@array[18..19];
	print OUT "$infstring1\t$infstring2\t$infstring3\n" if ($FMstring eq "NO" || $chr eq "chrmt");	

}

close I;
close OUT;



