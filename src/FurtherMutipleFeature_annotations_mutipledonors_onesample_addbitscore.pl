#!/usr/bin/perl
use strict;
use warnings;


#### This scripts is to add addtional features for mutiple blast


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","b:s","o:s","c:s","d:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{b} || !defined $opts{o}  || !defined $opts{c} || !defined $opts{d} ) {
	die "************************************************************************
	Usage: $0.pl -i Mutiple annotation but not contain mutiple feature
				-b Blast results for mutiple donors
				-c Annotation of different elements based on the chr,start
				-d chrosomes change file
				-o Output of multiple annotation with multiple features
************************************************************************\n";
}



my $multiple=$opts{i};

my $blast=$opts{b};

my $anno=$opts{c};
my $alteration=$opts{d};
my $output=$opts{o};



my %alter;

open CH,"$alteration" or die "cannot open file $alteration";
while (<CH>){
	chomp;
	s/\r//g;
	my ($i,$ch)=split /\t/,$_;

	$alter{$i}=$ch;

}
close CH;




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

my %element; my %inf;

while (<E>){
	chomp;
	my ($chra,$starta,$enda,$strand,$iden,$elementa,$dista)=(split/\t/,$_)[0,1,2,5,6,14,15];
	my $stringa=join ":",($elementa,$dista);
	
	my $indexa=join ":",($chra,$starta,$enda);
	
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
	

	
	if ($array[0] =~/\.(\S+)/){
		
		my @inf=split/\n/,$str{$array[0]}->{string};
		
		my $FMstring;
		
		if ($#inf<1){
			$FMstring="NO";
		}else{
			
			foreach my $i (1..$#inf){
			 		my ($Mchr,$Midentity,$Mmatch,$Mstart,$Mend,$score,$evlaue)=(split/\t/,$inf[$i])[1,2,3,9,10,12,13];
			#
			 		my $Mchr2=$alter{$Mchr};
			 		my $Mtype2=($Mstart<$Mend)?"+":"-";
			 		my $Mstart2=($Mstart<$Mend)?$Mstart:$Mend;
			 		my $Mend2=($Mstart<$Mend)?$Mend:$Mstart;
			#
			 		my $com=join ":",($Mchr2,$Mstart2,$Mend2);
			# 		##### if need to add the Melement ### we need to generate the file
			 		my $Melement=$element{$com};
			#
			#
			#
			 		my $Mstring=join":",($Mchr2,$Midentity,$score,$evlaue,$Mmatch,$Mtype2,$Mstart2,$Mend2,$Melement);
			#
			 		$FMstring.=$Mstring.";";
			
				}
		
		}
		my $infstring=join "\t",@array[0..16];
	
		my $infstring3=join "\t",@array[20..24];
		print OUT "$infstring\t$FMstring\tNO\tNO\t$infstring3\n";	
		
		
	}else{
		print OUT "$_\n";
	}
	
}

close I;
close OUT;



