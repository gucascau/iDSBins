#!/usr/bin/perl
use strict;
use warnings;


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s","g:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{o} || !defined $opts{g}) {
	die "************************************************************************
	Usage: extract_fastq.pl  -g All blast results(mutiple pieces) 
				-i input of unique assembled fasta (with coverage information) 
				-o Blast_output (Here divide them into mutiple types, including no-mapping (with MAT), single donor, two donors, three donors and four donors)
************************************************************************\n";
}

my $output=$opts{o};
my $input=$opts{i};
my $index=$opts{g};


#my %conf1; my %inf1; my %int1; my %confInsert1; my %potential1_two; my %potential1_chr; my %check_all1; my %basic1;
 my %num;  my %str;
 my $id;
my %hash; my %len; my %cov; my %sequence;
open I,"$input" or die $!;
while (<I>) {
	chomp;
	if (/^>(\S+)/){
		$id=$1;
		#my $coverage=(split/\|/,$1)[];
		#$cov{$id}=$3;
		$hash{$id}++;
	}else{
		$sequence{$id}.=$_;
		$len{$id}=length $_;
	}
}
close I;

###
my %mutiple; my %inf;
open G,"$index" or die $!;
open OUT,">$output.overall.txt" or die $!;

while (<G>){
	chomp;
	my ($qid,$pid,$identity,$match,$gap,$mismatch,$qstart,$qend,$qlength,$pstart,$pend,$plength,$score,$evlaue)=split/\s+/,$_;
	
	if (!exists $hash{$qid}){
#		print  ">$qid\n$sequence{$qid}\n";
		next;
	}
	# if (($qstart >55 || ($plength-$pend)>55) && ($pid eq "ref-NC_001135-" && $pstart >= 13650 && $pend<= 13850) || ($pid eq "ref-NC_001135-" && $pstart >= 200750 && $pend <= 201000) || ($pid eq "ref-NC_001135-" && $pstart >= 294300 && $pend <= 294500)){
	# 	print EX "$cov{$qid}\t$_\n";
	# 	next;
	# }
	next if ($pid eq "ref-NC_001135-" && $pstart >= 13650 && $pend<= 13850);
	next if ($pid eq "ref-NC_001135-" && $pstart >= 200750 && $pend <= 201000);
	
	
	#next if ($pid eq "ref-NC_001135-" && $pstart >= 294300 && $pend <= 294500);
	next if (($qstart<35 || abs($qlength-$qend)<20) && $match <60);
	my $information=$_;
	
	#$length{$qid}=$qlength;
	if (exists $mutiple{$qid}){
		
		my $min=($qlength<$inf{$qid}->{length})?$qlength:$inf{$qid}->{length};
		
		if ($qstart >= $inf{$qid}->{max} || ( ($inf{$qid}->{max} -$qstart <=20) && ($inf{$qid}->{max} -$qstart <= $min/2) )){

			$num{$qid}++;
			#$inf{$qid}->{string}.=$information."\n";
			print OUT "$information\n";

			$inf{$qid}->{max}=$qend;
			$inf{$qid}->{length}=$qlength;
			#push @array,($qstart,$qend);

			#$inf{$qid}->{pre}=$qstart;
			$str{$qid}.=$information."\n";

			my $mapped=join "\t",($pid,$qstart,$qend,$pstart,$pend,$match);
			$inf{$qid}->{information}.="\t".$mapped;


		}elsif ($qend <= $inf{$qid}->{min} || ( ($qend - $inf{$qid}->{min} <=20) && ($qend - $inf{$qid}->{min} <=$min/2) )){
			$num{$qid}++;
			#$inf{$qid}->{string}.=$information."\n";
			print OUT "$information\n";
			$inf{$qid}->{min}=$qstart;
			$str{$qid}.=$information."\n";
			$inf{$qid}->{length}=$qlength;

			my $mapped=join "\t",($pid,$qstart,$qend,$pstart,$pend,$match);
			$inf{$qid}->{information}.="\t".$mapped;

		}

	}else{
		$mutiple{$qid}++;
		$inf{$qid}->{min}=$qstart;
		$inf{$qid}->{max}=$qend;
		$num{$qid}++;
		#$inf{$qid}->{string}=$information."\n";
		print OUT "$information\n";
		$str{$qid}.=$information."\n";
		#print OUT "$inf{$qid}->{string}";

		my $mapped=join "\t",($pid,$qstart,$qend,$pstart,$pend,$match);
		$inf{$qid}->{information}=$mapped;

	}

}

#### create a information for single donor, two donor, four donor using blast results

open SINGLE, ">$output.single.txt" or die $!;
open TWO, ">$output.two.txt" or die $!;
open THREE, ">$output.three.txt" or die $!;
open FOUR, ">$output.four.txt" or die $!;
open EX,">$output.mat.fasta" or die $!;


##### create a information for two or more donor:
#### format: ID	firstmappedrefchr	firstmappedreadstart	firstmappedreadend	firstmappedrefstart	firstmappedreflength	secondmappedreadstart	secondmappedreadend	secondmappedrefstart	secondmappedreflength
#####
open TWOINF, ">$output.twoinf.txt" or die $!;
open THREEINF, ">$output.threeinf.txt" or die $!;
open FOURINF, ">$output.fourinf.txt" or die $!;


foreach my $i (keys %hash){
	if (!exists $num{$i}){
		print EX ">$i\n$sequence{$i}\n";
	}elsif ($num{$i}==1){
		print SINGLE "$str{$i}";
	}elsif($num{$i}==2){
		print TWO "$str{$i}";
		my @array=split/\t/,$inf{$i}->{information};

		my %start2;

		$start2{$array[1]}=join "\t",($array[0],$array[1],$array[2],$array[3],$array[4],$array[5]);

		$start2{$array[7]}=join "\t",($array[6],$array[7],$array[8],$array[9],$array[10],$array[11]);


		my $twoinf;
		foreach my $i ( sort {$a<=>$b} keys %start2){
			$twoinf.="\t".$start2{$i};
		}

		#my $twoinf=($startA{start}>$startB{start})?($inf{$i}->{information}):$twooption;
		print TWOINF "$len{$i}\t$i$twoinf\n";

	}elsif($num{$i}==3){
		print THREE "$str{$i}";
		my @array=split/\t/,$inf{$i}->{information};

		my %start3;
		$start3{$array[1]}=join "\t",($array[0],$array[1],$array[2],$array[3],$array[4],$array[5]);

		$start3{$array[7]}=join "\t",($array[6],$array[7],$array[8],$array[9],$array[10],$array[11]);


		$start3{$array[13]}=join "\t",($array[12],$array[13],$array[14],$array[15],$array[16],$array[17]);

		my $threeinf;
		foreach my $i ( sort {$a<=>$b} keys %start3){
			$threeinf.="\t".$start3{$i};
		}



		print THREEINF "$len{$i}\t$i$threeinf\n";
	}elsif($num{$i}==4){
		print FOUR "$str{$i}";

		my @array=split/\t/,$inf{$i}->{information};
		my %start4;


		$start4{$array[1]}=join "\t",($array[0],$array[1],$array[2],$array[3],$array[4],$array[5]);

		$start4{$array[7]}=join "\t",($array[6],$array[7],$array[8],$array[9],$array[10],$array[11]);

		$start4{$array[13]}=join "\t",($array[12],$array[13],$array[14],$array[15],$array[16],$array[17]);

		$start4{$array[19]}=join "\t",($array[18],$array[19],$array[20],$array[21],$array[22],$array[23]);

		my $fourinf;
		foreach my $i ( sort {$a<=>$b} keys %start4){
			$fourinf.="\t".$start4{$i};
		}
		print FOURINF "$len{$i}\t$i$fourinf\n";
	}

}
close SINGLE;
close TWO;
close THREE;
close FOUR;

close TWOINF;
close THREEINF;
close FOURINF;








