#!/usr/bin/perl
use strict;
use warnings;
my $version="1.0 version";
use Getopt::Long;

#Author:Xin Wang 
#email: xin.wang@childrens.harvard.edu
#PI: Kaifu Chen

### This script is to divided all the insertion into different categories:
### requirement: 1
#####1.	We didived them into one donor, two donor, three donor and without donors


my %opts;
GetOptions(\%opts,"i:s","o:s","g:s","c:s","t:s","e:s","m:s","h:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{g} ||!defined $opts{o}||!defined $opts{c} || defined $opts{h}) {
	die "************************************************************************
	Usage: extract_fastq.pl -g All blast results(mutiple pieces) -i Assembled fasta file -c Clustering Files With Coverage -o Blast_output 
	
	### Here divide them into mutiple types, including no-mapping (with MAT), single donor, two donors, three donors and four donors

	
	Request Parameters:
	-g BlastFile for Assembled Fasta
	-i Assembled Fasta file
	-c Clustering File with Coverage. This provide the coverage of each represent insertions generated before.
	-o Output substring of mutiple insertion types
	
	Optional Restriction Parameters:
	-t The maxmium microhomologous length between donors(default 20)
	-e The overlapping region between donor, less than half of min donor size (0.5)
	-m The mininus length of large insertions (default 10bp)
	-h Help
	
************************************************************************\n";
}



my $output=$opts{o};
my $input=$opts{i};
my $index=$opts{g};


#my %conf1; my %inf1; my %int1; my %confInsert1; my %potential1_two; my %potential1_chr; my %check_all1; my %basic1;
 my %num;  my %str;
 my $id;

 ### This index is for the sequence, length index;
my %hash; my %len;  my %sequence;
open I,"$input" or die $!;
while (<I>) {
	chomp;
	if (/^>(\S+)/){
		$id=$1;
		$hash{$id}++;
	}else{
		$sequence{$id}.=$_;
		$len{$id}=length $_;
	}
}
close I;


### This index is for the representive Read counts and reads identity, the clustering results generated from the first round of deduplication 

my $indexR=$opts{c}; 
my %cov; my %qual; my %identity;

open IR,"$indexR" or die $!;
while (<IR>) {
	chomp;
	my ($status,$RCount,$RId,$RQu,$RIden)=(split/\t/,$_)[0,2,3,4,5];
	$cov{$RId}=$RCount;
	$qual{$RId}=$RQu;
	$identity{$RId}=$RIden;
}
close IR;




###
my %mutiple; my %inf;
open G,"$index" or die $!;
open OUT,">$output.overall.txt" or die $!;

while (<G>){
	chomp;
	my ($qid,$pid,$identity,$match,$gap,$mismatch,$qstart,$qend,$qlength,$pstart,$pend,$plength,$score,$evlaue)=split/\s+/,$_;
	next if (!exists $hash{$qid});
	# if (($qstart >55 || ($plength-$pend)>55) && ($pid eq "chrIII" && $pstart >= 13650 && $pend<= 13850) || ($pid eq "chrIII" && $pstart >= 200750 && $pend <= 201000) || ($pid eq "chrIII" && $pstart >= 294300 && $pend <= 294500)){
	# 	print EX "$cov{$qid}\t$_\n";
	# 	next;
	# }
	next if ($pid eq "chrIII" && $pstart >= 13650 && $pend<= 13850);
	next if ($pid eq "chrIII" && $pstart >= 200750 && $pend <= 201000);
	next if ($pid eq "chrIII" && $pstart >= 294300 && $pend <= 294500);
	
	my $information=$_;
	#$length{$qid}=$qlength;
	if (exists $mutiple{$qid}){
		
		my $min=($qlength<$inf{$qid}->{length})?$qlength:$inf{$qid}->{length};
		
		
		### Here we require the microhomology less than 20bp and less than half of insertion fragments
		if ($qstart >= $inf{$qid}->{max} || ( ($inf{$qid}->{max} -$qstart <=20) && ($inf{$qid}->{max} -$qstart <= $min/2) )){
		#if ($qstart >= $inf{$qid}->{max} || ($qstart > $inf{$qid}->{min} && $qstart < $inf{$qid}->{max} && ($qend- $inf{$qid}->{max})>10 )){

			$num{$qid}++;
			#$inf{$qid}->{string}.=$information."\n";
			print OUT "$information\n";

			$inf{$qid}->{max}=$qend;
			$inf{$qid}->{length}=$qlength;

			#push @array,($qstart,$qend);

			#$inf{$qid}->{pre}=$qstart;
			$str{$qid}.=$information."\n";

			my $mapped=join "\t",($pid,$qstart,$qend,$pstart,$pend,$match,$identity);
			$inf{$qid}->{information}.="\t".$mapped;

		}elsif ($qend <= $inf{$qid}->{min} || ( ($qend - $inf{$qid}->{min} <=20) && ($qend - $inf{$qid}->{min} <=$min/2) )){
			#}elsif ($qend <= $inf{$qid}->{min} || ( ($inf{$qid}->{min}-$qstart) >10 && $qend>$inf{$qid}->{min} && $qend <$inf{$qid}->{max})){
			$num{$qid}++;
			#$inf{$qid}->{string}.=$information."\n";
			print OUT "$information\n";
			$inf{$qid}->{min}=$qstart;
			$str{$qid}.=$information."\n";

			my $mapped=join "\t",($pid,$qstart,$qend,$pstart,$pend,$match,$identity);
			$inf{$qid}->{information}.="\t".$mapped;
			$inf{$qid}->{length}=$qlength;


		}

	}else{
		$mutiple{$qid}++;
		$inf{$qid}->{min}=$qstart;
		$inf{$qid}->{max}=$qend;
		$inf{$qid}->{length}=$qlength;
		$num{$qid}++;
		#$inf{$qid}->{string}=$information."\n";
		print OUT "$information\n";
		$str{$qid}.=$information."\n";
		#print OUT "$inf{$qid}->{string}";

		my $mapped=join "\t",($pid,$qstart,$qend,$pstart,$pend,$match,$identity);
		$inf{$qid}->{information}=$mapped;

	}

}

#### create a information for single donor, two donor, four donor using blast results

open SINGLE, ">$output.single.txt" or die $!;
open TWO, ">$output.two.txt" or die $!;
open THREE, ">$output.three.txt" or die $!;
open FOUR, ">$output.four.txt" or die $!;
open EX,">$output.noaligned.fasta" or die $!;
open EXS,">$output.noaligned.txt" or die $!;

##### create a information for two or more donor:
#### format: ID	firstmappedrefchr	firstmappedreadstart	firstmappedreadend	firstmappedrefstart	firstmappedreflength	secondmappedreadstart	secondmappedreadend	secondmappedrefstart	secondmappedreflength
#####
open TWOINF, ">$output.twoinf.txt" or die $!;
open THREEINF, ">$output.threeinf.txt" or die $!;
open FOURINF, ">$output.fourinf.txt" or die $!;

my $fnum=0;
foreach my $i (keys %hash){
	$fnum++;
	if (!exists $num{$i}){
		print EX ">$i\t$cov{$i}\t$qual{$i}\n$sequence{$i}\n";
		print EXS "UnBlast\tUnBlast$fnum\t$cov{$i}\t$i\t$qual{$i}\t$identity{$i}\n";
	}elsif ($num{$i}==1){
		print SINGLE "$cov{$i}\t$str{$i}";
	}elsif($num{$i}==2){
		print TWO "$str{$i}";
		my @array=split/\t/,$inf{$i}->{information};

		my %start2;

		$start2{$array[1]}=join "\t",($array[0],$array[1],$array[2],$array[3],$array[4],$array[5],$array[6]);

		$start2{$array[8]}=join "\t",($array[7],$array[8],$array[9],$array[10],$array[11],$array[12],$array[13]);
		
		

		my $twoinf;
		foreach my $m ( sort {$a<=>$b} keys %start2){
			$twoinf.="\t".$start2{$m};
		}

		#my $twoinf=($startA{start}>$startB{start})?($inf{$i}->{information}):$twooption;
		print TWOINF "$cov{$i}\t$len{$i}\t$i$twoinf\n";

	}elsif($num{$i}==3){
		print THREE "$str{$i}";
		my @array=split/\t/,$inf{$i}->{information};

		my %start3;
		$start3{$array[1]}=join "\t",($array[0],$array[1],$array[2],$array[3],$array[4],$array[5],$array[6]);

		$start3{$array[8]}=join "\t",($array[7],$array[8],$array[9],$array[10],$array[11],$array[12],$array[13]);

		$start3{$array[15]}=join "\t",($array[14],$array[15],$array[16],$array[17],$array[18],$array[19],$array[20]);

		my $threeinf;
		foreach my $m ( sort {$a<=>$b} keys %start3){
			$threeinf.="\t".$start3{$m};
		}



		print THREEINF "$cov{$i}\t$len{$i}\t$i$threeinf\n";
	}elsif($num{$i}==4){
		print FOUR "$str{$i}";

		my @array=split/\t/,$inf{$i}->{information};
		my %start4;


		$start4{$array[1]}=join "\t",($array[0],$array[1],$array[2],$array[3],$array[4],$array[5],$array[6]);

		$start4{$array[8]}=join "\t",($array[7],$array[8],$array[9],$array[10],$array[11],$array[12],$array[13]);

		$start4{$array[15]}=join "\t",($array[14],$array[15],$array[16],$array[17],$array[18],$array[19],$array[20]);
		
		$start4{$array[22]}=join "\t",($array[21],$array[22],$array[23],$array[24],$array[25],$array[26],$array[27]);

		my $fourinf;
		foreach my $m ( sort {$a<=>$b} keys %start4){
			$fourinf.="\t".$start4{$m};
		}
		print FOURINF "$cov{$i}\t$len{$i}\t$i$fourinf\n";
	}

}
close SINGLE;
close TWO;
close THREE;
close FOUR;

close TWOINF;
close THREEINF;
close FOURINF;








