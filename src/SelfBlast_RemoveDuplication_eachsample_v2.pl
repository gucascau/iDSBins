#!/usr/bin/perl
use strict;
use warnings;
#use Data::Dump qw(dump);
#Author:Xin Wang 
#email: xin.wang@childrens.harvard.edu
#PI: Kaifu Chen

### function: This script is to extract the potential elimanated insertions and added the removed insertion coverage into previous insertions
### The script is based on the previous self treatment 




my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"f:s","r:s","o:s","i:s","b:s","h:s","t:s","g:s","m:s","c:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{f} || !defined $opts{r} ||!defined $opts{o} ||!defined $opts{i}||!defined $opts{b}|| defined $opts{h}) {
	die "************************************************************************
	Usage: $0.pl -f Ffastq -r Rfastq -i Read with evaluation -b SelfBlast -o Output of Quality Control Reads
	
	Request Parameters:
	-b Blast Results (Self Blast results))
	-i The evaluation of each inserted reads(quality, identity)
	-f Forward reads
	-r Reverse reads
	-o The final results strings of files, including the final forward/reverse deduplicated reads, final clustering file with statistical resutls and their representive read
	
	Optional Parameters:
	-t Identity of two reads (default 95)
	-g Gap size (default 2)
	-m mismatches (default 6)
	-c Coverage (Matched size/Full length, default 95%)
	-h Help
************************************************************************\n";
}



########################################################################################
#### The first step to generate a unique cluster results from the self blast results

#### Criterion : Identity more than 95%, less than 6 mismatches and 2 gapsize, two reads coverage more than 95%
########################################################################################

my $blast=$opts{b};

my $ident=(defined $opts{t})?$opts{t}:95;
my $cover=(defined $opts{c})?$opts{c}:0.95;
my $mismatches = (defined $opts{m})?$opts{m}:6;
my $gap=(defined $opts{g})?$opts{g}:2;
my $output=$opts{o};

my %remove; my %con; my %contain; my %hash; my %name; my $n=0;

open BLAST,"$blast" or die $!;
while (<BLAST>){
	chomp;
	my ($id1,$id2,$qidentity,$qmatches,$qmismatches,$qgapsize,$qlenth)=(split/\t/,$_)[0,1,2,3,4,5,8];
	next if ($id1 eq $id2);
	my $av=$qmatches/$qlenth;
	next unless ($av>=$cover && $qidentity >=$ident && $qmismatches <=6 && $qgapsize<=2);
	
	if (!exists $name{$id1} && !exists $name{$id2}){
		$n++;
		#push @{$name{$id1}},($id1,$id2);
		my $string=join "\t",($id1,$id2);
		$name{$id1}=$n;
		$name{$id2}=$n;
		$hash{$n}=$string;
	}elsif (exists $name{$id1} && !exists $name{$id2} ){
		my $num=$name{$id1};
		#push @{$name{$id1}},$id2;
		$name{$id2}=$num;
		$hash{$num}.="\t$id2";
	}elsif (exists $name{$id2} && !exists $name{$id1} ){
		my $num=$name{$id2};
		#push @{$name{$id2}},$id1;
		$hash{$num}.="\t$id1";
		$name{$id1}=$num;
	}elsif (exists $name{$id2} && exists $name{$id1} && $name{$id2} != $name{$id1}){
		print "$name{$id2}\t$name{$id1}\n";
		#identify the min value and, add the keys from max value to min value and minumium keys,  delete the max hash and value
		my ($max,$min)=($name{$id2} > $name{$id1})?($name{$id2},$name{$id1}):($name{$id1},$name{$id2});
		my @arrayR=split/\t/,$hash{$max};
		foreach my $q (@arrayR){
			$name{$q}=$min;
			$hash{$min}.="\t$q";
		}
		delete $hash{$max};
	}	
	
	
}


######################################################################################################################################
### make the index files: 
### 				Forward reads, Reverse reads, Quality file with the large insertion events
######################################################################################################################################

### read forward reads and reverse reads
my $Fread=$opts{f};
open FASTQ1,"$Fread" or die $!;

my $Fid;  my %hash1;  
while (<FASTQ1>) {
	chomp;
    # print "$_" ;

    if ($. % 4 == 1)  {
		$Fid=(split/\s+/,$_)[0];
		$Fid=~s/^@//;
		#print "$id1\n";
		$hash1{$Fid}=$_."\n";
	}else{
		$hash1{$Fid}.=$_."\n";
	}
	
}
close FASTQ1;

my $Rread=$opts{r};

open FASTQ2,"$Rread" or die $!;

my $Rid;  my %hash2; 
while (<FASTQ2>) {
	chomp;
    # print "$_" ;

    if ($. % 4 == 1)  {
		$Rid=(split/\s+/,$_)[0];
		$Rid=~s/^@//;
		#print "$id1\n";
		$hash2{$Rid}=$_."\n";
	}else{
		$hash2{$Rid}.=$_."\n";
	}
	
}
close FASTQ2;

#### Read the statistical file with the identify and quality 
my $rqual=$opts{i}; my %inf;

open RQU,"$rqual" or die $!;
while (<RQU>){
	chomp;
	my ($id,$qual,$ide)=(split/\t/,$_)[0,1,2];
	next if ($id eq "ID");
	$inf{$id}->{quality}=$qual;
	$inf{$id}->{identity}=$ide;
}

close RQU;


######################################################################################################################################
## Here we generated the clustering results for further eliminating the duplicates and determine the represent read for each cluster
## Criterion: 1. Ranking the clustering reads with the idendenty and quality. We used the rank first read to represent eah cluster
###
######################################################################################################################################


open CLS,">$output.cls" or die $!;

my %uniq; my %repeat;
foreach my $i (sort keys %hash){
	
	
	my @array=split /\t/,$hash{$i};
	my $clsN=@array;
	print CLS "$i\t$clsN\t";
	my @infQ; my @infI; my %qual=(); my %iden=();my %scoreI=();my %scoreQ=();
	foreach my $p (@array){
		$qual{$p}=$inf{$p}->{quality};
		$iden{$p}=$inf{$p}->{identity};
	}

	
	my @keyQ = sort { $qual{$b} <=> $qual{$a} or $a cmp $b } keys %qual;

	my ($prev,$rankQ);
	for my $k (@keyQ) {
	    $rankQ++ unless defined($prev) && $prev==$qual{$k};
		$scoreQ{$k}=$rankQ;
	    $prev = $qual{$k};
	}

	my @keyI = sort { $iden{$b} <=> $iden{$a} or $a cmp $b } keys %iden;

	my ($prevI,$rankI);
	for my $j (@keyI) {
	    $rankI++ unless defined($prevI) && $prevI==$iden{$j};
		$scoreI{$j}=$rankI;
	    $prev = $iden{$j};
	}

	### rank the clustering read using the ranking scores from identify, 
	### if you would like to combine the quality and identiy, please use the ranking combination
	my %finalscore=();
	foreach my $q (@array){
		### only consider the Identity ranking
		$finalscore{$q}=$scoreI{$q};
		
		### consider both identity and quality
		#$finalscore{$q}=$scoreQ{$q}+$scoreI{$q};
	}
	### print the final clustering according to the ranking scores, and push all the clustering reads into hash
	my $Fstring;
	
	my @keyF= sort { $finalscore{$a} <=> $finalscore{$b} } keys %finalscore;
	
	#### push the representive read to hash, and put the number of reads from each cluster into hash value
	$uniq{$keyF[0]}=$clsN;
	foreach my $m(sort { $finalscore{$a} <=> $finalscore{$b} } keys %finalscore){
		$Fstring.="$finalscore{$m}; ";
		$repeat{$m}++;
		print CLS "\t$m";
		
	}
	
	print CLS "\t$Fstring\n";
	

	
}
	

######################################################################################################################################
## Here we generated the unique fastq files for large insertion
######################################################################################################################################

open FOUT,">$output.Funiq.forward.fastq" or die $!;
open ROUT,">$output.Funiq.reverse.fastq" or die $!;
my $uncls=0;
foreach my $f (keys %inf){
	
	if (exists $uniq{$f}){
		print FOUT "$hash1{$f}";
		print ROUT "$hash2{$f}";
		
	}elsif(!exists $uniq{$f} && !exists $repeat{$f} ){
		$uncls++;
		print FOUT "$hash1{$f}";
		print ROUT "$hash2{$f}";
		print CLS "Undefined$uncls\t1\t$f\n";
	}
	
}

close CLS;
close FOUT;
close ROUT;

