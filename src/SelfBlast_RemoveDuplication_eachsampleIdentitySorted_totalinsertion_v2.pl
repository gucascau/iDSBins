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
GetOptions(\%opts,"f:s","o:s","i:s","b:s","h:s","t:s","g:s","m:s","c:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{f} ||!defined $opts{o} ||!defined $opts{i}||!defined $opts{b}|| defined $opts{h}) {
	die "************************************************************************
	Usage: $0.pl -f Insfasta -i Read with evaluation -b SelfBlast -o Output of Quality Control Reads
	
	Request Parameters:
	-b Blast Results (Self Blast results))
	-i The evaluation of each inserted reads(quality, identity)
	-f Assembled reads (After first round of deduplication base)
	-o The final results strings of files, including the final forward/reverse deduplicated reads, final clustering file with statistical resutls and their representive read
	
	Optional Parameters:
	-t Identity of two reads (default 95)
	-g Gap size (default 6)
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
my $gap=(defined $opts{g})?$opts{g}:6;
my $output=$opts{o};

my %remove; my %con; my %contain; my %hash; my %name; my $n=0;

open BLAST,"$blast" or die $!;
while (<BLAST>){
	chomp;
	my ($id1,$id2,$qidentity,$qmatches,$qmismatches,$qgapsize,$qlenth)=(split/\t/,$_)[0,1,2,3,4,5,8];
	next if ($id1 eq $id2);
	my $av=$qmatches/$qlenth;
	next unless ($av>=$cover && $qidentity >=$ident && $qmismatches <=$mismatches && $qgapsize<=$gap);
	
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
		#print "$name{$id2}\t$name{$id1}\n";
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
### 		asssembled fasta with the large insertion events
######################################################################################################################################

### read the assembled fasta
my $Assembled=$opts{f};
open FASTA,"$Assembled" or die $!;
#
my $Fid;  my %sequence; my %hashF;
while (<FASTA>) {
	chomp;
    # print "$_" ;
	if (/^>(\S+)/){
		$Fid=$1;
		$hashF{$Fid}++;
	}else{
		$sequence{$Fid}.=$_;
		#$len{$Fid}=length $_;
	}

}
close FASTA;


#### Read the statistical file with the identify and quality 
my $rqual=$opts{i}; my %inf;

### this file generated by the (yWH475-1_S10_detected_assemblysep.CoverageWithAssembly.txt)
open RQU,"$rqual" or die $!;
while (<RQU>){
	chomp;
	my ($Rcount,$id,$qual,$iden)=(split/\t/,$_)[2,3,4,5];
	next if ($id eq "ID");
	$inf{$id}->{quality}=$qual;
	$inf{$id}->{cov}=$Rcount;
	$inf{$id}->{identity}=$iden;
}

close RQU;


######################################################################################################################################
## Here we generated the clustering results for further eliminating the duplicates and determine the represent read for each cluster
## Criterion: 1. Ranking the clustering reads with the identiy, quality and Read Count. We used the rank first read to represent each cluster
###
######################################################################################################################################


open CLS,">$output.cls" or die $!;

print CLS "AsemStatus\tInsertionID\tReadNumber\tRepRead\tRepQual\tRepIden\tClsReads\tClsIden\tClsQual\tRankNum\n";

my %uniq; my %repeat;
foreach my $i (sort keys %hash){
	
	my @array=split /\t/,$hash{$i};
	my $clsN=@array;
	print CLS "Assembly\t$i\t";
	my @infQ; my @infI; my %qual=(); my %iden=();my %scoreI=();my %scoreQ=(); my %cov=();
	
	
	### here is to define the quality, identity of the clustering reads	
	foreach my $p (@array){
		
		$qual{$p}=$inf{$p}->{quality};
		$cov{$p}=$inf{$p}->{cov};
		### We set up the no-aligned identity as 0, 
		### which would be better for the identity comparison that we can put them into a pretty low priority.
		### If all the groups cannot have the identity, we then only consider the inserted quality.
		$iden{$p}=($inf{$p}->{identity} eq "NA")?0:$inf{$p}->{identity};
	}

	
	
#### Here we rank the reads first with quality and second with identity. In this case, we obtained the representive result with the best alignment and high quality	

	my @keyF = sort {$cov{$b} <=> $cov{$a} || $iden{$b} <=> $iden{$a} || $qual{$b} <=> $qual{$a} } keys %cov;

	#### push the representive read to hash, and put the number of reads from each cluster into hash value
	$uniq{$keyF[0]}=$clsN; 
	
	# # we sorted the read quality first
# 	my @keyQ = sort { $qual{$b} <=> $qual{$a} or $a cmp $b } keys %qual;

	my ($prev,$prevC, $rankQ);
	for my $k (@keyF) {
	    $rankQ++ unless defined($prev) && $prev==$qual{$k} && $prevC == $cov{$k};
		$scoreQ{$k}=$rankQ;
	    $prev = $qual{$k};
		$prevC =$cov{$k};
	}
	
	
	my $num=0; my $Fstring; my $Fidentity; my $FQual; my $SortID;my $FCounts=0; my $FCov;
	
	#### 
	
	foreach my $m(@keyF){
		$num++;
		$Fstring.="$num;";
		$Fidentity.="$inf{$m}->{identity};";
		$FQual.="$inf{$m}->{quality};";
		$FCov.="$inf{$m}->{cov};";
		$repeat{$m}++;
		$SortID.="$m;";
		$FCounts +=$inf{$m}->{cov};	
	}
	

	print CLS "$FCounts\t$keyF[0]\t$inf{$keyF[0]}->{quality}\t$inf{$keyF[0]}->{identity}\t$SortID\t$Fidentity\t$FQual\t$FCov\t$Fstring\n";
	

	
}
	

######################################################################################################################################
## Here we generated the unique fastq files for large insertion
######################################################################################################################################

open FOUT,">$output.FAssuniq.fasta" or die $!;
#open ROUT,">$output.Funiq.reverse.fastq" or die $!;
my $uncls=0;
foreach my $f (keys %hashF){
	
	if (exists $uniq{$f}){
		print FOUT ">$f\n$sequence{$f}\n";
		
	}elsif(!exists $uniq{$f} && !exists $repeat{$f} ){
		$uncls++;
		my $Fcoverage=$inf{$f}->{cov};
		print FOUT ">$f\n$sequence{$f}\n";
		print CLS "Assembly\tUndefined$uncls\t$Fcoverage\t$f\t$inf{$f}->{quality}\t$inf{$f}->{identity}\tNoClustered\t$inf{$f}->{identity}\t$inf{$f}->{quality}\tNoRank\n";
	}
	
}

close CLS;
close FOUT;


