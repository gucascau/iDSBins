#!/usr/bin/perl
use strict;
use warnings;

####
### This script is to extract Uniq blast results (Single donors) and statistic their information.
### requirement: 1
#####1

my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s","a:s","b:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{o}|| !defined $opts{a} || !defined $opts{b}   ) {
	die "************************************************************************
	Usage: $0.pl -i Two_insertion.txt -a gapsize of reads -b gapsize of mapping regions -o Blast_statatistic
************************************************************************\n";
}


### write a report for first round of insertion files
### number of reads in different chromsomes including duplciates reads
### number of uniq events in different chromosome using the uniq reads

### inserted length distribution
### uniq inserted length distribution
#
#
#
# ##### insertion #####
#12	500	 M06255:3:000000000-CGCML:1:2106:11859:4639	ref-NC_001144-	187	223	454078	454114	37	ref-NC_001144-	225	385	455024	454864	161
# 1	522	M06255:3:000000000-CGCML:1:2108:7400:6298	ref-NC_001148-	40	361	806767	807088	322	ref-NC_001148-	362	481	267805	267686	120
# 1	523	M06255:3:000000000-CGCML:1:1117:29285:13920	ref-NC_001148-	40	361	806767	807088	322	ref-NC_001148-	362	481	267805	267686	120
# 1	524	M06255:3:000000000-CGCML:1:1104:24754:8244	ref-NC_001148-	40	361	806767	807088	322	ref-NC_001148-	362	481	267805	267686	120
# 1	528	M06255:3:000000000-CGCML:1:1118:4168:15800	ref-NC_001148-	39	360	806767	807088	322	ref-NC_001148-	361	480	267805	267686	120
# 10224	529	M06255:3:000000000-CGCML:1:2107:24431:9059	ref-NC_001148-	40	361	806767	807088	322	ref-NC_001148-	362	481	267805	267686	120
# 2	519	M06255:3:000000000-CGCML:1:2106:9996:3535	ref-NC_001148-	40	361	806767	807088	322	ref-NC_001148-	362	481	267805	267686	120


my $input=$opts{i};
my $output=$opts{o};

my $gapA=$opts{a};
my $gapB=$opts{b};

my %num; my %match; my %hash;
open I,"$input" or die $!;
my $n=1;
while (<I>){
	chomp;
	s/\r//g;
	my ($cov,$len,$qid,$p1id,$p1qstart,$p1qend,$p1pstart,$p1pend,$p1match,$p2id,$p2qstart,$p2qend,$p2pstart,$p2pend,$p2match)=split/\s+/,$_;

	#$Fqstart,$Rqstart,$Fqend,$Rqend,$Fpstart,$Rpstart,$Fpend,$Rpend,$Flength,$Rlength)=split/\s+/,$_;
	my $information=$_;
	### split into two categories: one had no other insertions
	###  							another had
	next if ($qid eq "qID");
	#next if (exists $hash{$qid});
	# $pstart=($pstart<$pend)?$pstart:$pend;
# 	$pend=($pstart<$pend)?$pend:$pstart;
## my $nstart=$Fqstart;
# 	my $nend=$Rqend;
	my $m=0;
	my $q=$n-1;
	foreach my $i (1..$q){
		if ($p1id eq $hash{$i}->{p1id} && (abs($p1qstart-$hash{$i}->{p1qstart}) <=$gapA) && (abs($p1qend-$hash{$i}->{p1qend}) <=$gapA)  && (abs($p1pstart-$hash{$i}->{p1pstart}) <=$gapB) && (abs($p1pend-$hash{$i}->{p1pend}) <=$gapB) && $p2id eq $hash{$i}->{p2id} && (abs($p2qstart-$hash{$i}->{p2qstart}) <=$gapA) && (abs($p2qend-$hash{$i}->{p2qend}) <=$gapA)  && (abs($p2pstart-$hash{$i}->{p2pstart}) <=$gapB) && (abs($p2pend-$hash{$i}->{p2pend}) <=$gapB) ){

		#if ($pid eq $hash{$i}->{pid} && (abs($qstart-$hash{$i}->{qstart}) <=$gapA) && (abs($qend-$hash{$i}->{qend}) <=$gapA)  && (abs($pstart-$hash{$i}->{pstart}) <=$gapB) && (abs($pend-$hash{$i}->{pend}) <=$gapB) )

			$hash{$i}->{num} += $cov;
			$hash{$i}->{string}.="\t".$qid;
			$m++;

			$hash{$i}->{inf}=$information if ($cov > $hash{$i}->{cov});
			
			$hash{$i}->{cov}=$cov if ($cov > $hash{$i}->{cov});
			print "$qid\t$m\n";
		}

	}

	if ($m==0){

		#print "$n\n";
		$hash{$n}->{cov}=$cov;
		#### set the nonduplicated index ####
		### first block
		$hash{$n}->{p1id}=$p1id;
		$hash{$n}->{p1qstart}=$p1qstart;
		$hash{$n}->{p1qend}=$p1qend;
		$hash{$n}->{p1pstart}=$p1pstart;
		$hash{$n}->{p1pend}=$p1pend;

		### second block
		$hash{$n}->{p2id}=$p2id;
		$hash{$n}->{p2qstart}=$p2qstart;
		$hash{$n}->{p2qend}=$p2qend;
		$hash{$n}->{p2pstart}=$p2pstart;
		$hash{$n}->{p2pend}=$p2pend;

		### push the
		$hash{$n}->{inf}=$information;

		$hash{$n}->{num}=$cov;
		$hash{$n}->{string}=$qid;
		$n++;

	}
}

close I;

####
open OUT,">$output.uniq.txt" or die $!;
open CLT,">$output.cls.txt" or die $!;
foreach my $num(keys %hash){
	print  OUT "$hash{$num}->{inf}\t$hash{$num}->{num}\n";
	print CLT "$hash{$num}->{num}\t$hash{$num}->{string}\n";
}

close OUT;
close CLT;