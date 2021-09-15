#!/usr/bin/perl
use strict;
use warnings;

####
### This script is to extract Uniq blast results (two donors) and statistic their information.
### requirement: 1
#####1


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s","a:s","b:s","q:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{o} || !defined $opts{q} ) {
	die "************************************************************************
	Usage: $0.pl -i Single_insertion.txt -q quality file of each reads -o Unique insertions output string
	
	Request Parameters for two donor deduplication:
	-i Insertion with detected two donor
	-q The quality information of each read
	-o Output of unique insertion that only detected with two donor, this including the clusting file and unique insertion event file
	
	Optional Parameters:
	-a The difference allowed between two reads, start site, end start site on the inserted reads (default: 6)
	-b The difference allowed between two reads, start site, end start site on the chromosome locus of donor (default: 6)
	-h Help
	
	
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

my $gapA=(defined $opts{a})?$opts{a}:6;
my $gapB=(defined $opts{b})?$opts{b}:6;



### read the quality file
my $Quality=$opts{q};

my %qual;
open QUAL,"$Quality" or die $!;
while (<QUAL>){
	chomp;
	my ($id,$qu)=(split/\t/,$_)[0,5];
	$qual{$id}=$qu;
}
close QUAL;


my %num; my %match; my %hash; my %fiden; my %finfo; my %Rcount;  my %inf;
open I,"$input" or die $!;
my $n=1;
while (<I>){
	chomp;
	s/\r//g;
	my ($cov,$len,$qid,$p1id,$p1qstart,$p1qend,$p1pstart,$p1pend,$p1match,$p1ident,$p2id,$p2qstart,$p2qend,$p2pstart,$p2pend,$p2match,$p2ident)=split/\s+/,$_;
	$fiden{$qid}=($p1ident+$p2ident)/2;

	#$Fqstart,$Rqstart,$Fqend,$Rqend,$Fpstart,$Rpstart,$Fpend,$Rpend,$Flength,$Rlength)=split/\s+/,$_;
	my $information=$_;
	$Rcount{$qid}=$cov;

	#$Fqstart,$Rqstart,$Fqend,$Rqend,$Fpstart,$Rpstart,$Fpend,$Rpend,$Flength,$Rlength)=split/\s+/,$_;
	$inf{$qid}=$information;
	
	### split into two categories: one had no other insertions
	###  							another had
	next if ($qid eq "qID");
	$finfo{$qid}=$_;
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
			#print "$qid\t$m\n";
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
		$hash{$n}->{qual}=$qual{$qid};
		$n++;

	}
}

close I;

####
open OUT,">$output.uniq.txt" or die $!;
open CLT,">$output.cls.txt" or die $!;

my %Fnum;
foreach my $num(keys %hash){	
	my @array=split/\t/,$hash{$num}->{string};
	
	### Here we also considered the quality. If the coverage the same, we only considered the high quality one
	### identify the representive reads
	my $fqual; my $fiden; my $fcov; my $fid; my $FRcount;
	
	my %Lqual=(); my %Lcov=(); my %Liden=();
	foreach my $i (@array){
		$Lqual{$i}=$qual{$i};
		$Lcov{$i}=$Rcount{$i};
		$Liden{$i}=$fiden{$i};
		
		$FRcount += $Rcount{$i};
		$fid.=$i.";";
		$fcov.=$Rcount{$i}.";";
		$fiden.=$fiden{$i}.";";
		$fqual.= $qual{$i}.";";
	}
	
	my @keyF = sort {$Lcov{$b} <=> $Lcov{$a} || $Liden{$b} <=> $Liden{$a} || $Lqual{$b} <=> $Lqual{$a} } keys %Lcov;
	
	my $RepRead=$keyF[0];
	next if (exists $Fnum{$RepRead});
	$Fnum{$RepRead}++;
	
	print  OUT "$FRcount\t$RepRead\t$inf{$RepRead}\n";

	print CLT "SecondFinalTwo\tTID$num\t$FRcount\t$RepRead\t$qual{$RepRead}\t$fiden{$RepRead}\t$fid\t$fqual\t$fiden\t$fcov\n";
		
}

close OUT;
close CLT;