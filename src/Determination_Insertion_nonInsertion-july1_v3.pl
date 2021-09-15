#!/usr/bin/perl
use strict;
use warnings;
my $version="1.0 version";
use Getopt::Long;

#Author:Xin Wang 
#email: xin.wang@childrens.harvard.edu
#PI: Kaifu Chen

### This script is to extract the reads that potenital contain large insertion events:
### requirement: 1
#####1.	We considered the reads without insertion if they completely blasted against reference genome from start to end of the reads, ignoring whether they have small indels. 
#####    And as long as one read span the whole MATA region, we considered these reads without insertion events. 
#####2.	Eliminating the reads with non-insertion, We considered the read as insertion event following the below criterion:
#####	2.1 If either of the reads can mapped against other region, we confidently considered them as large insertion
#####	2.2 If there is no other alignment, We requsted the shortest insertion larger than 10bp (This size can also be adjustable).


my %opts;
GetOptions(\%opts,"i:s","o:s","g:s","c:s","t:s","e:s","m:s","f:s","r:s","q:s","j:s","n:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{g} ||!defined $opts{f}||!defined $opts{r} ||!defined $opts{q} ||!defined $opts{o}||!defined $opts{j}) {
	die "************************************************************************
	Usage: extract_fastq.pl -i R1_Blast -g R2_Blast -f R1_fastq -r R2_fastq  -o Blast_categories
	
	Request Parameters:
	-i Forward Blast file
	-g Reverse Blast file
	-f Forward fastq file
	-r Reverse fastq file
	-q Quality file generated before
	-j Add the name of categories (First, SecondAssembly, or SecondUnAssembly) based on each step. This result will generated for the final detected inserted reads.
	-o Output substring of quality control Reads
	
	Optional Parameters:
	-c The chromosme of MATA reference position (chrIII:294345-294428, default chrIII)
	-t The start site of MATA reference position (chrIII:294345-294428, default 294300)
	-e The end site of MATA reference position (chrIII:294345-294428, default 294500)
	-m The mininus length of large insertions (default 10bp)
	-n The size of MATA region (default 90bp)
	-h Help
	
************************************************************************\n";
}

my $output=$opts{o};
my $input=$opts{i};
my $index=$opts{g};
my $string=$opts{j};


my $chr=(defined $opts{c})?$opts{c}:"chrIII";
my $Rstart= (defined $opts{t})?$opts{t}:294300;
my $Rend= (defined $opts{e})?$opts{e}:294500;
my $Isize= (defined $opts{m})?$opts{m}:10;
my $Msize= (defined $opts{n})?$opts{n}:90;

my $MinRegionSize=$Isize+$Msize;

#my $yeast=$opts{m};

my $Fread=$opts{f};
open FASTQ1,"$Fread" or die $!;

my $id1;  my %hash1;  my %len1;
while (<FASTQ1>) {
	chomp;
    # print "$_" ;

    if ($. % 4 == 1)  {
		$id1=(split/\s+/,$_)[0];
		$id1=~s/^@//;
		#print "$id1\n";
		$hash1{$id1}=$_."\n";
	}elsif($.%4 ==2){
			$hash1{$id1}.=$_."\n";
			$len1{$id1}=length $_;
	}else{
		$hash1{$id1}.=$_."\n";
	}
	
}
close FASTQ1;

my $Rread=$opts{r};

open FASTQ2,"$Rread" or die $!;

my $id2;  my %hash2; my %len2;
while (<FASTQ2>) {
	chomp;
    # print "$_" ;

    if ($. % 4 == 1)  {
		$id2=(split/\s+/,$_)[0];
		$id2=~s/^@//;
		#print "$id1\n";
		$hash2{$id2}=$_."\n";
	}elsif($.%4 ==2){
		$hash2{$id2}.=$_."\n";
		$len2{$id2}=length $_;
	}else{
		$hash2{$id2}.=$_."\n";
	}
	
}
close FASTQ2;

my $Quality=$opts{q};

my %Rqual;

open QUAL,"$Quality" or die $!;
while (<QUAL>){
	chomp;
	my ($id,$qu)=(split/\t/,$_)[0,5];
	$Rqual{$id}=$qu;
}
close QUAL;


my %conf1; my %inf1; my %int1; my %confInsert1; my %potential1_two; my %potential1_chr; my %check_all1; my %basic1; my %unclassfied1;

open I,"$input" or die "cannot open file $input";
while (<I>) {
	chomp;
    # print "$_" ;
	my ($qid,$pid,$identity,$match,$gap,$mismatch,$qstart,$qend,$qlength,$pstart,$pend,$plength,$score,$evlaue)=split/\s+/,$_;
	my $information=$_;

	
	#### Here due to mutiple alignments of MAT right side HMR, HML, we ignore the homologous sequences.
	next if ($pid eq "chrIII" && $pstart >= 13650 && $pend<= 13850);
	next if ($pid eq "chrIII" && $pstart >= 200750 && $pend <= 201000);
	
	$basic1{$qid}=$information;
	$inf1{$qid}->{length}=$qlength;
	$inf1{$qid}->{Identity}=$identity;
	$inf1{$qid}->{Matches}=$match;
	### Usually the best matches were identified from the first one, so if we could identifed them as non-inserted, potential insertion, or confident insertion, we will ignore the following blast results.
	
	next if (exists $conf1{$qid});
	next if (exists $confInsert1{$qid});
	next if (exists $potential1_chr{$qid});
	
	my $cov=$match/$qlength;
	#print "$cov\n";
	

	#### Confident insertions from other chromosome:
	if ($pid ne "chrIII"){
		$confInsert1{$qid}=$information;	
		
	###### confident non-insertion, might contain some deletion events:
	}elsif($pstart>=$Rstart  && $pstart <=$Rend && $pend>=$Rstart && $pend<=$Rend && $match>=70){
		$conf1{$qid}=$information;
		
	##### Confident insertion from other location of chr3:
	}elsif($pstart>=$Rend  || $pstart<=$Rstart){
		$confInsert1{$qid}=$information;
		
		
	#### potential split the fragment with insertion (telemere repetitive elements or MAT elements), as long as they are longer than 100bp
	}elsif($qlength >$MinRegionSize ){
		$potential1_two{$qid}=$information;

	}else{
		$unclassfied1{$qid}=$information;
		
	}

}
close I;


my %conf2; my %inf2;my %int2; my %confInsert2;my %potential2_two; my %potential2_chr; my %check_all2; my %basic2; my %unclassfied2;
open G,"$index" or die "cannot open file $index";
while (<G>) {
	chomp;
    # print "$_" ;
	my ($qid,$pid,$identity,$match,$gap,$mismatch,$qstart,$qend,$qlength,$pstart,$pend,$plength,$score,$evlaue)=split/\s+/,$_;
	my $information=$_;

	
	#### Here due to mutiple alignments of MAT right side HMR, HML. 
	#### This region should be based on the homologous alignments of the target region
	next if ($pid eq "chrIII" && $pstart >= 13650 && $pend<= 13850 );
	next if ($pid eq "chrIII" && $pstart >= 200750 && $pend <= 201000);
	
	$basic2{$qid}=$information;
	$inf2{$qid}->{length}=$qlength;
	$inf2{$qid}->{Identity}=$identity;
	$inf2{$qid}->{Matches}=$match;
	### Usually the best matches were identified from the first one, so if we could identifed them as non-inserted, potential insertion, or confident insertion, we will ignore the following blast results.
	
	next if (exists $conf2{$qid});
	next if (exists $confInsert2{$qid});
	next if (exists $potential2_chr{$qid});

	my $cov=$match/$qlength;
	

		#### potential other chromosome insertion
	if ($pid ne "chrIII"){
		$confInsert2{$qid}=$information;
		###### confident non-insertion: we require at least 70bp map against the MATA region
	}elsif($pstart>=$Rstart  && $pstart <=$Rend && $pend>=$Rstart && $pend<=$Rend && $match>=70){
		$conf2{$qid}=$information;
		##### confident insertion from other location of chr3
	}elsif($pstart>=$Rend  || $pstart<=$Rstart){
		$confInsert2{$qid}=$information;
	

		#potential insertion that are longer than 100bp
	}elsif($qlength >$MinRegionSize){
		$potential2_two{$qid}=$information;
	}else{
		$unclassfied2{$qid}=$information;
		
	}


}
close G;

unless (-e $output){
 	system ("mkdir $output");
 }else{
	 print "$output have already existed, no need to create! Continue ...\n";
}



open OUT,">$output/$output.noninsert.confident.txt" or die $!;

### including the insertion from other chromosome and chromosome:
open CI,">$output/$output.insert.blast.txt" or die $!;
#### one with unknown insertion another with insertion from other location.

#### Potential insertion only showed two edge of insertion events 4-44 for two PE, and reads length is longer than 90bp

open UN,">$output/$output.unclassified.txt" or die $!;

open OUTF,">$output/$output.Ainsertion.R1.fastq" or die $!;
open OUTR,">$output/$output.Ainsertion.R2.fastq" or die $!;

open OUTS,">$output/$output.Ainsertion.quality.txt" or die $!;

open OUTU,">$output/$output.Anoinsertion.quality.txt" or die $!;

print OUTS "ID\tReadQuality\tRIdentity\tRMatches\tRLength\tType\n";




open OUT,">$output/$output.noninsert.confident.txt" or die $!;

open R1,">$output/$output.noninsert.R1.txt" or die $!;
open R2,">$output/$output.noninsert.R2.txt" or die $!;




#open LIST,"$yeast" or die $!;


### This part only focused on the potential inserted reads

foreach my $i (keys %hash1){

	##### print confident no-insertion. If any one of forward reads and reseverse reads spanned the whole region, we considered them as no-insertions
	##### we need to firstly eliminate all these non-inserted reads;

	if (exists $conf2{$i} && exists $conf1{$i}){
		
		my ($identity1,$match1,$qlength1)=(split/\t/,$conf1{$i})[2,3,8];
		my ($identity2,$match2,$qlength2)=(split/\t/,$conf2{$i})[2,3,8];
		
		my $Fidentity=($identity1+$identity2)/2;
		my $Fmatch= ($match1+$match2)/2;
		
		print OUT "$conf1{$i}\n$conf2{$i}\n";
		print OUTU "$i\t$Rqual{$i}\t$Fidentity\t$Fmatch\t$qlength1\t$string\n";
		
	}elsif(exists $conf1{$i} && !exists $conf2{$i}){
				
		my ($identity1,$match1,$qlength1)=(split/\t/,$conf1{$i})[2,3,8];
		print R1 "$conf1{$i}\n";
		print OUTU "$i\t$Rqual{$i}\t$identity1\t$match1\t$qlength1\t$string\n";
		
	}elsif(exists $conf2{$i} && !exists $conf1{$i}){
		
		my ($identity2,$match2,$qlength2)=(split/\t/,$conf2{$i})[2,3,8];
		print R2 "$conf2{$i}\n";
		print OUTU "$i\t$Rqual{$i}\t$identity2\t$match2\t$qlength2\t$string\n";
	
		
		
	#### print confident insertion (insertion from other chromomsome and other locus of chrIII)
	}elsif(exists $confInsert1{$i} && exists $confInsert2{$i}){
		
		my ($identity1,$match1,$qlength1)=(split/\t/,$confInsert1{$i})[2,3,8];
		my ($identity2,$match2,$qlength2)=(split/\t/,$confInsert2{$i})[2,3,8];
		
		my $Fidentity=($identity1+$identity2)/2;
		my $Fmatch= ($match1+$match2)/2;
		
		print CI "$confInsert1{$i}\n$confInsert2{$i}\n";
		print OUTF "$hash1{$i}";
		print OUTR "$hash2{$i}";
		print OUTS "$i\t$Rqual{$i}\t$Fidentity\t$Fmatch\t$qlength1\t$string\n";
		

		
	#### print insertion from other locus confirmed by one read, but another reads might not contain nondetected insertion information
	### requirement: both reads larger than 100bp
	}elsif(exists $confInsert1{$i} && exists $potential2_two{$i}){
		print CI "$confInsert1{$i}\n$potential2_two{$i}\n" ;
		print CI "$confInsert1{$i}\n";	
		
		my ($identity1,$match1,$qlength1)=(split/\t/,$confInsert1{$i})[2,3,8];
		
		print OUTF "$hash1{$i}";
		print OUTR "$hash2{$i}";
		
		print OUTS "$i\t$Rqual{$i}\t$identity1\t$match1\t$qlength1\t$string\n";
		
		
	}elsif(exists $confInsert2{$i} && exists $potential2_two{$i}){
		print CI "$potential2_two{$i}\n$confInsert2{$i}\n";
		
		my ($identity2,$match2,$qlength2)=(split/\t/,$confInsert2{$i})[2,3,8];
		
		
		print OUTF "$hash1{$i}";
		print OUTR "$hash2{$i}";
		
		print OUTS "$i\t$Rqual{$i}\t$identity2\t$match2\t$qlength2\t$string\n";
		
	##### print insertions from other nondetected insertion, here we request both reads longer than 100bp

	}elsif ($len1{$i} >= $MinRegionSize && $len2{$i} >= $MinRegionSize ){
		
		print CI "$i\tnondetectable\n";
		
		my $qlength=($len1{$i}+$len2{$i})/2;
		
		print OUTF "$hash1{$i}";
		print OUTR "$hash2{$i}";
		
		print OUTS "$i\t$Rqual{$i}\tNA\tNA\t$qlength\t$string\n";
		
	### This unclassified cannot be reads with large insertion, might be with short indels.	
	}else{
		my $qlength=($len1{$i}+$len2{$i})/2;
		print UN "$i\tunclassfied\n";
		print OUTU "$i\t$Rqual{$i}\tNA\tNA\t$qlength\t$string\n";
	}

}

close R1;
close R2;

close CI;
close UN;



close OUTF;
close OUTR;
close OUTS;
close OUTU;
