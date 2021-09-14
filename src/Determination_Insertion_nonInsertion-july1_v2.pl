#!/usr/bin/perl
use strict;
use warnings;
my $version="1.0 version";
use Getopt::Long;

### This script is to extract the reads that potenital contain large insertion events:
### requirement: 1
#####1.	We considered the reads without insertion if they completely blasted against reference genome from start to end of the reads, ignoring whether they have small indels. 
#####    And as long as one read span the whole MATA region, we considered these reads without insertion events. 
#####2.	Eliminating the reads with non-insertion, We considered the read as insertion event following the below criterion:
#####	2.1 If either of the reads can mapped against other region, we confidently considered them as large insertion
#####	2.2 If there is no other alignment, We requsted the shortest insertion larger than 10bp (This size can also be adjustable).


my %opts;
GetOptions(\%opts,"i:s","o:s","g:s","c:s","t:s","e:s","m:s","f:s","r:s","q:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{g} ||!defined $opts{f}||!defined $opts{r} ||!defined $opts{q} ||!defined $opts{o}) {
	die "************************************************************************
	Usage: extract_fastq.pl -i R1_Blast -g R2_Blast -f R1_fastq -r R2_fastq  -o Blast_categories
	
	Request Parameters:
	-i Forward Blast file
	-g Reverse Blast file
	-f Forward fastq file
	-r Reverse fastq file
	-q Quality file generated before
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

my $chr=(defined $opts{c})?$opts{c}:"chrIII";
my $Rstart= (defined $opts{t})?$opts{t}:294300;
my $Rend= (defined $opts{e})?$opts{e}:294500;
my $Isize= (defined $opts{m})?$opts{m}:10;
my $Msize= (defined $opts{n})?$opts{n}:90;

my $MinRegionSize=$Isize+$Msize;

#my $yeast=$opts{m};

my $Fread=$opts{f};
open FASTQ1,"$Fread" or die $!;

my $id1;  my %hash1;  
while (<FASTQ1>) {
	chomp;
    # print "$_" ;

    if ($. % 4 == 1)  {
		$id1=(split/\s+/,$_)[0];
		$id1=~s/^@//;
		#print "$id1\n";
		$hash1{$id1}=$_."\n";
	}else{
		$hash1{$id1}.=$_."\n";
	}
	
}
close FASTQ1;

my $Rread=$opts{r};

open FASTQ2,"$Rread" or die $!;

my $id2;  my %hash2; 
while (<FASTQ2>) {
	chomp;
    # print "$_" ;

    if ($. % 4 == 1)  {
		$id2=(split/\s+/,$_)[0];
		$id2=~s/^@//;
		#print "$id1\n";
		$hash2{$id2}=$_."\n";
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
	### Usually the best matches were identified from the first one, so if we could identifed them as non-inserted, potential insertion, or confident insertion, we will ignore the following blast results.
	
	next if (exists $conf1{$qid});
	next if (exists $confInsert1{$qid});
	next if (exists $potential1_chr{$qid});
	
	my $cov=$match/$qlength;
	print "$cov\n";
	

	#### Confident insertions from other chromosome:
	if ($pid ne "chrIII"){
		$confInsert1{$qid}=$information;	
		
	###### confident non-insertion, might contain some deletion events:
	}elsif($pstart>=$Rstart  && $pstart <=$Rend && $pend>=$Rstart && $pend<=$Rend && $match>=70){
		$conf1{$qid}=$information;
		
	##### Confident insertion from other location of chr3:
	}elsif($pstart>=$Rend  || $pstart<=$Rstart){
		$potential1_chr{$qid}=$information;
		
		
	#### potential split the fragment with insertion (telemere repetitive elements or MAT elements)
	}elsif($pstart>=$Rstart  && $pstart <=$Rend && $pend>=$Rstart && $pend<=$Rend && $match<=50 ){
		$potential1_two{$qid}=$information;
		####################################	
		$inf1{$qid}->{start}=$qstart;
		$inf1{$qid}->{match}=$match;
		$inf1{$qid}->{length}=$qlength;
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
	### Usually the best matches were identified from the first one, so if we could identifed them as non-inserted, potential insertion, or confident insertion, we will ignore the following blast results.
	
	next if (exists $conf2{$qid});
	next if (exists $confInsert2{$qid});
	next if (exists $potential2_chr{$qid});

	my $cov=$match/$qlength;
	

		#### potential other chromosome insertion
	if ($pid ne "chrIII"){
		$confInsert2{$qid}=$information;
		###### confident non-insertion:
	}elsif($pstart>=$Rstart  && $pstart <=$Rend && $pend>=$Rstart && $pend<=$Rend && $match>=70){
		$conf2{$qid}=$information;
		##### confident insertion from other location of chr3
	}elsif($pstart>=$Rend  || $pstart<=$Rstart){
		$potential2_chr{$qid}=$information;
		#### potential insertion or deletion but need to check original data  
		#### potential split the fragment with insertion (telemere repetitive elements or MAT, therefore these regions located with two sides of MATA)
	}elsif($pstart>=$Rstart  && $pstart <=$Rend && $pend>=$Rstart && $pend<=$Rend && $match<=50 ){
		$potential2_two{$qid}=$information;
	####################################	
		$inf2{$qid}->{start}=$qstart;
		$inf2{$qid}->{match}=$match;
		$inf2{$qid}->{length}=$qlength;
	}else{
		$unclassfied2{$qid}=$information;
		
	}


}
close G;







open OUT,">$output/$output.noninsert.confident.txt" or die $!;

open R1,">$output/$output.noninsert.R1.txt" or die $!;
open R2,">$output/$output.noninsert.R2.txt" or die $!;

### including the insertion from other chromosome and chromosome:
open CI,">$output/$output.insert.two_confident.txt" or die $!;
#### one with unknown insertion another with insertion from other location.
open OI,">$output/$output.insert.potential_otherchr.txt" or die $!;
open OIC,">$output/$output.insert.potential_chr3.txt" or die $!;
#### Potential insertion only showed two edge of insertion events 4-44 for two PE, and reads length is longer than 90bp
open GII, ">$output/$output.Pinsertion.txt" or die $!;


###### if we removed the bias location 
####next if ($pid eq "chrIII" && $pstart<13800 && $pstart>13600);####next if ( $pid eq "chrIII" && $pstart<200950 && $pstart>200800);

#open PII,">$output.Pinsertiondoublecheck.txt" or die $!;
#open PID,">$output.Pdeletiondoublecheck.txt" or die $!;
#open PIU,">$output.Pdeletion_uncertaindoublecheck.txt" or die $!;

open UN,">$output/$output.unclassified.txt" or die $!;


#open LIST,"$yeast" or die $!;


### This part only focused on the potential inserted reads

foreach my $i (keys %basic1){

	
	##### print confident no-insertion. If any one of forward reads and reseverse reads spanned the whole region, we considered them as no-insertions
	if (exists $conf2{$i} && exists $conf1{$i}){
		print OUT "$conf1{$i}\n$conf2{$i}\n";
	}elsif(exists $conf1{$i} && !exists $conf2{$i}){
		
		
		print R1 "$conf1{$i}\n$basic2{$i}\n" if (exists $basic2{$i});
		print R1 "$conf1{$i}\n" if (!exists $basic2{$i});
	}elsif(exists $conf2{$i} && !exists $conf1{$i}){
		print R2 "$basic1{$i}\n$conf2{$i}\n";
		
	
		
	#### print confident insertion (insertion from other chromomsome)
	}elsif(exists $confInsert1{$i} && exists $confInsert2{$i}){
		print CI "$confInsert1{$i}\n$confInsert2{$i}\n";
		
		#### print insertion from other chromosome confirmed by one read, but another reads might not contain insertion information
	}elsif(exists $confInsert1{$i} && !exists $confInsert2{$i}){
		print OI "$confInsert1{$i}\n$basic2{$i}\n" if (exists $basic2{$i});
		print OI "$confInsert1{$i}\n" if (!exists $basic2{$i});
		
	}elsif(exists $confInsert2{$i} && !exists $confInsert1{$i}){
		print OI "$basic1{$i}\n$confInsert2{$i}\n";
		
	#### print insertion from Chromosome 3 other location
	}elsif(exists $potential1_chr{$i} && exists $potential2_chr{$i}){
		print CI "$potential1_chr{$i}\n$potential2_chr{$i}\n";	
		##### print insertion from Chromosome 3 other location but other elements are also inside
	}elsif(exists $potential1_chr{$i} && !exists $potential2_chr{$i}){
		
		print OIC "$basic1{$i}\n$basic2{$i}\n" if (exists $basic2{$i});	
		print OIC "$basic1{$i}\n" if (!exists $basic2{$i});	
			
	}elsif(exists $potential2_chr{$i} && !exists $potential1_chr{$i}){
		print OIC "$basic1{$i}\n$basic2{$i}\n";
			
	#### print unconfirmed insertion or deletions of exactly 294345-294384 
	### normally this happened with small insertion or small deletion based on their correpsonding size
	### Here we perform and function in order to reconsider the best last that showed aligned to 13600- 13700
	### here we further divide the potential insertion and confirm their insertion events if the length larger than 90bp
	### The following are not quite confident due to the 
	
	}elsif(exists $potential1_two{$i} && exists $potential2_two{$i} ){
		#print GI "$basic1{$i}\n$basic2{$i}\n";
		
		### We considered both read length longer than 90bp containing potential insertion events
		print GII "$basic1{$i}\n$basic2{$i}\n" if ($inf1{$i}->{length} >= $MinRegionSize && $inf2{$i}->{length} >= $MinRegionSize);

		### print unknown information
	}elsif ($inf1{$i}->{length} >= $MinRegionSize ){
		print GII "$basic1{$i}\n$basic2{$i}\n" if (exists $basic2{$i});
		print GII "$basic1{$i}\n" if (!exists $basic2{$i});	
	}else{
		### This unclassified cannot be reads with large insertion
		print UN "$basic1{$i}\n$basic2{$i}\n" if (exists $basic2{$i});	
		print UN "$basic1{$i}\n"  if (!exists $basic2{$i});	
	}

}

close R1;
close R2;
close OI;

close CI;
close UN;

close GII;
close OIC;

# close PID;
# close PII;
# close PIU;
