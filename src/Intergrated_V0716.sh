#!/bin/bash
# Sample batchscript to run a simple job array to create 5 different files, and modify the files independently with one batchscript
 
#SBATCH --partition=bch-compute # queue to be used
#SBATCH --time=30:00:00 # Running time (in hours-minutes-seconds)
#SBATCH --job-name=WTAD16_LargeInsertion # Job name
#SBATCH --mail-type=BEGIN,END,FAIL # send and email when the job begins, ends or fails
#SBATCH --mail-user=xin.wang@childrens.harvard.edu # Email address to send the job status
#SBATCH --output=output_%A_%a.txt # Name of the output file
#SBATCH --nodes=1 # Number of compute nodes
#SBATCH --ntasks=10 # Number of cpu cores on one node
#SBATCH --mem=20G



#### Here is the last version of large insertion events

 
echo " Detection job array is started ..."
date
 
i=yYY517A-M-16d_S19

### Set up path file:


in=/lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast

mkdir ${i}
cd ${in}/${i}

###############################################################################################################
#### 1. Quality Control
###############################################################################################################

#### 1.1 we removed trimmed low quailty sequences (less than 10)
### also removed all reads that have a 31-mer match to PhiX, allowing two mismatches. 
####detemrined the indext information

### First round the quality control
mkdir QualityControl
cd ${in}/${i}/QualityControl

perl  /temp_work/ch220812/Project/script/DSBinsertion/Extract_fastq_withindex.pl -f ../../${i}_L001_R1_001.fastq -r ../../${i}_L001_R2_001.fastq -i TGG -g AGG -o ${i}_filter

#### 1.2 filter phix sequences

/home/ch220812/software/bbmap/bbduk.sh -in1=${i}_filter.R1.fastq in2=${i}_filter.R2.fastq out1=${i}_filter.unmaped.R1.fastq out2=${i}_filter.unmaped.R2.fastq outm1=${i}_filter.mappedphix.R1.fastq outm2=${i}_filter.mappedphix.R2.fastq ref=/home/ch220812/software/bbmap/resources/phix_adapters.fa.gz k=31 hdist=2 stats=stats.txt

#### 1.3 set the stringint quality of beginning of 25 and ending of 15, here we also generate the fasta file of each type of reads 

perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/Divide_high_lowQuality_v2.pl -f ${i}_filter.unmaped.R1.fastq -r ${i}_filter.unmaped.R2.fastq -u 25 -g 15 -o ${i}_QC


###############################################################################################################
# 2 .Detection
###############################################################################################################


cd ${in}/${i}/

mkdir Detection

cd Detection

#### 2.1  First round with the unmerged reads

### We perform the blast for each forward and reverse read

/home/ch220812/software/ncbi-blast-2.8.1+/bin/blastn -query ../QualityControl/${i}_QC.highquality.R1.fasta -out ${i}.forward.blast.tbl  -db /home/ch220812/database/YeastGenomeIndex/S288C/S288C_reference_genome_R64-2-1_20150113/modified/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa -num_threads 8 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -max_target_seqs 1

/home/ch220812/software/ncbi-blast-2.8.1+/bin/blastn -query ../QualityControl/${i}_QC.highquality.R2.fasta -out ${i}.reverse.blast.tbl -db /home/ch220812/database/YeastGenomeIndex/S288C/S288C_reference_genome_R64-2-1_20150113/modified/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa -num_threads 8 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -max_target_seqs 1


### Then detemine the insertion events based on the forward and reverse reads, we could set up different parameter accordingly. 
### We generated a folder that contains the reads with large insertion events and a statistical file with the information of read quality, identity, mapped size.
### This statistical file will be further used to detemine the unique large insertion event.
### Here we used all the reads: the reads information put to the fastq file

perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/Determination_Insertion_nonInsertion-july1_v3.pl -i ${i}.forward.blast.tbl -g ${i}.reverse.blast.tbl  -f ../QualityControl/${i}_QC.highquality.R1.fastq -r ../QualityControl/${i}_QC.highquality.R2.fastq -q ../QualityControl/${i}_QC.Allquality.stat -o ${i}_detected_first -j First

date


### 2.2 To further confirm the insertion event of each read, second round with the merge methods were applied.
## Assemble All the reads, in this case we would eliminate lots of sequence errors generated at the end of reads

# This will generate some assembled files and unassembled files
/home/ch220812/software/pear-0.9.11-linux-x86_64/bin/pear -f ../QualityControl/${i}_QC.highquality.R1.fastq -r ../QualityControl/${i}_QC.highquality.R2.fastq -o ${i}_merged

### tranfer into fastq
 perl /temp_work/ch220812/Project/script/DSBinsertion/fastq2fasta.pl -i ${i}_merged.assembled.fastq -o ${i}_merged.assembled.fasta

 perl /temp_work/ch220812/Project/script/DSBinsertion/fastq2fasta.pl -i ${i}_merged.unassembled.forward.fastq -o ${i}_merged.unassembled.forward.fasta
 perl /temp_work/ch220812/Project/script/DSBinsertion/fastq2fasta.pl -i ${i}_merged.unassembled.reverse.fastq -o ${i}_merged.unassembled.reverse.fasta

 

### perform the blast analyes
/home/ch220812/software/ncbi-blast-2.8.1+/bin/blastn -query ${i}_merged.assembled.fasta -out ${i}_merged.assembled.blast.tbl  -db  /home/ch220812/database/YeastGenomeIndex/S288C/S288C_reference_genome_R64-2-1_20150113/modified/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa -num_threads 8 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -max_target_seqs 1

/home/ch220812/software/ncbi-blast-2.8.1+/bin/blastn -query ${i}_merged.unassembled.forward.fasta -out ${i}_merged.unassembled.forward.blast.tbl  -db /home/ch220812/database/YeastGenomeIndex/S288C/S288C_reference_genome_R64-2-1_20150113/modified/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa -num_threads 8 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -max_target_seqs 1

/home/ch220812/software/ncbi-blast-2.8.1+/bin/blastn -query ${i}_merged.unassembled.reverse.fasta -out ${i}_merged.unassembled.reverse.blast.tbl -db /home/ch220812/database/YeastGenomeIndex/S288C/S288C_reference_genome_R64-2-1_20150113/modified/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa -num_threads 8 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -max_target_seqs 1


### Detect the insertion for the second methods

# for assmbled one
perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/Determination_Insertion_nonInsertion_single_blast_after_merge_v3.pl -i ${i}_merged.assembled.blast.tbl -q ../QualityControl/${i}_QC.Allquality.stat -f ${i}_merged.assembled.fasta -o ${i}_detected_second_assembled -j SecondAssembly


# for unassmbled ones
### Here we used unassembled the reads: the reads information put to the fastq file
perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/Determination_Insertion_nonInsertion-july1_v3.pl -i ${i}_merged.unassembled.forward.blast.tbl -g ${i}_merged.unassembled.reverse.blast.tbl  -f ${i}_merged.unassembled.forward.fastq -r ${i}_merged.unassembled.forward.fastq -q ../QualityControl/${i}_QC.Allquality.stat -o  ${i}_detected_second_unassembled -j SecondUnAssembly

## 2.3 We further combined two round of large insertion events and confirm whether each read contain insertion event consistently. 

## If reads identified in both method, we confirm these reads as inserted reads. If reads only confirmed only in one method, we will check their read read and determined based on the identity (90) and quality ()

mkdir ${i}_detected_final
cd ${i}_detected_final



cp ../${i}_detected_second_assembled/${i}_detected_second_assembled.assmebled.Ainsertion.quality.txt ../${i}_detected_second_unassembled/${i}_detected_second_unassembled.Ainsertion.quality.txt ../${i}_detected_first/${i}_detected_first.Ainsertion.quality.txt ./
cat ${i}_detected_second_assembled.assmebled.Ainsertion.quality.txt ${i}_detected_second_unassembled.Ainsertion.quality.txt >${i}_detected_second.Ainsertion.quality.txt
## check consistency
perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/CheckConsistentInsertion.pl -i ${i}_detected_first.Ainsertion.quality.txt -g ${i}_detected_second.Ainsertion.quality.txt -o ${i}_detected_combined
 
 
 
echo "First detection job array is finished ..."


 
 #### 
echo "Deduplication job array is started ..."

cd ${in}/${i}/
mkdir DeDuplication
cd DeDuplication
mkdir DeduplicationFirst
cd DeduplicationFirst

echo "First round of deduplication is started: junction based ..."


###########################################################################################################
#####  2. First deduplication:
#####.  cut the edge of two sides of each inserted reads and perform the first round of deduplication
#####.  criterion: the each junction sequence contains 19bp inserted sequence and 11bp MATA sequences
###########################################################################################################


cp ${in}/${i}/Detection/${i}_detected_final/${i}_detected_combined.highqual.txt ./

### 1.1  we combined the two junction sequences of all detected reads
perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/Cut_twoedge_sequences_forunassembled_v2.pl  -f ${in}/${i}/QualityControl/${i}_QC.highquality.R1.fastq -r ${in}/${i}/QualityControl/${i}_QC.highquality.R2.fastq -i ${i}_detected_combined.highqual.txt -o ${i}_detected_cut60.fasta

#### 1.2 Blast against itself
/home/ch220812/software/ncbi-blast-2.8.1+/bin/makeblastdb -dbtype nucl -in ${i}_detected_cut60.fasta
 
/home/ch220812/software/ncbi-blast-2.8.1+/bin/blastn -query ${i}_detected_cut60.fasta -db ${i}_detected_cut60.fasta -out ${i}_detected_cut60.blast -num_threads 10 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue'

### 1.3 We perform the first round of deduplications using the junction sequences
perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/SelfBlast_RemoveDuplication_eachsampleIdentitySorted_v2.pl  -b ${i}_detected_cut60.blast -i ${i}_detected_combined.highqual.txt -f ${in}/${i}/QualityControl/${i}_QC.highquality.R1.fastq -r ${in}/${i}/QualityControl/${i}_QC.highquality.R2.fastq -o ${i}_detected_FirstDep


echo "First round of deduplication is ended: junction based ..."

date

###########################################################################################################
#####  2. Second deduplication:
#####.  we used the referenced based method to perform the deduplication, allowing the frameshift 
#####.  criterion: 
###########################################################################################################

echo "Second round of deduplication is started: reference based ..."

cd ..
mkdir DeduplicationSecond
cd DeduplicationSecond
cp ../DeduplicationFirst/${i}_detected_FirstDep.cls ./

### 2.1 Prepair the required information of second round deduplication: the unassembled reads, assembled reads, and their blast alignments. 


###    Extract the representive reads detected from the first round, then seperate them into assembled and unassembled reads.

perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/SperatedInsertionbasedAssembly.pl -g ${in}/${i}/Detection/${i}_merged.assembled.fastq -f ${in}/${i}/Detection/${i}_merged.unassembled.forward.fastq -r ${in}/${i}/Detection/${i}_merged.unassembled.reverse.fastq -i ${i}_detected_FirstDep.cls -o ${i}_detected_assemblysep


#### for the assmbled sequences 
### 2.2 Divide the assembled file into different categories based on the number of insertions

### For assmebled we determine the number of large insertions

## double check the blastn results
/home/ch220812/software/ncbi-blast-2.8.1+/bin/blastn -query ${i}_detected_assemblysep.Insassembled.fasta -out ${i}_detected_assemblysep.Insassembled.blast_all_detected.tbl -db /home/ch220812/database/YeastGenomeIndex/S288C/S288C_reference_genome_R64-2-1_20150113/modified/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa -num_threads 8 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue'

perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/Determine_insertion_number_formutipleinserts_forallassembled_universal_updated_v3.pl -g ${i}_detected_assemblysep.Insassembled.blast_all_detected.tbl -i ${i}_detected_assemblysep.Insassembled.fasta  -c ${i}_detected_assemblysep.CoverageWithAssembly.txt -o ${i}_detected_number

### Due to many insertions cannot identify blast results, we set a BLAST method using blastn-short

/home/ch220812/software/ncbi-blast-2.8.1+/bin/blastn  -task blastn-short  -word_size 11 -evalue 0.001 -dust no -soft_masking false -gapopen 5 -penalty -1 -perc_identity 80 -query  ${i}_detected_number.noaligned.fasta -out ${i}_detected_number.noaligned.blast.tbl -db /home/ch220812/database/YeastGenomeIndex/S288C/S288C_reference_genome_R64-2-1_20150113/modified/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa -num_threads 8 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue'

perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/Determine_insertion_number_formutipleinserts_forallassembled_universal_updated_v3.pl  -g ${i}_detected_number.noaligned.blast.tbl -i ${i}_detected_number.noaligned.fasta -c ${i}_detected_number.noaligned.txt -o ${i}_detected_number.noaligned


### We integrated all these insertion with one donor from assembly and unassembly



#### 1 cat single donor insertion, here included some mutiple inseritons but only detected one donor, single donor and single estimated donor### we allow for 6bp frameshife of 

cat ${i}_detected_number.single.txt ${i}_detected_number.noaligned.single.txt >${i}_detected_second_singleassembly.txt
perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/Single_confident_donor_coverage_deduplication_v3.pl -i ${i}_detected_second_singleassembly.txt -a 6 -b 6  -q ${in}/${i}/QualityControl/${i}_QC.highquality.stat -o ${i}_detected_second_singleassembly


### 2. two donor insertion 

cat ${i}_detected_number.twoinf.txt ${i}_detected_number.noaligned.twoinf.txt >${i}_detected_second_twoassembly.txt

perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/Two_confident_donor_coverage_deduplication_universal_v3.pl -i ${i}_detected_second_twoassembly.txt -a 6 -b 6 -q ${in}/${i}/QualityControl/${i}_QC.Allquality.stat -o ${i}_detected_second_twoassembly

.
### 3. Three donor insertion
cat ${i}_detected_number.threeinf.txt ${i}_detected_number.noaligned.threeinf.txt >${i}_detected_second_threeassembly.txt
perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/Three_confident_donor_coverage_deduplication_universal_v3.pl -i ${i}_detected_second_threeassembly.txt -a 6 -b 6 -q ${in}/${i}/QualityControl/${i}_QC.Allquality.stat -o ${i}_detected_second_threeassembly


### 4. Four donor insertion

cat ${i}_detected_number.fourinf.txt ${i}_detected_number.noaligned.fourinf.txt >${i}_detected_second_fourassembly.txt

perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/Four_confident_donor_coverage_deduplication_universal_v3.pl -i ${i}_detected_second_fourassembly.txt -a 6 -b 6 -q ${in}/${i}/QualityControl/${i}_QC.Allquality.stat -o ${i}_detected_second_fourassembly


## 5. Final Undetected insertions, we will double check whether they came from MAT



mkdir MAT
cd MAT
perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/Seperate_unmappable_mat_uniq_v3.pl -f  ../${i}_detected_number.noaligned.noaligned.fasta -b ../../../../${i}_merged.assembled.blast.tbl -o ${i}_assembled_final_undetermined



##### 6 single estimated donor insertion

#### for the unassmbled sequences #### 
### for the unassembled sequence, we mainly estimated the large insertion if two reads blast against two side of simialr region within 5kb

### Double check the unassembled sequence and generated the confident blast file of large insertions that can identity the donor information in both reads.
### Here we used unassembled the reads: the reads information put to the fastq file
cd ../
perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/Determination_Insertion_nonInsertion-july1_v3.pl -i ${in}/${i}/Detection/${i}_merged.unassembled.forward.blast.tbl -g ${in}/${i}/Detection/${i}_merged.unassembled.reverse.blast.tbl  -f ${i}_detected_assemblysep.InsUnassembled.R1.fastq -r ${i}_detected_assemblysep.InsUnassembled.R2.fastq -q ${in}/${i}/QualityControl/${i}_QC.Allquality.stat -o ${i}_unassembledEstimate -j UnassEst

cd Test_unassembledEstimate

  ##### estimate the single insertion  #### here we can only determine the single insertion
perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/Yeast_overalapping_insertion_from_R1R2_v3.pl -i ${i}_unassembledEstimate.insert.blast.txt -c ../${i}_detected_assemblysep.CoverageWithAssembly.txt -g 5000 -o ${i}_unassembledEstimate

# Here we will not consider the insertions that cannot well aligned to the end of R1 and the begining of R2, we allow 5bp shift at each end 
perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/Single_confident_donor_coverage_deduplication_estimated_v3.pl -i ${i}_unassembledEstimate.oneinsert.txt -a 6 -b 6 -s 5 -q ${in}/${i}/QualityControl/${i}_QC.Allquality.stat -o ${i}_detected_second_estimated

echo "Referenced deduplication job array is finshed ..."

date

##############################################################################
### Part3.  extract the sequence and export the number of donors and annotation####
##################################################################################



#### Generate the final insertion events


echo "Perform annotation of each categories of large insertions"

## 3.1 Prepair all the poteintial feature information of each assembled reads
# Identify all blast of each assembled reads with default paramter:


### For assembled: Blast aginst all assembled (First step using the default blast)
cd ../

#/home/ch220812/software/ncbi-blast-2.8.1+/bin/blastn -query ${i}_detected_assemblysep.Insassembled.fasta -out ${i}_detected_assemblysep.Insassembled.blast_all_detected.tbl -db /home/ch220812/database/YeastGenomeIndex/S288C/S288C_reference_genome_R64-2-1_20150113/modified/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa -num_threads 8 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' 

### integrate the default blast and the unaligned blast (shorter parameter):
cat ${i}_detected_number.noaligned.blast.tbl ${i}_detected_assemblysep.Insassembled.blast_all_detected.tbl >${i}_detected_assemblysep.Insassembled.FinallyBlast.tbl

### Generate the bedfile index:
perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/Create_bedfile_forblast_determinemutiplelocus_v3.pl -i ${i}_detected_assemblysep.Insassembled.FinallyBlast.tbl -o ${i}_detected_assemblysep.Insassembled.FinallyBlast.bed

## sort and annotation
sort -k 1V,1 -k 2n,2 ${i}_detected_assemblysep.Insassembled.FinallyBlast.bed >${i}_detected_assemblysep.Insassembled.FinallyBlast.sorted.bed
/home/ch220812/software/bedtools2/bin/bedtools closest -D b -t first -b /home/ch220812/database/YeastGenomeIndex/S288C/S288C_reference_genome_R64-2-1_20150113/modified/saccharomyces_cerevisiae_R64-2-1_20150113_version6.bed -a ${i}_detected_assemblysep.Insassembled.FinallyBlast.sorted.bed >${i}_detected_assemblysep.Insassembled.FinallyBlast.annotation.bed

### Annotate the insertion with single donor, we need the blast result to confirm whether they have the alternative feature
perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/Annotation_withupdownstream_insertion_For1donor_v3.pl -a ${i}_detected_assemblysep.Insassembled.fasta -b ${in}/${i}/QualityControl/${i}_QC.Allquality.stat -c ${i}_detected_second_singleassembly.uniq.txt -d ${i}_detected_assemblysep.Insassembled.FinallyBlast.tbl -e ${i}_detected_assemblysep.Insassembled.FinallyBlast.annotation.bed -t ${i} -o ${i}_final_single


### Annotate the insertion for two donor
perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/Annotation_withupdownstream_insertion_For2donors_v3.pl -a ${i}_detected_assemblysep.Insassembled.fasta -b ${in}/${i}/QualityControl/${i}_QC.Allquality.stat -c ${i}_detected_second_twoassembly.uniq.txt -d ${i}_detected_assemblysep.Insassembled.FinallyBlast.tbl -e ${i}_detected_assemblysep.Insassembled.FinallyBlast.annotation.bed -t ${i} -o ${i}_final_two

### Annotate the insertion for three donor
perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/Annotation_withupdownstream_insertion_For3donors_v3.pl -a ${i}_detected_assemblysep.Insassembled.fasta -b ${in}/${i}/QualityControl/${i}_QC.Allquality.stat -c ${i}_detected_second_threeassembly.uniq.txt -d ${i}_detected_assemblysep.Insassembled.FinallyBlast.tbl -e ${i}_detected_assemblysep.Insassembled.FinallyBlast.annotation.bed -t ${i} -o ${i}_final_three

### Annotate the insertion for four donor

perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/Annotation_withupdownstream_insertion_For4donors_v3.pl -a ${i}_detected_assemblysep.Insassembled.fasta -b ${in}/${i}/QualityControl/${i}_QC.Allquality.stat -c ${i}_detected_second_fourassembly.uniq.txt -e ${i}_detected_assemblysep.Insassembled.FinallyBlast.annotation.bed -d ${i}_detected_assemblysep.Insassembled.FinallyBlast.tbl -t ${i} -o ${i}_final_four

### Annotate the estimated insertion for esimated donor

cd ${i}_unassembledEstimate

perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/Change_ID_fordeduplicatedestmateone_insertion_tobedsequence_versionOct_v3.pl -g /home/ch220812/database/YeastGenomeIndex/S288C/S288C_reference_genome_R64-2-1_20150113/modified/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa -m ${i}_unassembledEstimate.Ainsertion.R1.fastq -n ${i}_unassembledEstimate.Ainsertion.R2.fastq -f ${i}_detected_second_estimated.uniq.txt -t ${i} -o ${i}_estimated

### double check whether they do have alternative feature

/home/ch220812/software/ncbi-blast-2.8.1+/bin/blastn -query ${i}_estimated.wholeseq.fasta -out ${i}_estimated.wholeseq.blast -db /home/ch220812/database/YeastGenomeIndex/S288C/S288C_reference_genome_R64-2-1_20150113/modified/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa -num_threads 8 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue'


### Annotate the estimated insertion.
perl -alne '{my ($start,$end)=($F[10]>$F[9])?($F[9],$F[10]):($F[10],$F[9]); my $strand= ($F[10]>$F[9])?"+":"-"; my $string= join "\t", $F[1],$start,$end; next if (exists $hash{$string}) ; $hash{$string}++; print "$string\t$F[0]\t0\t$strand"}' ${i}_estimated.wholeseq.blast |sort -k 1,1V -k 2,2n >${i}_estimated.wholeseq.bed
/home/ch220812/software/bedtools2/bin/bedtools closest -D b -t first -b /home/ch220812/database/YeastGenomeIndex/S288C/S288C_reference_genome_R64-2-1_20150113/modified/saccharomyces_cerevisiae_R64-2-1_20150113_version6.bed -a ${i}_estimated.wholeseq.bed >${i}_estimated.wholeseq.annotation.bed

### Add the alterative feature for estimated insertion
perl /lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/20210712_Age/Finished_blast/scripts/FurtherEstimatedFeature_annotations_v3.pl -i ${i}_estimated.txt -b ${i}_estimated.wholeseq.blast -c ${i}_estimated.wholeseq.annotation.bed -o ${i}_estimated.final.txt


############################################
### Generate final insertion events
############################################

cd ../../
mkdir DeduplicationFinal
cd DeduplicationFinal

cp ../DeduplicationSecond/${i}_final_single.Ftable ../DeduplicationSecond/${i}_final_two.Ftable ../DeduplicationSecond/${i}_final_three.Ftable ../DeduplicationSecond/${i}_final_four.Ftable ../DeduplicationSecond/${i}_detected_number.noaligned.noaligned.fasta  ../DeduplicationSecond/${i}_unassembledEstimate/${i}_estimated.final.txt ./
