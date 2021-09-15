#!/bin/sh

#### you need to set up the software environment, ENV is where you install your all softwares blast, pear, bedtools and iDSBins
ENV="/scratch/ch220812/Software/"
sh ${ENV}/iDSBins/src/Final_intergrated.sh -a Test -b $PWD -c CAG -d ATC -f Test_L001_R1_001.fastq -r Test_L001_R2_001.fastq -p ${ENV}