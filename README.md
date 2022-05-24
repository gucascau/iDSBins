
# iDSBins
**High-throughput identification of large and complex DNA insertions at DNA double strand breaks.**

Insertions of mobile elements, mitochondrial DNA and fragments of nuclear chromosomes at DNA double-strand breaks sites (DSBs) threaten genome integrity and are common in cancer. Despite extensive efforts, our knowledge of these insertions still remains unknown. These large insertions were previously profiled along with quantitatively Sanger sequencing, but have not yet been combined with the massively parallel approaches to tackle the complexity of insertion events. Here, we introduced a high-throughput sequencing Break-Ins (DNA double strand break insertion sequencing), a method that can detect large and complicated insertion/deletion events at DNA double strand break sites. In parallel with the novel technology, we also developed a simplified, standardized and fully automated data analysis software toolkit, iDSBins, which enables routine scoring and interpretation of large-scale insertions at targeted double strand break site. 

iDSBins consists of four parts: iDSBquality, iDSBdetection, iDSBdeduplication and iDSBannotation.

	1. iDSBquality is used to eliminate the error index reads, phix reads and customed low quality reads.
	2. iDSBdetection is used to detect the reads that contained large insertion events.
	3. iDSBdeduplication is used to eliminate the duplicated large insertion events that caused by sequence errors, identify the representive read sequence and quality, and measure the read counts of final insertion events.
	4. iDSBannotation is used to annoate the number of donors, the genetic feature of donors.

This repository contain a set of simple scripts that carry out the key anaylses for identification of large and complex insertions, as well as small insertion and deletions at DNA double strand breaks in Yeast genome. These toolkits also contain scripts to define the microhomologies and small insertion/deletions at junction sequence for large insertion events, functionally annotate the insertion events, calculate the inserted coverage of rDNA, LTR, evaluate the relationship between large insertions and Rloop, ARS, hotspots, G4, Rereplication ect.  Similar strategies can also apply to the detection of insertion events at double strand break in human genome.

![GithubPipeline](https://user-images.githubusercontent.com/23031126/170120136-310bb443-5f29-4520-869f-939312b0610c.png)

# Availability 
1. Detect the large and complex large DNA insertion evens at DSB sites.
2. Measure the read counts and qualities of unique large insertion events.
3. Determine the donor number of large insertion events.
4. Functionally annotate the inserted donors.
5. Identify the hotspots of large insertions.
6. Measure the distance to the closest genomic feature.
7. Extract the junction flanking sequence (50bp) of donor.


# Dependencies

Perl is used to run the scripts. The following softwares are also required:

. blast (ncbi-blast-2.8.1+)(https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

. pear (pear-0.9.11) (https://cme.h-its.org/exelixis/web/software/pear/doc.html)

. bedtools (bedtools-2.25.0) (https://bedtools.readthedocs.io/en/latest/)

# Install

```
    cd ~
    git clone https://github.com/gucascau/iDSBins.git
```   

# Usage
```
Usage: sh Final_intergrated.sh [-a SampleID] [-b Working Directory] [-c Forward Index] [-d Reverse Index] [-f forward Reads] [-r Reverse Reads] [-p Software Installed Directory] [options]
		 
Request Parameters:
	-a Sample Id (Example: Test)
	-b The working directory, where the raw read stored and we perform the analyses of the large insertion
	-c Index of forward reads(Example: CTC)
	-d Index of reverse reads (Example: ACC)
	-f Forward reads (Example: Test_L001_R1_001.fastq)
	-r Reverse reads (Example: Test_L001_R2_001.fastq)
	-p Software installed Path, Attation: This required to install DSBins, BLAST, PEAR, Bedtools in the same folder (Default:"")

Optional Parameters:

Optional Parameters -- overall requirements:
	-n Number of threads (Default: 15)
	-gs Genome sequence (Must corrected with chromosome ID)
	-ga Genome annotation (Default BED file format)

Optional Parameters of each step:
iDSBquality - Optional paramter for the quality control of reads:
	-mq The minimum quality requirment of MAT region (Default: 25)
	-iq The minimum quality requirment of Inserted region (Default: 15)

iDSBdetection - Define the MATA information:
	-il minimum length of large insertion (Default 10bp)
	-ms Total size of whole MATA region (Default 90)
	-mc Mapped chromosme of MATA reference position (Default chrIII)
	-ms Mapped start site of MATA reference position (Default 294300)
	-me Mapped end site of MATA reference position (Default 294500)

iDSBdeduplication - parameters for the duplicated reads:

	Parameters for first cuting-edge deduplication:
	-cs cut-off of upstream and downstream, the cut-off we collected mata 11bp and inserted 19bp sequences for forward and reverse reads(default 30bp, 11+19bp)
	-cf The cutting start site of Forward reads (default 33)
	-cr The cutting start site of reverse reads (default 39)

	Parameters for reference-based deuplication:
	-di The minimum identity for two deduplicated reads(Default: 95)
	-do The minimum overlapping ratio for two deduplicated reads(Default: 0.95)
	-dd The maximum number of indels for two deduplicated reads, we set up the same values for different processes of deduplcations(Default: 6)
	-dm The maximum number of mismatches for two deduplicated reads, we set up the same values for different processes of deduplcations(Default: 6)
	-dl The maximum aligned gapsize of forward and reverse read that (Default: 30000)
	-dc The minimum coverage supports of each insertion events (Default:1)
	-dq The minimum quality supports of each insertion events (Default:1)

iDSBannotation: parameters to detemine the number of donor and annotation of donors
	-hm The maxmium microhomologous length between donors (default 33)
	-ho Overlapping region between donor, less than half of min donor size (default 0.5)

	-h help

```


# Output
The software will generate two seperate files for large insertion events, including single donor, and multiple donors. The SampleID_final.Single.txt contains insertion events that only have one donor.  SampleID_final.Multiple.txt contains quite complex insertion events, including 2 donors, 3 donors, 4 donors, 1orMore donors, 2orMore donors, 3orMore donors, and unaligned donors. 

| column | explaination |
| ------| ------|
| 1st | Case ID |
| 2nd | Representative Read |
| 3rd | Sample ID |
| 4th | Whole sequence |
| 5th | Inserted Sequence |
| 6th | Insertion Size (bp)|
| 7th | Donor Number |
| 8th | Representive Read Identity |
| 9th | Total Read Count |
| 10th | Representive Read Quality |
| 11th | Donor  Chromosome |
| 12th | Donor Start |
| 13th | Donor End |
| 14th | Donor Strand |
| 15th | Donor Feature |
| 16th | Distance To Feature |
| 17th | Distance Description |
| 18th | Alternative Feature |
| 19th | Junction sequence at 5' |
| 20th | Junction sequence at 3'|
| 21th | Whethe estimated insertion based on gap Of reads status (3kb) |
| 22th | Start distance between neighbor insertions|
| 23th | End distance between neighbor insertions |

For each insertion event, we selected the most reliable representative read that showed the highest read quality, highest donor identity and the highest read count support. 

For more detail information, please feel free to contact: xin.wang@childrens.harvard.edu
