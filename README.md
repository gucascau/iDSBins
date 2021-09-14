# iDSBins
High-throughput identification of large and complex DNA insertions at DNA double strand breaks.

Insertions of mobile elements, mitochondrial DNA and fragements of nuclear chromosomes at DNA double-strand breaks(DSB) threaten genome instability. This genomic instability can rebalance the effects of specific knockout, i.e. DNA2 mutants.

This repository contain a set of simple scripts that carry out the key anaylses for identification of large and complex insertions, as well as small insertion and deletions at DNA double strand breaks in Yeast genome. These toolkits also contain scripts to define the microhomologies and small insertion/deletions at junction sequence for large insertion events, functionally annotate the insertion events, calculate the inserted coverage of rDNA, LTR, evaluate the relationship between large insertions and Rloop, ARS, hotspots, G4, Rereplication ect.  Similar strategies can also apply to the detection of insertion events at double strand break in human genome.

iDSBins consists of four parts: iDSBquality, iDSBdetection, iDSBdeduplication and iDSBannotation.

1. iDSBquality is used to eliminate the error index reads, phix reads and customed low quality reads.
2. iDSBdetection is used to detect the reads that contained large insertion events.
3. iDSBdeduplication is used to eliminate the duplicated large insertion events that caused by sequence errors, identify the representive read sequence and quality, and measure the read counts of final insertion events.
4. iDSBannotation is used to annoate the number of donors, the genetic feature of donors.

# Availability 
1. Detect the large and complex DNA insertions at DSB sites.
2. Detect the small deletion/insertion at DSB sites.
3. Perform the deduplication and calculate the coverage of each insertion events.
4. Detect the microhomologies and junctions for large insertion events.
5. Functionally annotate the inserted elements.
6. Identify the hotspots of large insertions.
7. Requirements


# Dependencies

Perl is used to run the scripts. The following softwares are also required:

. blast (https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

. bbduk (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/)

. pear (https://cme.h-its.org/exelixis/web/software/pear/doc.html)

. bedtools (https://bedtools.readthedocs.io/en/latest/)


# Usage
```
Usage: sh Final_intergrated.sh -a Sample ID -b Work Directory -c Forward Index -d Reverse Index -f forward reads -r reverse reads -p Software installed Directory

Request Parameters:
	-a Sample Id (Example: yYY398-B_S10)
	-b The working directory, where the raw read stored and we perform the analyses of the large insertion
	-c Index of forward reads(Example: CTC)
	-d Index of reverse reads (Example: ACC)
	-f forward reads (Example: SampleID_L001_R1_001.fastq)
	-r reverse reads (Example: SampleID_L001_R2_001.fastq)
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
For more detail information, please feel free to contact: xin.wang@childrens.harvard.edu
