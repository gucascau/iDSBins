# iDSBins
High-throughput identification of large and complex DNA insertions at DNA double strand breaks.

Insertions of mobile elements, mitochondrial DNA and fragements of nuclear chromosomes at DNA double-strand breaks(DSB) threaten genome instability. This genomic instability can rebalance the effects of specific knockout, i.e. DNA2 mutants.

This repository contain a set of simple scripts that carry out the key anaylses for identification of large and complex insertions, as well as small insertion and deletions at DNA double strand breaks in Yeast genome. These toolkits also contain scripts to define the microhomologies and small insertion/deletions at junction sequence for large insertion events, functionally annotate the insertion events, calculate the inserted coverage of rDNA, LTR, evaluate the relationship between large insertions and Rloop, ARS, hotspots, G4, Rereplication ect.  Similar strategies can also apply to the detection of insertion events at double strand break in human genome.

# Availability 
1. Detect the large and complex DNA insertions at DSB sites.
2. Detect the small deletion/insertion at DSB sites.
3. Perform the deduplication and calculate the coverage of each insertion events.
4. Detect the microhomologies and junctions for large insertion events.
5. Functionally annotate the inserted elements.
6. Identify the hotspots of large insertions.

# Dependencies

Perl is used to run the scripts. The following softwares are also required:

. blast (https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

. bbduk (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/)

. pear (https://cme.h-its.org/exelixis/web/software/pear/doc.html)

. bedtools (https://bedtools.readthedocs.io/en/latest/)


# Usage

Usage: sh $0 -a Sample ID -b Work Directory -c Forward Index -d Reverse Index -f forward reads -r reverse reads -p Software installed Directory

