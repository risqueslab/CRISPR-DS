# CRISPR-DS
Data processing pipeline for CRISPR-Duplex Sequencing data. This pipeline uses
a bash shell script, DS_PE_Unified.3.03.sh to process raw FastQ files from Illumina
Sequencing platforms to single strand and double strand FastQ files and bam files. Briefly:



# Consensus Making
This data processing pipeline relies on UnifiedConsensusMaker.py to create single strand and double strand consensus sequences using raw FastQ files from Illumina sequencing platforms. Documentation on the usage and dependencies of this script can be found here: https://github.com/loeblab/Duplex-Sequencing

# Dependencies for CRISPR
* Picard 2.2.1
* samtools 1.3.1
* htslib 1.3.1
* BWA
* X-code (text editor)
* Matplotlib
* GATK
* fgbio


Read-me file incompleted, pending completion date: September 12, 2017
