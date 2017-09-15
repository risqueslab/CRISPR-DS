# CRISPR-DS
Data processing pipeline for CRISPR-Duplex Sequencing data. This pipeline uses
a bash shell script, DS_PE_Unified.3.03.sh, to process raw FastQ files from Illumina
Sequencing platforms to calling variants from double strand consensus sequences. 
Note: Current pipeline is customized for sequencing gene, TP53, however this can
easily be customized to any genomic target of interest. 
An overview of the process:

![here](https://github.com/risqueslab/CRISPR-DS/blob/master/CRISPR-DS_data_processing.png)



# Consensus Making
This data processing pipeline relies on UnifiedConsensusMaker.py to create single strand and double strand consensus sequences using raw FastQ files from Illumina sequencing platforms. Documentation on the usage and dependencies of this script can be found here: https://github.com/loeblab/Duplex-Sequencing

# Dependencies for DS_PE_Unified.3.03.sh
* All dependencies for UnifiedConsensusMaker.py, see link above.
* bwa-0.7.15
* GATK 3.6
* Picard 2.5.0
* fgbio 0.2.0

# Input:
Raw, de-multiplexed, paired-end FastQ files from Illumina platform. (i.e. Sample.seq1.fastq, Sample.seq2.fastq)

# Output:
* Single consensus sequence (SSCS) FastQ 
* Duplex consensus sequence (DCS) FastQ
* SSCS aligned BAM file
* DCS aligned BAM file 
* Mutpos file:  
**Example format:**    
Chromosome  Ref_Base  Ref_Position  #DCS_Reads  #A_Muts #T_Muts #C_Muts #G_Muts #Ns
* Several other intermediate and supplementary QC files

# Usage:
The shell script has 3 parts that can be customized for each analysis run. 
(1) Run variables
(2) File paths
(3) Samples to analyze  
Each script then becomes a record of how the data was processed.  
  
  (1) Set Run Variables:  
  minMem=3            # Minimum number of reads to reach consensus  
  maxMem=200          # Maximum number of reads to reach consesnsus  
  cutOff=0.7          # % of nucleotides at a position in read that must be identical in order   
                      for consensus at that position  
  nCutOff=1           # Maximum fraction of Ns allowed in a consensus  
  tagLen=10           # Adapter sequence length  
  spacerLen=1         # Spacer sequence length  
  readLen=300         # Sequencer read length  
  clipBegin=7         # Number of bases to clip off of beginning of reads  
  clipEnd=30          # Number of bases to clip off of end of reads  
  
  (2) Set File Locations and Paths
  DS_PATH=/Users/RRisques/Desktop/Duplex_Sequencing
  PICARD_PATH=/Applications/Utilities/Seq_Analysis_Tools/picard-tools-2.2.1
  ALIGN_REF=/Users/RRisques/Desktop/Duplex_Sequencing/Reference/Human_Genome/human_g1k_v37.fasta
  GATK_PATH=/Applications/Utilities/Seq_Analysis_Tools/GenomeAnalysisTK-3.6
  FGBIO_PATH=/Applications/Utilities/Seq_Analysis_Tools/fgbio/target/scala-2.12
  REGION_BED=/Users/RRisques/Desktop/Duplex_Sequencing/Reference/P53/p53region.cutsites.bed  # Bed file with genomic region  
  REF_PATH=/Users/RRisques/Desktop/Duplex_Sequencing/Reference/P53  
    
  (3) Set Samples to Analyze
  Folder names separated by spaces containing R1 and R2 Fastq.gz files named: samplename.seq1.fastq.gz and  
  samplename.seq2.fastq.gz  
  NOTE: This script can use compressed Fastq files, no need to unzip files before running this script  
  
Run Script:
  From the terminal-  
  >> cd into the directory containing your SAMPLE FOLDERS and copy of this script  
  >> bash -x DS_PE_Unified.3.0.3.sh 2> DS_PE_Unified_Record.se   

  
Read-me file incompleted, pending completion date: September 12, 2017  
