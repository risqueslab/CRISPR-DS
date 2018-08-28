# CRISPR-DS
Data processing pipeline for CRISPR-Duplex Sequencing data. This pipeline uses
a bash shell script, DS_PE_Unified.3.03.sh, to process raw FastQ files from Illumina
Sequencing platforms to call variants from double strand consensus sequences.


**Note:** Current pipeline is customized for sequencing gene, *TP53*, however this can
easily be customized to any genomic target of interest.

An overview of the process:

![here](https://github.com/risqueslab/CRISPR-DS/blob/master/media/CRISPR-DS_data_processing.png)



# Consensus Making
This data processing pipeline relies on UnifiedConsensusMaker.py to create single strand and double strand consensus sequences using raw FastQ files from Illumina sequencing platforms.

Further documentation on the algorithm of UnifiedConsensusMaker.py can be found here: https://github.com/loeblab/Duplex-Sequencing

# Dependencies for DS_PE_Unified.3.03.sh

|Package| Version written with|
|:---:|:---|
|bwa     |0.7.15|
|GATK    |3.6|
|picard  |2.2.1|
|fgbio   |0.2.0|
|python  |2.7.X|
|samtools|>=1.3.1|
|pysam	 |>=0.9.0|
|matplotlib |	>=1.5.1 (optional)|


# Input:
Raw, de-multiplexed, paired-end FastQ files from Illumina platform. (i.e. sample.seq1.fastq, sample.seq2.fastq)

# Output:
* Single strand consensus sequence (SSCS) FastQ
* Double strand consensus sequence (DCS) FastQ
* SSCS aligned BAM file
* DCS aligned BAM file
* Mutpos file:  


> **Example format:**

> |Chr  |Ref  |Pos     |#DCS_Reads|#Muts|#T_Muts|#C_Muts|#G_Muts|#A_muts|#Ins |#Del |#Ns  |
|:---:|:---:|:---:   |:---:     |:---:|:---:  |:---:  |:---:  |:---:  |:---:|:---:|:---:|
|chr1 |A    |17000105|30000     |1    |0      |1      |0      |0      |0    |0    |100  |

* Several other intermediate and supplementary QC files

# Usage:
The shell script has 3 parts that can be customized for each analysis run.
1. Run variables
2. File paths
3. Samples to analyze  

Each script then becomes a record of how the data was processed.  

---
1. Set Run Variables:

```
  minMem=3            # Minimum number of reads to reach consensus      
  maxMem=200          # Maximum number of reads to reach consensus  
  cutOff=0.7          # Fraction of nucleotides at a position in read that must be identical in order   
                      for consensus at that position  
  nCutOff=1           # Maximum fraction of Ns allowed in a consensus  
  tagLen=10           # Adapter sequence length  
  spacerLen=1         # Spacer sequence length  
  readLen=300         # Sequencer read length  
  clipBegin=7         # Number of bases to clip off of beginning of reads  
  clipEnd=30          # Number of bases to clip off of end of reads  
```


  (2) Set File Locations and Paths
  ```
  DS_PATH=~/Duplex_Sequencing
  PICARD_PATH=~/picard-tools-2.2.1
  ALIGN_REF=~/human_g1k_v37.fasta
  GATK_PATH=~/GenomeAnalysisTK-3.6
  FGBIO_PATH=~/fgbio/target/scala-2.12
  REGION_BED=~/p53region.cutsites.bed   
```

  (3) Set Samples to Analyze
  Folder names separated by spaces containing read 1 and read 2 Fastq.gz files named:
  ```
  samplename.seq1.fastq.gz
  samplename.seq2.fastq.gz
  ```  
**Note:** This script can use compressed Fastq files, no need to unzip files before running this script  

Run Script:
```  
cd into the directory containing sample folders and copy of this script  

bash -x DS_PE_Unified.3.0.3.sh 2> DS_PE_Unified_Record.se   
```
