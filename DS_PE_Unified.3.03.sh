#!/bin/bash

# DS_Unified_PE.sh
#	V3.03
#	Last updated July 2017 by Daniela Nachmanson.
#
#	Historical DS analysis scripts:
#	V1 - Single End
#	V2 - Paired End
#
# Duplex Sequencing analysis pipeline using paired end reads and Scott Kennedy's unified consensus maker.
#	Consensus maker: 	   			 UnifiedConsensusMaker.py
#	Alignment algorithm used:  BWA MEM
#	Realignment:	             GATK IndelRealigner
#	End clipping: 		   			 GATK ClipReads
# Clipping overlap:	   			 fgbio ClipOverlappingReads
#	Pileup creator:		   			 samtools mpileup
#	Filter pileup:		   			 filter_pileup.py
# Create mutpos:	           mut-position.1.33.py (this version prints with indel information) ***
#	Create countMuts file:	   countMuts.py
#	Plot DCS depth by pos:	   plot_depth_by_position.py
# Plot i-Size histogram:	   picard CollectInsertSizeMetrics
#	Plot error rate per cycle: GATK ErrorRatePerCycle and plot_error_by_cycle.py
#	Print read stats:          get-flagstats-unified.sh
#
#	In this version:
#	*** NEW ORDER OF REALIGNMENT AND THEN CLIPPING***
# *** NEW TOOL FOR CLIPPING OVERLAPPING REGIONS ***
#
# Stop on any error inside or outside pipeline or on an unassigned variable.
set -e
set -o pipefail
set -u

# 1. SET RUN VARIABLES
minMem=3            # Minimum number of reads to reach consensus
maxMem=200          # Maximum number of reads to reach consesnsus
cutOff=0.7          # % of nucleotides at a position in read that must be identical in order for consensus at that position
nCutOff=1           # Maximum fraction of Ns allowed in a consensus
tagLen=10           # Adapter sequence length
spacerLen=1         # Spacer sequence length (bases between adapter and target sequence)
readLen=300         # Sequencer read length
clipBegin=7         # Number of bases to clip off of beginning of reads
clipEnd=30          # Number of bases to clip off of end of reads

# 2. SET FILE LOCATIONS AND PATHS
DS_PATH=~/Desktop/Duplex_Sequencing
PACKAGES_PATH=/Applications/Utilities/Seq_Analysis_Tools
PICARD_PATH=${PACKAGES_PATH}/picard-tools-2.2.1
GATK_PATH=${PACKAGES_PATH}GenomeAnalysisTK-3.6
FGBIO_PATH=${PACKAGES_PATH}fgbio/target/scala-2.12
ALIGN_REF=~/Desktop/Duplex_Sequencing/Reference/Human_Genome/human_g1k_v37.fasta
REGION_BED=~/Desktop/Duplex_Sequencing/Reference/P53/p53region.cutsites.bed # Bed file with genomic region
REF_PATH=~/Desktop/Duplex_Sequencing/Reference/P53


# Adjusting variables for script use
finalReadLen=$((readLen-tagLen-spacerLen))
endTrimStart=$(($finalReadLen-$clipEnd+1))

# 3. SET SAMPLES TO ANALYZE:
# Folder names separated by spaces. Each folder contains R1 and R2 Fastq.gz files named: samplename.seq1.fastq.gz and samplename.seq2.fastq.gz
# NOTE: This script can use compressed Fastq files, no need to unzip files before running this script

folderList='sample1 sample2 sample3 sample4 sample5 '

# 4. RUN SCRIPT:
# From the terminal-
# >> cd into the directory containing your SAMPLE FOLDERS and copy of this script
# >> bash -x DS_PE_Unified.3.0.3.sh 2> DS_PE_Unified_Record.se

# Automated work flow follows:

for elmt in $folderList
do
	cd ${elmt}

	# Consensus Maker
  java -jar $PICARD_PATH/picard.jar FastqToSam F1=${elmt}.seq1.fastq.gz F2=${elmt}.seq2.fastq.gz O=/dev/stdout SM=${elmt}|python $DS_PATH/Programs/UnifiedConsensusMaker.py --input /dev/stdin --taglen ${tagLen} --spacerlen ${spacerLen} --write-sscs --prefix ${elmt} --tagstats

	# Align forward and reverse SSCS reads using BWA algorithm MEM
	bwa mem $ALIGN_REF ${elmt}_read1_sscs.fq.gz ${elmt}_read2_sscs.fq.gz|samtools sort -o ${elmt}_mem.sscs.sort.bam -

	# Align forward and reverse DCS reads using BWA algorithm MEM
	bwa mem $ALIGN_REF ${elmt}_read1_dcs.fq.gz ${elmt}_read2_dcs.fq.gz|samtools sort -o ${elmt}_mem.dcs.sort.bam -

	# Index both SSCS and DCS files
	samtools index ${elmt}_mem.sscs.sort.bam
	samtools index ${elmt}_mem.dcs.sort.bam

	#***Currently processing just DCS reads for further analysis***

	# Filter out unmapped reads
	samtools view -F4 ${elmt}_mem.dcs.sort.bam | samtools view -Sb -T $ALIGN_REF - > ${elmt}.dcs.filt.bam
	samtools index ${elmt}.dcs.filt.bam

   # Put read groups on bam file, which is necessary for later tools
	java -jar -Xmx2g $PICARD_PATH/picard.jar AddOrReplaceReadGroups INPUT=${elmt}.dcs.filt.bam OUTPUT=${elmt}.dcs.filt.readgroups.bam RGLB=Eeny RGPL=Meeny RGPU=Miny RGSM=Moe
	samtools index ${elmt}.dcs.filt.readgroups.bam

	# Realign reads to help correct mapping errors with indels
	java -Xmx4g  -jar $GATK_PATH/GenomeAnalysisTK.jar -T RealignerTargetCreator -dfrac 1 -R $ALIGN_REF -I ${elmt}.dcs.filt.readgroups.bam -o ${elmt}.dcs.filt.bam.readgroups.intervals
	java -Xmx4g -jar $GATK_PATH/GenomeAnalysisTK.jar -T IndelRealigner -dfrac 1 -R $ALIGN_REF -I ${elmt}.dcs.filt.readgroups.bam -targetIntervals ${elmt}.dcs.filt.bam.readgroups.intervals -o ${elmt}.dcs.filt.realign.bam

	# End clipping DCS reads.
	# No need to re-index, GATK does it here.
	java -jar -Xmx4g $GATK_PATH/GenomeAnalysisTK.jar -T ClipReads -I ${elmt}.dcs.filt.realign.bam -o ${elmt}.dcs.filt.realign.clipped.bam -R $ALIGN_REF --cyclesToTrim "$endTrimStart"-"$finalReadLen,1"-"$clipBegin" --clipRepresentation SOFTCLIP_BASES

	# Clip overlapping nucleotides in read pairs
	java -jar $FGBIO_PATH/fgbio-0.2.0-SNAPSHOT.jar ClipOverlappingReads -i ${elmt}.dcs.filt.realign.clipped.bam -o ${elmt}.dcs.filt.no_overlap.bam -r $ALIGN_REF

 	# Create pileup. (Note: -Q0 filter is apparently critical for making mpileup work on softclipped reads)
 	samtools mpileup -Q0 -B -A -d 500000 -f $ALIGN_REF ${elmt}.dcs.filt.no_overlap.bam > ${elmt}.dcs.clipped.no_overlap.pileup

	# Create pileup with JUST genomic coordinates of interest
	python $DS_PATH/Programs/filter_pileup.py $REGION_BED ${elmt}.dcs.clipped.no_overlap.pileup ${elmt}.dcs.clipped.no_overlap.region.pileup N

	# By default, no clonality filter, minimum depth of 1 read to be reported. Minimum number of muts to be reported here is 1.
	cat ${elmt}.dcs.clipped.no_overlap.region.pileup | python $DS_PATH/Programs/mut-position.1.33.py -n 1 > ${elmt}.30clip.DCS-muts.txt

	# By default, no clonality filter, minimum depth of 1 read to be reported. Minimum number of muts to be reported here is 0.
	cat ${elmt}.dcs.clipped.no_overlap.region.pileup| python $DS_PATH/Programs/mut-position.1.33.py -d 1 > ${elmt}.DCS.pileup.mutpos

	# Generate statistics files:

	# Plot DCS depth by genomic coordinate
	python $DS_PATH/Programs/Plot/plot_depth_by_position.py ${elmt}.dcs.clipped.no_overlap.region.pileup

  # Plot insert-size histogram (using unfiltered and unclipped data)
	java -jar $PICARD_PATH/picard.jar CollectInsertSizeMetrics I=${elmt}_mem.dcs.sort.bam O=${elmt}.iSize_Metrics.txt H=${elmt}.iSize_Histogram.pdf M=0.5

	# Plot mutations by read cycle
	java -jar $GATK_PATH/GenomeAnalysisTK.jar -T ErrorRatePerCycle -R $ALIGN_REF -dfrac 1 -I ${elmt}.dcs.filt.readgroups.bam -o ${elmt}.ErrorRatePerCycle.txt

	python $DS_PATH/Programs/Plot/plot_error_by_cycle.py ${elmt}.ErrorRatePerCycle.txt

	# Print reads statistics
	bash $DS_PATH/Programs/get-flagstats-unified-CRISPR.sh ${elmt} > ${elmt}.flagstats.stats

	cd ..
done
