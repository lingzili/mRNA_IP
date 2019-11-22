#!/bin/bash

printf "\n***Description: This script takes trimmed fastq files of RNA-Seq data in the current directory, runs STAR alignment and FastQC, and outputs a counts matrix for all.***\n"

# Exit if an error occurs
set -e

# To debug, display the command before executing it by "set -x"; to turn this option off by "set +x"

# Decide on threads
printf "\n=========================================================================================\n"
read -p "***How many threads do you want to use?*** Answer: " nThread
echo "***Your input thread is $nThread***"
printf "=========================================================================================\n"

# Create directories
# -p option will make sure that mkdir will create the directory only if it does not exist, and it won’t throw an error if it does exist
mkdir -p starOutputTrimmed countOutputTrimmed fastqcOutputTrimmed

printf "=========================================================================================\n"
echo "STAR"
printf "=========================================================================================\n"

# STAR alignment
STARindex=/data/Genomes/mouse/mm10/star/index101bp/

ls *.fq.gz | cut -d"_" -f 1,2 | sort -u | while read id
do 
	echo "***The sample ID is $id***"

	# Set up output filenames and locations
	sortedByCoordBAM=./starOutputTrimmed/${id}_Aligned.sortedByCoord.out.bam
	sortedByNameBAM=./starOutputTrimmed/${id}_Aligned.out.bam
	indexedBAM=./starOutputTrimmed/${id}_Aligned.sortedByCoord.out.bam.bai
	statBAM=./starOutputTrimmed/${id}.stat.txt

	read1=${id}_R1_001_val_1.fq.gz
	read2=${id}_R2_001_val_2.fq.gz

	echo "***STAR alignment of ${read1} and ${read2}***"

	/data/home/lingzili/.conda/envs/lingzili/bin/STAR --genomeDir $STARindex \
	--runThreadN $nThread \
	--readFilesIn $read1 $read2 \
	--readFilesCommand zcat \
	--outFileNamePrefix ./starOutputTrimmed/${id}_ \
	--outBAMsortingThreadN $nThread \
	--outSAMtype BAM Unsorted SortedByCoordinate
	echo "***${id} is aligned***"
	echo "***Note: The unsorted BAM file from STAR is actually sorted by name***"

	echo "***Start index BAM file of ${id}***"
	samtools index -@ $nThread $sortedByCoordBAM $indexedBAM

	echo "***Start extract statistics from BAM file of ${id}***"
	samtools flagstat -@ $nThread $sortedByCoordBAM >$statBAM
done

printf "=========================================================================================\n"
echo "featureCounts"
printf "=========================================================================================\n"

geneAnnotation=/data/Genomes/mouse/mm10/mm10.refGene.gtf

featureCounts -T $nThread --verbose -p -a $geneAnnotation -o ./countOutputTrimmed/featureCount ./starOutputTrimmed/*.sortedByCoord.out.bam

# Quality control with FastQC and MultiQC
printf "=========================================================================================\n"
echo "FastQC"
printf "=========================================================================================\n"

# Disable pop-up display
unset DISPLAY

# Run FastQC on fastq and BAM files
fq=./*.fq.gz
BAM=./starOutputTrimmed/*.sortedByCoord.out.bam

fastqc -t $nThread -o fastqcOutputTrimmed ${fq} ${BAM}

# multiqc to summarize all stats
printf "=========================================================================================\n"
echo "MultiQC"
printf "=========================================================================================\n"

/usr/local/bin/multiqc -o fastqcOutputTrimmed .

# Extract information from count table
echo "***Please remember to use cut to extract information from count table***"