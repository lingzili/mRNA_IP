#!/bin/bash

printf "\n***Description: This script takes fastq files of RNA-Seq data in VDS2 directory and trims them.***\n"

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
mkdir -p trimOutput

printf "=========================================================================================\n"
echo "trim_galore"
printf "=========================================================================================\n"

# Define location of fastq files
vds2Dir=/data/home/lingzili/VDS2/NGS_Runs/190807_A00316_0089_AHCM7JDRXX/Data/Intensities/BaseCalls/Demultiplex/LL_SM

cat GcgcreTRAP_samples.txt | cut -d"_" -f 1,2 | sort -u | while read id
do 
	echo "***The sample ID is $id***"
 
  # Setup output filenames and locations
  read1=/data/home/lingzili/VDS2/NGS_Runs/190807_A00316_0089_AHCM7JDRXX/Data/Intensities/BaseCalls/Demultiplex/LL_SM/${id}_R1_001.fastq.gz
	read2=/data/home/lingzili/VDS2/NGS_Runs/190807_A00316_0089_AHCM7JDRXX/Data/Intensities/BaseCalls/Demultiplex/LL_SM/${id}_R2_001.fastq.gz
  
  echo "***Trimming starts on ${read1} and ${read2}***"
  
  trim_galore --cores $nThread --paired -o trimOutput ${read1} ${read2}
done
