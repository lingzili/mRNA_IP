#!/bin/bash

mkdir rRNA_count

# Count rRNA
split_bam="/data/home/lingzili/genomeTools/RSeQC-3.0.0/scripts/split_bam.py"
rRNA_bed="/data/home/lingzili/genomeTools/RSeQC-3.0.0/genomeBED/mm10_rRNA.bed"

# BAM file should be sorted and indexed
# Parallel count rRNA
ls *.sortedByCoord.out.bam | parallel -j 16 python $split_bam -i {} -r $rRNA_bed -o ./rRNA_count/{}.rRNA ">" ./rRNA_count/{}.rRNA.txt

# Extract information
cd rRNA_count

# Extract rRNA counts
for j in *.rRNA.txt
do 
baseFilename=`basename $j _Aligned.sortedByCoord.out.bam.rRNA.txt`
RRNAREADS=$(head ${j} -n 2 | tail -n 1 | awk -F ":" '{ print $2 }')
echo "$baseFilename RRNAREADS $RRNAREADS" >> rRNA.count.out
done

# Extract total counts
for j in *.rRNA.txt
do 
baseFilename=`basename $j _Aligned.sortedByCoord.out.bam.rRNA.txt`
TOTALREADS=$(head ${j} -n 1 | awk '{ print $3 }')
echo "$baseFilename TOTALREADS $TOTALREADS" >> total.count.out
done
