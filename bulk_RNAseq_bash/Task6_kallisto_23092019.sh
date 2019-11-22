#!/bin/bash

# Align pair-end reads with kallisto index from mouse GRCm38.p6 cDNA
# Quantify abundances of the transcripts with kallisto (default running mode is paired-end)
ls *.fq.gz | cut -d"_" -f 1,2 | sort -u | while read id
do 
echo "***The sample ID is $id***"
	read1=${id}_R1_001_val_1.fq.gz
	read2=${id}_R2_001_val_2.fq.gz
echo "***Read 1 is ${read1} and read 2 is ${read2}***"
kallisto quant -t 32 \
-i /data/home/lingzili/mm10_genome/Mus_musculus/UCSC/GRCm38_kallisto/GRCm38.cdna.all.idx \
--gtf /data/home/lingzili/mm10_genome/Mus_musculus/UCSC/GRCm38_kallisto/Mus_musculus.GRCm38.96.gtf.gz \
-o ${id} \
-b 100 \
$read1 $read2
echo "***${id} is done.***"
done

# Change directory name
for i in *_S*
do baseFilename=`cut -d"_" -f 1 $i`
echo "New directory name is ${baseFilename}"
mv $i $baseFilename
done 

for i in *_S*
do baseFilename=`cut -d"_" -f 1 $i`
echo "New directory name is ${baseFilename}"
done 



ls *_S* | while read id
do
echo "The old directory name is $id"

done