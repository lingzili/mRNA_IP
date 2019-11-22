#!/bin/bash
mkdir readOrigin

# Read genomic distribution
read_distribution="/data/home/lingzili/genomeTools/RSeQC-3.0.0/scripts/read_distribution.py"
RefBED="/data/home/lingzili/genomeTools/RSeQC-3.0.0/genomeBED/mm10_RefSeq.bed"

ls *.sortedByCoord.out.bam | parallel -j 20 python $read_distribution -i {} -r $RefBED ">" ./readOrigin/{}.readDist.txt
