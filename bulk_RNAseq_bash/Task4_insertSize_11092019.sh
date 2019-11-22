#!/bin/bash
mkdir insertSize

# Insert size distribution of paired-end library
picardPath="/data/home/lingzili/genomeTools/picard.jar"

ls *.sortedByCoord.out.bam | parallel -j 16 java -jar $picardPath CollectInsertSizeMetrics \
      I={} \
      O=./insertSize/{}.size.out \
      H=./insertSize/{}_histogram.pdf \
      INCLUDE_DUPLICATES=true \
      HISTOGRAM_WIDTH=800 \
      M=0.1