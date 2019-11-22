#!/bin/bash

# Gene body coverage
mkdir geneBodyCoverage

geneBody_coverage="/data/home/lingzili/genomeTools/RSeQC-3.0.0/scripts/geneBody_coverage.py"
housekeepBED="/data/home/lingzili/genomeTools/RSeQC-3.0.0/genomeBED/mm10.HouseKeepingGenes.bed"
prefix_output="SM4588-5170"

python $geneBody_coverage -r $housekeepBED -i /data/home/lingzili/Pulldown_Analysis/trimmedFastqFiles/starOutputTrimmed/ -o ./geneBodyCoverage/$prefix_output