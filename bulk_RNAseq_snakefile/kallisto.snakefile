# Command for dry run: snakemake --snakefile kallisto.snakefile -n -p -j 20
# Command for dag: snakemake --snakefile kallisto.snakefile --dag | dot -Tpdf > dag.pdf
# Command to execute: snakemake --snakefile kallisto.snakefile -j 20
# Command for run report: snakemake --snakefile kallisto.snakefile --report report.html

trimDir = "/data/home/lingzili/Pulldown_Analysis/19112019/trimOutput/"

# Extract sample ID from pair-end fastq files
import os

extractList = []

for fileName in os.listdir(trimDir):
	extractList.append('_'.join(fileName.split("_")[0:2]))

# Get unique sample ID
ID = list(set(extractList))

rule all:
    input:
	    "kallistoOutput/multiqc_report.html"

################################################################
## kallisto
################################################################
rule kallisto:
	input:
		read1 = "trimOutput/{sample}_R1_001_val_1.fq.gz", 
		read2 = "trimOutput/{sample}_R2_001_val_2.fq.gz"
	params:
		cDNA = "/data/home/lingzili/mm10_genome/Mus_musculus/UCSC/GRCm38_kallisto/GRCm38.cdna.all.idx",
		gtf = "/data/home/lingzili/mm10_genome/Mus_musculus/UCSC/GRCm38_kallisto/Mus_musculus.GRCm38.96.gtf.gz",
		outprefix = "kallistoOutput/{sample}"
	output:
		"kallistoOutput/{sample}.log"
	message: "kallisto for quantifying abundances of transcripts on {wildcards.sample}"
	threads: 10
	shell:
		"kallisto quant -t {threads} -i {params.cDNA} --gtf {params.gtf} -o {params.outprefix} -b 100 {input.read1} {input.read2} 2>{output}"

################################################################
## MultiQC
################################################################
rule multiqc: 
    input: 
	    expand(["kallistoOutput/{sample}.log"], sample = ID)
    output:
        "kallistoOutput/multiqc_report.html"
    message: "MultiQC"
	shell:
		"/usr/local/bin/multiqc -o kallistoOutput {input}"