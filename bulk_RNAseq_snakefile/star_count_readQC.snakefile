# Pipeline for trimming, STAR alignment, count and QC for fastq and BAM files

# Command for dry run: snakemake --snakefile star_count_readQC.snakefile -n -p -j 20
# Command for dag: snakemake --snakefile star_count_readQC.snakefile --dag | dot -Tpdf > dag.pdf
# Command to execute: snakemake --snakefile star_count_readQC.snakefile -j 10
# Command for run report: snakemake --snakefile star_count_readQC.snakefile --report report.html

os.system("mkdir -p trimOutput")
os.system("mkdir -p starOutput")
os.system("mkdir -p fastqcOutput")
os.system("mkdir -p rseqcOutput")

vds2Dir = "/data/home/lingzili/VDS2/NGS_Runs/191119_A00316_0110_BHHGTYDRXX/Data/Intensities/BaseCalls/Demultiplex/Kari_SM"

# Extract sample ID from pair-end fastq files
import os

extractList = []

for fileName in os.listdir(vds2Dir):
	extractList.append('_'.join(fileName.split("_")[0:2]))

# Get unique sample ID
ID = list(set(extractList))

rule all:
    input:
	    expand("starOutput/{sample}_Aligned.sortedByCoord.out.bam.bai", sample = ID),
	    expand("rseqcOutput/{sample}.rRNA.txt", sample = ID),
	    "fastqcOutput/multiqc_report.html",
	    "starOutput/featureCount_19112019"
	      
################################################################
## Read trimming
################################################################
rule trim:
	input:
		vds2Dir + "/{sample}_R1_001.fastq.gz", 
		vds2Dir + "/{sample}_R2_001.fastq.gz"
	output:
		"trimOutput/{sample}_R1_001_val_1.fq.gz",
		"trimOutput/{sample}_R2_001_val_2.fq.gz",
		"trimOutput/{sample}_R1_001.fastq.gz_trimming_report.txt",
		"trimOutput/{sample}_R2_001.fastq.gz_trimming_report.txt",
		"trimOutput/{sample}_R1_001_val_1_fastqc.zip",
		"trimOutput/{sample}_R2_001_val_2_fastqc.zip"
	message:"Trimming reads on sample {wildcards.sample}"
	threads: 6
	shell:
		"trim_galore --cores {threads} --paired --fastqc -o trimOutput {input}"

################################################################
## STAR alignment
################################################################
rule align_STAR:
	input:
		read1 = "trimOutput/{sample}_R1_001_val_1.fq.gz", 
		read2 = "trimOutput/{sample}_R2_001_val_2.fq.gz",
		STARindex = "/data/Genomes/mouse/mm10/star/index101bp/"
	output:
		"starOutput/{sample}_Aligned.sortedByCoord.out.bam",
		"starOutput/{sample}_Log.final.out"
	params:
		outprefix = "{sample}" + "_"
	message: "Running STAR Alignment on {wildcards.sample}"
	threads: 10
	shell:
		"/data/home/lingzili/.conda/envs/lingzili/bin/STAR --runThreadN {threads} \
		--genomeDir {input.STARindex} \
		--readFilesIn {input.read1} {input.read2} \
		--readFilesCommand zcat \
		--outFileNamePrefix ./starOutput/{params.outprefix} \
		--outSAMtype BAM SortedByCoordinate"

rule index_bam:
	input:
		"starOutput/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		"starOutput/{sample}_Aligned.sortedByCoord.out.bam.bai"
	message: "Indexing BAM files of {wildcards.sample}"
	threads: 10
	shell:
		"samtools index -@ {threads} {input} {output}"

################################################################
## featureCounts
################################################################
rule sort_bam_byName:
	input:
		"starOutput/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		"starOutput/{sample}.bam"
	message: "Sorting BAM files of {wildcards.sample} by name"
	threads: 10
	shell:
		"samtools sort -n -@ {threads} -o {output} {input}"

rule featureCounts:
	input:
		expand("starOutput/{sample}.bam", sample = ID)
	output:
		"starOutput/featureCount_19112019"
	params:
		"/data/Genomes/mouse/mm10/mm10.refGene.gtf"
	message: "featureCounts"
	threads: 10
	shell:
		"featureCounts -T {threads} --verbose -p -a {params} -o {output} {input}"

################################################################
## rRNA count
################################################################
rule rRNA_count:
	input:
		bam = "starOutput/{sample}_Aligned.sortedByCoord.out.bam",
		indexed_bam = "starOutput/{sample}_Aligned.sortedByCoord.out.bam.bai"
	output:
		"rseqcOutput/{sample}.rRNA.txt"
	params:
		split_bam = "/data/home/lingzili/genomeTools/RSeQC-3.0.0/scripts/split_bam.py",
		rRNA_bed = "/data/home/lingzili/genomeTools/RSeQC-3.0.0/genomeBED/mm10_rRNA.bed",
		outprefix = "{sample}" + ".rRNA"
	message:
		"Count ribosomal RNA on {wildcards.sample}"
	shell:
		"python {params.split_bam} -i {input.bam} -r {params.rRNA_bed} -o {params.outprefix} 1>{output}"

################################################################
## Read distribution
################################################################
rule read_distribution:
    input:
        "starOutput/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        "rseqcOutput/{sample}.readDist.txt"
    params:
	    read_distribution = "/data/home/lingzili/genomeTools/RSeQC-3.0.0/scripts/read_distribution.py",
	    RefBED="/data/home/lingzili/genomeTools/RSeQC-3.0.0/genomeBED/mm10_RefSeq.bed"
    message: 
	    "Running RseQC read distribution on {wildcards.sample}"
    shell:
        "python {params.read_distribution} -i {input} -r {params.RefBED} 1>{output}"

################################################################
## Insert size
################################################################
rule collect_insert_size:
    input:
        "starOutput/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        txt = "rseqcOutput/{sample}_size_metrics.txt",
        pdf = "rseqcOutput/{sample}_histogram.pdf"
    message: "Collecting insert size for {wildcards.sample}"
    shell:
        "java -jar /data/home/lingzili/miniconda3/share/picard-2.21.1-0/picard.jar CollectInsertSizeMetrics \
		I={input} O={output.txt} H={output.pdf} INCLUDE_DUPLICATES=true HISTOGRAM_WIDTH=800 M=0.1"

################################################################
## Gene body coverage
################################################################
rule geneBody_coverage: 
    input: 
	    bam = "starOutput/",
	    indexed_bam = expand("starOutput/{sample}_Aligned.sortedByCoord.out.bam.bai", sample = ID)
    output:
        "rseqcOutput/pulldown_19112019.geneBodyCoverage.txt"
    params:
	    geneBody_coverage = "/data/home/lingzili/genomeTools/RSeQC-3.0.0/scripts/geneBody_coverage.py",
	    housekeepBED = "/data/home/lingzili/genomeTools/RSeQC-3.0.0/genomeBED/mm10.HouseKeepingGenes.bed",
	    outprefix = "rseqcOutput/pulldown_19112019"
    message: "Gene body coverage"
	shell:
		"python {params.geneBody_coverage} -r {params.housekeepBED} -i {input.bam} -o {params.outprefix}"

################################################################
## MultiQC
################################################################
rule multiqc: 
    input:
	    "rseqcOutput/pulldown_19112019.geneBodyCoverage.txt",
	    expand(["trimOutput/{sample}_R1_001.fastq.gz_trimming_report.txt", "trimOutput/{sample}_R2_001.fastq.gz_trimming_report.txt", 
		"trimOutput/{sample}_R1_001_val_1_fastqc.zip", "trimOutput/{sample}_R2_001_val_2_fastqc.zip",
		"starOutput/{sample}_Log.final.out", "rseqcOutput/{sample}.readDist.txt", "rseqcOutput/{sample}_size_metrics.txt"], sample = ID)
    output:
        "fastqcOutput/multiqc_report.html"
    message: "MultiQC"
	shell:
		"/usr/local/bin/multiqc -o fastqcOutput {input}"
