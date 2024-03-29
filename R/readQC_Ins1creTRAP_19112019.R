library(tidyverse)
library(here)

# Load read QC data -------------------------------------------------------
RNAseq <- read.csv("~/mRNA_IP/count/readQC_19112019.csv")
View(RNAseq)
colnames(RNAseq)[1] <- "sampleName"

# Rank samples by ID
RNAseq$sampleName <- factor(RNAseq$sampleName, levels = RNAseq$sampleName[order(RNAseq$id)])
RNAseq$sampleName

# Alignment STAR ----------------------------------------------------------
# Gather columns into key-value pairs
align_stat <- RNAseq %>%
  gather(key = "Mapped", value = "Reads", "Uniquely.mapped", "Mapped.to.multiple.loci", "Mapped.to.too.many.loci", "Unmapped.too.short", "Unmapped.other")

# Rank the mapped labels
align_stat$Mapped <- factor(align_stat$Mapped, levels = c("Unmapped.other", "Unmapped.too.short", "Mapped.to.too.many.loci", "Mapped.to.multiple.loci", "Uniquely.mapped"))

# Barplot for alignment reads
# Set standard theme for barplot
standard_theme_barplot <- theme(
  axis.line = element_line(colour = "black"),
  axis.text.x = element_text(color = "black", size = 20, face = "bold", angle = 45, hjust = 1),
  axis.text.y = element_text(color = "black", size = 20, face = "bold"),
  axis.title.x = element_text(color = "black", size = 20, face = "bold"),
  axis.title.y = element_text(color = "black", size = 20, face = "bold"),
  strip.text.x = element_text(color = "black", size = 20, face = "bold"),
  strip.background = element_rect(fill = "white"),
  legend.title = element_blank(),
  legend.text = element_text(color = "black", size = 20, face = "bold"),
  legend.key = element_rect(fill = "white"), # Remove grey background of the legend
  legend.position = "bottom",
  legend.direction = "vertical",
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, size = 2),
  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
  plot.title = element_text(color = "black", size = 18, face = "bold")
)

star_plot <- align_stat %>%
  ggplot(aes(x = sampleName, y = Reads, fill = Mapped)) +
  geom_bar(stat = "identity") +
  labs(title = "Alignment to reference genome", x = NULL, y = "Reads") +
  standard_theme_barplot

star_plot

ggsave(here::here("graph/Ins1creTRAP_19112019/Alignment_star_12122019.png"), star_plot, width = 10, height = 7.5)

# Pseudoalignment kallisto ------------------------------------------------
pseudoalign_p1 <- RNAseq %>%
  gather(key = "Mapped", value = "Reads", "Pseudoaligned", "Not.aligned") %>%
  ggplot(aes(x = sampleName, y = Reads, fill = Mapped))

pseudoalign_p2 <- pseudoalign_p1 +
  geom_bar(stat = "identity") +
  labs(title = "Pseudoalignment with kallisto", x = NULL, y = "Reads") +
  standard_theme_barplot

pseudoalign_p2

ggsave(here::here("graph/Ins1creTRAP_19112019/pseudoalign_15122019.png"), pseudoalign_p2, width = 8, height = 6)

# Mark Duplicates ---------------------------------------------------------
mkDup_p1 <- RNAseq %>%
  gather(key = "Mapped", value = "Reads", "Read.Pair.Unique", "Unpaired.Read.Unique", "Read.Pair.Not.Optical.Duplicates", "Read.Pair.Optical.Duplicates", "Unpaired.Read.Duplicates") %>%
  ggplot(aes(x = sampleName, y = Reads, fill = Mapped))

mkDup_p2 <- mkDup_p1 +
  geom_bar(stat = "identity") +
  labs(title = "MarkDuplicates", x = NULL, y = "% Total reads") +
  standard_theme_barplot

mkDup_p2

ggsave(here::here("graph/Ins1creTRAP_19112019/MarkDuplicates_12122019.png"), mkDup_p2, width = 12, height = 8)

# Load duplicates data after removal of rRNA
RNAseq <- read.csv("~/mRNA_IP/count/rRNA_ex_deduplication.csv")
View(RNAseq)
colnames(RNAseq)[1] <- "sampleName"

# Rank samples by ID
RNAseq$sampleName <- factor(RNAseq$sampleName, levels = RNAseq$sampleName[order(RNAseq$id)])
RNAseq$sampleName

rRNA_ex_Dup <- RNAseq %>%
  gather(key = "Mapped", value = "Reads", "Read.Pair.Unique", "Unpaired.Read.Unique", "Read.Pair.Not.Optical.Duplicates", "Read.Pair.Optical.Duplicates", "Unpaired.Read.Duplicates") %>%
  ggplot(aes(x = sampleName, y = Reads, fill = Mapped)) +
  geom_bar(stat = "identity") +
  labs(title = "MarkDuplicates after removal of rRNA", x = NULL, y = "% Total reads") +
  standard_theme_barplot

ggsave(here::here("graph/Ins1creTRAP_19112019/rRNA_ex_MarkDuplicates_15122019.png"), rRNA_ex_Dup, width = 12, height = 8)

# Ribosomal RNA percentage ------------------------------------------------
rRNA_p1 <- RNAseq %>%
  ggplot(aes(x = sampleName, y = rRNA_percent, fill = sampleType))

rRNA_p2 <- rRNA_p1 +
  geom_bar(stat = "identity") +
  labs(title = "Alignment to ribosomal RNA", x = NULL, y = "% Reads") +
  ylim(0, 40) +
  standard_theme_barplot

rRNA_p2

# Compare rRNA percentage before and after removal of duplicates
rm_rRNA <- RNAseq %>%
  gather(key = "Dup", value = "Percent", "rRNA_percent", "rmDup_rRNA_percent")

# Rank labels
rm_rRNA$Dup <- factor(rm_rRNA$Dup,
  levels = c("rRNA_percent", "rmDup_rRNA_percent"),
  labels = c("Before MarkDuplicates", "After MarkDuplicates")
)
rm_rRNA$Dup

rm_rRNA_plot <- rm_rRNA %>%
  ggplot(aes(x = sampleName, y = Percent, fill = Dup)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Ribosomal RNA percentage", x = NULL, y = "% Total reads") +
  standard_theme_barplot

rm_rRNA_plot

ggsave(here::here("graph/Ins1creTRAP_19112019/rRNA_12122019.png"), rm_rRNA_plot, width = 12, height = 8)

# Read counts after MarkDuplicates and rRNA removal -----------------------
counts_filtered <- RNAseq %>%
  ggplot(aes(x = sampleName, y = rmDup_ex_count, fill = sampleType)) +
  geom_bar(stat = "identity") +
  labs(title = "Read counts after removal of duplicates and rRNA", x = NULL, y = "Read counts") +
  standard_theme_barplot

counts_filtered

ggsave(here::here("graph/Ins1creTRAP_19112019/counts_filtered_12122019.png"), counts_filtered, width = 8, height = 6)

# featureCounts -----------------------------------------------------------
fc <- RNAseq %>%
  gather(key = "Feature", value = "Fragments", "Assigned", "Unassigned_MultiMapping", "Unassigned_NoFeatures", "Unassigned_Ambiguity")

# Rank labels
fc$Feature <- factor(fc$Feature, levels = c("Unassigned_MultiMapping", "Unassigned_NoFeatures", "Unassigned_Ambiguity", "Assigned"))

fc_p2 <- fc %>%
  ggplot(aes(x = sampleName, y = Fragments, fill = Feature)) +
  geom_bar(stat = "identity") +
  labs(title = "featureCounts on primary alignments", x = NULL, y = "Alignments") +
  standard_theme_barplot

fc_p2

ggsave(here::here("graph/Ins1creTRAP_19112019/featureCounts_primary_152019.png"), fc_p2, width = 8, height = 8)

# genebody coverage -------------------------------------------------------
# Load genebody coverage data
genebody <- read.csv("~/mRNA_IP/count/rseqc_gene_body_coverage_plot.csv")
View(genebody)

# Rename column names
colnames(genebody) <- c(
  "Percentile", "TRAP input", "Ins1creTRAP input 1", "Ins1creTRAP input 2",
  "TRAP IP", "Ins1creTRAP IP 1", "Ins1creTRAP IP 2",
  "TRAP supernatant", "Ins1creTRAP supernatant 1", "Ins1creTRAP supernatant 2"
)

# Convert to long format
genebody_long <- genebody %>%
  gather(Sample, Coverage, 2:10)

# Set standard theme for line graph
standard_theme_line <- theme(
  axis.line = element_line(colour = "black"),
  axis.text.x = element_text(color = "black", size = 16, face = "bold"),
  axis.text.y = element_text(color = "black", size = 16, face = "bold"),
  axis.title.x = element_text(color = "black", size = 18, face = "bold"),
  axis.title.y = element_text(color = "black", size = 18, face = "bold"),
  strip.text.x = element_text(color = "black", size = 18, face = "bold"),
  strip.text.y = element_text(color = "black", size = 18, face = "bold"),
  strip.background = element_rect(fill = "white"),
  legend.title = element_blank(),
  legend.text = element_text(color = "black", size = 16, face = "bold"),
  legend.key = element_rect(fill = "white"), # Remove grey background of the legend
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, size = 2),
  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
  plot.title = element_text(color = "black", size = 16, face = "bold")
)

# Line graph for genebody coverage
genebody_p1 <- genebody_long %>%
  ggplot(aes(x = Percentile, y = Coverage, color = Sample))

genebody_p2 <- genebody_p1 +
  geom_line(size = 1) +
  labs(title = "Genebody coverage", x = "Gene body percentile (5'->3')", y = "Coverage") +
  standard_theme_line

genebody_p2

ggsave(here::here("graph/Ins1creTRAP_19112019/genebody_23112019.png"), genebody_p2)

# Read genomic origin -----------------------------------------------------
Origin <- read.csv("~/mRNA_IP/count/rseqc_read_distribution_plot.csv")
View(Origin)

# Rename column names
colnames(Origin)[1] <- "Sample"
colnames(Origin)[3] <- "5' UTR"
colnames(Origin)[4] <- "3' UTR"

# Convert to long format
Origin_long <- Origin %>%
  gather(Type, Percent, 2:12)

# Rank the mapped labels and samples
Origin_long$Type <- factor(Origin_long$Type, levels = c(
  "TSS_up_1kb", "TSS_up_5kb", "TSS_up_10kb",
  "TES_down_1kb", "TES_down_5kb", "TES_down_10kb",
  "Other_intergenic", "Introns", "5' UTR", "3' UTR", "Exons_CDS"
))

Origin_long$Sample <- factor(Origin_long$Sample, levels = c(
  "TRAP input", "Ins1creTRAP input 1", "Ins1creTRAP input 2",
  "TRAP IP", "Ins1creTRAP IP 1", "Ins1creTRAP IP 2",
  "TRAP supernatant", "Ins1creTRAP supernatant 1", "Ins1creTRAP supernatant 2"
))

# Barplot for read genomic origin
origin_p1 <- Origin_long %>%
  ggplot(aes(x = Sample, y = Percent, fill = Type))

origin_p2 <- origin_p1 +
  geom_bar(stat = "identity") +
  labs(title = "Read genomic origins", x = NULL, y = "% Tags") +
  scale_fill_brewer(palette = "Spectral") +
  standard_theme_barplot

origin_p2

ggsave(here::here("graph/Ins1creTRAP_19112019/readOrigin_23112019.png"), origin_p2)

# EstimateLibraryComplexity -----------------------------------------------
# Load data in while loop
sampleList <- c("5657_S18", "5659_S20", "5653_S14", "5656_S17", "5654_S15", "5660_S21", "5652_S13", "5655_S16", "5658_S19")

i <- 1

while (i <= 9) {
  assign(paste0("x", sampleList[i]), read.delim(paste0("~/mRNA_IP/count/", sampleList[i], ".lib.metrics.txt")))
  i <- i + 1
}

# Merge multiple data frames
LibraryComplexity <- merge_recurse(list(x5652_S13, x5653_S14, x5654_S15, x5655_S16, x5656_S17, x5657_S18, x5658_S19, x5659_S20, x5660_S21),
  by = "duplication_group_count"
)

comp_p1 <- LibraryComplexity %>%
  gather(id, dupCount, 2:10) %>%
  filter(duplication_group_count <= 19) %>%
  ggplot(aes(x = duplication_group_count, y = dupCount, color = id)) +
  geom_point(size = 2, alpha = .7)

comp_p2 <- comp_p1 +
  labs(title = "Read duplication counts", x = "Duplication group count", y = "Reads") +
  scale_color_hue(labels = c(
    "TRAP input", "Ins1creTRAP input 1", "Ins1creTRAP input 2",
    "TRAP IP", "Ins1creTRAP IP 1", "Ins1creTRAP IP 2",
    "TRAP supernatant", "Ins1creTRAP supernatant 1", "Ins1creTRAP supernatant 2"
  )) +
  scale_x_continuous(breaks = seq(1, 19, 2)) +
  standard_theme_line

comp_p2
