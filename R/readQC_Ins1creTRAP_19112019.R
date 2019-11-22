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
  axis.text.x = element_text(color = "black", size = 12, face = "bold", angle = 45, hjust = 1),
  axis.text.y = element_text(color = "black", size = 16, face = "bold"),
  axis.title.x = element_text(color = "black", size = 18, face = "bold"),
  axis.title.y = element_text(color = "black", size = 18, face = "bold"),
  strip.text.x = element_text(color = "black", size = 18, face = "bold"),
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

star_plot <- align_stat %>%
  ggplot(aes(x = sampleName, y = Reads, fill = Mapped)) +
  geom_bar(stat = "identity") +
  labs(title = "Alignment to reference genome", x = NULL, y = "Reads") +
  standard_theme_barplot

star_plot

ggsave(here::here("graph/Ins1creTRAP_19112019/Alignment_star_23112019.png"), star_plot)

# Pseudoalignment kallisto ------------------------------------------------
pseudoalign_p1 <- RNAseq %>%
  gather(key = "Mapped", value = "Reads", "Pseudoaligned", "Not.aligned") %>%
  ggplot(aes(x = sampleName, y = Reads, fill = Mapped))

pseudoalign_p2 <- pseudoalign_p1 +
  geom_bar(stat = "identity") +
  labs(title = "Pseudoalignment with kallisto", x = NULL, y = "Reads") +
  standard_theme_barplot

pseudoalign_p2

ggsave(here::here("graph/Ins1creTRAP_19112019/pseudoalign.png"), pseudoalign_p2, width = 7.86, height = 5.52)

# Ribosomal RNA percentage ------------------------------------------------
rRNA_p1 <- RNAseq %>%
  ggplot(aes(x = sampleName, y = rRNA_percent, fill = sampleType))

rRNA_p2 <- rRNA_p1 +
  geom_bar(stat = "identity") +
  labs(title = "Alignment to ribosomal RNA", x = NULL, y = "% Reads") +
  ylim(0, 40) +
  standard_theme_barplot

rRNA_p2

ggsave(here::here("graph/Ins1creTRAP_19112019/rRNA_23112019.png"), rRNA_p2)

# genebody coverage -------------------------------------------------------
##### Start from HERE


# Convert to long format
genebody_long <- genebody %>%
  gather(Loci, Coverage, 1:100)

# Change Loci to numeric
genebody_long$Loci <- gsub("V", "", genebody_long$Loci)

genebody_long$Loci <- as.numeric(genebody_long$Loci)

# Rank the sample
genebody_long$V101 <- factor(genebody_long$V101, levels = c("SM5169", "SM5165", "SM5170", "SM5166", "SM4996", "SM4998"))

# Set standard theme for line graph
standard_theme_facet_line <- theme(
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
  ggplot(aes(x = Loci, y = Coverage, color = V101, group = V101))

genebody_p2 <- genebody_p1 +
  geom_line(size = 1) +
  labs(title = "Genebody coverage", x = "Gene body percentile (5'->3')", y = "Coverage") +
  scale_color_manual(values = c("#FF9999", "#00AFBB", "#FF9999", "#00AFBB", "#FF9999", "#00AFBB"), labels = c("SM4996" = "GcgcreTRAP 2 Input", "SM4998" = "GcgcreTRAP 2 Pulldown", "SM5169" = "TRAP Input", "SM5165" = "TRAP Pulldown", "SM5170" = "GcgcreTRAP 1 Input", "SM5166" = "GcgcreTRAP 1 Pulldown")) +
  standard_theme_facet_line

genebody_p2

ggsave(here::here("graph/genebody_09092019.png"), genebody_p2, width = 8, height = 6)

# Read genomic origin -----------------------------------------------------
Origin <- read_excel("data/Read_Distribution.xlsx")
View(Origin)

# Set standard theme for barplot
standard_theme_barplot <- theme(
  axis.line = element_line(colour = "black"),
  axis.text.x = element_text(color = "black", size = 16, face = "bold", angle = 45, hjust = 1),
  axis.text.y = element_text(color = "black", size = 16, face = "bold"),
  axis.title.x = element_text(color = "black", size = 18, face = "bold"),
  axis.title.y = element_text(color = "black", size = 18, face = "bold"),
  strip.text.x = element_text(color = "black", size = 18, face = "bold"),
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

# Rank the sample
Origin$Sample <- factor(Origin$Sample, levels = c("SM5169", "SM5165", "SM5170", "SM5166", "SM4996", "SM4998", "SM4997"))

# Barplot for read genomic origin
origin_p1 <- Origin %>%
  ggplot(aes(x = Sample, y = Tag_Percent, fill = Origin))

origin_p2 <- origin_p1 +
  geom_bar(stat = "identity") +
  labs(title = "Read genomic origins", x = NULL, y = "% Tags") +
  scale_x_discrete(labels = c("SM4996" = "GcgcreTRAP 2 Input", "SM4998" = "GcgcreTRAP 2 Pulldown", "SM5169" = "TRAP Input", "SM5165" = "TRAP Pulldown", "SM5170" = "GcgcreTRAP 1 Input", "SM5166" = "GcgcreTRAP 1 Pulldown")) +
  ylim(0, 100) +
  scale_fill_manual(values = c("red", "#CC6666", "pink", "blue", "green", "darkgreen", "lightgreen", "yellow", "grey", "black")) +
  standard_theme_barplot

origin_p2

ggsave(here::here("graph/readOrigin_09092019.png"), origin_p2, width = 8, height = 6)
