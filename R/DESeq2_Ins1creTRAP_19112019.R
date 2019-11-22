# DESeq2 analysis on mRNA IP sequenced on 19.11.2019
library(tidyverse)
library(here)
library(DESeq2)
library(pheatmap)

# Load data ---------------------------------------------------------------
# Load metadata
metaData <- read.csv("~/RNA_Pulldown/data/metadata_19112019.csv", row.names = 1)
View(metaData)

# Load count table
data <- read.csv("~/RNA_Pulldown/data/featureCount_19112019_extraction.txt", row.names = 1, sep = "")

# Order by column names
data <- data %>% select(sort(names(.)))

# Across all columns, remove "X"
names(data) <- gsub("X", "", names(data))
View(data)

# Create dds --------------------------------------------------------------
# Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = data, colData = metaData, design = ~sampleType)

# Normalize samples
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# rlog (regularized-logarithm) transformation
rld <- rlog(dds, blind = TRUE)

# PCA plot ----------------------------------------------------------------
# Define standard plot themne
standard_theme <- theme(
  axis.line = element_line(colour = "black"),
  axis.text.x = element_text(color = "black", size = 16, face = "bold"),
  axis.text.y = element_text(color = "black", size = 16, face = "bold"),
  axis.title.x = element_text(color = "black", size = 18, face = "bold"),
  axis.title.y = element_text(color = "black", size = 18, face = "bold"),
  legend.title = element_blank(),
  legend.text = element_text(color = "black", size = 12, face = "bold"),
  legend.key = element_rect(fill = "white"), # Remove grey background of the legend
  strip.text.x = element_blank(),
  strip.background = element_rect(fill = "white"),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, size = 2),
  plot.title = element_text(color = "black", size = 20, face = "bold")
)

pca <- plotPCA(rld, intgroup = "sampleName") + standard_theme

pca

ggsave(here::here("graph/Ins1creTRAP_PCA_19112019.png"), pca)

# Extract the rlog matrix from the object
rld_mat <- assay(rld)

# Compute pairwise correlation values
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap(rld_cor, labels_row = metadata$sampleName, labels_col = metadata$sampleName)

dds <- DESeq(dds, parallel = TRUE)

# Inspecting the results table
res <- results(dds)

# Order by p value
res <- res[order(res$padj), ]
summary(res)

png("Ins1_22112019.png")
par(mar = c(8, 8, 6, 6), las = 2, cex.lab = 1.5, cex.axis = .8, cex.main = 1.5, cex.sub = .8)
plotCounts(dds, gene = "Ins1", intgroup = "sampleName", xlab = "")
dev.off()

png("Mafa_22112019.png")
par(mar = c(8, 8, 6, 6), las = 2, cex.lab = 1.5, cex.axis = .8, cex.main = 1.5, cex.sub = .8)
plotCounts(dds, gene = "Mafa", intgroup = "sampleName", xlab = "")
dev.off()

png("Gcg_22112019.png")
par(mar = c(8, 8, 6, 6), las = 2, cex.lab = 1.5, cex.axis = .8, cex.main = 1.5, cex.sub = .8)
plotCounts(dds, gene = "Gcg", intgroup = "sampleName", xlab = "")
dev.off()

png("Pnlip_22112019.png")
par(mar = c(8, 8, 6, 6), las = 2, cex.lab = 1.5, cex.axis = .8, cex.main = 1.5, cex.sub = .8)
plotCounts(dds, gene = "Pnlip", intgroup = "sampleName", xlab = "")
dev.off()