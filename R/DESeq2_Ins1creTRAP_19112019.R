# DESeq2 analysis on mRNA IP sequenced on 19.11.2019
library(tidyverse)
library(egg)
library(grid)
library(here)
library(DESeq2)
library(pheatmap)
library(org.Mm.eg.db)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(AnnotationHub)
library(ensembldb)

# Load data ---------------------------------------------------------------
# Load metadata
metaData <- read.csv("~/mRNA_IP/count/metadata_19112019.csv", row.names = 1)
View(metaData)

# Load count table
data <- read.csv("~/mRNA_IP/count/featureCount_19112019_extraction.txt", row.names = 1, sep = "")

# Order by column names
data <- data %>% dplyr::select(sort(names(.)))

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

# Marker genes ------------------------------------------------------------
# By R plot
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

# By ggplot2
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

# Rename column names
colnames(normalized_counts) <- c(
  "TRAP input", "Ins1creTRAP input 1", "Ins1creTRAP input 2",
  "TRAP IP", "Ins1creTRAP IP 1", "Ins1creTRAP IP 2",
  "TRAP supernatant", "Ins1creTRAP supernatant 1", "Ins1creTRAP supernatant 2"
)

counts <- data.frame(gene = row.names(normalized_counts), normalized_counts)

# Convert to long format
counts_long <- counts %>%
  gather(Sample, Count, 2:10)

counts_long$Sample <- factor(counts_long$Sample, levels = c(
  "Ins1creTRAP.input.1", "Ins1creTRAP.input.2", "TRAP.input",
  "Ins1creTRAP.IP.1", "Ins1creTRAP.IP.2", "TRAP.IP",
  "Ins1creTRAP.supernatant.1", "Ins1creTRAP.supernatant.2", "TRAP.supernatant"
))

# Ins1
Ins1_p1 <- counts_long %>%
  filter(gene == "Ins1") %>%
  ggplot(aes(x = Sample, y = Count, fill = Sample))

Ins1_p2 <- Ins1_p1 +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Ins1", x = NULL, y = "DESeq2-normalized counts") +
  scale_fill_brewer(palette = "Spectral") +
  standard_theme_barplot

Ins1_p2

# Mafa
Mafa_plot <- counts_long %>%
  filter(gene == "Mafa") %>%
  ggplot(aes(x = Sample, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Mafa", x = NULL, y = "DESeq2-normalized counts") +
  scale_fill_brewer(palette = "Spectral") +
  standard_theme_barplot

Mafa_plot

# Gcg
Gcg_plot <- counts_long %>%
  filter(gene == "Gcg") %>%
  ggplot(aes(x = Sample, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Gcg", x = NULL, y = "DESeq2-normalized counts") +
  scale_fill_brewer(palette = "Spectral") +
  standard_theme_barplot

Gcg_plot

# Pnlip
Pnlip_plot <- counts_long %>%
  filter(gene == "Pnlip") %>%
  ggplot(aes(x = Sample, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Pnlip", x = NULL, y = "DESeq2-normalized counts") +
  scale_fill_brewer(palette = "Spectral") +
  standard_theme_barplot

Pnlip_plot

# Count data transformations ----------------------------------------------
# rlog (regularized-logarithm) transformation
rld <- rlog(dds, blind = TRUE)

# Extract the rlog matrix from the object
rld_mat <- assay(rld)

# Compute pairwise correlation values
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap(rld_cor,
  labels_row = metaData$sampleName, labels_col = metaData$sampleName, angle_col = 45,
  legend = FALSE, cellwidth = 25, cellheight = 25
)

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

pca_fixed <- set_panel_size(pca, width = unit(10, "cm"), height = unit(4, "in"))

grid.newpage()
grid.draw(pca_fixed)

# Differential expression analysis ----------------------------------------
# Run DESeq2 differential expression analysis
dds <- DESeq(dds, parallel = TRUE)

# Check the fit of the dispersion estimates
plotDispEsts(dds)

# Inspecting the results table (adjusted p-value < 0.05)
res <- results(dds, contrast = c("sampleType", "IP", "Input"), alpha = 0.05)

# Order by p value
res <- res[order(res$padj), ]
summary(res)

# GO over-representation analysis -----------------------------------------
# Add Ensembl and Entrez IDs
res$entrez <- mapIds(org.Mm.eg.db, keys = row.names(res), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

res$ensembl <- mapIds(org.Mm.eg.db, keys = row.names(res), column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")

head(res)

# Create background dataset for hypergeometric testing using all genes tested for significance in the results
allOE_genes <- as.character(res$ensembl)

# Extract significant results
sigOE <- subset(res, padj < 0.05)

sigOE_genes <- as.character(sigOE$ensembl)

# Run GO enrichment analysis
## BP: Biological Process, MF: Molecular Function, CC: Cellular Component, or “ALL” for all three
ego <- enrichGO(
  gene = sigOE_genes, universe = allOE_genes, keyType = "ENSEMBL", OrgDb = org.Mm.eg.db,
  ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

png("graph/Ins1creTRAP_19112019/GO_BP_03120219.png", width = 600, height = 600, units = "px")
par(mar = c(2, 2, 2, 2))
dotplot(ego, showCategory = 25)
dev.off()

png("graph/Ins1creTRAP_19112019/GO_Clusters_03120219.png", width = 800, height = 800, units = "px")
par(mar = c(2, 2, 2, 2))
emapplot(ego, showCategory = 25)
dev.off()
