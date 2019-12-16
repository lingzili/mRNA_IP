library(edgeR)
library(limma)
library(tidyverse)
library(ggfortify)
library(RColorBrewer)

# Load featureCounts data -------------------------------------------------
counts <- read.delim("~/mRNA_IP/count/primary_featureCount_ext_19112019.txt", row.names = 1, comment.char = "#")
View(counts)

# Rename column names
colnames(counts) <- c("5652_S13", "5653_S14", "5654_S15", "5655_S16", "5656_S17", "5657_S18", "5658_S19", "5659_S20", "5660_S21")

# Create DGEList object ---------------------------------------------------
d0 <- DGEList(counts)
dim(d0)

# Rename samples
colnames(d0) <- c(
  "TRAP input", "Ins1creTRAP input 1", "Ins1creTRAP input 2",
  "TRAP IP", "Ins1creTRAP IP 1", "Ins1creTRAP IP 2",
  "TRAP supernatant", "Ins1creTRAP supernatant 1", "Ins1creTRAP supernatant 2"
)

# Add sample group names
d0$samples$group <- as.factor(c("Input", "Input", "Input", "IP", "IP", "IP", "Supernatant", "Supernatant", "Supernatant"))

d0$samples

# Filter low-expressed genes ----------------------------------------------
# How many genes have zero counts across all nine samples?
table(rowSums(d0$counts == 0) == 9)

# Cut off CPM value less than 1
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop, ]
dim(d) # number of genes left

# Calculate the log2-CPM --------------------------------------------------
cpm_log <- cpm(d, log = TRUE)

# Expression distribution
par(mar = c(12, 4, 2, 2))
fillColor <- brewer.pal(ncol(cpm_log), "Set1")
boxplot(cpm_log, las = 2, col = fillColor, pch = 20, pt.cex = 3, cex = 1)
title(ylab = "log2 CPM")

# Hierarchical clustering with heatmaps -----------------------------------
heatmap(cor(cpm_log), cexRow = 1, cexCol = 1, margins = c(12, 12))

# PCA ---------------------------------------------------------------------
# Define standard plot themne
standard_theme <- theme(
  axis.line = element_line(colour = "black"),
  axis.text.x = element_text(color = "black", size = 16, face = "bold"),
  axis.text.y = element_text(color = "black", size = 16, face = "bold"),
  axis.title.x = element_text(color = "black", size = 18, face = "bold"),
  axis.title.y = element_text(color = "black", size = 18, face = "bold"),
  legend.title = element_blank(),
  legend.text = element_text(color = "black", size = 16, face = "bold"),
  legend.key = element_rect(fill = "white"), # Remove grey background of the legend
  strip.text.x = element_blank(),
  strip.background = element_rect(fill = "white"),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, size = 2),
  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
  plot.title = element_text(color = "black", size = 20, face = "bold")
)

pca <- prcomp(t(cpm_log))

autoplot(pca, data = t(cpm_log)) +
  geom_point(aes(colour = rownames(t(cpm_log))), size = 3, alpha = 2) +
  standard_theme +
  scale_color_brewer(palette = "Set1")

# Plot individual genes ---------------------------------------------------
df_cpm_log <- data.frame(gene = row.names(cpm_log), cpm_log) %>%
  gather(Sample, Count, 2:10)

df_cpm_log$Sample <- factor(df_cpm_log$Sample, levels = c(
  "Ins1creTRAP.input.1", "Ins1creTRAP.input.2", "TRAP.input",
  "Ins1creTRAP.IP.1", "Ins1creTRAP.IP.2", "TRAP.IP",
  "Ins1creTRAP.supernatant.1", "Ins1creTRAP.supernatant.2", "TRAP.supernatant"
))

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
  legend.position = "none",
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, size = 2),
  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
  plot.title = element_text(color = "black", size = 18, face = "bold")
)

# Ins1
df_cpm_log %>%
  filter(gene == "Ins1") %>%
  ggplot(aes(x = Sample, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Ins1", x = NULL, y = "log2 CPM") +
  scale_fill_brewer(palette = "Spectral") +
  standard_theme_barplot

# Gcg
df_cpm_log %>%
  filter(gene == "Gcg") %>%
  ggplot(aes(x = Sample, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Gcg", x = NULL, y = "log2 CPM") +
  scale_fill_brewer(palette = "Spectral") +
  standard_theme_barplot

# Pnlip
df_cpm_log %>%
  filter(gene == "Pnlip") %>%
  ggplot(aes(x = Sample, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Pnlip", x = NULL, y = "log2 CPM") +
  scale_fill_brewer(palette = "Spectral") +
  standard_theme_barplot

# Sample size for RNAseq studies ------------------------------------------
library(RNASeqPower)
