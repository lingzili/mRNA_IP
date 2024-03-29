# Load library
library(tidyverse)
library(sleuth)
library(biomaRt)
library(here)

# Prepare design table ----------------------------------------------------
# Define file paths for kallisto directories
getwd()

sample_id <- c("5652_S13", "5653_S14", "5654_S15", "5655_S16", "5656_S17", "5657_S18", "5658_S19", "5659_S20", "5660_S21")

paths <- list(
  "C:/Users/lingzili/Documents/mRNA_IP/kallisto/Ins1creTRAP_19112019/5652_S13",
  "C:/Users/lingzili/Documents/mRNA_IP/kallisto/Ins1creTRAP_19112019/5653_S14",
  "C:/Users/lingzili/Documents/mRNA_IP/kallisto/Ins1creTRAP_19112019/5654_S15",
  "C:/Users/lingzili/Documents/mRNA_IP/kallisto/Ins1creTRAP_19112019/5655_S16",
  "C:/Users/lingzili/Documents/mRNA_IP/kallisto/Ins1creTRAP_19112019/5656_S17",
  "C:/Users/lingzili/Documents/mRNA_IP/kallisto/Ins1creTRAP_19112019/5657_S18",
  "C:/Users/lingzili/Documents/mRNA_IP/kallisto/Ins1creTRAP_19112019/5658_S19",
  "C:/Users/lingzili/Documents/mRNA_IP/kallisto/Ins1creTRAP_19112019/5659_S20",
  "C:/Users/lingzili/Documents/mRNA_IP/kallisto/Ins1creTRAP_19112019/5660_S21"
)

# Add sample names to file paths
names(paths) <- sample_id

# Load experimental design
s2c <- read.csv("C:/Users/lingzili/Documents/mRNA_IP/kallisto/Ins1creTRAP_19112019/sample_info.csv")
s2c

# Correct the name of the first column
names(s2c)[1] <- "sample"
names(s2c)

# Add file path to experimental design
s2c <- mutate(s2c, path = paths)
s2c[] <- lapply(s2c, as.character)
s2c

# Add gene names ----------------------------------------------------------
# Load mouse gene names from Ensembl
mart <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)

# Rename the columns
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

head(t2g)

# Add to the sleuth object
so <- sleuth_prep(s2c, target_mapping = t2g, extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)

# kallisto abundance table ------------------------------------------------
# Get kallisto abundance table from sleuth object
count_table <- kallisto_table(so, use_filtered = FALSE, normalized = FALSE, include_covariates = TRUE)
head(count_table)

# Export kallisto abundance table
write.csv(count_table, file = "kallisto_count_19112019.csv", row.names = FALSE)

# kallisto count distribution
# Choose the sample with most uniquely mapped read:5659_S20
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

geneCount_p1 <- count_table %>%
  filter(id == "5659_S20") %>%
  ggplot(aes(x = est_counts))

geneCount_p2 <- geneCount_p1 +
  geom_histogram(stat = "bin", bins = 100) +
  xlim(-5, 300) +
  labs(title = "Estimated count distribution", x = "Estimated transcript counts from kallisto", y = "Number of transcripts") +
  standard_theme_barplot

geneCount_p2

ggsave(here::here("graph/Ins1creTRAP_kallisto_distribution_22112019.png"), geneCount_p2)

# PCA ---------------------------------------------------------------------
# Define standard plot themne
standard_theme <- theme(
  axis.line = element_line(colour = "black"),
  axis.text.x = element_text(color = "black", size = 16, face = "bold"),
  axis.text.y = element_text(color = "black", size = 16, face = "bold"),
  axis.title.x = element_text(color = "black", size = 18, face = "bold"),
  axis.title.y = element_text(color = "black", size = 18, face = "bold"),
  legend.title = element_blank(),
  legend.text = element_text(color = "black", size = 18, face = "bold"),
  legend.key = element_rect(fill = "white"), # Remove grey background of the legend
  strip.text.x = element_blank(),
  strip.background = element_rect(fill = "white"),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, size = 2),
  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
  plot.title = element_text(color = "black", size = 20, face = "bold")
)

# Calculate PC variance
pc_variance <- plot_pc_variance(so)

list_variance <- pc_variance$data$var

# PCA plot
pca_p1 <- plot_pca(so, color_by = "sample", text_labels = FALSE, point_size = 4, point_alpha = 0.7)

pca_p2 <- pca_p1 +
  standard_theme +
  xlab(paste0("PC1: ", format(round(list_variance[1], 2), nsmall = 2), "% variance")) +
  ylab(paste0("PC2: ", format(round(list_variance[2], 2), nsmall = 2), "% variance")) +
  scale_color_brewer(palette = "Set1")

pca_p2

ggsave(here::here("graph/Ins1creTRAP_PCA_22112019.png"), pca_p2)

# Sample heatmap ----------------------------------------------------------
tpm_table <- count_table %>%
  dplyr::select(target_id, tpm, sample) %>%
  pivot_wider(names_from = sample, values_from = tpm)

colnames(tpm_table) <- factor(colnames(tpm_table), levels = c(
  "target_id", "Ins1creTRAP input 1", "Ins1creTRAP input 2", "TRAP input",
  "Ins1creTRAP IP 1", "Ins1creTRAP IP 2", "TRAP IP",
  "Ins1creTRAP supernatant 1", "Ins1creTRAP supernatant 2", "TRAP supernatant"
))

# Remove rows with zero values
filter_tpm_table <- tpm_table[rowSums(tpm_table[,] == 0) == 0, ]

png("graph/Ins1creTRAP_19112019/sleuth_Cor.png")
heatmap(cor(filter_tpm_table[, -1]), cexRow = 1.5, cexCol = 1.5, margins = c(16, 16))
dev.off()

# Expression distribution
png("graph/Ins1creTRAP_19112019/sleuth_distribution.png")
par(mar = c(12, 4, 2, 2))
fillColor <- brewer.pal(ncol(filter_tpm_table[, -1]), "Set1")
boxplot(filter_tpm_table[, -1], las = 2, col = fillColor, pch = 20, pt.cex = 3, cex = 1)
title(ylab = "TPM")
dev.off()

# Group density -----------------------------------------------------------
# It is assumed that all samples should have a similar range and distribution of expression values.
pgd <- plot_group_density(so, use_filtered = TRUE, trans = "log", grouping = "condition", offset = 1) +
  theme(
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 16, face = "bold"),
    axis.text.y = element_text(color = "black", size = 16, face = "bold"),
    axis.title.x = element_text(color = "black", size = 18, face = "bold"),
    axis.title.y = element_text(color = "black", size = 18, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(color = "black", size = 12, face = "bold"),
    legend.key = element_rect(fill = "white"), # Remove grey background of the legend
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 2),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), legend.position = "top"
  ) +
  labs(x = "log(mean(estimated counts)+1)", y = "Density")

pgd

ggsave(here::here("graph/Ins1creTRAP_Density_22112019.png"), pgd)

# Construct reduced and full model ----------------------------------------
so <- sleuth_fit(so, ~genotype, "reduced")

so <- sleuth_fit(so, ~ genotype + condition, "full")

# Likelihood ratio test (lrt)
so <- sleuth_lrt(so, "reduced", "full")

# Gene-level differential expression results
sleuth_table_gene <- sleuth_results(so, "reduced:full", "lrt", show_all = FALSE)
sleuth_table_gene <- dplyr::filter(sleuth_table_gene, qval <= 0.05)
head(sleuth_table_gene, 20)

# Consistent transcript-level differential expression results
sleuth_table_tx <- sleuth_results(so, "reduced:full", "lrt", show_all = FALSE, pval_aggregate = FALSE)
sleuth_table_tx <- dplyr::filter(sleuth_table_tx, qval <= 0.05)
head(sleuth_table_tx, 20)

## Transcript: Ins1
png("graph/InscreTRAP_Ins1_22110219.png", width = 800, height = 400, units = "px")

plot_bootstrap(so, "ENSMUST00000039652.5", units = "tpm", color_by = "condition") +
  standard_theme +
  labs(title = "Ins1", x = NULL, y = "TPM")

dev.off()

## Transcript: Gcg
png("graph/InscreTRAP_Gcg_22110219.png", width = 800, height = 400, units = "px")

plot_bootstrap(so, "ENSMUST00000102733.9", units = "tpm", color_by = "condition") +
  standard_theme +
  labs(title = "Gcg", x = NULL, y = "TPM")

dev.off()

# Transcript: Mafa
png("graph/InscreTRAP_Mafa_22110219.png", width = 800, height = 400, units = "px")

plot_bootstrap(so, "ENSMUST00000062002.5", units = "tpm", color_by = "condition") +
  standard_theme +
  labs(title = "Mafa", x = NULL, y = "TPM")

dev.off()

## Transcript: Pnlip
png("graph/InscreTRAP_Pnlip_22110219.png", width = 800, height = 400, units = "px")

plot_bootstrap(so, "ENSMUST00000057270.8", units = "tpm", color_by = "condition") +
  standard_theme +
  labs(title = "Pnlip", x = NULL, y = "TPM")

dev.off()

# Fold change -------------------------------------------------------------
so <- sleuth_wt(so, "conditionIP")

wt_results <- sleuth_results(so, "conditionIP", "wt", show_all = TRUE)

# Create a table of significantly differential genes
table(wt_results[, "qval"] < 0.05)

wt_significant <- dplyr::filter(wt_results, qval <= 0.05)

plot_volcano(so, "conditionIP", test_type = "wt", which_model = "full", sig_level = 0.05, point_alpha = 0.2) +
  standard_theme +
  labs(x = "Fold change")

# Launch Shinny app
sleuth_live(so)
