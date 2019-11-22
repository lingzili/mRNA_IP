# Load library
library(tidyverse)
library(sleuth)
library(biomaRt)
library(here)

# Prepare design table ----------------------------------------------------
# Define file paths for kallisto directories
# Remove SM4998 becauase of its low quality
sample_id <- c("SM4996", "SM5165", "SM5166", "SM5169", "SM5170")

paths <- list(
  "C:/Users/lingzili/Documents/RNA_Pulldown/data/kallisto/SM4996",
  "C:/Users/lingzili/Documents/RNA_Pulldown/data/kallisto/SM5165",
  "C:/Users/lingzili/Documents/RNA_Pulldown/data/kallisto/SM5166",
  "C:/Users/lingzili/Documents/RNA_Pulldown/data/kallisto/SM5169",
  "C:/Users/lingzili/Documents/RNA_Pulldown/data/kallisto/SM5170"
)

# Add sample names to file paths
names(paths) <- sample_id

# Load experimental design
s2c <- read.csv("C:/Users/lingzili/Documents/RNA_Pulldown/data/kallisto/sample_info.csv")

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
  scale_color_manual(values = c("darkred", "#E7B800", "#FC4E07", "#4E84C4", "green"))

pca_p2

ggsave(here::here("graph/GcgTRAP_PCA_26092019.png"), pca_p2)

# Group density -----------------------------------------------------------
# Plot group density
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

ggsave(here::here("graph/GcgTRAP_Density_24092019.png"), pgd)

# Wald test ---------------------------------------------------------------
# wt generates the beta statistic, which approximates to the fold change in expression between the 2 condition tested
so <- sleuth_fit(so, ~condition, "full")

so <- sleuth_fit(so, formula = ~1, fit_name = "reduced")

so <- sleuth_lrt(so, "reduced", "full")

so <- sleuth_wt(so, "conditionPulldown")

wt_results <- sleuth_results(so, "conditionPulldown", "wt", show_all = TRUE)

# Create a table of significantly differential genes
table(wt_results[, "qval"] < 0.05)

wt_significant <- dplyr::filter(wt_results, qval <= 0.05)

wt_significant

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

# MA plot
ma <- plot_ma(so, "conditionPulldown", test_type = "wt", which_model = "full", sig_level = 0.05, point_alpha = 0.2, sig_color = "red", highlight = NULL, highlight_color = "green") +
  standard_theme +
  labs(x = "mean(log(counts+0.5))", y = "beta value (Pulldown vs. Input)") +
  theme(legend.position = "top")

ma

ggsave(here::here("graph/MA_plpt_26092019.png"), ma)

# Plot sample heatmap
png("graph/heatmap.png")
plot_sample_heatmap(so,
  use_filtered = TRUE, color_high = "white",
  color_low = "dodgerblue", x_axis_angle = 50,
  annotation_cols = setdiff(colnames(so$sample_to_covariates), "sample"),
  cluster_bool = TRUE
)
dev.off()

# Boxplots
## Transcript: Mafb
png("graph/Mafb.png", width = 800, height = 400, units = "px")

plot_bootstrap(so, "ENSMUST00000099126.4", units = "tpm", color_by = "condition") +
  standard_theme +
  labs(title = "Mafb", x = NULL, y = "TPM")

dev.off()

## Transcript: Ins1
png("graph/Ins1.png", width = 800, height = 400, units = "px")

plot_bootstrap(so, "ENSMUST00000039652.5", units = "tpm", color_by = "condition") +
  standard_theme +
  labs(title = "Ins1", x = NULL, y = "TPM")

dev.off()

# Transcript: Prss2
png("graph/Prss2.png", width = 800, height = 400, units = "px")

plot_bootstrap(so, "ENSMUST00000070380.4", units = "tpm", color_by = "condition") +
  standard_theme +
  labs(title = "Prss2", x = NULL, y = "TPM")

dev.off()

# Transcript: Try10
png("graph/Try10.png", width = 800, height = 400, units = "px")

plot_bootstrap(so, "ENSMUST00000072103.6", units = "tpm", color_by = "condition") +
  standard_theme +
  labs(title = "Try10", x = NULL, y = "TPM")

dev.off()

# Sleuth Shiny app
sleuth_live(so)

# Heatmap with top 20 genes from TRAP input
# Load kallisto table downloaded from Shiny app
kallisto_table <- read.csv("~/RNA_Pulldown/data/kallisto_table_5samples.csv", row.names = 1)
View(kallisto_table)

# Extract gene names
gene_info <- wt_results[, c("ext_gene", "target_id")]

# Merge
merge_table <- merge(kallisto_table, gene_info, by = "target_id", all.x = TRUE)
View(merge_table)

# Rank the sample
merge_table$sample <- factor(merge_table$sample, levels = c("TRAP Input", "TRAP Pulldown", "GcgcreTRAP 1 Input", "GcgcreTRAP 1 Pulldown", "GcgcreTRAP 2 Input"))

# Prepare data for pheatmap
tpm_heatmap <- merge_table %>%
  dplyr::select(c(sample, tpm, ext_gene, target_id)) %>%
  spread(sample, tpm)

# Order the data by TRAP input
tpm_heatmap <- tpm_heatmap[order(-tpm_heatmap$`TRAP Input`), ]

# Remove rows with duplicated gene names
tpm_heatmap <- tpm_heatmap %>% distinct(ext_gene, .keep_all = TRUE)

library(pheatmap)

# Acinar markers: Try5, Prss2, Pnlip, Spink1
# Alpha cell markers: Gc, Plce1, Mafb
# Beta cell markers: Ins1
tpm_marker <- tpm_heatmap %>% filter(ext_gene %in% c("Try5", "Prss2", "Pnlip", "Amy2a4", "Gc", "Plce1", "Mafb", "Ins1"))

rownames(tpm_marker) <- tpm_marker$ext_gene

pheatmap(tpm_marker[, 3:7], scale = "row", fontsize = 14, cluster_cols = FALSE)

# Subset top 20
tpm_heatmap_top20 <- tpm_heatmap[1:20, ]

rownames(tpm_heatmap_top20) <- tpm_heatmap_top20$ext_gene
