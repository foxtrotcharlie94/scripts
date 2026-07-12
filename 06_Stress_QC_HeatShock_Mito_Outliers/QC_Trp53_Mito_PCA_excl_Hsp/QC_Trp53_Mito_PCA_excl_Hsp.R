library(ggplot2)
library(dplyr)
library(patchwork)
library(ggrepel)

# ── 1. Load original counts ───────────────────────────────────────────────────
counts_hsc <- read.table("LT-HSCs_gene_nUMI_per_cell_for_Foivos.txt",
                         header = TRUE, row.names = 1, sep = "\t")

cat("Dimensions:", dim(counts_hsc), "\n")

# ── 2. Basic metadata ─────────────────────────────────────────────────────────
sample_ids  <- sub("_.*", "", colnames(counts_hsc))
total_umis  <- colSums(counts_hsc)

# ── 3. Identify gene sets ─────────────────────────────────────────────────────
# Trp53
trp53_gene <- grep("^Trp53$", rownames(counts_hsc), value = TRUE, ignore.case = TRUE)
cat("Trp53 gene found:", trp53_gene, "\n")

# Heat shock genes
hsp_genes <- grep("^Hsp", rownames(counts_hsc), value = TRUE, ignore.case = TRUE)
cat("Heat shock genes found (", length(hsp_genes), "):\n")
print(hsp_genes)

# Mitochondrial genes
mt_genes <- grep("^mt-", rownames(counts_hsc), value = TRUE, ignore.case = TRUE)
cat("Mitochondrial genes found (", length(mt_genes), "):\n")
print(mt_genes)

# ── 4. Compute per-cell ratios ────────────────────────────────────────────────
trp53_ratio <- as.numeric(counts_hsc[trp53_gene, ]) / total_umis
hsp_ratio   <- colSums(counts_hsc[hsp_genes, ])     / total_umis
mt_ratio    <- colSums(counts_hsc[mt_genes, ])      / total_umis

df <- data.frame(
  cell        = colnames(counts_hsc),
  sample      = sample_ids,
  exposure    = sub("[0-9]+$", "", sample_ids),
  total_umis  = total_umis,
  trp53_ratio = trp53_ratio,
  hsp_ratio   = hsp_ratio,
  mt_ratio    = mt_ratio
)

# ── 5. Plot helper: violin + boxplot per sample ───────────────────────────────
plot_ratio <- function(df, y_var, y_label, title) {
  ggplot(df, aes_string(x = "sample", y = y_var, fill = "sample")) +
    geom_violin(scale = "width", alpha = 0.7, trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9, fill = "white") +
    scale_fill_brewer(palette = "Set2") +
    scale_y_log10() +
    labs(title = title, x = NULL, y = y_label) +
    theme_classic() +
    theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 9),
          legend.position = "none",
          plot.title    = element_text(face = "bold", size = 11))
}

# ── 6. Trp53 plot ─────────────────────────────────────────────────────────────
p_trp53 <- plot_ratio(df, "trp53_ratio",
                      y_label = "Trp53 UMIs / total UMIs",
                      title   = "Trp53 expression per cell per sample - LT-HSC")

ggsave("hsc_qc_trp53_per_sample.pdf", plot = p_trp53, width = 10, height = 6)

# ── 7. Heat shock plot ────────────────────────────────────────────────────────
p_hsp <- plot_ratio(df, "hsp_ratio",
                    y_label = "Heat shock UMIs / total UMIs",
                    title   = paste0("Heat shock gene expression per cell per sample - LT-HSC\n(",
                                     length(hsp_genes), " genes: ",
                                     paste(hsp_genes[1:min(5, length(hsp_genes))], collapse = ", "),
                                     ifelse(length(hsp_genes) > 5, "...", ""), ")"))

ggsave("hsc_qc_hsp_per_sample.pdf", plot = p_hsp, width = 10, height = 6)

# ── 8. Mitochondrial plot ─────────────────────────────────────────────────────
p_mt <- plot_ratio(df, "mt_ratio",
                   y_label = "Mitochondrial UMIs / total UMIs",
                   title   = paste0("Mitochondrial gene expression per cell per sample - LT-HSC\n(",
                                    length(mt_genes), " genes)"))

ggsave("hsc_qc_mt_per_sample.pdf", plot = p_mt, width = 10, height = 6)

# ── 9. Combined 3-panel ───────────────────────────────────────────────────────
p_combined_qc <- p_trp53 / p_hsp / p_mt
ggsave("hsc_qc_combined.pdf", plot = p_combined_qc, width = 10, height = 16)

# ── 10. Pseudobulk PCA ────────────────────────────────────────────────────────
cat("Computing pseudobulk aggregation...\n")

pseudobulk <- t(sapply(unique(sample_ids), function(s) {
  cells <- which(sample_ids == s)
  rowSums(counts_hsc[, cells, drop = FALSE])
}))

# Normalize to CPM and log transform
pseudobulk_cpm <- pseudobulk / rowSums(pseudobulk) * 1e6
pseudobulk_log <- log1p(pseudobulk_cpm)

# PCA metadata
pca_meta <- data.frame(
  sample   = rownames(pseudobulk_log),
  exposure = sub("[0-9]+$", "", rownames(pseudobulk_log))
)

# ── PCA 1: All genes ──────────────────────────────────────────────────────────
pca_all     <- prcomp(pseudobulk_log, scale. = TRUE, center = TRUE)
pca_all_df  <- as.data.frame(pca_all$x[, 1:2])
pca_all_df$sample   <- rownames(pca_all_df)
pca_all_df$exposure <- sub("[0-9]+$", "", pca_all_df$sample)
var_all <- round(summary(pca_all)$importance[2, 1:2] * 100, 1)

p_pca_all <- ggplot(pca_all_df, aes(x = PC1, y = PC2, color = exposure, label = sample)) +
  geom_point(size = 5) +
  geom_text_repel(size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = c("LB" = "#2196F3", "HB" = "#F44336")) +
  labs(title = "Pseudobulk PCA - All genes - LT-HSC",
       x = paste0("PC1 (", var_all[1], "%)"),
       y = paste0("PC2 (", var_all[2], "%)"),
       color = "Exposure") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"))

# ── PCA 2: Excluding HSP + MT genes ──────────────────────────────────────────
exclude_genes        <- c(hsp_genes, mt_genes)
genes_keep           <- setdiff(colnames(pseudobulk_log), exclude_genes)
pseudobulk_log_filt  <- pseudobulk_log[, genes_keep]

pca_filt    <- prcomp(pseudobulk_log_filt, scale. = TRUE, center = TRUE)
pca_filt_df <- as.data.frame(pca_filt$x[, 1:2])
pca_filt_df$sample   <- rownames(pca_filt_df)
pca_filt_df$exposure <- sub("[0-9]+$", "", pca_filt_df$sample)
var_filt <- round(summary(pca_filt)$importance[2, 1:2] * 100, 1)

p_pca_filt <- ggplot(pca_filt_df, aes(x = PC1, y = PC2, color = exposure, label = sample)) +
  geom_point(size = 5) +
  geom_text_repel(size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = c("LB" = "#2196F3", "HB" = "#F44336")) +
  labs(title = paste0("Pseudobulk PCA - Excluding HSP + MT genes - LT-HSC\n(",
                      length(exclude_genes), " genes removed)"),
       x = paste0("PC1 (", var_filt[1], "%)"),
       y = paste0("PC2 (", var_filt[2], "%)"),
       color = "Exposure") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"))

# ── 11. Save PCA plots ────────────────────────────────────────────────────────
ggsave("hsc_pseudobulk_pca_all_genes.pdf",  plot = p_pca_all,           width = 7,  height = 6)
ggsave("hsc_pseudobulk_pca_no_hsp_mt.pdf",  plot = p_pca_filt,          width = 7,  height = 6)
ggsave("hsc_pseudobulk_pca_combined.pdf",   plot = p_pca_all | p_pca_filt, width = 14, height = 6)

# ── 12. Save ratio table ──────────────────────────────────────────────────────
write.csv(df, "hsc_qc_ratios_per_cell.csv", row.names = FALSE)

cat("All QC and PCA plots saved!\n")

