library(ggplot2)
library(dplyr)
library(patchwork)
getwd()
setwd("C:/Users/fc809/Downloads")
# ── 1. Load original counts ───────────────────────────────────────────────────
counts_hsc  <- read.table("LT-HSCs_gene_nUMI_per_cell_for_Foivos.txt",
                          header = TRUE, row.names = 1, sep = "\t")
counts_mono <- read.table("classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt",
                          header = TRUE, row.names = 1, sep = "\t")

# ── 2. Function to compute Lenti+ frequency per sample ───────────────────────
lenti_freq <- function(counts, cell_type_label) {
  total_umis  <- colSums(counts)
  lenti_umis  <- as.numeric(counts["Lenti", ])
  lenti_ratio <- lenti_umis / total_umis
  
  lenti_status <- ifelse(lenti_ratio >= 4e-5, "Lenti+",
                         ifelse(lenti_ratio <= 9e-7, "Lenti-", "Unclassified"))
  
  sample_ids <- sub("_.*", "", colnames(counts))
  exposure   <- sub("[0-9]+$", "", sample_ids)
  
  df <- data.frame(sample      = sample_ids,
                   exposure    = exposure,
                   lenti_status = lenti_status)
  
  # Frequency per sample
  freq <- df %>%
    group_by(sample, exposure) %>%
    summarise(
      n_total    = n(),
      n_lenti_pos = sum(lenti_status == "Lenti+"),
      n_lenti_neg = sum(lenti_status == "Lenti-"),
      n_unclass   = sum(lenti_status == "Unclassified"),
      pct_lenti_pos = n_lenti_pos / n_total * 100,
      pct_lenti_neg = n_lenti_neg / n_total * 100,
      pct_unclass   = n_unclass   / n_total * 100,
      .groups = "drop"
    ) %>%
    mutate(cell_type = cell_type_label,
           sample = factor(sample, levels = sort(unique(sample))))
  
  return(freq)
}

freq_hsc  <- lenti_freq(counts_hsc,  "LT-HSC")
freq_mono <- lenti_freq(counts_mono, "Classical Monocytes")
freq_all  <- bind_rows(freq_hsc, freq_mono)

# Print summary
cat("LT-HSC Lenti+ frequency per sample:\n")
print(freq_hsc[, c("sample", "exposure", "n_total", "n_lenti_pos", "pct_lenti_pos")])
cat("\nClassical Monocyte Lenti+ frequency per sample:\n")
print(freq_mono[, c("sample", "exposure", "n_total", "n_lenti_pos", "pct_lenti_pos")])

# ── 3. Plot: % Lenti+ per sample ──────────────────────────────────────────────
plot_lenti_freq <- function(freq_df, title) {
  ggplot(freq_df, aes(x = sample, y = pct_lenti_pos, fill = exposure)) +
    geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.3) +
    geom_text(aes(label = paste0(round(pct_lenti_pos, 1), "%")),
              vjust = -0.4, size = 3, fontface = "bold") +
    scale_fill_manual(values = c("HB" = "#F44336", "LB" = "#2196F3")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = title,
         x = NULL, y = "% Lenti+ cells",
         fill = "Exposure",
         caption = paste0("Lenti+ threshold: UMI ratio ≥ 4×10⁻⁵\n",
                          "n shown above bars = % of total cells per sample")) +
    theme_classic() +
    theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 10),
          plot.title    = element_text(face = "bold", size = 12),
          plot.caption  = element_text(size = 8, color = "grey50"),
          legend.position = "top")
}

p_hsc  <- plot_lenti_freq(freq_hsc,  "Lenti+ frequency per sample - LT-HSC")
p_mono <- plot_lenti_freq(freq_mono, "Lenti+ frequency per sample - Classical Monocytes")

# ── 4. Stacked bar showing all 3 categories ───────────────────────────────────
plot_lenti_stacked <- function(freq_df, title) {
  freq_long <- freq_df %>%
    select(sample, exposure, pct_lenti_pos, pct_lenti_neg, pct_unclass) %>%
    tidyr::pivot_longer(cols = starts_with("pct_"),
                        names_to = "status", values_to = "pct") %>%
    mutate(status = recode(status,
                           "pct_lenti_pos" = "Lenti+",
                           "pct_lenti_neg" = "Lenti-",
                           "pct_unclass"   = "Unclassified"),
           status = factor(status, levels = c("Lenti+", "Unclassified", "Lenti-")))
  
  ggplot(freq_long, aes(x = sample, y = pct, fill = status)) +
    geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = c("Lenti+"       = "#2CA02C",
                                 "Lenti-"       = "#9467BD",
                                 "Unclassified" = "#CCCCCC")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    facet_grid(~ exposure, scales = "free_x", space = "free_x") +
    labs(title = title,
         x = NULL, y = "% of cells", fill = "Lenti status",
         caption = "Lenti+: ratio ≥ 4×10⁻⁵  |  Lenti-: ratio ≤ 9×10⁻⁷  |  Grey: unclassified") +
    theme_classic() +
    theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 10),
          plot.title     = element_text(face = "bold", size = 12),
          plot.caption   = element_text(size = 8, color = "grey50"),
          strip.text     = element_text(face = "bold"),
          legend.position = "top")
}

p_hsc_stack  <- plot_lenti_stacked(freq_hsc,  "Lenti classification per sample - LT-HSC")
p_mono_stack <- plot_lenti_stacked(freq_mono, "Lenti classification per sample - Classical Monocytes")

# ── 5. Combined plot: both cell types side by side ────────────────────────────
p_combined_bar   <- p_hsc   / p_mono
p_combined_stack <- p_hsc_stack / p_mono_stack

# ── 6. Save ───────────────────────────────────────────────────────────────────
ggsave("lenti_freq_hsc_barplot.pdf",         plot = p_hsc,            width = 10, height = 6)
ggsave("lenti_freq_mono_barplot.pdf",        plot = p_mono,           width = 10, height = 6)
ggsave("lenti_freq_hsc_stacked.pdf",         plot = p_hsc_stack,      width = 10, height = 6)
ggsave("lenti_freq_mono_stacked.pdf",        plot = p_mono_stack,     width = 10, height = 6)
ggsave("lenti_freq_combined_barplot.pdf",    plot = p_combined_bar,   width = 10, height = 12)
ggsave("lenti_freq_combined_stacked.pdf",    plot = p_combined_stack, width = 10, height = 12)

# ── 7. Save frequency table ───────────────────────────────────────────────────
write.csv(freq_all, "lenti_frequency_per_sample.csv", row.names = FALSE)

cat("All Lenti frequency plots saved!\n")

library(tidyr)

# ── Ratio: Lenti+% Mono / Lenti+% HSC per sample ─────────────────────────────
ratio_df <- freq_hsc %>%
  select(sample, exposure, pct_lenti_pos) %>%
  rename(pct_hsc = pct_lenti_pos) %>%
  inner_join(
    freq_mono %>% select(sample, pct_lenti_pos) %>% rename(pct_mono = pct_lenti_pos),
    by = "sample"
  ) %>%
  mutate(ratio = pct_mono / pct_hsc)

cat("Lenti+% Mono / Lenti+% HSC per sample:\n")
print(ratio_df)

p_ratio <- ggplot(ratio_df, aes(x = sample, y = ratio, fill = exposure)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = round(ratio, 2)),
            vjust = -0.4, size = 3.5, fontface = "bold") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  scale_fill_manual(values = c("HB" = "#F44336", "LB" = "#2196F3")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Ratio of Lenti+ frequency: Classical Monocytes / LT-HSC",
       subtitle = "Dashed line = ratio of 1 (equal Lenti+ frequency in both cell types)",
       x = NULL, y = "Lenti+% Mono / Lenti+% HSC",
       fill = "Exposure",
       caption = "Values >1 = higher Lenti+ in monocytes than HSCs | Values <1 = higher Lenti+ in HSCs") +
  theme_classic() +
  theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 10),
        plot.title     = element_text(face = "bold", size = 12),
        plot.subtitle  = element_text(size = 9, color = "grey40"),
        plot.caption   = element_text(size = 8, color = "grey50"),
        legend.position = "top")

ggsave("lenti_freq_ratio_mono_over_hsc.pdf", plot = p_ratio, width = 10, height = 6)

# Save updated table with ratio
write.csv(ratio_df, "lenti_freq_ratio_per_sample.csv", row.names = FALSE)

cat("Ratio plot saved!\n")

library(Seurat)

# ── Load Seurat objects ───────────────────────────────────────────────────────
seu_hsc  <- readRDS("LT_HSC_seurat_final.rds")
seu_mono <- readRDS("classical_mono_seurat_final.rds")

# ── Cluster labels ────────────────────────────────────────────────────────────
hsc_labels <- c(
  "0" = "0: Quiescent LT-HSC",
  "1" = "1: Metabolically active HSC",
  "2" = "2: Lymphoid/MPP-biased",
  "3" = "3: IFN-activated HSC",
  "4" = "4: Cycling HSC (S-phase)",
  "5" = "5: Mk-biased HSC"
)

mono_labels <- c(
  "0" = "0: Mature classical mono",
  "1" = "1: Ly6C-high mono",
  "2" = "2: IFN-stimulated mono",
  "3" = "3: Proliferating mono",
  "4" = "4: Immature/emergency mono",
  "5" = "5: MHC-II high mono",
  "6" = "6: Anti-inflammatory mono",
  "7" = "7: Transitional mono"
)

# ── Helper: Lenti+ frequency per cluster per sample ───────────────────────────
lenti_cluster_freq <- function(seu, cluster_labels, cell_type_label) {
  meta <- seu@meta.data %>%
    filter(lenti_status %in% c("Lenti+", "Lenti-")) %>%
    mutate(cluster_label = cluster_labels[as.character(seurat_clusters)])
  
  freq <- meta %>%
    group_by(sample, exposure, cluster_label) %>%
    summarise(
      n_total     = n(),
      n_lenti_pos = sum(lenti_status == "Lenti+"),
      pct_lenti_pos = n_lenti_pos / n_total * 100,
      .groups = "drop"
    ) %>%
    mutate(cell_type = cell_type_label)
  
  return(freq)
}

freq_cluster_hsc  <- lenti_cluster_freq(seu_hsc,  hsc_labels,  "LT-HSC")
freq_cluster_mono <- lenti_cluster_freq(seu_mono, mono_labels, "Classical Monocytes")

# ── Helper: dotplot with mean ± SEM per cluster ───────────────────────────────
plot_cluster_lenti <- function(freq_df, cell_type_label) {
  
  # Mean ± SEM per cluster per exposure
  summary_df <- freq_df %>%
    group_by(exposure, cluster_label) %>%
    summarise(
      mean_pct = mean(pct_lenti_pos),
      sem_pct  = sd(pct_lenti_pos) / sqrt(n()),
      .groups  = "drop"
    )
  
  ggplot() +
    # Individual sample dots
    geom_point(data = freq_df,
               aes(x = cluster_label, y = pct_lenti_pos, color = exposure),
               position = position_dodge(width = 0.6),
               size = 2.5, alpha = 0.7) +
    # Mean ± SEM
    geom_point(data = summary_df,
               aes(x = cluster_label, y = mean_pct, color = exposure),
               position = position_dodge(width = 0.6),
               size = 5, shape = 18) +
    geom_errorbar(data = summary_df,
                  aes(x = cluster_label,
                      ymin = mean_pct - sem_pct,
                      ymax = mean_pct + sem_pct,
                      color = exposure),
                  position = position_dodge(width = 0.6),
                  width = 0.25, linewidth = 0.8) +
    scale_color_manual(values = c("HB" = "#F44336", "LB" = "#2196F3")) +
    labs(title = paste0("Lenti+ frequency per cluster - ", cell_type_label),
         subtitle = "Diamond = mean  |  Lines = ±SEM  |  Dots = individual samples",
         x = NULL, y = "% Lenti+ cells", color = "Exposure") +
    theme_classic() +
    theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 9),
          plot.title     = element_text(face = "bold", size = 12),
          plot.subtitle  = element_text(size = 9, color = "grey40"),
          legend.position = "top")
}

p_dot_hsc  <- plot_cluster_lenti(freq_cluster_hsc,  "LT-HSC")
p_dot_mono <- plot_cluster_lenti(freq_cluster_mono, "Classical Monocytes")

# ── Combined ──────────────────────────────────────────────────────────────────
p_dot_combined <- p_dot_hsc / p_dot_mono

# ── Save ──────────────────────────────────────────────────────────────────────
ggsave("lenti_cluster_dotplot_hsc.pdf",      plot = p_dot_hsc,      width = 12, height = 6)
ggsave("lenti_cluster_dotplot_mono.pdf",     plot = p_dot_mono,     width = 14, height = 6)
ggsave("lenti_cluster_dotplot_combined.pdf", plot = p_dot_combined, width = 14, height = 12)

write.csv(freq_cluster_hsc,  "lenti_cluster_freq_hsc.csv",  row.names = FALSE)
write.csv(freq_cluster_mono, "lenti_cluster_freq_mono.csv", row.names = FALSE)

cat("Cluster Lenti+ dotplots saved!\n")


##########check##########

library(Seurat)

# ── 1. Load objects and original counts ───────────────────────────────────────
seu_hsc  <- readRDS("LT_HSC_seurat_final.rds")
seu_mono <- readRDS("classical_mono_seurat_final.rds")

counts_hsc  <- read.table("LT-HSCs_gene_nUMI_per_cell_for_Foivos.txt",
                          header = TRUE, row.names = 1, sep = "\t")
counts_mono <- read.table("classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt",
                          header = TRUE, row.names = 1, sep = "\t")

# ── 2. Recompute Lenti classification from scratch ────────────────────────────
classify_lenti <- function(counts) {
  total_umis  <- colSums(counts)
  lenti_umis  <- as.numeric(counts["Lenti", ])
  lenti_ratio <- lenti_umis / total_umis
  status <- ifelse(lenti_ratio >= 4e-5, "Lenti+",
                   ifelse(lenti_ratio <= 9e-7, "Lenti-", "Unclassified"))
  return(data.frame(
    cell         = colnames(counts),
    sample       = sub("_.*", "", colnames(counts)),
    lenti_ratio  = lenti_ratio,
    lenti_status = status
  ))
}

lenti_hsc  <- classify_lenti(counts_hsc)
lenti_mono <- classify_lenti(counts_mono)

# ── 3. Compare with what's stored in Seurat metadata ─────────────────────────
cat("=== LT-HSC ===\n")
cat("Cells in Seurat object:", ncol(seu_hsc), "\n")
cat("Cells in original counts:", ncol(counts_hsc), "\n")

# Check if cell names match
hsc_cells_in_seurat <- colnames(seu_hsc)
hsc_lenti_seurat    <- seu_hsc$lenti_status
hsc_lenti_recomputed <- lenti_hsc$lenti_status[match(hsc_cells_in_seurat, lenti_hsc$cell)]

cat("\nLenti status in Seurat object:\n")
print(table(hsc_lenti_seurat))
cat("\nLenti status recomputed from original counts (for Seurat cells only):\n")
print(table(hsc_lenti_recomputed))
cat("\nAre they identical?\n")
print(all(hsc_lenti_seurat == hsc_lenti_recomputed, na.rm = TRUE))
cat("Number of mismatches:", sum(hsc_lenti_seurat != hsc_lenti_recomputed, na.rm = TRUE), "\n")

cat("\n=== Classical Monocytes ===\n")
cat("Cells in Seurat object:", ncol(seu_mono), "\n")
cat("Cells in original counts:", ncol(counts_mono), "\n")

mono_cells_in_seurat  <- colnames(seu_mono)
mono_lenti_seurat     <- seu_mono$lenti_status
mono_lenti_recomputed <- lenti_mono$lenti_status[match(mono_cells_in_seurat, lenti_mono$cell)]

cat("\nLenti status in Seurat object:\n")
print(table(mono_lenti_seurat))
cat("\nLenti status recomputed from original counts (for Seurat cells only):\n")
print(table(mono_lenti_recomputed))
cat("\nAre they identical?\n")
print(all(mono_lenti_seurat == mono_lenti_recomputed, na.rm = TRUE))
cat("Number of mismatches:", sum(mono_lenti_seurat != mono_lenti_recomputed, na.rm = TRUE), "\n")

# ── 4. Check Lenti+ frequency per sample from Seurat vs from raw counts ───────
cat("\n=== Lenti+ % per sample - HSC ===\n")
hsc_freq_seurat <- table(seu_hsc$sample, seu_hsc$lenti_status)
hsc_freq_pct    <- round(prop.table(hsc_freq_seurat, margin = 1) * 100, 1)
print(hsc_freq_pct)

cat("\n=== Lenti+ % per sample - Mono ===\n")
mono_freq_seurat <- table(seu_mono$sample, seu_mono$lenti_status)
mono_freq_pct    <- round(prop.table(mono_freq_seurat, margin = 1) * 100, 1)
print(mono_freq_pct)

# ── 5. Check Lenti ratio distribution per sample ──────────────────────────────
cat("\n=== Lenti ratio summary per sample - HSC ===\n")
lenti_hsc$sample <- sub("_.*", "", lenti_hsc$cell)
print(tapply(lenti_hsc$lenti_ratio, lenti_hsc$sample, function(x)
  c(min=min(x), median=median(x), max=max(x), pct_above_threshold=mean(x >= 4e-5)*100)))

cat("\n=== Lenti ratio summary per sample - Mono ===\n")
lenti_mono$sample <- sub("_.*", "", lenti_mono$cell)
print(tapply(lenti_mono$lenti_ratio, lenti_mono$sample, function(x)
  c(min=min(x), median=median(x), max=max(x), pct_above_threshold=mean(x >= 4e-5)*100)))

#####lenti_number_v_fract##########

library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)

# ── 1. Load original counts ───────────────────────────────────────────────────
counts_hsc  <- read.table("LT-HSCs_gene_nUMI_per_cell_for_Foivos.txt",
                          header = TRUE, row.names = 1, sep = "\t")
counts_mono <- read.table("classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt",
                          header = TRUE, row.names = 1, sep = "\t")

# ── 2. Helper: compute both classifications ───────────────────────────────────
classify_lenti_both <- function(counts, cell_type_label) {
  total_umis   <- colSums(counts)
  lenti_umis   <- as.numeric(counts["Lenti", ])
  lenti_ratio  <- lenti_umis / total_umis
  sample_ids   <- sub("_.*", "", colnames(counts))
  exposure     <- sub("[0-9]+$", "", sample_ids)
  
  df <- data.frame(
    cell         = colnames(counts),
    sample       = sample_ids,
    exposure     = exposure,
    lenti_umis   = lenti_umis,
    total_umis   = total_umis,
    lenti_ratio  = lenti_ratio,
    # Classification 1: current cutoff (ratio >= 4e-5)
    lenti_ratio_cutoff = ifelse(lenti_ratio >= 4e-5, "Lenti+", "Lenti-"),
    # Classification 2: at least 1 UMI
    lenti_umi_cutoff   = ifelse(lenti_umis >= 1, "Lenti+", "Lenti-"),
    cell_type    = cell_type_label
  )
  return(df)
}

df_hsc  <- classify_lenti_both(counts_hsc,  "LT-HSC")
df_mono <- classify_lenti_both(counts_mono, "Classical Monocytes")
df_all  <- bind_rows(df_hsc, df_mono)

# ── 3. Compute frequency per sample for both cutoffs ──────────────────────────
compute_freq <- function(df, cutoff_col, cutoff_label) {
  df %>%
    group_by(cell_type, sample, exposure) %>%
    summarise(
      pct_lenti_pos = mean(.data[[cutoff_col]] == "Lenti+") * 100,
      .groups = "drop"
    ) %>%
    mutate(cutoff = cutoff_label)
}

freq_ratio_hsc  <- compute_freq(df_hsc,  "lenti_ratio_cutoff", "Ratio ≥ 4×10⁻⁵")
freq_umi_hsc    <- compute_freq(df_hsc,  "lenti_umi_cutoff",   "≥1 UMI")
freq_ratio_mono <- compute_freq(df_mono, "lenti_ratio_cutoff", "Ratio ≥ 4×10⁻⁵")
freq_umi_mono   <- compute_freq(df_mono, "lenti_umi_cutoff",   "≥1 UMI")

freq_hsc_both  <- bind_rows(freq_ratio_hsc,  freq_umi_hsc)
freq_mono_both <- bind_rows(freq_ratio_mono, freq_umi_mono)

# ── 4. Barplot comparison function ────────────────────────────────────────────
plot_cutoff_comparison <- function(freq_df, cell_type_label) {
  freq_df <- freq_df %>%
    mutate(sample = factor(sample, levels = sort(unique(sample))),
           cutoff = factor(cutoff, levels = c("Ratio ≥ 4×10⁻⁵", "≥1 UMI")))
  
  ggplot(freq_df, aes(x = sample, y = pct_lenti_pos, fill = cutoff)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7, color = "white", linewidth = 0.3) +
    geom_text(aes(label = paste0(round(pct_lenti_pos, 1), "%")),
              position = position_dodge(width = 0.8),
              vjust = -0.4, size = 2.8, fontface = "bold") +
    scale_fill_manual(values = c("Ratio ≥ 4×10⁻⁵" = "#2CA02C",
                                 "≥1 UMI"          = "#FF7F0E")) +
    facet_grid(~ exposure, scales = "free_x", space = "free_x") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = paste0("Lenti+ frequency comparison - ", cell_type_label),
         subtitle = "Green = current ratio cutoff  |  Orange = ≥1 UMI cutoff",
         x = NULL, y = "% Lenti+ cells", fill = "Cutoff") +
    theme_classic() +
    theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 9),
          plot.title     = element_text(face = "bold", size = 12),
          plot.subtitle  = element_text(size = 9, color = "grey40"),
          legend.position = "top",
          strip.text     = element_text(face = "bold"))
}

p_bar_hsc  <- plot_cutoff_comparison(freq_hsc_both,  "LT-HSC")
p_bar_mono <- plot_cutoff_comparison(freq_mono_both, "Classical Monocytes")
p_bar_combined <- p_bar_hsc / p_bar_mono

ggsave("lenti_cutoff_comparison_hsc.pdf",      plot = p_bar_hsc,      width = 10, height = 6)
ggsave("lenti_cutoff_comparison_mono.pdf",     plot = p_bar_mono,     width = 10, height = 6)
ggsave("lenti_cutoff_comparison_combined.pdf", plot = p_bar_combined, width = 10, height = 12)

# ── 5. Distribution plots: ratio vs absolute UMI per cell per sample ──────────
plot_distributions <- function(df, cell_type_label, file_prefix) {
  
  samples_order <- sort(unique(df$sample))
  
  # ── Plot A: Lenti UMI ratio per cell (violin) ─────────────────────────────
  # Remove zero-ratio cells for log scale visualization
  df_nonzero <- df %>% filter(lenti_ratio > 0)
  
  p_ratio <- ggplot(df_nonzero, aes(x = sample, y = lenti_ratio, fill = exposure)) +
    geom_violin(scale = "width", alpha = 0.7, trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.9) +
    # Mark the ratio cutoff threshold
    geom_hline(yintercept = 4e-5, linetype = "dashed",
               color = "#2CA02C", linewidth = 0.8) +
    geom_hline(yintercept = 9e-7, linetype = "dashed",
               color = "#9467BD", linewidth = 0.8) +
    scale_y_log10(labels = scales::scientific) +
    scale_fill_manual(values = c("HB" = "#F44336", "LB" = "#2196F3")) +
    scale_x_discrete(limits = samples_order) +
    labs(title = paste0("Lenti UMIs / Total UMIs per cell - ", cell_type_label),
         subtitle = "Green dashed = Lenti+ threshold (4×10⁻⁵) | Purple dashed = Lenti- threshold (9×10⁻⁷)\nOnly cells with >0 Lenti UMIs shown | Y-axis log scale",
         x = NULL, y = "Lenti UMIs / Total UMIs", fill = "Exposure") +
    theme_classic() +
    theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 9),
          plot.title     = element_text(face = "bold", size = 11),
          plot.subtitle  = element_text(size = 8, color = "grey40"),
          legend.position = "top")
  
  # ── Plot B: Absolute Lenti UMIs per cell - Lenti+ cells only ─────────────
  df_lenti_pos <- df %>% filter(lenti_ratio >= 4e-5)  # only Lenti+ cells
  
  p_abs <- ggplot(df_lenti_pos, aes(x = sample, y = lenti_umis, fill = exposure)) +
    geom_violin(scale = "width", alpha = 0.7, trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.9) +
    geom_hline(yintercept = 1, linetype = "dashed",
               color = "#FF7F0E", linewidth = 0.8) +
    scale_y_log10() +
    scale_fill_manual(values = c("HB" = "#F44336", "LB" = "#2196F3")) +
    scale_x_discrete(limits = samples_order) +
    labs(title = paste0("Absolute Lenti UMIs per cell (Lenti+ cells only) - ", cell_type_label),
         subtitle = "Only cells classified as Lenti+ (ratio >= 4e-5⁵) | Orange dashed = 1 UMI | Y-axis log scale",
         x = NULL, y = "Lenti UMI count", fill = "Exposure") +
    theme_classic() +
    theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 9),
          plot.title     = element_text(face = "bold", size = 11),
          plot.subtitle  = element_text(size = 8, color = "grey40"),
          legend.position = "top")
  
  # ── Plot C: % cells with >0 Lenti UMIs per sample ─────────────────────────
  pct_nonzero <- df %>%
    group_by(sample, exposure) %>%
    summarise(pct_any_lenti = mean(lenti_umis > 0) * 100, .groups = "drop") %>%
    mutate(sample = factor(sample, levels = samples_order))
  
  p_pct_nonzero <- ggplot(pct_nonzero, aes(x = sample, y = pct_any_lenti, fill = exposure)) +
    geom_bar(stat = "identity", width = 0.7, color = "white") +
    geom_text(aes(label = paste0(round(pct_any_lenti, 1), "%")),
              vjust = -0.4, size = 3, fontface = "bold") +
    scale_fill_manual(values = c("HB" = "#F44336", "LB" = "#2196F3")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = paste0("% cells with at least 1 Lenti UMI - ", cell_type_label),
         x = NULL, y = "% cells with ≥1 Lenti UMI", fill = "Exposure") +
    theme_classic() +
    theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 9),
          plot.title     = element_text(face = "bold", size = 11),
          legend.position = "top")
  
  # ── Combine and save ───────────────────────────────────────────────────────
  p_dist_combined <- p_ratio / p_abs / p_pct_nonzero
  
  ggsave(paste0(file_prefix, "_lenti_ratio_distribution.pdf"),    plot = p_ratio,        width = 11, height = 6)
  ggsave(paste0(file_prefix, "_lenti_abs_distribution.pdf"),      plot = p_abs,          width = 11, height = 6)
  ggsave(paste0(file_prefix, "_lenti_pct_nonzero.pdf"),           plot = p_pct_nonzero,  width = 11, height = 6)
  ggsave(paste0(file_prefix, "_lenti_distributions_combined.pdf"),plot = p_dist_combined,width = 11, height = 16)
  
  cat("Distribution plots saved for", cell_type_label, "\n")
}

plot_distributions(df_hsc,  "LT-HSC",                "hsc")
plot_distributions(df_mono, "Classical Monocytes",    "mono")

cat("All plots saved!\n")

#####subsampling_vs_Lentipositives########

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ── 2. Downsampling function ──────────────────────────────────────────────────
downsample_cell <- function(cell_counts, target) {
  total <- sum(cell_counts)
  if (total <= target) return(cell_counts)
  sampled <- sample(rep(seq_along(cell_counts), cell_counts), target)
  tabulate(sampled, nbins = length(cell_counts))
}

# ── 3. Compute Lenti+ % at different UMI depths ───────────────────────────────
lenti_at_depth <- function(counts, target_umis, cell_type_label, seed = 42) {
  set.seed(seed)
  sample_ids <- sub("_.*", "", colnames(counts))
  
  # Remove cells with fewer UMIs than target
  total_umis <- colSums(counts)
  keep_cells <- total_umis >= target_umis
  counts_sub <- counts[, keep_cells]
  sids_sub   <- sample_ids[keep_cells]
  
  # Downsample
  counts_ds <- apply(counts_sub, 2, downsample_cell, target = target_umis)
  rownames(counts_ds) <- rownames(counts_sub)
  
  # Classify Lenti+
  lenti_umis  <- as.numeric(counts_ds["Lenti", ])
  total_ds    <- colSums(counts_ds)
  lenti_ratio <- lenti_umis / total_ds
  lenti_pos   <- lenti_ratio >= 4e-5
  
  # Frequency per sample
  data.frame(sample       = sids_sub,
             exposure     = sub("[0-9]+$", "", sids_sub),
             lenti_pos    = lenti_pos,
             target_depth = target_umis,
             cell_type    = cell_type_label) %>%
    group_by(sample, exposure, target_depth, cell_type) %>%
    summarise(n_cells     = n(),
              n_lenti_pos = sum(lenti_pos),
              pct_lenti   = mean(lenti_pos) * 100,
              .groups = "drop")
}

# ── 4. Run across a range of UMI depths ──────────────────────────────────────
# Use "original" as no-subsampling baseline
original_lenti <- function(counts, cell_type_label) {
  sample_ids  <- sub("_.*", "", colnames(counts))
  total_umis  <- colSums(counts)
  lenti_umis  <- as.numeric(counts["Lenti", ])
  lenti_ratio <- lenti_umis / total_umis
  lenti_pos   <- lenti_ratio >= 4e-5
  
  data.frame(sample    = sample_ids,
             exposure  = sub("[0-9]+$", "", sample_ids),
             lenti_pos = lenti_pos,
             cell_type = cell_type_label) %>%
    group_by(sample, exposure, cell_type) %>%
    summarise(n_cells     = n(),
              n_lenti_pos = sum(lenti_pos),
              pct_lenti   = mean(lenti_pos) * 100,
              .groups = "drop") %>%
    mutate(target_depth = Inf)  # marks as original
}

# Define depths to test
depths <- c(1000, 2000, 3000, 4000, 5000, 7500, 10000, 15000)

cat("Running subsampling across depths for LT-HSC...\n")
hsc_results <- bind_rows(
  original_lenti(counts_hsc, "LT-HSC"),
  bind_rows(lapply(depths, function(d) {
    cat("  Depth:", d, "\n")
    lenti_at_depth(counts_hsc, d, "LT-HSC")
  }))
)

cat("Running subsampling across depths for Classical Monocytes...\n")
mono_results <- bind_rows(
  original_lenti(counts_mono, "Classical Monocytes"),
  bind_rows(lapply(depths, function(d) {
    cat("  Depth:", d, "\n")
    lenti_at_depth(counts_mono, d, "Classical Monocytes")
  }))
)

all_results <- bind_rows(hsc_results, mono_results)

# Convert depth to factor for plotting (with "Original" label)
all_results <- all_results %>%
  mutate(depth_label = ifelse(is.infinite(target_depth), "Original",
                              as.character(target_depth)),
         depth_label = factor(depth_label,
                              levels = c(as.character(depths), "Original")))
# ── Updated subsampling with minimum cell filter ──────────────────────────────
min_cells <- 50  # minimum cells required to include a sample at a given depth

lenti_at_depth_safe <- function(counts, target_umis, cell_type_label, seed = 42) {
  set.seed(seed)
  sample_ids <- sub("_.*", "", colnames(counts))
  total_umis <- colSums(counts)
  
  # Check per sample how many cells survive
  keep_cells <- total_umis >= target_umis
  n_per_sample <- table(sample_ids[keep_cells])
  
  # Only keep samples with enough cells
  valid_samples <- names(n_per_sample[n_per_sample >= min_cells])
  
  if (length(valid_samples) == 0) {
    cat("  Depth", target_umis, "- no samples with >=", min_cells, "cells, skipping\n")
    return(NULL)
  }
  
  cat("  Depth:", target_umis, "- valid samples:", paste(valid_samples, collapse=", "), "\n")
  
  keep_cells_valid <- keep_cells & (sample_ids %in% valid_samples)
  counts_sub <- counts[, keep_cells_valid]
  sids_sub   <- sample_ids[keep_cells_valid]
  
  # Downsample
  counts_ds <- apply(counts_sub, 2, downsample_cell, target = target_umis)
  rownames(counts_ds) <- rownames(counts_sub)
  
  lenti_umis  <- as.numeric(counts_ds["Lenti", ])
  total_ds    <- colSums(counts_ds)
  lenti_ratio <- lenti_umis / total_ds
  lenti_pos   <- lenti_ratio >= 4e-5
  
  data.frame(sample       = sids_sub,
             exposure     = sub("[0-9]+$", "", sids_sub),
             lenti_pos    = lenti_pos,
             target_depth = target_umis,
             cell_type    = cell_type_label) %>%
    group_by(sample, exposure, target_depth, cell_type) %>%
    summarise(n_cells     = n(),
              n_lenti_pos = sum(lenti_pos),
              pct_lenti   = mean(lenti_pos) * 100,
              .groups = "drop")
}

depths <- c(1000, 2000, 3000, 4000, 5000, 7500, 10000, 15000)

cat("LT-HSC:\n")
hsc_results_safe <- bind_rows(
  original_lenti(counts_hsc, "LT-HSC"),
  bind_rows(Filter(Negate(is.null),
                   lapply(depths, function(d) lenti_at_depth_safe(counts_hsc, d, "LT-HSC"))))
)

cat("\nClassical Monocytes:\n")
mono_results_safe <- bind_rows(
  original_lenti(counts_mono, "Classical Monocytes"),
  bind_rows(Filter(Negate(is.null),
                   lapply(depths, function(d) lenti_at_depth_safe(counts_mono, d, "Classical Monocytes"))))
)

all_results_safe <- bind_rows(hsc_results_safe, mono_results_safe) %>%
  mutate(depth_label = ifelse(is.infinite(target_depth), "Original",
                              as.character(target_depth)),
         depth_label = factor(depth_label,
                              levels = c(as.character(depths), "Original")))

# ── Heatmap with n_cells annotated ───────────────────────────────────────────
plot_heatmap_safe <- function(df, cell_type_label) {
  df_ct <- df %>%
    filter(cell_type == cell_type_label) %>%
    mutate(label = paste0(round(pct_lenti, 1), "%\n(n=", n_cells, ")"))
  
  ggplot(df_ct, aes(x = depth_label, y = sample, fill = pct_lenti)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = label), size = 2.5) +
    scale_fill_gradient2(low = "white", mid = "#AED6F1", high = "#1A5276",
                         midpoint = median(df_ct$pct_lenti, na.rm = TRUE),
                         name = "% Lenti+") +
    facet_grid(exposure ~ ., scales = "free_y", space = "free_y") +
    labs(title = paste0("Lenti+ % heatmap across UMI depths - ", cell_type_label),
         subtitle = paste0("Samples excluded at depths where <", min_cells,
                           " cells survive | n = cells used at each depth"),
         x = "UMI subsampling depth", y = NULL) +
    theme_classic() +
    theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 9),
          plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          strip.text    = element_text(face = "bold"))
}

# ── Line plot with safe results ───────────────────────────────────────────────
plot_line <- function(df, cell_type_label) {
  df_ct <- df %>%
    filter(cell_type == cell_type_label,
           !is.infinite(target_depth))  # exclude Original from line plot
  
  summary_df <- df_ct %>%
    group_by(exposure, depth_label) %>%
    summarise(mean_pct = mean(pct_lenti),
              sem_pct  = sd(pct_lenti) / sqrt(n()),
              .groups  = "drop")
  
  ggplot(summary_df, aes(x = depth_label, y = mean_pct,
                         color = exposure, group = exposure)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_pct - sem_pct,
                      ymax = mean_pct + sem_pct),
                  width = 0.2, linewidth = 0.7) +
    scale_color_manual(values = c("HB" = "#F44336", "LB" = "#2196F3")) +
    labs(title = paste0("Impact of UMI subsampling on Lenti+ detection - ", cell_type_label),
         subtitle = "Mean +/- SEM across samples | Ratio cutoff >= 4e-5 | Samples excluded where <50 cells survive",
         x = "UMI subsampling depth", y = "% Lenti+ cells",
         color = "Exposure") +
    theme_classic() +
    theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 9),
          plot.title     = element_text(face = "bold", size = 11),
          plot.subtitle  = element_text(size = 9, color = "grey40"),
          legend.position = "top")
}

# ── Per sample plot with disconnected Original point ─────────────────────────
plot_per_sample <- function(df, cell_type_label) {
  df_ct <- df %>%
    filter(cell_type == cell_type_label) %>%
    mutate(is_original = is.infinite(target_depth))
  
  df_lines <- df_ct %>% filter(!is_original)
  df_orig  <- df_ct %>% filter(is_original)
  
  ggplot() +
    # Lines + points for subsampled depths only
    geom_line(data = df_lines,
              aes(x = depth_label, y = pct_lenti,
                  color = exposure, group = sample),
              linewidth = 0.8, alpha = 0.8) +
    geom_point(data = df_lines,
               aes(x = depth_label, y = pct_lenti,
                   color = exposure, group = sample),
               size = 2, alpha = 0.9) +
    # Original as isolated diamond, not connected
    geom_point(data = df_orig,
               aes(x = depth_label, y = pct_lenti,
                   color = exposure, group = sample),
               size = 4, shape = 18, alpha = 0.9) +
    scale_color_manual(values = c("HB" = "#F44336", "LB" = "#2196F3")) +
    facet_wrap(~ exposure, nrow = 1) +
    labs(title = paste0("Lenti+ % per sample across UMI depths - ", cell_type_label),
         subtitle = "Lines connect subsampled depths (>=50 cells) | Diamond = original unsubsampled | Lines stop at last valid depth per sample",
         x = "UMI subsampling depth", y = "% Lenti+ cells",
         color = "Exposure") +
    theme_classic() +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 9),
          plot.title      = element_text(face = "bold", size = 11),
          plot.subtitle   = element_text(size = 9, color = "grey40"),
          legend.position = "none",
          strip.text      = element_text(face = "bold"))
}

# ── Generate all plots ────────────────────────────────────────────────────────
p_heat_hsc_safe  <- plot_heatmap_safe(all_results_safe, "LT-HSC")
p_heat_mono_safe <- plot_heatmap_safe(all_results_safe, "Classical Monocytes")

p_line_hsc_safe  <- plot_line(all_results_safe, "LT-HSC")
p_line_mono_safe <- plot_line(all_results_safe, "Classical Monocytes")

p_sample_hsc_safe  <- plot_per_sample(all_results_safe, "LT-HSC")
p_sample_mono_safe <- plot_per_sample(all_results_safe, "Classical Monocytes")

# ── Save ──────────────────────────────────────────────────────────────────────
ggsave("lenti_subsampling_heatmap_hsc_safe.pdf",
       plot = p_heat_hsc_safe,  width = 13, height = 7)
ggsave("lenti_subsampling_heatmap_mono_safe.pdf",
       plot = p_heat_mono_safe, width = 13, height = 7)
ggsave("lenti_subsampling_line_hsc_safe.pdf",
       plot = p_line_hsc_safe,  width = 10, height = 6)
ggsave("lenti_subsampling_line_mono_safe.pdf",
       plot = p_line_mono_safe, width = 10, height = 6)
ggsave("lenti_subsampling_line_combined_safe.pdf",
       plot = p_line_hsc_safe / p_line_mono_safe, width = 10, height = 12)
ggsave("lenti_subsampling_persample_hsc_safe.pdf",
       plot = p_sample_hsc_safe,  width = 12, height = 6)
ggsave("lenti_subsampling_persample_mono_safe.pdf",
       plot = p_sample_mono_safe, width = 12, height = 6)

write.csv(all_results_safe, "lenti_subsampling_results_safe.csv", row.names = FALSE)

cat("All safe subsampling plots saved!\n")

library(ggplot2)
library(dplyr)

# ── 1. Load original counts ───────────────────────────────────────────────────
counts_hsc  <- read.table("LT-HSCs_gene_nUMI_per_cell_for_Foivos.txt",
                          header = TRUE, row.names = 1, sep = "\t")
counts_mono <- read.table("classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt",
                          header = TRUE, row.names = 1, sep = "\t")

# ── 2. Compute Lenti+ % per sample ───────────────────────────────────────────
compute_lenti_freq <- function(counts, cell_type_label) {
  total_umis  <- colSums(counts)
  lenti_umis  <- as.numeric(counts["Lenti", ])
  lenti_ratio <- lenti_umis / total_umis
  lenti_pos   <- lenti_ratio >= 4e-5
  sample_ids  <- sub("_.*", "", colnames(counts))
  exposure    <- sub("[0-9]+$", "", sample_ids)
  
  data.frame(cell      = colnames(counts),
             sample    = sample_ids,
             exposure  = exposure,
             lenti_pos = lenti_pos,
             cell_type = cell_type_label) %>%
    group_by(sample, exposure, cell_type) %>%
    summarise(pct_lenti = mean(lenti_pos) * 100, .groups = "drop")
}

freq_hsc  <- compute_lenti_freq(counts_hsc,  "LT-HSC")
freq_mono <- compute_lenti_freq(counts_mono, "Classical Monocytes")
freq_all  <- bind_rows(freq_hsc, freq_mono) %>%
  mutate(cell_type = factor(cell_type, levels = c("LT-HSC", "Classical Monocytes")))

# ── 3. Statistical tests ──────────────────────────────────────────────────────
run_tests <- function(freq_df, cell_type_label) {
  df <- freq_df %>% filter(cell_type == cell_type_label)
  hb <- df$pct_lenti[df$exposure == "HB"]
  lb <- df$pct_lenti[df$exposure == "LB"]
  
  mw <- wilcox.test(hb, lb, exact = FALSE)
  wt <- t.test(hb, lb, var.equal = FALSE)
  tt <- t.test(hb, lb, var.equal = TRUE)
  
  cat("\n===", cell_type_label, "- HB vs LB ===\n")
  cat("Mann-Whitney U p-value:  ", signif(mw$p.value, 3), "\n")
  cat("Welch's t-test p-value:  ", signif(wt$p.value, 3), "\n")
  cat("Student's t-test p-value:", signif(tt$p.value, 3), "\n")
  
  data.frame(cell_type = cell_type_label,
             p_mw      = mw$p.value,
             p_welch   = wt$p.value,
             p_student = tt$p.value)
}

stats_hsc  <- run_tests(freq_all, "LT-HSC")
stats_mono <- run_tests(freq_all, "Classical Monocytes")
stats_all  <- bind_rows(stats_hsc, stats_mono)

# ── 4. Summary stats ──────────────────────────────────────────────────────────
summary_all <- freq_all %>%
  group_by(cell_type, exposure) %>%
  summarise(mean_pct = mean(pct_lenti),
            sem_pct  = sd(pct_lenti) / sqrt(n()),
            .groups  = "drop")

# ── 5. Format p-values ────────────────────────────────────────────────────────
fmt_p <- function(p) {
  if (p < 0.001) return("p<0.001")
  paste0("p=", signif(p, 2))
}

make_label <- function(ct) {
  s <- stats_all %>% filter(cell_type == ct)
  paste0(
    "Mann-Whitney: ", fmt_p(s$p_mw),      "\n",
    "Welch's t:    ", fmt_p(s$p_welch),   "\n",
    "Student's t:  ", fmt_p(s$p_student)
  )
}

# ── 6. Bracket and annotation positions ───────────────────────────────────────
bracket_df <- freq_all %>%
  group_by(cell_type) %>%
  summarise(y_bracket = max(pct_lenti) * 1.12, .groups = "drop")

annot_df <- bracket_df %>%
  mutate(label = sapply(cell_type, make_label),
         y_pos = y_bracket * 1.03)

# ── 7. Plot ───────────────────────────────────────────────────────────────────
p_final <- ggplot() +
  # Bars
  geom_bar(data = summary_all,
           aes(x = exposure, y = mean_pct, fill = exposure),
           stat = "identity", width = 0.55,
           color = "white", linewidth = 0.3, alpha = 0.85) +
  # Error bars
  geom_errorbar(data = summary_all,
                aes(x = exposure,
                    ymin = mean_pct - sem_pct,
                    ymax = mean_pct + sem_pct),
                width = 0.12, linewidth = 0.9) +
  # Individual sample points
  geom_point(data = freq_all,
             aes(x = exposure, y = pct_lenti, fill = exposure),
             shape = 21, size = 3, color = "white",
             position = position_jitter(width = 0.07, seed = 42)) +
  # Bracket horizontal line
  geom_segment(data = bracket_df,
               aes(x = 1, xend = 2,
                   y = y_bracket, yend = y_bracket),
               inherit.aes = FALSE,
               linewidth = 0.6, color = "grey30") +
  # Bracket left tick
  geom_segment(data = bracket_df,
               aes(x = 1, xend = 1,
                   y = y_bracket * 0.97, yend = y_bracket),
               inherit.aes = FALSE,
               linewidth = 0.6, color = "grey30") +
  # Bracket right tick
  geom_segment(data = bracket_df,
               aes(x = 2, xend = 2,
                   y = y_bracket * 0.97, yend = y_bracket),
               inherit.aes = FALSE,
               linewidth = 0.6, color = "grey30") +
  # Stats text
  geom_text(data = annot_df,
            aes(x = 1.5, y = y_pos, label = label),
            size = 2.8, hjust = 0.5, vjust = 0,
            family = "mono", color = "grey20") +
  scale_fill_manual(values = c("HB" = "#F44336", "LB" = "#2196F3")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.5))) +
  facet_wrap(~ cell_type, scales = "free_y") +
  labs(title    = "Lenti+ frequency in LT-HSCs and Classical Monocytes",
       subtitle = "Mean +/- SEM | Individual samples shown | Cutoff: ratio >= 4e-5",
       x = NULL, y = "% Lenti+ cells") +
  theme_classic() +
  theme(legend.position = "none",
        strip.text      = element_text(face = "bold", size = 11),
        plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(size = 9, color = "grey40"),
        axis.text.x     = element_text(size = 11, face = "bold"),
        axis.ticks.x    = element_blank(),
        panel.spacing   = unit(2, "lines"))

ggsave("lenti_freq_HBvsLB_stats.pdf", plot = p_final, width = 10, height = 7)
write.csv(stats_all, "lenti_freq_stats.csv", row.names = FALSE)

cat("Plot and stats saved!\n")


library(dplyr)

# ── UMI depth distribution per cell per sample ────────────────────────────────
make_umi_df <- function(counts, cell_type_label) {
  data.frame(
    sample    = sub("_.*", "", colnames(counts)),
    exposure  = sub("[0-9]+$", "", sub("_.*", "", colnames(counts))),
    total_umis = colSums(counts),
    cell_type = cell_type_label
  )
}

umi_df <- bind_rows(
  make_umi_df(counts_hsc,  "LT-HSC"),
  make_umi_df(counts_mono, "Classical Monocytes")
) %>%
  mutate(cell_type = factor(cell_type, levels = c("LT-HSC", "Classical Monocytes")),
         sample    = factor(sample, levels = sort(unique(sample))))

# ── Plot 1: Violin per sample faceted by cell type ────────────────────────────
p_umi_violin <- ggplot(umi_df, aes(x = sample, y = total_umis, fill = exposure)) +
  geom_violin(scale = "width", alpha = 0.75, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.9) +
  scale_fill_manual(values = c("HB" = "#F44336", "LB" = "#2196F3")) +
  scale_y_log10(labels = scales::comma) +
  facet_wrap(~ cell_type, scales = "free_x", nrow = 1) +
  labs(title    = "Total UMIs per cell per sample",
       subtitle = "Y-axis log scale | White boxplot shows median and IQR",
       x = NULL, y = "Total UMIs per cell",
       fill = "Exposure") +
  theme_classic() +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 9),
        plot.title      = element_text(face = "bold", size = 12),
        plot.subtitle   = element_text(size = 9, color = "grey40"),
        legend.position = "top",
        strip.text      = element_text(face = "bold", size = 11))

# ── Plot 2: Median UMI per sample as dotplot ──────────────────────────────────
median_umi <- umi_df %>%
  group_by(cell_type, sample, exposure) %>%
  summarise(median_umis = median(total_umis),
            mean_umis   = mean(total_umis),
            .groups = "drop")

p_umi_median <- ggplot(median_umi, aes(x = sample, y = median_umis, fill = exposure)) +
  geom_bar(stat = "identity", width = 0.7,
           color = "white", linewidth = 0.3, alpha = 0.85) +
  geom_text(aes(label = scales::comma(round(median_umis))),
            vjust = -0.4, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("HB" = "#F44336", "LB" = "#2196F3")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                     labels = scales::comma) +
  facet_wrap(~ cell_type, scales = "free_x", nrow = 1) +
  labs(title    = "Median UMIs per cell per sample",
       subtitle = "Values shown above bars",
       x = NULL, y = "Median UMIs per cell",
       fill = "Exposure") +
  theme_classic() +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 9),
        plot.title      = element_text(face = "bold", size = 12),
        plot.subtitle   = element_text(size = 9, color = "grey40"),
        legend.position = "top",
        strip.text      = element_text(face = "bold", size = 11))

# ── Combined ──────────────────────────────────────────────────────────────────
p_umi_combined <- p_umi_violin / p_umi_median

ggsave("umi_depth_violin_per_sample.pdf",  plot = p_umi_violin,   width = 14, height = 6)
ggsave("umi_depth_median_per_sample.pdf",  plot = p_umi_median,   width = 14, height = 6)
ggsave("umi_depth_combined.pdf",           plot = p_umi_combined, width = 14, height = 12)

cat("UMI depth plots saved!\n")
