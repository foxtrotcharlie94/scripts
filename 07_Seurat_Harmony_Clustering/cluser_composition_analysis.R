library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# ── Set working directory ─────────────────────────────────────────────────────
setwd("C:/Users/fc809/Downloads")

# ── 1. Load objects ───────────────────────────────────────────────────────────
seu_hsc  <- readRDS("LT_HSC_seurat_final.rds")
seu_mono <- readRDS("classical_mono_seurat_final.rds")

# ── 2. Find markers ───────────────────────────────────────────────────────────
markers_hsc  <- FindAllMarkers(seu_hsc,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_mono <- FindAllMarkers(seu_mono, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# ── 3. Save top 50 per cluster as TXT files ───────────────────────────────────
save_cluster_txt <- function(markers_df, prefix) {
  clusters <- unique(markers_df$cluster)
  for (cl in clusters) {
    top50 <- markers_df %>%
      filter(cluster == cl) %>%
      slice_max(order_by = avg_log2FC, n = 50) %>%
      select(gene, avg_log2FC, pct.1, pct.2, p_val_adj)
    filename <- paste0(prefix, "_cluster", cl, "_top50_markers.txt")
    write.table(top50, filename, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("Saved:", filename, "\n")
  }
}

save_cluster_txt(markers_hsc,  "hsc")
save_cluster_txt(markers_mono, "mono")

# ── 4. Top 10 per cluster for heatmaps ───────────────────────────────────────
top10_hsc <- markers_hsc %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup()

top10_mono <- markers_mono %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup()

# ── 5. Heatmap - LT-HSC ──────────────────────────────────────────────────────
p_heat_hsc <- DoHeatmap(seu_hsc,
                        features = unique(top10_hsc$gene),
                        group.by = "seurat_clusters",
                        size = 3) +
  ggtitle("Top 10 markers per cluster - LT-HSC") +
  theme(axis.text.y = element_text(size = 7))

ggsave("hsc_13_cluster_markers_heatmap.pdf", plot = p_heat_hsc, width = 14, height = 12)

# ── 6. Heatmap - Classical Monocytes ─────────────────────────────────────────
p_heat_mono <- DoHeatmap(seu_mono,
                         features = unique(top10_mono$gene),
                         group.by = "seurat_clusters",
                         size = 3) +
  ggtitle("Top 10 markers per cluster - Classical Monocytes") +
  theme(axis.text.y = element_text(size = 7))

ggsave("mono_13_cluster_markers_heatmap.pdf", plot = p_heat_mono, width = 14, height = 12)

# ── 7. Combined heatmap ───────────────────────────────────────────────────────
p_combined_heatmap <- p_heat_hsc / p_heat_mono +
  plot_annotation(title = "Cluster markers - LT-HSC (top) and Classical Monocytes (bottom)")

ggsave("combined_cluster_markers_heatmap.pdf",
       plot = p_combined_heatmap,
       width = 14, height = 24)

cat("All marker files and heatmaps saved!\n")

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# ── 1. Cluster annotations ────────────────────────────────────────────────────
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

# ── 2. Helper function ────────────────────────────────────────────────────────
plot_cluster_frequency <- function(seurat_obj, cluster_labels, title_prefix, file_prefix) {
  
  # Build metadata table, keep only Lenti+ and Lenti-
  meta <- seurat_obj@meta.data %>%
    filter(lenti_status %in% c("Lenti+", "Lenti-")) %>%
    mutate(
      group = paste0(exposure, "-", lenti_status),
      cluster_label = cluster_labels[as.character(seurat_clusters)]
    )
  
  # Count cells per sample per cluster per group
  freq_sample <- meta %>%
    group_by(sample, group, exposure, lenti_status, cluster_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(sample, group) %>%
    mutate(freq = n / sum(n) * 100) %>%
    ungroup()
  
  # Mean + SEM per group per cluster
  freq_summary <- freq_sample %>%
    group_by(group, cluster_label) %>%
    summarise(
      mean_freq = mean(freq),
      sem_freq  = sd(freq) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(group = factor(group, levels = c("LB-Lenti+", "LB-Lenti-", "HB-Lenti+", "HB-Lenti-")))
  
  # ── Plot 1: Grouped bar chart (mean ± SEM) ──────────────────────────────────
  p_bar <- ggplot(freq_summary, aes(x = cluster_label, y = mean_freq, fill = group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = mean_freq - sem_freq, ymax = mean_freq + sem_freq),
                  position = position_dodge(width = 0.8), width = 0.25, linewidth = 0.5) +
    scale_fill_manual(values = c(
      "LB-Lenti+" = "#2196F3",
      "LB-Lenti-" = "#90CAF9",
      "HB-Lenti+" = "#F44336",
      "HB-Lenti-" = "#FFCDD2"
    )) +
    labs(title = paste0(title_prefix, " - Cluster frequency by group"),
         subtitle = "Mean ± SEM across samples",
         x = NULL, y = "% of cells", fill = "Group") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          legend.position = "top")
  
  ggsave(paste0(file_prefix, "_cluster_freq_barplot.pdf"),
         plot = p_bar, width = 14, height = 7)
  
  # ── Plot 2: Dot plot with individual samples ──────────────────────────────────
  freq_sample <- freq_sample %>%
    mutate(group = factor(group, levels = c("LB-Lenti+", "LB-Lenti-", "HB-Lenti+", "HB-Lenti-")))
  
  p_dot <- ggplot(freq_summary, aes(x = cluster_label, y = mean_freq, color = group)) +
    geom_bar(aes(fill = group), stat = "identity",
             position = position_dodge(width = 0.8), width = 0.7, alpha = 0.4) +
    geom_errorbar(aes(ymin = mean_freq - sem_freq, ymax = mean_freq + sem_freq),
                  position = position_dodge(width = 0.8), width = 0.25, linewidth = 0.5) +
    geom_point(data = freq_sample, aes(x = cluster_label, y = freq, color = group),
               position = position_dodge(width = 0.8), size = 1.5, alpha = 0.8) +
    scale_fill_manual(values = c(
      "LB-Lenti+" = "#2196F3",
      "LB-Lenti-" = "#90CAF9",
      "HB-Lenti+" = "#F44336",
      "HB-Lenti-" = "#FFCDD2"
    )) +
    scale_color_manual(values = c(
      "LB-Lenti+" = "#1565C0",
      "LB-Lenti-" = "#42A5F5",
      "HB-Lenti+" = "#B71C1C",
      "HB-Lenti-" = "#EF9A9A"
    )) +
    labs(title = paste0(title_prefix, " - Cluster frequency by group"),
         subtitle = "Mean ± SEM with individual samples overlaid",
         x = NULL, y = "% of cells", fill = "Group", color = "Group") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          legend.position = "top")
  
  ggsave(paste0(file_prefix, "_cluster_freq_dotplot.pdf"),
         plot = p_dot, width = 14, height = 7)
  
  # ── Plot 3: Stacked bar per sample ───────────────────────────────────────────
  p_stack <- ggplot(freq_sample, aes(x = sample, y = freq, fill = cluster_label)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ group, scales = "free_x", nrow = 1) +
    labs(title = paste0(title_prefix, " - Cluster composition per sample"),
         x = NULL, y = "% of cells", fill = "Cluster") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "bottom",
          legend.text = element_text(size = 8),
          strip.text = element_text(face = "bold"))
  
  ggsave(paste0(file_prefix, "_cluster_freq_stacked.pdf"),
         plot = p_stack, width = 16, height = 8)
  
  # ── Save frequency table ──────────────────────────────────────────────────────
  write.csv(freq_sample,   paste0(file_prefix, "_cluster_freq_per_sample.csv"),  row.names = FALSE)
  write.csv(freq_summary,  paste0(file_prefix, "_cluster_freq_summary.csv"),     row.names = FALSE)
  
  cat("Frequency plots saved for", title_prefix, "\n")
}

# ── 3. Run for LT-HSC ────────────────────────────────────────────────────────
seu_hsc  <- readRDS("LT_HSC_seurat_final.rds")
plot_cluster_frequency(seu_hsc, hsc_labels, "LT-HSC", "hsc")

# ── 4. Run for Classical Monocytes ───────────────────────────────────────────
seu_mono <- readRDS("classical_mono_seurat_final.rds")
plot_cluster_frequency(seu_mono, mono_labels, "Classical Monocytes", "mono")

cat("All frequency plots saved!\n")

