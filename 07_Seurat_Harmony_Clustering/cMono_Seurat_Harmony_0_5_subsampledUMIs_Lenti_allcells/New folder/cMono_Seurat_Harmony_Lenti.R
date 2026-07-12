# ── Clean up ──────────────────────────────────────────────────────────────────
rm(list = setdiff(ls(), "counts_mono"))
gc()

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(MASS)

# ── 1. Classify Lenti+/Lenti- BEFORE downsampling ─────────────────────────────
lenti_umis  <- as.numeric(counts_mono["Lenti", ])
total_umis  <- colSums(counts_mono)
lenti_ratio <- lenti_umis / total_umis

lenti_status <- ifelse(lenti_ratio >= 4e-5, "Lenti+",
                       ifelse(lenti_ratio <= 9e-7, "Lenti-", "Unclassified"))

cat("Lenti classification (before downsampling):\n")
print(table(lenti_status))
print(round(prop.table(table(lenti_status)) * 100, 1))

# Per sample breakdown
cat("\nLenti status per sample:\n")
print(table(sub("_.*", "", colnames(counts_mono)), lenti_status))

# ── 2. Downsample to 4149 UMIs ────────────────────────────────────────────────
target_umi <- 4149

downsample_cell <- function(cell_counts, target) {
  total <- sum(cell_counts)
  if (total <= target) return(cell_counts)
  sampled <- sample(rep(seq_along(cell_counts), cell_counts), target)
  result  <- tabulate(sampled, nbins = length(cell_counts))
  return(result)
}

set.seed(42)
cat("Downsampling all cells to", target_umi, "UMIs...\n")
counts_ds_mono <- apply(counts_mono, 2, downsample_cell, target = target_umi)
rownames(counts_ds_mono) <- rownames(counts_mono)
colnames(counts_ds_mono) <- colnames(counts_mono)

cat("Median UMIs per sample after downsampling:\n")
print(tapply(colSums(counts_ds_mono), sub("_.*", "", colnames(counts_ds_mono)), median))

# ── 3. Create Seurat object ───────────────────────────────────────────────────
seu_mono <- CreateSeuratObject(counts = counts_ds_mono,
                               min.cells = 0,
                               min.features = 0,
                               project = "Classical_Mono")

seu_mono$sample       <- sub("_.*", "", colnames(seu_mono))
seu_mono$exposure     <- sub("[0-9]+$", "", seu_mono$sample)
seu_mono$batch        <- ifelse(seu_mono$sample %in% c("LB1", "LB2", "HB1", "HB2"), "batch1",
                                ifelse(seu_mono$sample %in% c("LB3", "LB4", "HB3", "HB4"), "batch2",
                                       ifelse(seu_mono$sample == "HB5", "batch3", NA)))
seu_mono$lenti_status <- lenti_status  # add Lenti classification

cat("Cells per sample:\n")
print(table(seu_mono$sample))
cat("Lenti status in Seurat object:\n")
print(table(seu_mono$lenti_status))

# ── 4. Normalize, variable features, scale, PCA ───────────────────────────────
seu_mono <- NormalizeData(seu_mono)
seu_mono <- FindVariableFeatures(seu_mono, selection.method = "vst", nfeatures = 2000)
seu_mono <- ScaleData(seu_mono, features = rownames(seu_mono))
seu_mono <- RunPCA(seu_mono, npcs = 50, verbose = FALSE)

# ── 5. Harmony ────────────────────────────────────────────────────────────────
seu_mono <- RunHarmony(seu_mono, group.by.vars = "batch", dims.use = 1:30)

# ── 6. Cluster and UMAP ───────────────────────────────────────────────────────
seu_mono <- FindNeighbors(seu_mono, reduction = "harmony", dims = 1:30)
seu_mono <- FindClusters(seu_mono, resolution = 0.5)
seu_mono <- RunUMAP(seu_mono, reduction = "harmony", dims = 1:30)

# ── 7. Main UMAP plots ────────────────────────────────────────────────────────
p1 <- DimPlot(seu_mono, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Clusters")
p2 <- DimPlot(seu_mono, reduction = "umap", group.by = "sample") +
  ggtitle("By Sample")
p3 <- DimPlot(seu_mono, reduction = "umap", group.by = "exposure",
              cols = c("LB" = "#2196F3", "HB" = "#F44336")) +
  ggtitle("LB vs HB")
p4 <- DimPlot(seu_mono, reduction = "umap", group.by = "seurat_clusters",
              split.by = "exposure", label = TRUE) +
  ggtitle("Split by Exposure")

# Lenti plots
p5_lenti <- DimPlot(seu_mono, reduction = "umap", group.by = "lenti_status",
                    cols = c("Lenti+" = "#2CA02C", "Lenti-" = "#9467BD",
                             "Unclassified" = "#CCCCCC"),
                    order = c("Lenti+", "Lenti-", "Unclassified")) +
  ggtitle("Lenti+ vs Lenti-")

p6_lenti <- DimPlot(seu_mono, reduction = "umap", group.by = "lenti_status",
                    split.by = "exposure",
                    cols = c("Lenti+" = "#2CA02C", "Lenti-" = "#9467BD",
                             "Unclassified" = "#CCCCCC"),
                    order = c("Lenti+", "Lenti-", "Unclassified")) +
  ggtitle("Lenti status split by Exposure")

p7_lenti <- DimPlot(seu_mono, reduction = "umap", group.by = "lenti_status",
                    split.by = "sample",
                    cols = c("Lenti+" = "#2CA02C", "Lenti-" = "#9467BD",
                             "Unclassified" = "#CCCCCC"),
                    order = c("Lenti+", "Lenti-", "Unclassified"),
                    ncol = 5) +
  ggtitle("Lenti status per sample")

# Combined overview
(p1 | p3) / (p2 | p4)
p5_lenti
(p5_lenti | p6_lenti)

# ── 8. Save main plots ────────────────────────────────────────────────────────
p_combined <- (p1 | p3) / (p2 | p4)
ggsave("mono_00_combined.pdf",            plot = p_combined,         width = 16, height = 12)
ggsave("mono_01_clusters.pdf",            plot = p1,                 width = 8,  height = 6)
ggsave("mono_02_by_sample.pdf",           plot = p2,                 width = 8,  height = 6)
ggsave("mono_03_LB_vs_HB.pdf",            plot = p3,                 width = 8,  height = 6)
ggsave("mono_04_split_by_exposure.pdf",   plot = p4,                 width = 12, height = 6)
ggsave("mono_10_lenti_status.pdf",        plot = p5_lenti,           width = 8,  height = 6)
ggsave("mono_11_lenti_split_exposure.pdf",plot = p6_lenti,           width = 12, height = 6)
ggsave("mono_12_lenti_per_sample.pdf",    plot = p7_lenti,           width = 25, height = 10)

# ── 9. Per-sample highlight plots (HB and LB) ─────────────────────────────────
plot_sample_clean <- function(seurat_obj, sample_name) {
  umap_df <- as.data.frame(Embeddings(seurat_obj, "umap"))
  umap_df$highlight <- ifelse(rownames(umap_df) %in%
                                colnames(seurat_obj)[seurat_obj$sample == sample_name],
                              "yes", "no")
  umap_df <- umap_df[order(umap_df$highlight), ]
  ggplot(umap_df, aes(x = umap_1, y = umap_2, color = highlight,
                      size = highlight, alpha = highlight)) +
    geom_point() +
    scale_color_manual(values = c("no" = "#CCCCCC", "yes" = "#E63946")) +
    scale_size_manual(values  = c("no" = 0.1,       "yes" = 0.8)) +
    scale_alpha_manual(values = c("no" = 0.3,       "yes" = 1.0)) +
    ggtitle(sample_name) +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(size = 12, face = "bold"))
}

hb_samples <- sort(unique(seu_mono$sample[seu_mono$exposure == "HB"]))
lb_samples <- sort(unique(seu_mono$sample[seu_mono$exposure == "LB"]))

p_HB <- wrap_plots(lapply(hb_samples, plot_sample_clean, seurat_obj = seu_mono), nrow = 1) +
  plot_annotation(title = "HB samples - individual projection")
p_LB <- wrap_plots(lapply(lb_samples, plot_sample_clean, seurat_obj = seu_mono), nrow = 1) +
  plot_annotation(title = "LB samples - individual projection")

ggsave("mono_05_HB_per_sample.pdf", plot = p_HB, width = 5 * length(hb_samples), height = 5)
ggsave("mono_06_LB_per_sample.pdf", plot = p_LB, width = 5 * length(lb_samples), height = 5)

# ── 10. Helper: density difference plots ──────────────────────────────────────
density_diff_plots <- function(seurat_obj, group_col, group_a, group_b,
                               label_a, label_b, color_a, color_b, prefix) {
  umap_coords <- as.data.frame(Embeddings(seurat_obj, "umap"))
  umap_coords$group <- seurat_obj[[group_col]][, 1]
  
  cells_a <- umap_coords[umap_coords$group == group_a, ]
  cells_b <- umap_coords[umap_coords$group == group_b, ]
  
  x_range <- range(umap_coords$umap_1)
  y_range <- range(umap_coords$umap_2)
  n_grid  <- 200
  epsilon <- 1e-10
  
  dens_a <- kde2d(cells_a$umap_1, cells_a$umap_2, n = n_grid, lims = c(x_range, y_range))
  dens_b <- kde2d(cells_b$umap_1, cells_b$umap_2, n = n_grid, lims = c(x_range, y_range))
  
  dens_a$z <- dens_a$z / sum(dens_a$z)
  dens_b$z <- dens_b$z / sum(dens_b$z)
  
  n_a <- nrow(cells_a)
  n_b <- nrow(cells_b)
  counts_a <- dens_a$z * n_a
  counts_b <- dens_b$z * n_b
  total    <- counts_a + counts_b
  
  df <- expand.grid(umap_1 = dens_a$x, umap_2 = dens_a$y)
  df$diff         <- as.vector(dens_b$z - dens_a$z)
  df$log2fc       <- as.vector(log2((dens_b$z + epsilon) / (dens_a$z + epsilon)))
  df$log2fc_capped <- pmax(pmin(df$log2fc, 3), -3)
  df$pct_b        <- as.vector(counts_b / (total + epsilon) * 100)
  
  subtitle <- paste0("Blue = ", label_b, " enriched  |  Red = ", label_a, " enriched")
  
  # Plot 1: normalized density difference
  p1 <- ggplot(df, aes(x = umap_1, y = umap_2, fill = diff)) +
    geom_tile() +
    scale_fill_gradient2(low = color_a, mid = "white", high = color_b,
                         midpoint = 0, name = paste0(label_b, " - ", label_a, "\ndensity")) +
    geom_point(data = umap_coords, aes(x = umap_1, y = umap_2),
               inherit.aes = FALSE, size = 0.05, alpha = 0.1, color = "black") +
    labs(title = paste0("Differential density: ", label_b, " vs ", label_a),
         subtitle = subtitle,
         caption = "Units: difference in normalized probability density (each group sums to 1).\nPositive values = higher relative probability of blue group; negative = higher relative probability of red group.") +
    theme_classic() +
    theme(plot.caption = element_text(size = 8, color = "grey40", hjust = 0))
  ggsave(paste0(prefix, "_07_differential_density.pdf"), plot = p1, width = 8, height = 6)
  
  # Plot 2: log2 fold-change
  p2 <- ggplot(df, aes(x = umap_1, y = umap_2, fill = log2fc_capped)) +
    geom_tile() +
    scale_fill_gradient2(low = color_a, mid = "white", high = color_b,
                         midpoint = 0, limits = c(-3, 3),
                         name = paste0("log2(", label_b, "/", label_a, ")\n(capped ±3)")) +
    geom_point(data = umap_coords, aes(x = umap_1, y = umap_2),
               inherit.aes = FALSE, size = 0.05, alpha = 0.1, color = "black") +
    labs(title = paste0("Log2 fold-change density: ", label_b, " vs ", label_a),
         subtitle = subtitle) +
    theme_classic()
  ggsave(paste0(prefix, "_08_log2fc_density.pdf"), plot = p2, width = 8, height = 6)
  
  # Plot 3: % of group_b cells
  p3 <- ggplot(df, aes(x = umap_1, y = umap_2, fill = pct_b)) +
    geom_tile() +
    scale_fill_gradient2(low = color_a, mid = "white", high = color_b,
                         midpoint = 50, limits = c(0, 100),
                         name = paste0("% ", label_b, " cells")) +
    geom_point(data = umap_coords, aes(x = umap_1, y = umap_2),
               inherit.aes = FALSE, size = 0.05, alpha = 0.1, color = "black") +
    labs(title = paste0("% ", label_b, " cells per UMAP region"),
         subtitle = subtitle) +
    theme_classic()
  ggsave(paste0(prefix, "_09_pct_density.pdf"), plot = p3, width = 8, height = 6)
  
  cat("Density plots saved for", label_a, "vs", label_b, "\n")
}

# ── 11. Run density plots for HB vs LB ────────────────────────────────────────
density_diff_plots(seu_mono,
                   group_col = "exposure",
                   group_a = "HB", group_b = "LB",
                   label_a = "HB", label_b = "LB",
                   color_a = "#F44336", color_b = "#2196F3",
                   prefix = "mono")

# ── 12. Run density plots for Lenti- vs Lenti+ ────────────────────────────────
# Only use classified cells (exclude Unclassified)
seu_lenti <- subset(seu_mono, subset = lenti_status %in% c("Lenti+", "Lenti-"))
cat("Cells used for Lenti density plots:", ncol(seu_lenti), "\n")

density_diff_plots(seu_lenti,
                   group_col = "lenti_status",
                   group_a = "Lenti-", group_b = "Lenti+",
                   label_a = "Lenti-", label_b = "Lenti+",
                   color_a = "#9467BD", color_b = "#2CA02C",
                   prefix = "mono_lenti")

# ── 13. Save Seurat object ────────────────────────────────────────────────────
saveRDS(seu_mono, "classical_mono_seurat_final.rds")

cat("All plots saved!\n")
