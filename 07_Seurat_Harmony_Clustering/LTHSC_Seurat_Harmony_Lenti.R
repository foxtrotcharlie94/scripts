# ── Clean up ──────────────────────────────────────────────────────────────────
rm(list = ls())
gc()

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(MASS)

# ── 1. Load data ──────────────────────────────────────────────────────────────
counts_hsc <- read.table("LT-HSCs_gene_nUMI_per_cell_for_Foivos.txt",
                         header = TRUE, row.names = 1, sep = "\t")
cat("Dimensions:", dim(counts_hsc), "\n")

# ── 2. Classify Lenti+/Lenti- BEFORE downsampling ─────────────────────────────
lenti_umis  <- as.numeric(counts_hsc["Lenti", ])
total_umis  <- colSums(counts_hsc)
lenti_ratio <- lenti_umis / total_umis

lenti_status <- ifelse(lenti_ratio >= 4e-5, "Lenti+",
                       ifelse(lenti_ratio <= 9e-7, "Lenti-", "Unclassified"))

cat("Lenti classification (before downsampling):\n")
print(table(lenti_status))
print(round(prop.table(table(lenti_status)) * 100, 1))

cat("\nLenti status per sample:\n")
print(table(sub("_.*", "", colnames(counts_hsc)), lenti_status))

# ── 3. Downsample to 5000 UMIs ────────────────────────────────────────────────
target_umi <- 5000

downsample_cell <- function(cell_counts, target) {
  total <- sum(cell_counts)
  if (total <= target) return(cell_counts)
  sampled <- sample(rep(seq_along(cell_counts), cell_counts), target)
  result  <- tabulate(sampled, nbins = length(cell_counts))
  return(result)
}

set.seed(42)
cat("Downsampling all cells to", target_umi, "UMIs...\n")
counts_ds_hsc <- apply(counts_hsc, 2, downsample_cell, target = target_umi)
rownames(counts_ds_hsc) <- rownames(counts_hsc)
colnames(counts_ds_hsc) <- colnames(counts_hsc)

cat("Median UMIs per sample after downsampling:\n")
print(tapply(colSums(counts_ds_hsc), sub("_.*", "", colnames(counts_ds_hsc)), median))

# ── 4. Create Seurat object ───────────────────────────────────────────────────
seu_hsc <- CreateSeuratObject(counts = counts_ds_hsc,
                              min.cells = 0,
                              min.features = 0,
                              project = "LT-HSC")

seu_hsc$sample       <- sub("_.*", "", colnames(seu_hsc))
seu_hsc$exposure     <- sub("[0-9]+$", "", seu_hsc$sample)
seu_hsc$batch        <- ifelse(seu_hsc$sample %in% c("LB1", "LB2", "HB1", "HB2"), "batch1",
                               ifelse(seu_hsc$sample %in% c("LB3", "LB4", "HB3", "HB4"), "batch2",
                                      ifelse(seu_hsc$sample == "HB5", "batch3", NA)))
seu_hsc$lenti_status <- lenti_status

cat("Cells per sample:\n")
print(table(seu_hsc$sample))
cat("Lenti status in Seurat object:\n")
print(table(seu_hsc$lenti_status))

# ── 5. Normalize, variable features, scale, PCA ───────────────────────────────
seu_hsc <- NormalizeData(seu_hsc)
seu_hsc <- FindVariableFeatures(seu_hsc, selection.method = "vst", nfeatures = 2000)
seu_hsc <- ScaleData(seu_hsc, features = rownames(seu_hsc))
seu_hsc <- RunPCA(seu_hsc, npcs = 50, verbose = FALSE)

# ── 6. Harmony ────────────────────────────────────────────────────────────────
seu_hsc <- RunHarmony(seu_hsc, group.by.vars = "batch", dims.use = 1:30)

# ── 7. Cluster and UMAP ───────────────────────────────────────────────────────
seu_hsc <- FindNeighbors(seu_hsc, reduction = "harmony", dims = 1:30)
seu_hsc <- FindClusters(seu_hsc, resolution = 0.5)
seu_hsc <- RunUMAP(seu_hsc, reduction = "harmony", dims = 1:30)

# ── 8. Main UMAP plots ────────────────────────────────────────────────────────
p1 <- DimPlot(seu_hsc, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Clusters")
p2 <- DimPlot(seu_hsc, reduction = "umap", group.by = "sample") +
  ggtitle("By Sample")
p3 <- DimPlot(seu_hsc, reduction = "umap", group.by = "exposure",
              cols = c("LB" = "#2196F3", "HB" = "#F44336")) +
  ggtitle("LB vs HB")
p4 <- DimPlot(seu_hsc, reduction = "umap", group.by = "seurat_clusters",
              split.by = "exposure", label = TRUE) +
  ggtitle("Split by Exposure")

p5_lenti <- DimPlot(seu_hsc, reduction = "umap", group.by = "lenti_status",
                    cols = c("Lenti+" = "#2CA02C", "Lenti-" = "#9467BD",
                             "Unclassified" = "#CCCCCC"),
                    order = c("Lenti+", "Lenti-", "Unclassified")) +
  ggtitle("Lenti+ vs Lenti-")

p6_lenti <- DimPlot(seu_hsc, reduction = "umap", group.by = "lenti_status",
                    split.by = "exposure",
                    cols = c("Lenti+" = "#2CA02C", "Lenti-" = "#9467BD",
                             "Unclassified" = "#CCCCCC"),
                    order = c("Lenti+", "Lenti-", "Unclassified")) +
  ggtitle("Lenti status split by Exposure")

p7_lenti <- DimPlot(seu_hsc, reduction = "umap", group.by = "lenti_status",
                    split.by = "sample",
                    cols = c("Lenti+" = "#2CA02C", "Lenti-" = "#9467BD",
                             "Unclassified" = "#CCCCCC"),
                    order = c("Lenti+", "Lenti-", "Unclassified"),
                    ncol = 5) +
  ggtitle("Lenti status per sample")

(p1 | p3) / (p2 | p4)
p5_lenti
(p5_lenti | p6_lenti)

# ── 9. Cluster proportions ────────────────────────────────────────────────────
cluster_sample_table <- prop.table(table(seu_hsc$sample, seu_hsc$seurat_clusters), margin = 1) * 100
round(cluster_sample_table, 1)

# ── 10. Save main plots ───────────────────────────────────────────────────────
p_combined <- (p1 | p3) / (p2 | p4)
ggsave("hsc_00_combined.pdf",             plot = p_combined, width = 16, height = 12)
ggsave("hsc_01_clusters.pdf",             plot = p1,         width = 8,  height = 6)
ggsave("hsc_02_by_sample.pdf",            plot = p2,         width = 8,  height = 6)
ggsave("hsc_03_LB_vs_HB.pdf",             plot = p3,         width = 8,  height = 6)
ggsave("hsc_04_split_by_exposure.pdf",    plot = p4,         width = 12, height = 6)
ggsave("hsc_10_lenti_status.pdf",         plot = p5_lenti,   width = 8,  height = 6)
ggsave("hsc_11_lenti_split_exposure.pdf", plot = p6_lenti,   width = 12, height = 6)
ggsave("hsc_12_lenti_per_sample.pdf",     plot = p7_lenti,   width = 25, height = 10)

# ── 11. Per-sample highlight plots ────────────────────────────────────────────
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

hb_samples <- sort(unique(seu_hsc$sample[seu_hsc$exposure == "HB"]))
lb_samples <- sort(unique(seu_hsc$sample[seu_hsc$exposure == "LB"]))

p_HB <- wrap_plots(lapply(hb_samples, plot_sample_clean, seurat_obj = seu_hsc), nrow = 1) +
  plot_annotation(title = "HB samples - individual projection")
p_LB <- wrap_plots(lapply(lb_samples, plot_sample_clean, seurat_obj = seu_hsc), nrow = 1) +
  plot_annotation(title = "LB samples - individual projection")

ggsave("hsc_05_HB_per_sample.pdf", plot = p_HB, width = 5 * length(hb_samples), height = 5)
ggsave("hsc_06_LB_per_sample.pdf", plot = p_LB, width = 5 * length(lb_samples), height = 5)

# ── 12. Density difference helper function ────────────────────────────────────
density_diff_plots <- function(seurat_obj, group_col, group_a, group_b,
                               label_a, label_b, color_a, color_b, 
                               color_name_a = "Red", color_name_b = "Blue",  # NEW
                               prefix) {
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
  
  n_a      <- nrow(cells_a)
  n_b      <- nrow(cells_b)
  counts_a <- dens_a$z * n_a
  counts_b <- dens_b$z * n_b
  total    <- counts_a + counts_b
  
  df <- expand.grid(umap_1 = dens_a$x, umap_2 = dens_a$y)
  df$diff          <- as.vector(dens_b$z - dens_a$z)
  df$log2fc        <- as.vector(log2((dens_b$z + epsilon) / (dens_a$z + epsilon)))
  df$log2fc_capped <- pmax(pmin(df$log2fc, 3), -3)
  df$pct_b         <- as.vector(counts_b / (total + epsilon) * 100)
  
  # Now uses the actual color names
  subtitle <- paste0(color_name_b, " = ", label_b, " enriched  |  ", 
                     color_name_a, " = ", label_a, " enriched")
  
  p1 <- ggplot(df, aes(x = umap_1, y = umap_2, fill = diff)) +
    geom_tile() +
    scale_fill_gradient2(low = color_a, mid = "white", high = color_b,
                         midpoint = 0, name = paste0(label_b, " - ", label_a, "\ndensity")) +
    geom_point(data = umap_coords, aes(x = umap_1, y = umap_2),
               inherit.aes = FALSE, size = 0.05, alpha = 0.1, color = "black") +
    labs(title = paste0("Differential density: ", label_b, " vs ", label_a),
         subtitle = subtitle,
         caption = "Units: difference in normalized probability density (each group sums to 1).\nPositive values = higher relative probability of enriched group; negative = higher relative probability of other group.") +
    theme_classic() +
    theme(plot.caption = element_text(size = 8, color = "grey40", hjust = 0))
  ggsave(paste0(prefix, "_07_differential_density.pdf"), plot = p1, width = 8, height = 6)
  
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

# ── HB vs LB (Blue and Red) ───────────────────────────────────────────────────
density_diff_plots(seu_hsc,
                   group_col = "exposure",
                   group_a = "HB", group_b = "LB",
                   label_a = "HB", label_b = "LB",
                   color_a = "#F44336", color_b = "#2196F3",
                   color_name_a = "Red", color_name_b = "Blue",
                   prefix = "hsc")

# ── Lenti+/Lenti- (Green and Purple) ─────────────────────────────────────────
seu_lenti <- subset(seu_hsc, subset = lenti_status %in% c("Lenti+", "Lenti-"))

density_diff_plots(seu_lenti,
                   group_col = "lenti_status",
                   group_a = "Lenti-", group_b = "Lenti+",
                   label_a = "Lenti-", label_b = "Lenti+",
                   color_a = "#9467BD", color_b = "#2CA02C",
                   color_name_a = "Purple", color_name_b = "Green",
                   prefix = "hsc_lenti")

