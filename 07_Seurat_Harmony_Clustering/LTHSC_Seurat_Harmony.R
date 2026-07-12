# ── 0. Clean up ───────────────────────────────────────────────────────────────
rm(list = ls())
gc()
plan("sequential")

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(MASS)

# ── 1. Load full data ─────────────────────────────────────────────────────────
counts <- read.table("LT-HSCs_gene_nUMI_per_cell_for_Foivos.txt",
                     header = TRUE, row.names = 1, sep = "\t")
cat("Full matrix dimensions:", dim(counts), "\n")

# ── 2. Check UMI depth per sample ────────────────────────────────────────────
umi_per_cell   <- colSums(counts)
sample_ids_all <- sub("_.*", "", colnames(counts))
cat("Median UMIs per sample (before downsampling):\n")
print(tapply(umi_per_cell, sample_ids_all, median))

# ── 3. Downsample all cells to 5000 UMIs ─────────────────────────────────────
target_umi <- 5000

downsample_cell <- function(cell_counts, target) {
  total <- sum(cell_counts)
  if (total <= target) return(cell_counts)
  sampled <- sample(rep(seq_along(cell_counts), cell_counts), target)
  result  <- tabulate(sampled, nbins = length(cell_counts))
  return(result)
}

set.seed(42)
cat("Downsampling all cells to", target_umi, "UMIs... (this may take a while)\n")
counts_ds <- apply(counts, 2, downsample_cell, target = target_umi)
rownames(counts_ds) <- rownames(counts)
colnames(counts_ds) <- colnames(counts)

# Verify
cat("Median UMIs per sample after downsampling:\n")
print(tapply(colSums(counts_ds), sample_ids_all, median))

# ── 4. Create Seurat object ───────────────────────────────────────────────────
seu <- CreateSeuratObject(counts = counts_ds,
                          min.cells = 0,
                          min.features = 0,
                          project = "LT-HSC")

seu$sample   <- sub("_.*", "", colnames(seu))
seu$exposure <- sub("[0-9]+$", "", seu$sample)
seu$batch    <- ifelse(seu$sample %in% c("LB1", "LB2", "HB1", "HB2"), "batch1",
                       ifelse(seu$sample %in% c("LB3", "LB4", "HB3", "HB4"), "batch2",
                              ifelse(seu$sample == "HB5", "batch3", NA)))

cat("Cells per sample:\n")
print(table(seu$sample))
cat("Cells per batch:\n")
print(table(seu$batch))

# ── 5. Normalize, variable features, scale, PCA ───────────────────────────────
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)

ElbowPlot(seu, ndims = 50)
n_pcs <- 30

# ── Harmony ───────────────────────────────────────────────────────────────────
seu <- RunHarmony(seu, group.by.vars = "batch", dims.use = 1:n_pcs)

# ── Cluster and UMAP ──────────────────────────────────────────────────────────
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:n_pcs)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:n_pcs)

# ── Plots ─────────────────────────────────────────────────────────────────────
p1 <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("Clusters")
p2 <- DimPlot(seu, reduction = "umap", group.by = "sample") + ggtitle("By Sample")
p3 <- DimPlot(seu, reduction = "umap", group.by = "exposure",
              cols = c("LB" = "#2196F3", "HB" = "#F44336")) + ggtitle("LB vs HB")
p4 <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters",
              split.by = "exposure", label = TRUE) + ggtitle("Split by Exposure")

(p1 | p3) / (p2 | p4)

# ── Cluster proportions ───────────────────────────────────────────────────────
cluster_sample_table <- prop.table(table(seu$sample, seu$seurat_clusters), margin = 1) * 100
round(cluster_sample_table, 1)

# ── 1. Save combined 4-plot figure ────────────────────────────────────────────
p_combined <- (p1 | p3) / (p2 | p4)
ggsave("00_combined.pdf",          plot = p_combined, width = 16, height = 12)

# ── 2. Save each plot individually ────────────────────────────────────────────
ggsave("01_clusters.pdf",          plot = p1, width = 8, height = 6)
ggsave("02_by_sample.pdf",         plot = p2, width = 8, height = 6)
ggsave("03_LB_vs_HB.pdf",          plot = p3, width = 8, height = 6)
ggsave("04_split_by_exposure.pdf", plot = p4, width = 12, height = 6)

# ── 3. Per-sample highlight plots ─────────────────────────────────────────────
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

hb_samples <- sort(unique(seu$sample[seu$exposure == "HB"]))
lb_samples <- sort(unique(seu$sample[seu$exposure == "LB"]))

hb_plots <- lapply(hb_samples, plot_sample_clean, seurat_obj = seu)
lb_plots <- lapply(lb_samples, plot_sample_clean, seurat_obj = seu)

p_HB <- wrap_plots(hb_plots, nrow = 1) + plot_annotation(title = "HB samples - individual projection")
p_LB <- wrap_plots(lb_plots, nrow = 1) + plot_annotation(title = "LB samples - individual projection")

ggsave("05_HB_per_sample.pdf", plot = p_HB, width = 5 * length(hb_samples), height = 5)
ggsave("06_LB_per_sample.pdf", plot = p_LB, width = 5 * length(lb_samples), height = 5)

# ── 4. Differential density: LB vs HB ─────────────────────────────────────────
umap_coords <- as.data.frame(Embeddings(seu, "umap"))
umap_coords$exposure <- seu$exposure

hb_cells <- umap_coords[umap_coords$exposure == "HB", ]
lb_cells  <- umap_coords[umap_coords$exposure == "LB", ]

x_range <- range(umap_coords$umap_1)
y_range <- range(umap_coords$umap_2)
n_grid  <- 200

dens_HB <- kde2d(hb_cells$umap_1, hb_cells$umap_2, n = n_grid, lims = c(x_range, y_range))
dens_LB <- kde2d(lb_cells$umap_1, lb_cells$umap_2, n = n_grid, lims = c(x_range, y_range))

dens_HB$z <- dens_HB$z / sum(dens_HB$z)
dens_LB$z <- dens_LB$z / sum(dens_LB$z)

diff_z  <- dens_LB$z - dens_HB$z
diff_df <- expand.grid(umap_1 = dens_HB$x, umap_2 = dens_HB$y)
diff_df$diff <- as.vector(diff_z)

p_diff <- ggplot(diff_df, aes(x = umap_1, y = umap_2, fill = diff)) +
  geom_tile() +
  scale_fill_gradient2(low = "#F44336", mid = "white", high = "#2196F3",
                       midpoint = 0, name = "LB - HB\ndensity") +
  geom_point(data = umap_coords, aes(x = umap_1, y = umap_2),
             inherit.aes = FALSE, size = 0.05, alpha = 0.1, color = "black") +
  labs(title = "Differential density: LB vs HB",
       subtitle = "Blue = LB enriched  |  Red = HB enriched",
       caption = "Units: difference in normalized probability density (each group sums to 1).\nPositive values = higher relative probability of LB cells; negative = higher relative probability of HB cells.\nValues close to 0 indicate equal representation of both groups.") +
  theme_classic() +
  theme(plot.caption = element_text(size = 8, color = "grey40", hjust = 0))

ggsave("07_differential_density.pdf", plot = p_diff, width = 8, height = 6)

# ── 5. Save Seurat object ─────────────────────────────────────────────────────
saveRDS(seu, "LT-HSC_seurat_final.rds")

cat("All plots saved!\n")

# ── Shared density setup ───────────────────────────────────────────────────────
umap_coords <- as.data.frame(Embeddings(seu, "umap"))
umap_coords$exposure <- seu$exposure

hb_cells <- umap_coords[umap_coords$exposure == "HB", ]
lb_cells  <- umap_coords[umap_coords$exposure == "LB", ]

x_range <- range(umap_coords$umap_1)
y_range <- range(umap_coords$umap_2)
n_grid  <- 200

dens_HB <- kde2d(hb_cells$umap_1, hb_cells$umap_2, n = n_grid, lims = c(x_range, y_range))
dens_LB <- kde2d(lb_cells$umap_1, lb_cells$umap_2, n = n_grid, lims = c(x_range, y_range))

# Normalize to sum to 1
dens_HB$z <- dens_HB$z / sum(dens_HB$z)
dens_LB$z <- dens_LB$z / sum(dens_LB$z)

# Base data frame
diff_df <- expand.grid(umap_1 = dens_HB$x, umap_2 = dens_HB$y)

# ── 1. Original: normalized density difference (LB - HB) ──────────────────────
diff_df$diff <- as.vector(dens_LB$z - dens_HB$z)

p_diff <- ggplot(diff_df, aes(x = umap_1, y = umap_2, fill = diff)) +
  geom_tile() +
  scale_fill_gradient2(low = "#F44336", mid = "white", high = "#2196F3",
                       midpoint = 0, name = "LB - HB\ndensity") +
  geom_point(data = umap_coords, aes(x = umap_1, y = umap_2),
             inherit.aes = FALSE, size = 0.05, alpha = 0.1, color = "black") +
  labs(title = "Differential density: LB vs HB",
       subtitle = "Blue = LB enriched  |  Red = HB enriched") +
  theme_classic()

ggsave("07_differential_density.pdf", plot = p_diff, width = 8, height = 6)

# ── 2. Log2 fold-change: log2(LB / HB) ────────────────────────────────────────
epsilon <- 1e-10
n_HB <- nrow(hb_cells)
n_LB <- nrow(lb_cells)

hb_counts <- dens_HB$z * n_HB
lb_counts  <- dens_LB$z * n_LB
total      <- hb_counts + lb_counts

diff_df$log2fc        <- as.vector(log2((dens_LB$z + epsilon) / (dens_HB$z + epsilon)))
diff_df$log2fc_capped <- pmax(pmin(diff_df$log2fc, 3), -3)

p_fc <- ggplot(diff_df, aes(x = umap_1, y = umap_2, fill = log2fc_capped)) +
  geom_tile() +
  scale_fill_gradient2(low = "#F44336", mid = "white", high = "#2196F3",
                       midpoint = 0,
                       limits = c(-3, 3),
                       name = "log2(LB/HB)\n(capped ±3)") +
  geom_point(data = umap_coords, aes(x = umap_1, y = umap_2),
             inherit.aes = FALSE, size = 0.05, alpha = 0.1, color = "black") +
  labs(title = "Log2 fold-change density: LB vs HB",
       subtitle = "Blue = LB enriched  |  Red = HB enriched  |  White = equal") +
  theme_classic()

ggsave("08_log2fc_density.pdf", plot = p_fc, width = 8, height = 6)

# ── 3. Percentage of cells per grid cell ──────────────────────────────────────
diff_df$pct_LB <- as.vector(lb_counts / (total + epsilon) * 100)

p_pct <- ggplot(diff_df, aes(x = umap_1, y = umap_2, fill = pct_LB)) +
  geom_tile() +
  scale_fill_gradient2(low = "#F44336", mid = "white", high = "#2196F3",
                       midpoint = 50,
                       limits = c(0, 100),
                       name = "% LB cells") +
  geom_point(data = umap_coords, aes(x = umap_1, y = umap_2),
             inherit.aes = FALSE, size = 0.05, alpha = 0.1, color = "black") +
  labs(title = "% LB cells per UMAP region",
       subtitle = "Blue = LB enriched  |  Red = HB enriched  |  White = 50/50") +
  theme_classic()

ggsave("09_pct_LB_density.pdf", plot = p_pct, width = 8, height = 6)

cat("All three density plots saved!\n")

