library(edgeR)
library(GSVA)
library(BiocParallel)
library(msigdbr)
library(ggplot2)
library(dplyr)

# ── 0. Setup ──────────────────────────────────────────────────────────────────
workers <- min(6, parallel::detectCores() - 1)
BPPARAM <- SnowParam(workers, type = "SOCK", progressbar = TRUE)
cat("Using", workers, "workers\n")

OUTDIR <- "C:/Users/fc809/Downloads/classical_monocytes_pathway_PCA_Hallmark"
dir.create(OUTDIR, showWarnings = FALSE)

# ── 1. Load original count matrix ────────────────────────────────────────────
cat("Reading count matrix...\n")
counts_df <- read.delim(
  "C:/Users/fc809/Downloads/classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt",
  stringsAsFactors = FALSE,
  row.names = 1
)
cat("Count matrix:", nrow(counts_df), "genes x", ncol(counts_df), "cells\n")

# ── 2. Assign lenti status and build metadata ─────────────────────────────────
cat("Checking for Lenti gene...\n")
if (!"Lenti" %in% rownames(counts_df)) stop("'Lenti' gene not found in rownames!")

cell_names   <- colnames(counts_df)
sample_ids   <- sub("_[A-Z]+\\..*$", "", cell_names)
exposure     <- ifelse(grepl("^HB", sample_ids), "HB", "LB")
lenti_status <- ifelse(counts_df["Lenti", ] >= 1, "Lenti+", "Lenti-")

cat("\nLenti status breakdown:\n")
print(table(lenti_status))
cat("\nLenti status per sample:\n")
print(table(sample_ids, lenti_status))

cell_meta <- data.frame(
  cell         = cell_names,
  sample       = sample_ids,
  exposure     = exposure,
  lenti_status = as.character(lenti_status),
  stringsAsFactors = FALSE
)
cell_meta$pb_group <- paste(cell_meta$sample, cell_meta$lenti_status, sep = "_")

cat("\nCells per pseudo-bulk group:\n")
print(sort(table(cell_meta$pb_group)))

# ── 2b. Filter genes: present in ≥5% of cells in every sample ───────────────
cat("\nFiltering genes by per-sample detection rate...\n")

samples_unique <- unique(cell_meta$sample)

detection <- do.call(cbind, lapply(samples_unique, function(s) {
  cells <- cell_meta$cell[cell_meta$sample == s]
  rowMeans(counts_df[, cells, drop = FALSE] > 0)
}))
colnames(detection) <- samples_unique

pass <- rowSums(detection >= 0.05) == length(samples_unique)
cat("Genes before filter:", nrow(counts_df), "\n")
cat("Genes passing >=20% detection in all samples:", sum(pass), "\n")

counts_filt <- as.matrix(counts_df[pass, ])
cat("Dimensions after filter:", nrow(counts_filt), "genes x",
    ncol(counts_filt), "cells\n")

# ── 3. Pseudo-bulk aggregation ────────────────────────────────────────────────
cat("\nAggregating counts per pseudo-bulk group...\n")
pb_groups <- unique(cell_meta$pb_group)

pb_mat <- do.call(cbind, lapply(pb_groups, function(g) {
  cells <- cell_meta$cell[cell_meta$pb_group == g]
  rowSums(counts_filt[, cells, drop = FALSE])
}))
colnames(pb_mat) <- pb_groups

cat("NAs in pb_mat:", sum(is.na(pb_mat)), "\n")

# Remove Lenti — vector not endogenous gene
pb_mat <- pb_mat[rownames(pb_mat) != "Lenti", ]

cat("Pseudo-bulk matrix:", nrow(pb_mat), "genes x", ncol(pb_mat), "samples\n")
print(sort(colnames(pb_mat)))

# ── 4. TMM-normalize → logCPM ─────────────────────────────────────────────────
y    <- DGEList(counts = pb_mat)
y    <- calcNormFactors(y, method = "TMM")
keep <- rowSums(cpm(y) > 1) >= 2
y    <- y[keep, , keep.lib.sizes = FALSE]
cat("Genes after TMM filtering:", nrow(y), "\n")
logcpm <- cpm(y, log = TRUE)

# ── 5. Gene sets: Hallmark only ───────────────────────────────────────────────
msig_h <- msigdbr(species = "Mus musculus", category = "H") |>
  dplyr::select(gs_name, gene_symbol)

gene_sets <- split(msig_h$gene_symbol, msig_h$gs_name)
cat("Total gene sets:", length(gene_sets), "\n")

# ── 6. GSVA ───────────────────────────────────────────────────────────────────
gsva_param <- gsvaParam(
  exprData = logcpm,
  geneSets = gene_sets,
  minSize  = 10,
  maxSize  = 500
)

gsva_scores <- gsva(gsva_param, verbose = TRUE, BPPARAM = BPPARAM)
cat("GSVA score matrix:", nrow(gsva_scores), "pathways x",
    ncol(gsva_scores), "samples\n")

# ── 7. Sample metadata for plotting ──────────────────────────────────────────
sample_meta <- cell_meta |>
  dplyr::distinct(pb_group, .keep_all = TRUE) |>
  dplyr::arrange(pb_group)

sample_meta <- sample_meta[match(colnames(gsva_scores),
                                 sample_meta$pb_group), ]
rownames(sample_meta) <- sample_meta$pb_group

# ── 8. PCA ────────────────────────────────────────────────────────────────────
pca_res <- prcomp(t(gsva_scores), scale. = TRUE, center = TRUE)
pct_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pca_df <- as.data.frame(pca_res$x[, 1:4]) |>
  tibble::rownames_to_column("pb_group") |>
  dplyr::left_join(sample_meta, by = "pb_group")

# ── 9. PCA plot ───────────────────────────────────────────────────────────────
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2,
                            color = exposure,
                            shape = lenti_status,
                            label = pb_group)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text(nudge_y = 0.4, size = 2.8, show.legend = FALSE) +
  scale_color_manual(values = c("HB" = "#B71C1C", "LB" = "#1A237E")) +
  scale_shape_manual(values = c("Lenti+" = 16, "Lenti-" = 1)) +
  labs(x     = paste0("PC1 (", pct_var[1], "%)"),
       y     = paste0("PC2 (", pct_var[2], "%)"),
       title = "Pathway PCA — LT-HSCs (Hallmark only, TMM only)",
       color = "Exposure", shape = "Lenti status") +
  theme_classic(base_size = 12)

print(p_pca)

# ── 10. Loading plots ─────────────────────────────────────────────────────────
loadings <- as.data.frame(pca_res$rotation[, 1:2]) |>
  tibble::rownames_to_column("pathway") |>
  dplyr::mutate(
    label = gsub("^HALLMARK_", "", pathway) |>
      tolower() |>
      gsub("_", " ", x = _)
  )

plot_loadings <- function(df, pc, n = 15, title) {
  df |>
    dplyr::arrange(desc(abs(.data[[pc]]))) |>
    head(n) |>
    ggplot(aes(x    = .data[[pc]],
               y    = reorder(label, .data[[pc]]),
               fill = .data[[pc]] > 0)) +
    geom_col() +
    scale_fill_manual(values = c("TRUE" = "#B71C1C", "FALSE" = "#1A237E"),
                      guide = "none") +
    labs(x = paste("Loading on", pc), y = NULL, title = title) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_text(size = 8))
}

p_load1 <- plot_loadings(loadings, "PC1", 15,
                         paste0("PC1 (", pct_var[1], "%) — top loadings"))
p_load2 <- plot_loadings(loadings, "PC2", 15,
                         paste0("PC2 (", pct_var[2], "%) — top loadings"))

print(p_load1)
print(p_load2)

# ── 11. Save ──────────────────────────────────────────────────────────────────
ggsave(file.path(OUTDIR, "PCA_samples.png"),  p_pca,   width = 7, height = 5, dpi = 150)
ggsave(file.path(OUTDIR, "Loadings_PC1.png"), p_load1, width = 8, height = 6, dpi = 150)
ggsave(file.path(OUTDIR, "Loadings_PC2.png"), p_load2, width = 8, height = 6, dpi = 150)

saveRDS(gsva_scores, file.path(OUTDIR, "gsva_scores.rds"))
saveRDS(pca_res,     file.path(OUTDIR, "pca_res.rds"))
write.csv(t(gsva_scores),
          file.path(OUTDIR, "gsva_scores_samplesXpathways.csv"),
          row.names = TRUE)
write.csv(pca_df,
          file.path(OUTDIR, "pca_coords.csv"),
          row.names = FALSE)

message("\nDone. Outputs in: ", OUTDIR)

##########PERMUTATION_TEST##############


library(ggplot2)
library(dplyr)

# 1. Define Paths
input_file <- "C:/Users/fc809/Downloads/LTHSC_pathway_PCA_noHKnorm/pca_coords.csv"
output_dir <- "C:/Users/fc809/Downloads/Permutation_Analysis_Results"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 2. Load Data
df <- read.csv(input_file)

# 3. Define Function to Calculate Centroid Distance
calc_dist <- function(data, labels) {
  centroids <- data %>%
    mutate(group = labels) %>%
    group_by(group) %>%
    summarise(mean_PC1 = mean(PC1), mean_PC2 = mean(PC2))
  
  # Euclidean Distance between HB and LB centroids
  hb <- centroids[centroids$group == "HB", c("mean_PC1", "mean_PC2")]
  lb <- centroids[centroids$group == "LB", c("mean_PC1", "mean_PC2")]
  
  dist_val <- sqrt((hb$mean_PC1 - lb$mean_PC1)^2 + (hb$mean_PC2 - lb$mean_PC2)^2)
  return(dist_val)
}

# 4. Calculate Observed Distance
observed_dist <- calc_dist(df, df$exposure)

# 5. Run Permutation Test
set.seed(42) # For reproducibility
n_permutations <- 10000
perm_distances <- numeric(n_permutations)

for (i in 1:n_permutations) {
  shuffled_labels <- sample(df$exposure)
  perm_distances[i] <- calc_dist(df, shuffled_labels)
}

# Calculate p-value
p_val <- sum(perm_distances >= observed_dist) / n_permutations

# 6. Visualization
# Plot 1: PCA Plot
p1 <- ggplot(df, aes(x = PC1, y = PC2, color = exposure, shape = exposure)) +
  geom_point(size = 4, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA: HB vs LB Separation",
       subtitle = paste("Observed Centroid Distance:", round(observed_dist, 3))) +
  scale_color_manual(values = c("HB" = "#E41A1C", "LB" = "#377EB8"))

# Plot 2: Null Distribution Histogram
perm_df <- data.frame(dist = perm_distances)
p2 <- ggplot(perm_df, aes(x = dist)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "white") +
  geom_vline(xintercept = observed_dist, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = observed_dist, y = n_permutations/20, 
           label = paste("Observed p =", p_val), color = "red", vjust = -1, angle = 90) +
  theme_minimal() +
  labs(title = "Permutation Test Distribution",
       x = "Centroid Distance", y = "Frequency")

# 7. Save Outputs
ggsave(file.path(output_dir, "PCA_separation.png"), plot = p1, width = 6, height = 5)
ggsave(file.path(output_dir, "Permutation_distribution.png"), plot = p2, width = 6, height = 5)

# Save summary stats to text file
results_text <- paste(
  "Permutation Test Results",
  "=======================",
  paste("Observed Distance:", observed_dist),
  paste("Number of Permutations:", n_permutations),
  paste("P-value:", p_val),
  sep = "\n"
)
writeLines(results_text, file.path(output_dir, "stats_results.txt"))

print(paste("Analysis complete. Results saved in:", output_dir))



######################################

if (!require("vegan")) install.packages("vegan")
library(vegan)

# 1. Load data
df <- read.csv("C:/Users/fc809/Downloads/LTHSC_pathway_PCA_noKHNorm/pca_coords.csv")

# 2. Define the permutation scheme (The "Paired" part)
# This ensures we shuffle 'exposure' at the level of 'sample'
# but do not shuffle lenti status within a sample.
perm_scheme <- how(blocks = df$sample, nperm = 9999)

# Note: For a between-subject factor like Exposure where pairs move together,
# it is often cleaner to aggregate by sample first:
df_agg <- df %>%
  group_by(sample, exposure) %>%
  summarise(across(starts_with("PC"), mean))

# 3. Run PERMANOVA on aggregated data (Recommended for HB vs LB)
permanova_paired <- adonis2(df_agg[, c("PC1", "PC2", "PC3", "PC4")] ~ exposure, 
                            data = df_agg, 
                            method = "euclidean", 
                            permutations = 9999)

# 4. Alternatively, use the 'strata' approach if you want to keep all 18 points:
# We use 'strata' to tell the model that rows within the same sample are related.
# However, adonis2 'strata' is mostly for within-block effects. 
# For your case, the aggregated method above is the statistically rigorous standard.

print("--- PERMANOVA Results (Aggregated per Sample) ---")
print(permanova_paired)

# Save result
capture.output(permanova_paired, file = "C:/Users/fc809/Downloads/permanova_paired_results.txt")

