# ══════════════════════════════════════════════════════════════════════════════
#  Heat Shock Module Score Analysis — Classical Monocytess
#  Computes per-cell heat shock score as mean z-score across signature genes
#  Expression is normalized to total UMIs per cell before z-scoring, so the
#  score reflects transcriptomic fraction rather than raw count depth
#  (equal gene contribution, sequencing-depth independent)
#  Plots: KDE ridges, violin, CDF, and HB4 threshold rescue analysis
# ══════════════════════════════════════════════════════════════════════════════
#
#  Required packages:
#  install.packages(c("tidyverse", "ggridges", "patchwork"))
#
# ══════════════════════════════════════════════════════════════════════════════

rm(list = ls())   # clear stale objects from any previous run

library(tidyverse)
library(ggridges)
library(patchwork)

setwd("C:/Users/fc809/Downloads")

base_dir  <- "Classical_Monocytes (1)"
hs_file   <- "Classical_Monocytes (1)/nUMIs_in_each_cells_alldata_res0.6_cluster0.csv"
tot_file  <- "Classical_Monocytes (1)/total_nUMIs_in_each_cell_alldata_res0.6_cluster0.csv"
out_dir   <- "Classical_Monocytes (1)/heatshock_analysis"
dir.create(out_dir, showWarnings = FALSE)

# ── Step 1: Load stress-gene UMI counts (genes × cells) ──────────────────────
raw <- read_csv(hs_file, show_col_types = FALSE) %>%
  column_to_rownames(colnames(.)[1])   # first col = gene names

# rows = genes, cols = cells → transpose to cells × genes
mat   <- t(raw)
genes <- colnames(mat)

# ── Step 2: Load total UMIs per cell and normalize ───────────────────────────
tot_raw <- read_csv(tot_file, show_col_types = FALSE)

# FIX: use quoted strings to avoid namespace conflict with cell() function
total_umis <- tot_raw %>%
  dplyr::rename(cell_id = "cell", total = "total_nUMI") %>%
  dplyr::select(cell_id, total)

cat("Total-UMI file loaded:", nrow(total_umis), "cells\n")
cat("Stress-gene matrix:   ", nrow(mat), "cells ×", ncol(mat), "genes\n")

# Keep only cells present in both files
common_cells <- intersect(rownames(mat), total_umis$cell_id)
cat("Cells in common:      ", length(common_cells), "\n")

if (length(common_cells) == 0)
  stop("No cell IDs matched between the two files. Check ID format.")

mat        <- mat[common_cells, , drop = FALSE]
total_umis <- total_umis %>% filter(cell_id %in% common_cells) %>%
  arrange(match(cell_id, common_cells))

# Normalize: UMIs per gene / total UMIs → fraction of transcriptome
# Multiply by 1e4 (similar to CP10k) so values are on a readable scale
total_vec  <- total_umis$total
mat_norm   <- sweep(mat, 1, total_vec, FUN = "/") * 1e4

cat("\nNormalization check — first 3 cells, row sums of fraction (should equal 1e4 if all genes summed):\n")
# (will be << 1e4 because mat_norm only covers the small stress-gene subset)
print(round(rowSums(mat_norm)[1:3], 4))

# ── Step 3: Build cell × gene data frame ─────────────────────────────────────
df <- as_tibble(mat_norm, rownames = "cell_id") %>%
  mutate(sample = str_split_fixed(cell_id, "_", 2)[, 1])

cat("\nGenes in signature:", paste(genes, collapse = ", "), "\n")
cat("Total cells after matching:", nrow(df), "\n")
cat("Cells per sample:\n")
print(count(df, sample) %>% arrange(sample))

# ── Step 4: Z-score each gene across all cells, then average ─────────────────
# Z-scoring on normalized fractions: equal gene contribution,
# independent of both sequencing depth (handled by normalization) and
# expression level (handled by z-scoring)
gene_mat <- df %>% dplyr::select(all_of(genes)) %>% as.matrix()
z_mat    <- scale(gene_mat, center = TRUE, scale = TRUE)   # z-score per gene
df$hs_score <- rowMeans(z_mat)

cat("\nPer-sample module score summary:\n")
df %>%
  group_by(sample) %>%
  summarise(
    n             = n(),
    mean          = round(mean(hs_score), 3),
    median        = round(median(hs_score), 3),
    sd            = round(sd(hs_score), 3),
    min           = round(min(hs_score), 3),
    max           = round(max(hs_score), 3),
    pct_above_1SD = round(mean(hs_score > 1) * 100, 1)
  ) %>%
  print()

# ── Colour scheme ─────────────────────────────────────────────────────────────
sample_order <- c("HB1","HB2","HB3","HB4","HB5","LB1","LB2","LB3","LB4")

sample_colors <- c(
  HB1 = "#E74C3C", HB2 = "#E74C3C",
  HB3 = "#FF6B00", HB4 = "#FF6B00",   # orange = outliers
  HB5 = "#E74C3C",
  LB1 = "#3498DB", LB2 = "#3498DB",
  LB3 = "#3498DB", LB4 = "#3498DB"
)

df$sample <- factor(df$sample, levels = sample_order)

# ── Plot 1: KDE ridge plot ────────────────────────────────────────────────────
p_ridge <- ggplot(df, aes(x = hs_score, y = sample, fill = sample, color = sample)) +
  geom_density_ridges(
    alpha = 0.55, scale = 0.9,
    quantile_lines = TRUE, quantiles = 2,
    bandwidth = 0.06
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  coord_cartesian(xlim = c(NA, 6)) +
  scale_fill_manual(values  = sample_colors, guide = "none") +
  scale_color_manual(values = sample_colors, guide = "none") +
  scale_y_discrete(limits = rev(sample_order)) +
  scale_x_continuous(breaks = function(lims) sort(unique(c(scales::breaks_pretty()(lims), 1)))) +
  annotate("text", x = 1, y = 0.6, label = "1", size = 3.5, color = "grey40",
           hjust = -0.3, fontface = "bold") +
  labs(
    title    = "Classical Monocytes Heat Shock Module Score — Distribution per Sample",
    subtitle = "Mean z-score across signature genes | CP10k-normalized",
    x = "Module score (mean z-score of CP10k-normalized counts)",
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.text     = element_text(size = 12, color = "black"),
    axis.title    = element_text(size = 13, color = "black"),
    axis.text.y   = element_text(color = rev(sample_colors[sample_order]))
  )

ggsave(file.path(out_dir, "heatshock_ridge.pdf"), p_ridge, width = 9, height = 9)
ggsave(file.path(out_dir, "heatshock_ridge.png"), p_ridge, width = 9, height = 9, dpi = 160)

# ── Plot 2: Violin + jitter ───────────────────────────────────────────────────
annot <- df %>%
  group_by(sample) %>%
  summarise(
    median   = median(hs_score),
    pct_high = round(mean(hs_score > 1) * 100, 1)
  )

p_violin <- ggplot(df, aes(x = sample, y = hs_score, fill = sample, color = sample)) +
  geom_violin(alpha = 0.55, scale = "width", trim = FALSE) +
  geom_jitter(width = 0.15, alpha = 0.12, size = 0.6) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.4, color = "black", linewidth = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_text(data = annot,
            aes(y = max(df$hs_score) * 0.92,
                label = paste0(pct_high, "%")),
            size = 3.5, fontface = "bold", color = "black") +
  scale_fill_manual(values  = sample_colors, guide = "none") +
  scale_color_manual(values = sample_colors, guide = "none") +
  labs(
    title    = "Classical Monocytes Heat Shock Module Score per Cell",
    subtitle = "Crossbar = median  |  % = cells with score > 1 SD  |  Dashed = z = 1",
    x = NULL,
    y = "Module score (mean z-score of CP10k-normalized counts)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.text     = element_text(size = 12, color = "black"),
    axis.title    = element_text(size = 13, color = "black"),
    axis.text.x   = element_text(color = sample_colors[sample_order])
  )

ggsave(file.path(out_dir, "heatshock_violin.pdf"), p_violin, width = 11, height = 6)
ggsave(file.path(out_dir, "heatshock_violin.png"), p_violin, width = 11, height = 6, dpi = 160)

# ── Plot 3: CDF per sample ────────────────────────────────────────────────────
cdf_df <- df %>%
  group_by(sample) %>%
  arrange(hs_score, .by_group = TRUE) %>%
  mutate(cdf = seq_along(hs_score) / n())

p_cdf_all <- ggplot(cdf_df, aes(x = hs_score, y = cdf, color = sample,
                                linewidth = sample %in% c("HB3","HB4"),
                                linetype  = str_detect(sample, "LB"))) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50", linewidth = 0.5) +
  scale_color_manual(values = sample_colors, name = "Sample") +
  scale_linewidth_manual(values = c(`TRUE` = 1.8, `FALSE` = 1.1), guide = "none") +
  scale_linetype_manual(values = c(`TRUE` = "dashed", `FALSE` = "solid"), guide = "none") +
  labs(
    title = "Classical Monocytes Heat Shock Module Score — Cumulative Distribution",
    x = "Module score", y = "Cumulative fraction of cells"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title  = element_text(face = "bold", size = 14, color = "black"),
    axis.text   = element_text(size = 12, color = "black"),
    axis.title  = element_text(size = 13, color = "black"),
    legend.text = element_text(size = 11)
  )

p_cdf_zoom <- p_cdf_all +
  coord_cartesian(xlim = c(0.2, quantile(df$hs_score, 0.999)),
                  ylim = c(0.4, 1.0)) +
  labs(title = "Zoomed: upper half")

p_cdf <- p_cdf_all + p_cdf_zoom +
  plot_annotation(
    title = "Classical Monocytes Heat Shock Module Score — Cumulative Distribution",
    theme = theme(plot.title = element_text(face = "bold", size = 15, color = "black"))
  ) &
  theme(plot.title = element_text(face = "bold", size = 13))

ggsave(file.path(out_dir, "heatshock_cdf.pdf"), p_cdf, width = 14, height = 6)
ggsave(file.path(out_dir, "heatshock_cdf.png"), p_cdf, width = 14, height = 6, dpi = 160)

# ── Plot 4: HB4 threshold rescue analysis ─────────────────────────────────────
hb4        <- df %>% filter(sample == "HB4")
thresholds <- seq(-1, 2, by = 0.1)

rescue_df <- tibble(
  threshold   = thresholds,
  n_rescued   = map_int(thresholds, ~ sum(hb4$hs_score < .x)),
  pct_rescued = map_dbl(thresholds, ~ mean(hb4$hs_score < .x) * 100)
)

p_rescue <- ggplot(rescue_df, aes(x = threshold, y = pct_rescued)) +
  geom_line(color = "#FF6B00", linewidth = 1.5) +
  geom_point(color = "#FF6B00", size = 2.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 50, linetype = "dotted", color = "grey60") +
  annotate("text", x = 0.05, y = 2,
           label = "z = 0\n(global mean)", hjust = 0, size = 3.5, color = "grey40") +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Classical Monocytes HB4: Fraction of Cells Rescued by Score Threshold",
    subtitle = "Cells with score < threshold considered 'not heat shocked'",
    x = "Score threshold (z-score)",
    y = "% of HB4 cells rescued"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.text     = element_text(size = 12, color = "black"),
    axis.title    = element_text(size = 13, color = "black")
  )

pct_below <- round(mean(hb4$hs_score < 0) * 100, 1)
dens_hb4  <- density(hb4$hs_score, bw = 0.06)
ymax_hb4  <- max(dens_hb4$y) * 0.85

p_hb4_density <- ggplot(hb4, aes(x = hs_score)) +
  geom_density(fill = "#FF6B00", alpha = 0.4, color = "#FF6B00",
               linewidth = 1.5, bw = 0.06) +
  annotate("rect", xmin = -Inf, xmax = 0,
           ymin = 0, ymax = Inf, fill = "green4", alpha = 0.08) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  annotate("text", x = -0.2, y = ymax_hb4,
           label = paste0(pct_below, "% of cells\nbelow z = 0"),
           hjust = 0.5, vjust = 1, size = 4.5, color = "green4", fontface = "bold") +
  labs(
    title    = "Classical Monocytes HB4: Heat Shock Score Distribution",
    subtitle = "Green shading = cells below z = 0 (potential rescue fraction)",
    x = "Module score (mean z-score)", y = "Density"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.text     = element_text(size = 12, color = "black"),
    axis.title    = element_text(size = 13, color = "black")
  )

p_hb4 <- p_hb4_density + p_rescue +
  plot_annotation(
    title = "Classical Monocytes HB4: Cell-Level Rescue Analysis",
    theme = theme(plot.title = element_text(face = "bold", size = 15, color = "black"))
  )

ggsave(file.path(out_dir, "heatshock_HB4_rescue.pdf"), p_hb4, width = 14, height = 6)
ggsave(file.path(out_dir, "heatshock_HB4_rescue.png"), p_hb4, width = 14, height = 6, dpi = 160)

# ── Plot 5: Rescue curve for ALL samples ──────────────────────────────────────
rescue_all <- map_dfr(sample_order, function(s) {
  sdata <- df %>% filter(sample == s) %>% pull(hs_score)
  tibble(
    sample      = s,
    threshold   = thresholds,
    pct_rescued = map_dbl(thresholds, ~ mean(sdata < .x) * 100)
  )
}) %>%
  mutate(sample = factor(sample, levels = sample_order))

rescue_colors <- c(
  HB1 = "#FF4500", HB2 = "#FF8C00", HB3 = "#FFC300", HB4 = "#B8860B", HB5 = "#CD3700",
  LB1 = "#1A6EBF", LB2 = "#4AABDB", LB3 = "#0A3D6B", LB4 = "#72C8E8"
)

rescue_linetypes <- setNames(rep("solid", 9), sample_order)
rescue_widths    <- setNames(rep(1.4,    9), sample_order)

p_rescue_all <- ggplot(rescue_all,
                       aes(x = threshold, y = pct_rescued,
                           color     = sample,
                           linetype  = sample,
                           linewidth = sample)) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_hline(yintercept = 50, linetype = "dotted", color = "grey60", linewidth = 0.5) +
  annotate("text", x = 0.05, y = 2, label = "z = 0", hjust = 0, size = 3.5, color = "grey40") +
  scale_color_manual(values = rescue_colors, name = "Sample") +
  scale_linetype_manual(values = rescue_linetypes, name = "Sample") +
  scale_linewidth_manual(values = rescue_widths, name = "Sample") +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Classical Monocytes Heat Shock — Fraction of Cells Below Score Threshold",
    subtitle = "Orange = HB samples  |  Blue = LB samples  |  CP10k-normalized",
    x = "Score threshold (z-score)",
    y = "% of cells below threshold"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.text     = element_text(size = 12, color = "black"),
    axis.title    = element_text(size = 13, color = "black"),
    legend.text   = element_text(size = 11)
  )

ggsave(file.path(out_dir, "heatshock_rescue_allsamples.pdf"), p_rescue_all, width = 9, height = 6)
ggsave(file.path(out_dir, "heatshock_rescue_allsamples.png"), p_rescue_all, width = 9, height = 6, dpi = 160)

# ── Save per-cell scores ──────────────────────────────────────────────────────
df %>%
  dplyr::select(cell_id, sample, hs_score) %>%
  write_csv(file.path(out_dir, "heatshock_scores_per_cell.csv"))

cat("\n✓ All outputs saved to:", out_dir, "\n")

