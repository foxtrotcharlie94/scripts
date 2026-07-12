# ══════════════════════════════════════════════════════════════════════════════
#  Heat Shock Transcriptomic Fraction Analysis — Classical Monocytes
#  Metric: sum of HS gene UMIs / total UMIs per cell × 100
#  = % of each cell's transcriptome devoted to heat shock genes
#  No z-scoring, no normalization assumptions — fully transparent
#  Plots: KDE ridges, violin, CDF, and HB4 threshold rescue analysis
# ══════════════════════════════════════════════════════════════════════════════
#
#  Required packages:
#  install.packages(c("tidyverse", "ggridges", "patchwork"))
#
# ══════════════════════════════════════════════════════════════════════════════

rm(list = ls())

library(tidyverse)
library(ggridges)
library(patchwork)

setwd("C:/Users/fc809/Downloads")

base_dir  <- "Classical_Monocytes (1)"
hs_file   <- "Classical_Monocytes (1)/nUMIs_in_each_cells_alldata_res0.6_cluster0.csv"
tot_file  <- "Classical_Monocytes (1)/total_nUMIs_in_each_cell_alldata_res0.6_cluster0.csv"
out_dir   <- "Classical_Monocytes (1)/heatshock_fraction_analysis"
dir.create(out_dir, showWarnings = FALSE)

# ── Step 1: Load stress-gene UMI counts (genes × cells) ──────────────────────
raw <- read_csv(hs_file, show_col_types = FALSE) %>%
  column_to_rownames(colnames(.)[1])

mat   <- t(raw)   # cells × genes
genes <- colnames(mat)

# ── Step 2: Load total UMIs per cell ─────────────────────────────────────────
tot_raw <- read_csv(tot_file, show_col_types = FALSE)

total_umis <- tot_raw %>%
  dplyr::rename(cell_id = "cell", total = "total_nUMI") %>%
  dplyr::select(cell_id, total)

cat("Total-UMI file loaded:", nrow(total_umis), "cells\n")
cat("Stress-gene matrix:   ", nrow(mat), "cells ×", ncol(mat), "genes\n")

common_cells <- intersect(rownames(mat), total_umis$cell_id)
cat("Cells in common:      ", length(common_cells), "\n")

if (length(common_cells) == 0)
  stop("No cell IDs matched between the two files. Check ID format.")

mat        <- mat[common_cells, , drop = FALSE]
total_umis <- total_umis %>% filter(cell_id %in% common_cells) %>%
  arrange(match(cell_id, common_cells))

# ── Step 3: Compute HS fraction ───────────────────────────────────────────────
# hs_fraction = sum of HS gene UMIs / total UMIs per cell × 100
# Expressed as % of transcriptome — no z-scoring, no CP10k
total_vec   <- total_umis$total
hs_umis     <- rowSums(mat)                        # total HS UMIs per cell
hs_fraction <- (hs_umis / total_vec) * 100         # % of transcriptome

df <- tibble(
  cell_id     = common_cells,
  hs_fraction = hs_fraction,
  hs_umis     = hs_umis,
  total_umis  = total_vec,
  sample      = str_split_fixed(common_cells, "_", 2)[, 1]
)

cat("\nGenes in signature:", paste(genes, collapse = ", "), "\n")
cat("Total cells:", nrow(df), "\n")
cat("Cells per sample:\n")
print(count(df, sample) %>% arrange(sample))

cat("\nPer-sample HS fraction summary (% of transcriptome):\n")
df %>%
  group_by(sample) %>%
  summarise(
    n              = n(),
    mean           = round(mean(hs_fraction), 4),
    median         = round(median(hs_fraction), 4),
    sd             = round(sd(hs_fraction), 4),
    min            = round(min(hs_fraction), 4),
    max            = round(max(hs_fraction), 4),
    pct_above_1pct = round(mean(hs_fraction > 1) * 100, 1)
  ) %>%
  print()

# ── Colour scheme ─────────────────────────────────────────────────────────────
sample_order <- c("HB1","HB2","HB3","HB4","HB5","LB1","LB2","LB3","LB4")

sample_colors <- c(
  HB1 = "#E74C3C", HB2 = "#E74C3C",
  HB3 = "#FF6B00", HB4 = "#FF6B00",
  HB5 = "#E74C3C",
  LB1 = "#3498DB", LB2 = "#3498DB",
  LB3 = "#3498DB", LB4 = "#3498DB"
)

df$sample <- factor(df$sample, levels = sample_order)

# Reference line: max 90th percentile across samples excluding HB3 and HB4
# = most conservative "normal" threshold: 90% of cells in each clean sample
# fall below this value
p90_per_sample <- df %>%
  filter(!sample %in% c("HB3", "HB4")) %>%
  group_by(sample) %>%
  summarise(p90 = quantile(hs_fraction, 0.90))

threshold_line <- max(p90_per_sample$p90)
cat("\n90th percentile per sample (excl. HB3/HB4):\n")
print(p90_per_sample)
cat("Threshold line (max p90):", round(threshold_line, 4), "%\n")

# ── Plot 1: KDE ridge plot ────────────────────────────────────────────────────
p_ridge <- ggplot(df, aes(x = hs_fraction, y = sample, fill = sample, color = sample)) +
  geom_density_ridges(
    alpha = 0.55, scale = 0.9,
    quantile_lines = TRUE, quantiles = 2
  ) +
  geom_vline(xintercept = threshold_line, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  annotate("text", x = threshold_line, y = 0.6,
           label = paste0("p90 threshold\n", round(threshold_line, 3), "%"),
           size = 3.2, color = "grey40", hjust = -0.15, fontface = "bold") +
  scale_fill_manual(values  = sample_colors, guide = "none") +
  scale_color_manual(values = sample_colors, guide = "none") +
  scale_y_discrete(limits = rev(sample_order)) +
  labs(
    title    = "Classical Monocytes — Heat Shock Transcriptomic Fraction",
    subtitle = "% of total UMIs per cell mapping to HS signature genes  |  Dashed = max p90 across samples excl. HB3/HB4",
    x = "HS gene UMIs / total UMIs × 100 (%)",
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

ggsave(file.path(out_dir, "heatshock_fraction_ridge.pdf"), p_ridge, width = 9, height = 9)
ggsave(file.path(out_dir, "heatshock_fraction_ridge.png"), p_ridge, width = 9, height = 9, dpi = 160)

# ── Plot 2: Violin + jitter ───────────────────────────────────────────────────
annot <- df %>%
  group_by(sample) %>%
  summarise(
    median   = median(hs_fraction),
    pct_high = round(mean(hs_fraction > threshold_line) * 100, 1)
  )

p_violin <- ggplot(df, aes(x = sample, y = hs_fraction, fill = sample, color = sample)) +
  geom_violin(alpha = 0.55, scale = "width", trim = FALSE) +
  geom_jitter(width = 0.15, alpha = 0.12, size = 0.6) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.4, color = "black", linewidth = 0.6) +
  geom_hline(yintercept = threshold_line, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  geom_text(data = annot,
            aes(y = max(df$hs_fraction) * 0.92,
                label = paste0(pct_high, "%")),
            size = 3.5, fontface = "bold", color = "black") +
  scale_fill_manual(values  = sample_colors, guide = "none") +
  scale_color_manual(values = sample_colors, guide = "none") +
  labs(
    title    = "Classical Monocytes — Heat Shock Transcriptomic Fraction per Cell",
    subtitle = paste0("Crossbar = median  |  % = cells above p90 threshold (", round(threshold_line, 3), "%)  |  Dashed = p90 threshold"),
    x = NULL,
    y = "HS gene UMIs / total UMIs × 100 (%)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.text     = element_text(size = 12, color = "black"),
    axis.title    = element_text(size = 13, color = "black"),
    axis.text.x   = element_text(color = sample_colors[sample_order])
  )

ggsave(file.path(out_dir, "heatshock_fraction_violin.pdf"), p_violin, width = 11, height = 6)
ggsave(file.path(out_dir, "heatshock_fraction_violin.png"), p_violin, width = 11, height = 6, dpi = 160)

# ── Plot 3: CDF per sample ────────────────────────────────────────────────────
cdf_df <- df %>%
  group_by(sample) %>%
  arrange(hs_fraction, .by_group = TRUE) %>%
  mutate(cdf = seq_along(hs_fraction) / n())

p_cdf_all <- ggplot(cdf_df, aes(x = hs_fraction, y = cdf, color = sample,
                                linewidth = sample %in% c("HB3","HB4"),
                                linetype  = str_detect(sample, "LB"))) +
  geom_line() +
  geom_vline(xintercept = threshold_line, linetype = "dotted",
             color = "grey50", linewidth = 0.5) +
  scale_color_manual(values = sample_colors, name = "Sample") +
  scale_linewidth_manual(values = c(`TRUE` = 1.8, `FALSE` = 1.1), guide = "none") +
  scale_linetype_manual(values = c(`TRUE` = "dashed", `FALSE` = "solid"), guide = "none") +
  labs(
    title = "Classical Monocytes — Heat Shock Fraction: Cumulative Distribution",
    x = "HS gene UMIs / total UMIs × 100 (%)",
    y = "Cumulative fraction of cells"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title  = element_text(face = "bold", size = 14, color = "black"),
    axis.text   = element_text(size = 12, color = "black"),
    axis.title  = element_text(size = 13, color = "black"),
    legend.text = element_text(size = 11)
  )

p_cdf_zoom <- p_cdf_all +
  coord_cartesian(xlim = c(threshold_line, quantile(df$hs_fraction, 0.999)),
                  ylim = c(0.4, 1.0)) +
  labs(title = "Zoomed: upper half")

p_cdf <- p_cdf_all + p_cdf_zoom +
  plot_annotation(
    title = "Classical Monocytes — Heat Shock Fraction: Cumulative Distribution",
    theme = theme(plot.title = element_text(face = "bold", size = 15, color = "black"))
  ) &
  theme(plot.title = element_text(face = "bold", size = 13))

ggsave(file.path(out_dir, "heatshock_fraction_cdf.pdf"), p_cdf, width = 14, height = 6)
ggsave(file.path(out_dir, "heatshock_fraction_cdf.png"), p_cdf, width = 14, height = 6, dpi = 160)

# ── Plot 4: HB4 threshold rescue analysis ─────────────────────────────────────
hb4        <- df %>% filter(sample == "HB4")
frac_max   <- quantile(df$hs_fraction, 0.99)
thresholds <- seq(0, frac_max, length.out = 50)

rescue_df <- tibble(
  threshold   = thresholds,
  n_rescued   = map_int(thresholds, ~ sum(hb4$hs_fraction < .x)),
  pct_rescued = map_dbl(thresholds, ~ mean(hb4$hs_fraction < .x) * 100)
)

p_rescue <- ggplot(rescue_df, aes(x = threshold, y = pct_rescued)) +
  geom_line(color = "#FF6B00", linewidth = 1.5) +
  geom_point(color = "#FF6B00", size = 2.5) +
  geom_vline(xintercept = threshold_line, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 50, linetype = "dotted", color = "grey60") +
  annotate("text", x = threshold_line * 1.05, y = 2,
           label = paste0("p90 threshold\n", round(threshold_line, 3), "%"),
           hjust = 0, size = 3.5, color = "grey40") +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Classical Monocytes HB4: Fraction of Cells Rescued by Threshold",
    subtitle = "Cells with HS fraction < threshold considered 'not heat shocked'",
    x = "HS fraction threshold (%)",
    y = "% of HB4 cells rescued"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.text     = element_text(size = 12, color = "black"),
    axis.title    = element_text(size = 13, color = "black")
  )

pct_below <- round(mean(hb4$hs_fraction < threshold_line) * 100, 1)
dens_hb4  <- density(hb4$hs_fraction)
ymax_hb4  <- max(dens_hb4$y) * 0.85

p_hb4_density <- ggplot(hb4, aes(x = hs_fraction)) +
  geom_density(fill = "#FF6B00", alpha = 0.4, color = "#FF6B00", linewidth = 1.5) +
  annotate("rect", xmin = -Inf, xmax = threshold_line,
           ymin = 0, ymax = Inf, fill = "green4", alpha = 0.08) +
  geom_vline(xintercept = threshold_line, linetype = "dashed",
             color = "grey30", linewidth = 0.8) +
  annotate("text", x = threshold_line * 0.6, y = ymax_hb4,
           label = paste0(pct_below, "% of cells\nbelow p90 threshold"),
           hjust = 0.5, vjust = 1, size = 4.5, color = "green4", fontface = "bold") +
  labs(
    title    = "Classical Monocytes HB4: Heat Shock Fraction Distribution",
    subtitle = "Green shading = cells below p90 threshold (potential rescue fraction)",
    x = "HS gene UMIs / total UMIs × 100 (%)", y = "Density"
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

ggsave(file.path(out_dir, "heatshock_fraction_HB4_rescue.pdf"), p_hb4, width = 14, height = 6)
ggsave(file.path(out_dir, "heatshock_fraction_HB4_rescue.png"), p_hb4, width = 14, height = 6, dpi = 160)

# ── Plot 5: Rescue curve for ALL samples ──────────────────────────────────────
rescue_all <- map_dfr(sample_order, function(s) {
  sdata <- df %>% filter(sample == s) %>% pull(hs_fraction)
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
  geom_vline(xintercept = threshold_line, linetype = "dashed",
             color = "grey40", linewidth = 0.6) +
  geom_hline(yintercept = 50, linetype = "dotted", color = "grey60", linewidth = 0.5) +
  annotate("text", x = threshold_line * 1.05, y = 2,
           label = "p90\nthreshold", hjust = 0, size = 3.5, color = "grey40") +
  scale_color_manual(values = rescue_colors, name = "Sample") +
  scale_linetype_manual(values = rescue_linetypes, name = "Sample") +
  scale_linewidth_manual(values = rescue_widths, name = "Sample") +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Classical Monocytes — Fraction of Cells Below HS Fraction Threshold",
    subtitle = "Orange = HB samples  |  Blue = LB samples  |  Dashed = p90 threshold (excl. HB3/HB4)",
    x = "HS fraction threshold (%)",
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

ggsave(file.path(out_dir, "heatshock_fraction_rescue_allsamples.pdf"), p_rescue_all, width = 9, height = 6)
ggsave(file.path(out_dir, "heatshock_fraction_rescue_allsamples.png"), p_rescue_all, width = 9, height = 6, dpi = 160)

# ── Save per-cell fractions ───────────────────────────────────────────────────
df %>%
  dplyr::select(cell_id, sample, hs_fraction, hs_umis, total_umis) %>%
  write_csv(file.path(out_dir, "heatshock_fraction_per_cell.csv"))

cat("\n✓ All outputs saved to:", out_dir, "\n")

# ══════════════════════════════════════════════════════════════════════════════
#  Heat Shock Transcriptomic Fraction Analysis — LT-HSCs
#  Metric: sum of HS gene UMIs / total UMIs per cell × 100
#  = % of each cell's transcriptome devoted to heat shock genes
#  No z-scoring, no normalization assumptions — fully transparent
#  Plots: KDE ridges, violin, CDF, and HB4 threshold rescue analysis
# ══════════════════════════════════════════════════════════════════════════════
#
#  Required packages:
#  install.packages(c("tidyverse", "ggridges", "patchwork"))
#
# ══════════════════════════════════════════════════════════════════════════════

rm(list = ls())

library(tidyverse)
library(ggridges)
library(patchwork)

setwd("C:/Users/fc809/Downloads")

base_dir  <- "LT-HSCs (1)"
hs_file   <- "LT-HSCs (1)/nUMIs_in_each_cells_HSPCsubset_res0.3_cluster0.csv"
tot_file  <- "LT-HSCs (1)/total_nUMIs_in_each_cell_HSPCsubset_res0.3_cluster0.csv"
out_dir   <- "LT-HSCs (1)/heatshock_fraction_analysis"
dir.create(out_dir, showWarnings = FALSE)

# ── Step 1: Load stress-gene UMI counts (genes × cells) ──────────────────────
raw <- read_csv(hs_file, show_col_types = FALSE) %>%
  column_to_rownames(colnames(.)[1])

mat   <- t(raw)   # cells × genes
genes <- colnames(mat)

# ── Step 2: Load total UMIs per cell ─────────────────────────────────────────
tot_raw <- read_csv(tot_file, show_col_types = FALSE)

total_umis <- tot_raw %>%
  dplyr::rename(cell_id = "cell", total = "total_nUMI") %>%
  dplyr::select(cell_id, total)

cat("Total-UMI file loaded:", nrow(total_umis), "cells\n")
cat("Stress-gene matrix:   ", nrow(mat), "cells ×", ncol(mat), "genes\n")

common_cells <- intersect(rownames(mat), total_umis$cell_id)
cat("Cells in common:      ", length(common_cells), "\n")

if (length(common_cells) == 0)
  stop("No cell IDs matched between the two files. Check ID format.")

mat        <- mat[common_cells, , drop = FALSE]
total_umis <- total_umis %>% filter(cell_id %in% common_cells) %>%
  arrange(match(cell_id, common_cells))

# ── Step 3: Compute HS fraction ───────────────────────────────────────────────
# hs_fraction = sum of HS gene UMIs / total UMIs per cell × 100
# Expressed as % of transcriptome — no z-scoring, no CP10k
total_vec   <- total_umis$total
hs_umis     <- rowSums(mat)                        # total HS UMIs per cell
hs_fraction <- (hs_umis / total_vec) * 100         # % of transcriptome

df <- tibble(
  cell_id     = common_cells,
  hs_fraction = hs_fraction,
  hs_umis     = hs_umis,
  total_umis  = total_vec,
  sample      = str_split_fixed(common_cells, "_", 2)[, 1]
)

cat("\nGenes in signature:", paste(genes, collapse = ", "), "\n")
cat("Total cells:", nrow(df), "\n")
cat("Cells per sample:\n")
print(count(df, sample) %>% arrange(sample))

cat("\nPer-sample HS fraction summary (% of transcriptome):\n")
df %>%
  group_by(sample) %>%
  summarise(
    n              = n(),
    mean           = round(mean(hs_fraction), 4),
    median         = round(median(hs_fraction), 4),
    sd             = round(sd(hs_fraction), 4),
    min            = round(min(hs_fraction), 4),
    max            = round(max(hs_fraction), 4),
    pct_above_1pct = round(mean(hs_fraction > 1) * 100, 1)
  ) %>%
  print()

# ── Colour scheme ─────────────────────────────────────────────────────────────
sample_order <- c("HB1","HB2","HB3","HB4","HB5","LB1","LB2","LB3","LB4")

sample_colors <- c(
  HB1 = "#E74C3C", HB2 = "#E74C3C",
  HB3 = "#FF6B00", HB4 = "#FF6B00",
  HB5 = "#E74C3C",
  LB1 = "#3498DB", LB2 = "#3498DB",
  LB3 = "#3498DB", LB4 = "#3498DB"
)

df$sample <- factor(df$sample, levels = sample_order)

# Reference line: max 90th percentile across samples excluding HB3 and HB4
# = most conservative "normal" threshold: 90% of cells in each clean sample
# fall below this value
p90_per_sample <- df %>%
  filter(!sample %in% c("HB3", "HB4")) %>%
  group_by(sample) %>%
  summarise(p90 = quantile(hs_fraction, 0.90))

threshold_line <- max(p90_per_sample$p90)
cat("\n90th percentile per sample (excl. HB3/HB4):\n")
print(p90_per_sample)
cat("Threshold line (max p90):", round(threshold_line, 4), "%\n")

# ── Plot 1: KDE ridge plot ────────────────────────────────────────────────────
p_ridge <- ggplot(df, aes(x = hs_fraction, y = sample, fill = sample, color = sample)) +
  geom_density_ridges(
    alpha = 0.55, scale = 0.9,
    quantile_lines = TRUE, quantiles = 2
  ) +
  geom_vline(xintercept = threshold_line, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  annotate("text", x = threshold_line, y = 0.6,
           label = paste0("p90 threshold\n", round(threshold_line, 3), "%"),
           size = 3.2, color = "grey40", hjust = -0.15, fontface = "bold") +
  scale_fill_manual(values  = sample_colors, guide = "none") +
  scale_color_manual(values = sample_colors, guide = "none") +
  scale_y_discrete(limits = rev(sample_order)) +
  labs(
    title    = "LT-HSCs — Heat Shock Transcriptomic Fraction",
    subtitle = "% of total UMIs per cell mapping to HS signature genes  |  Dashed = max p90 across samples excl. HB3/HB4",
    x = "HS gene UMIs / total UMIs × 100 (%)",
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

ggsave(file.path(out_dir, "heatshock_fraction_ridge.pdf"), p_ridge, width = 9, height = 9)
ggsave(file.path(out_dir, "heatshock_fraction_ridge.png"), p_ridge, width = 9, height = 9, dpi = 160)

# ── Plot 2: Violin + jitter ───────────────────────────────────────────────────
annot <- df %>%
  group_by(sample) %>%
  summarise(
    median   = median(hs_fraction),
    pct_high = round(mean(hs_fraction > threshold_line) * 100, 1)
  )

p_violin <- ggplot(df, aes(x = sample, y = hs_fraction, fill = sample, color = sample)) +
  geom_violin(alpha = 0.55, scale = "width", trim = FALSE) +
  geom_jitter(width = 0.15, alpha = 0.12, size = 0.6) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.4, color = "black", linewidth = 0.6) +
  geom_hline(yintercept = threshold_line, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  geom_text(data = annot,
            aes(y = max(df$hs_fraction) * 0.92,
                label = paste0(pct_high, "%")),
            size = 3.5, fontface = "bold", color = "black") +
  scale_fill_manual(values  = sample_colors, guide = "none") +
  scale_color_manual(values = sample_colors, guide = "none") +
  labs(
    title    = "LT-HSCs — Heat Shock Transcriptomic Fraction per Cell",
    subtitle = paste0("Crossbar = median  |  % = cells above p90 threshold (", round(threshold_line, 3), "%)  |  Dashed = p90 threshold"),
    x = NULL,
    y = "HS gene UMIs / total UMIs × 100 (%)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.text     = element_text(size = 12, color = "black"),
    axis.title    = element_text(size = 13, color = "black"),
    axis.text.x   = element_text(color = sample_colors[sample_order])
  )

ggsave(file.path(out_dir, "heatshock_fraction_violin.pdf"), p_violin, width = 11, height = 6)
ggsave(file.path(out_dir, "heatshock_fraction_violin.png"), p_violin, width = 11, height = 6, dpi = 160)

# ── Plot 3: CDF per sample ────────────────────────────────────────────────────
cdf_df <- df %>%
  group_by(sample) %>%
  arrange(hs_fraction, .by_group = TRUE) %>%
  mutate(cdf = seq_along(hs_fraction) / n())

p_cdf_all <- ggplot(cdf_df, aes(x = hs_fraction, y = cdf, color = sample,
                                linewidth = sample %in% c("HB3","HB4"),
                                linetype  = str_detect(sample, "LB"))) +
  geom_line() +
  geom_vline(xintercept = threshold_line, linetype = "dotted",
             color = "grey50", linewidth = 0.5) +
  scale_color_manual(values = sample_colors, name = "Sample") +
  scale_linewidth_manual(values = c(`TRUE` = 1.8, `FALSE` = 1.1), guide = "none") +
  scale_linetype_manual(values = c(`TRUE` = "dashed", `FALSE` = "solid"), guide = "none") +
  labs(
    title = "LT-HSCs — Heat Shock Fraction: Cumulative Distribution",
    x = "HS gene UMIs / total UMIs × 100 (%)",
    y = "Cumulative fraction of cells"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title  = element_text(face = "bold", size = 14, color = "black"),
    axis.text   = element_text(size = 12, color = "black"),
    axis.title  = element_text(size = 13, color = "black"),
    legend.text = element_text(size = 11)
  )

p_cdf_zoom <- p_cdf_all +
  coord_cartesian(xlim = c(threshold_line, quantile(df$hs_fraction, 0.999)),
                  ylim = c(0.4, 1.0)) +
  labs(title = "Zoomed: upper half")

p_cdf <- p_cdf_all + p_cdf_zoom +
  plot_annotation(
    title = "LT-HSCs — Heat Shock Fraction: Cumulative Distribution",
    theme = theme(plot.title = element_text(face = "bold", size = 15, color = "black"))
  ) &
  theme(plot.title = element_text(face = "bold", size = 13))

ggsave(file.path(out_dir, "heatshock_fraction_cdf.pdf"), p_cdf, width = 14, height = 6)
ggsave(file.path(out_dir, "heatshock_fraction_cdf.png"), p_cdf, width = 14, height = 6, dpi = 160)

# ── Plot 4: HB4 threshold rescue analysis ─────────────────────────────────────
hb4        <- df %>% filter(sample == "HB4")
frac_max   <- quantile(df$hs_fraction, 0.99)
thresholds <- seq(0, frac_max, length.out = 50)

rescue_df <- tibble(
  threshold   = thresholds,
  n_rescued   = map_int(thresholds, ~ sum(hb4$hs_fraction < .x)),
  pct_rescued = map_dbl(thresholds, ~ mean(hb4$hs_fraction < .x) * 100)
)

p_rescue <- ggplot(rescue_df, aes(x = threshold, y = pct_rescued)) +
  geom_line(color = "#FF6B00", linewidth = 1.5) +
  geom_point(color = "#FF6B00", size = 2.5) +
  geom_vline(xintercept = threshold_line, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 50, linetype = "dotted", color = "grey60") +
  annotate("text", x = threshold_line * 1.05, y = 2,
           label = paste0("p90 threshold\n", round(threshold_line, 3), "%"),
           hjust = 0, size = 3.5, color = "grey40") +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "LT-HSCs HB4: Fraction of Cells Rescued by Threshold",
    subtitle = "Cells with HS fraction < threshold considered 'not heat shocked'",
    x = "HS fraction threshold (%)",
    y = "% of HB4 cells rescued"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.text     = element_text(size = 12, color = "black"),
    axis.title    = element_text(size = 13, color = "black")
  )

pct_below <- round(mean(hb4$hs_fraction < threshold_line) * 100, 1)
dens_hb4  <- density(hb4$hs_fraction)
ymax_hb4  <- max(dens_hb4$y) * 0.85

p_hb4_density <- ggplot(hb4, aes(x = hs_fraction)) +
  geom_density(fill = "#FF6B00", alpha = 0.4, color = "#FF6B00", linewidth = 1.5) +
  annotate("rect", xmin = -Inf, xmax = threshold_line,
           ymin = 0, ymax = Inf, fill = "green4", alpha = 0.08) +
  geom_vline(xintercept = threshold_line, linetype = "dashed",
             color = "grey30", linewidth = 0.8) +
  annotate("text", x = threshold_line * 0.6, y = ymax_hb4,
           label = paste0(pct_below, "% of cells\nbelow p90 threshold"),
           hjust = 0.5, vjust = 1, size = 4.5, color = "green4", fontface = "bold") +
  labs(
    title    = "LT-HSCs HB4: Heat Shock Fraction Distribution",
    subtitle = "Green shading = cells below p90 threshold (potential rescue fraction)",
    x = "HS gene UMIs / total UMIs × 100 (%)", y = "Density"
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
    title = "LT-HSCs HB4: Cell-Level Rescue Analysis",
    theme = theme(plot.title = element_text(face = "bold", size = 15, color = "black"))
  )

ggsave(file.path(out_dir, "heatshock_fraction_HB4_rescue.pdf"), p_hb4, width = 14, height = 6)
ggsave(file.path(out_dir, "heatshock_fraction_HB4_rescue.png"), p_hb4, width = 14, height = 6, dpi = 160)

# ── Plot 5: Rescue curve for ALL samples ──────────────────────────────────────
rescue_all <- map_dfr(sample_order, function(s) {
  sdata <- df %>% filter(sample == s) %>% pull(hs_fraction)
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
  geom_vline(xintercept = threshold_line, linetype = "dashed",
             color = "grey40", linewidth = 0.6) +
  geom_hline(yintercept = 50, linetype = "dotted", color = "grey60", linewidth = 0.5) +
  annotate("text", x = threshold_line * 1.05, y = 2,
           label = "p90\nthreshold", hjust = 0, size = 3.5, color = "grey40") +
  scale_color_manual(values = rescue_colors, name = "Sample") +
  scale_linetype_manual(values = rescue_linetypes, name = "Sample") +
  scale_linewidth_manual(values = rescue_widths, name = "Sample") +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "LT-HSCs — Fraction of Cells Below HS Fraction Threshold",
    subtitle = "Orange = HB samples  |  Blue = LB samples  |  Dashed = p90 threshold (excl. HB3/HB4)",
    x = "HS fraction threshold (%)",
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

ggsave(file.path(out_dir, "heatshock_fraction_rescue_allsamples.pdf"), p_rescue_all, width = 9, height = 6)
ggsave(file.path(out_dir, "heatshock_fraction_rescue_allsamples.png"), p_rescue_all, width = 9, height = 6, dpi = 160)

# ── Save per-cell fractions ───────────────────────────────────────────────────
df %>%
  dplyr::select(cell_id, sample, hs_fraction, hs_umis, total_umis) %>%
  write_csv(file.path(out_dir, "heatshock_fraction_per_cell.csv"))

cat("\n✓ All outputs saved to:", out_dir, "\n")

