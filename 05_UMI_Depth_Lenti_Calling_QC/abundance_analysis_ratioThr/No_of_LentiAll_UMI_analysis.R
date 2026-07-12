# ============================================================
# QC for LentiPos / LentiNeg cutoff choice
# Per cluster (cells pooled across HB+LB), examine the
# relationship between sequencing depth (total UMIs) and
# the LentiPos / LentiNeg / ambiguous classification.
# ============================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(tidyr); library(scales)
})
data.table::setDTthreads(10)

THR_POS <- 5e-5
THR_NEG <- 1e-5

ratio_csv <- "C:/Users/fc809/Downloads/lenti_ratio_per_cluster/per_cell_lenti_ratio.csv"
out_dir   <- "C:/Users/fc809/Downloads/lenti_cutoff_QC/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cells <- fread(ratio_csv)
cells[, total := as.numeric(total)]
cells[, lenti := as.numeric(lenti)]
cells[, genotype := fcase(
  ratio >= THR_POS, "LentiPos",
  ratio <= THR_NEG, "LentiNeg",
  default = "ambiguous"
)]
cells[, genotype := factor(genotype, levels = c("LentiNeg", "ambiguous", "LentiPos"))]
suppressWarnings({
  cl_num <- as.numeric(unique(cells$cluster))
  cells[, cluster := factor(cluster,
                            levels = if (any(is.na(cl_num))) sort(unique(cluster))
                            else as.character(sort(cl_num)))]
})
geno_colors <- c("LentiNeg" = "grey60", "ambiguous" = "orange", "LentiPos" = "steelblue")

# ── Per-cluster numeric summary ────────────────────────────────────────────────
summary_tbl <- cells[, .(
  n_cells           = .N,
  pct_pos           = round(100 * mean(genotype == "LentiPos"), 2),
  pct_amb           = round(100 * mean(genotype == "ambiguous"), 2),
  pct_neg           = round(100 * mean(genotype == "LentiNeg"), 2),
  med_total_pos     = median(total[genotype == "LentiPos"]),
  med_total_amb     = median(total[genotype == "ambiguous"]),
  med_total_neg     = median(total[genotype == "LentiNeg"]),
  med_lenti_pos     = median(lenti[genotype == "LentiPos"]),
  wilcox_p_PvN_depth = suppressWarnings(tryCatch(
    wilcox.test(total[genotype == "LentiPos"], total[genotype == "LentiNeg"])$p.value,
    error = function(e) NA_real_))
), by = cluster][order(suppressWarnings(as.numeric(as.character(cluster))))]
print(summary_tbl)
fwrite(summary_tbl, file.path(out_dir, "per_cluster_lentiStatus_QC_summary.csv"))

# ── Plot 1: Scatter — log10(total) vs log10(LentiAll+1), per cluster ──────────
# Subsample for plotting density: cap at 4000 cells per cluster
set.seed(1)
plot_pts <- cells[, .SD[sample(.N, min(.N, 4000))], by = cluster]

x_range <- range(log10(cells$total))
line_x  <- seq(x_range[1], x_range[2], length.out = 300)
cutoff_lines <- data.frame(
  log_total = rep(line_x, 2),
  log_lenti = c(log10(THR_POS * 10^line_x + 1),
                log10(THR_NEG * 10^line_x + 1)),
  cutoff    = rep(c(sprintf("LentiPos (ratio>=%g)", THR_POS),
                    sprintf("LentiNeg (ratio<=%g)", THR_NEG)),
                  each = length(line_x))
)

p_scatter <- ggplot(plot_pts, aes(log10(total), log10(lenti + 1))) +
  geom_point(aes(color = genotype), alpha = 0.25, size = 0.3) +
  geom_line(data = cutoff_lines,
            aes(x = log_total, y = log_lenti, linetype = cutoff),
            color = "black", linewidth = 0.5, inherit.aes = FALSE) +
  scale_color_manual(values = geno_colors) +
  scale_linetype_manual(values = setNames(c("solid","dashed"),
                                          unique(cutoff_lines$cutoff))) +
  facet_wrap(~ cluster) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)),
         linetype = guide_legend(title = NULL)) +
  labs(x = "log10(total UMIs per cell)",
       y = "log10(LentiAll UMIs + 1)",
       title = "LentiAll vs total UMIs by cluster (subsampled to 4k/cluster)",
       color = NULL) +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "1_scatter_lenti_vs_total_per_cluster.pdf"),
       p_scatter, width = 15, height = 11)

# ── Plot 2: Depth distribution by genotype, per cluster ────────────────────────
p_depth <- ggplot(cells, aes(genotype, total, fill = genotype)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  scale_y_log10(labels = scales::comma) +
  scale_fill_manual(values = geno_colors) +
  facet_wrap(~ cluster) +
  labs(x = NULL, y = "total UMIs per cell (log10)",
       title = "Total UMI depth by Lenti status, per cluster",
       subtitle = "If LentiPos sits higher than LentiNeg within a cluster, depth is biasing the call") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(out_dir, "2_depth_by_genotype_per_cluster.pdf"),
       p_depth, width = 15, height = 11)

# ── Plot 3: LentiPos / Neg / ambiguous fraction vs depth bin, per cluster ─────
n_bins <- 8
breaks <- quantile(cells$total, probs = seq(0, 1, length.out = n_bins + 1))
breaks[length(breaks)] <- breaks[length(breaks)] + 1
breaks <- unique(breaks)
cells[, depth_bin := cut(total, breaks = breaks, include.lowest = TRUE, labels = FALSE)]

bin_summary <- cells[, .(
  n_bin        = .N,
  median_depth = median(total),
  pct_neg      = 100 * mean(genotype == "LentiNeg"),
  pct_amb      = 100 * mean(genotype == "ambiguous"),
  pct_pos      = 100 * mean(genotype == "LentiPos")
), by = .(cluster, depth_bin)]
bin_summary <- bin_summary[n_bin >= 30]  # drop bins with too few cells

bin_long <- bin_summary |>
  pivot_longer(cols = c(pct_neg, pct_amb, pct_pos),
               names_to = "category", values_to = "pct") |>
  mutate(category = recode(category,
                           "pct_neg" = "LentiNeg",
                           "pct_amb" = "ambiguous",
                           "pct_pos" = "LentiPos"),
         category = factor(category, levels = c("LentiNeg","ambiguous","LentiPos")))

p_binned <- ggplot(bin_long, aes(median_depth, pct,
                                 color = category, group = category)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.4) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = geno_colors) +
  facet_wrap(~ cluster, scales = "free_y") +
  labs(x = "median total UMIs in bin (log10)",
       y = "% of cells",
       title = "Lenti status vs sequencing depth, per cluster",
       subtitle = "Bin cells by total-UMI deciles; a rising LentiPos% with depth = depth-driven calls",
       color = NULL) +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "3_lentiStatus_by_depthBin_per_cluster.pdf"),
       p_binned, width = 15, height = 11)

# ── Plot 4: What LentiAll count is required at each depth? (reference curves) ─
ref_df <- data.frame(total = 10^seq(2.5, 5.5, length.out = 300)) |>
  mutate(min_pos = THR_POS * total,
         max_neg = THR_NEG * total)

p_ref <- ggplot(ref_df, aes(total)) +
  geom_line(aes(y = min_pos, color = sprintf("LentiPos: min LentiAll UMIs (ratio>=%g)", THR_POS)),
            linewidth = 0.9) +
  geom_line(aes(y = max_neg, color = sprintf("LentiNeg: max LentiAll UMIs (ratio<=%g)", THR_NEG)),
            linewidth = 0.9) +
  geom_hline(yintercept = c(0.5, 1, 2, 5, 10), linetype = "dotted", color = "grey60") +
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 2, 3, 5, 10, 20, 50, 100)) +
  scale_color_manual(values = setNames(c("steelblue","grey30"),
                                       c(sprintf("LentiPos: min LentiAll UMIs (ratio>=%g)", THR_POS),
                                         sprintf("LentiNeg: max LentiAll UMIs (ratio<=%g)", THR_NEG)))) +
  labs(x = "Cell total UMIs (log10)",
       y = "LentiAll UMI count",
       title = "What LentiAll count does each threshold correspond to?",
       subtitle = "Cells above the blue line are called LentiPos; cells below the grey line are LentiNeg",
       color = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top", legend.text = element_text(size = 10))
ggsave(file.path(out_dir, "4_lenti_count_thresholds_vs_depth.pdf"),
       p_ref, width = 9, height = 6)

# ── Plot 5: LentiAll UMI distribution in LentiPos cells, per cluster ──────────
pos_cells <- cells[genotype == "LentiPos"]
x_max <- as.numeric(quantile(pos_cells$lenti, 0.99))

p_pos_dist <- ggplot(pos_cells, aes(lenti)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white",
                 boundary = 0.5) +
  coord_cartesian(xlim = c(0, x_max)) +
  facet_wrap(~ cluster, scales = "free_y") +
  labs(x = "LentiAll UMIs per cell (LentiPos cells only)",
       y = "cells",
       title = "Distribution of LentiAll UMIs in LentiPos cells, per cluster",
       subtitle = sprintf("Integer-bin histograms; x capped at p99 (%.0f UMIs)", x_max)) +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "5_LentiAll_distribution_LentiPos_per_cluster.pdf"),
       p_pos_dist, width = 15, height = 11)

cat("\nDone. Outputs in:", out_dir, "\n")

# ── Plot 6: Downsampling robustness of LentiPos calls, per cluster ────────────
# For each cluster, take LentiPos cells. Iterate N = 4K, 6K, 8K, ...:
# keep cells with original total >= N, hypergeometric-downsample each cell's
# UMIs to k=N, recompute ratio, ask how many remain LentiPos. Per-cluster
# stop rule: drop the cluster at any N where fewer than 500 cells qualify.
set.seed(42)

N_start <- 500L
N_step  <- 500L
min_cells_per_cluster <- 100L

N_levels <- seq(N_start, max(pos_cells$total), by = N_step)

downsample_one_N <- function(sub, n) {
  if (nrow(sub) < min_cells_per_cluster) return(NULL)
  new_lenti <- rhyper(nrow(sub),
                      m = sub$lenti, nn = sub$total - sub$lenti, k = n)
  new_ratio <- new_lenti / n
  data.table(N              = n,
             n_cells        = nrow(sub),
             n_still_pos    = sum(new_ratio >= THR_POS),
             frac_still_pos = mean(new_ratio >= THR_POS))
}

results <- pos_cells[, {
  dat <- copy(.SD)
  rbindlist(lapply(N_levels, function(n) downsample_one_N(dat[total >= n], n)))
}, by = cluster, .SDcols = c("lenti", "total")]

fwrite(results, file.path(out_dir, "6_lentiPos_downsampling_robustness.csv"))

p_down <- ggplot(results, aes(N, 100 * frac_still_pos,
                              color = cluster, group = cluster)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.1) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(limits = c(NA, 100),
                     breaks = seq(0, 100, by = 20)) +
  labs(x = "Downsampling target N (total UMIs per cell)",
       y = "% of LentiPos cells still called LentiPos after downsampling",
       title = "Downsampling robustness of LentiPos calls, per cluster",
       subtitle = sprintf(
         "Hypergeometric downsample to k=N; only cells with original total >= N kept (>= %d cells/cluster)",
         min_cells_per_cluster),
       color = "cluster") +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "6_lentiPos_downsampling_robustness.pdf"),
       p_down, width = 11, height = 7)

cat("\nDone. Outputs in:", out_dir, "\n")

# ── Plot 7: Same downsampling analysis but pooled across all clusters ─────────
pooled_results <- rbindlist(lapply(N_levels, function(n) {
  downsample_one_N(pos_cells[total >= n, .(lenti, total)], n)
}))
fwrite(pooled_results, file.path(out_dir, "7_lentiPos_downsampling_robustness_pooled.csv"))

p_down_pooled <- ggplot(pooled_results, aes(N, 100 * frac_still_pos)) +
  geom_line(linewidth = 0.9, color = "steelblue") +
  geom_point(size = 2, color = "steelblue") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(limits = c(NA, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "Downsampling target N (total UMIs per cell)",
       y = "% of LentiPos cells still called LentiPos after downsampling",
       title = "Downsampling robustness of LentiPos calls (all clusters pooled)",
       subtitle = sprintf(
         "Hypergeometric downsample to k=N; only cells with original total >= N kept (>= %d cells)",
         min_cells_per_cluster)) +
  theme_bw(base_size = 13)
ggsave(file.path(out_dir, "7_lentiPos_downsampling_robustness_pooled.pdf"),
       p_down_pooled, width = 9, height = 6)

cat("\nDone. Outputs in:", out_dir, "\n")

# ── Plot 8: Cells surviving total-UMI cutoffs, per cluster ────────────────────
cutoffs <- c(0, 5000, 10000, 15000, 20000)

surv <- rbindlist(lapply(cutoffs, function(co) {
  cells[total >= co, .(cutoff = co, n_cells = .N), by = cluster]
}))
surv[, cutoff_label := factor(
  ifelse(cutoff == 0, "all cells",
         paste0(">= ", format(cutoff, big.mark = ","))),
  levels = c("all cells",
             paste0(">= ", format(cutoffs[cutoffs > 0], big.mark = ","))))]

# Wide CSV for easy reading
fwrite(dcast(surv, cluster ~ cutoff_label, value.var = "n_cells"),
       file.path(out_dir, "8_cells_per_cluster_at_UMI_cutoffs.csv"))

p_surv <- ggplot(surv, aes(cluster, n_cells, fill = cutoff_label)) +
  geom_col(position = position_dodge2(preserve = "single", padding = 0.1),
           width = 0.85) +
  scale_fill_brewer(palette = "YlGnBu", name = "total UMIs") +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Cluster",
       y = "Number of cells",
       title = "Cells per cluster surviving total-UMI cutoffs",
       subtitle = "Counts include all cells (any Lenti status)") +
  theme_bw(base_size = 12) +
  theme(axis.text.x      = element_text(size = 11),
        panel.grid.major.x = element_blank())
ggsave(file.path(out_dir, "8_cells_per_cluster_at_UMI_cutoffs.pdf"),
       p_surv, width = 13, height = 6)

cat("\nDone. Outputs in:", out_dir, "\n")

