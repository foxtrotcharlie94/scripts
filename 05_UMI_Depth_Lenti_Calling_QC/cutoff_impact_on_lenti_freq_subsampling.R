getwd()
setwd("C:/Users/fc809/Downloads")
library(ggplot2)
library(dplyr)
library(patchwork)

setwd("C:/Users/fc809/Downloads")

# ── Load counts ───────────────────────────────────────────────────────────────
counts_hsc  <- read.table("LT-HSCs_gene_nUMI_per_cell_for_Foivos.txt",
                          header = TRUE, row.names = 1, sep = "\t")
counts_mono <- read.table("classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt",
                          header = TRUE, row.names = 1, sep = "\t")

# ── Downsample function ───────────────────────────────────────────────────────
downsample_cell <- function(cell_counts, target) {
  total <- sum(cell_counts)
  if (total <= target) return(cell_counts)
  sampled <- sample(rep(seq_along(cell_counts), cell_counts), target)
  tabulate(sampled, nbins = length(cell_counts))
}

min_cells <- 50

# ── Original (no subsampling) ─────────────────────────────────────────────────
original_lenti <- function(counts, cell_type_label, cutoff) {
  sample_ids  <- sub("_.*", "", colnames(counts))
  total_umis  <- colSums(counts)
  lenti_umis  <- as.numeric(counts["Lenti", ])
  lenti_ratio <- lenti_umis / total_umis
  lenti_pos   <- lenti_ratio >= cutoff
  
  data.frame(sample    = sample_ids,
             exposure  = sub("[0-9]+$", "", sample_ids),
             lenti_pos = lenti_pos,
             cell_type = cell_type_label) %>%
    group_by(sample, exposure, cell_type) %>%
    summarise(n_cells     = n(),
              n_lenti_pos = sum(lenti_pos),
              pct_lenti   = mean(lenti_pos) * 100,
              .groups     = "drop") %>%
    mutate(target_depth = Inf)
}

# ── Subsampling with min cell filter ─────────────────────────────────────────
lenti_at_depth <- function(counts, target_umis, cell_type_label,
                           cutoff = 1e-4, seed = 42) {
  set.seed(seed)
  sample_ids <- sub("_.*", "", colnames(counts))
  total_umis <- colSums(counts)
  
  keep_cells    <- total_umis >= target_umis
  n_per_sample  <- table(sample_ids[keep_cells])
  valid_samples <- names(n_per_sample[n_per_sample >= min_cells])
  
  if (length(valid_samples) == 0) {
    cat("  Depth", target_umis, "- no valid samples, skipping\n")
    return(NULL)
  }
  
  cat("  Depth:", target_umis, "- valid samples:",
      paste(valid_samples, collapse = ", "), "\n")
  
  keep_cells_valid <- keep_cells & (sample_ids %in% valid_samples)
  counts_sub <- counts[, keep_cells_valid]
  sids_sub   <- sample_ids[keep_cells_valid]
  
  counts_ds <- apply(counts_sub, 2, downsample_cell, target = target_umis)
  rownames(counts_ds) <- rownames(counts_sub)
  
  lenti_umis  <- as.numeric(counts_ds["Lenti", ])
  total_ds    <- colSums(counts_ds)
  lenti_ratio <- lenti_umis / total_ds
  lenti_pos   <- lenti_ratio >= cutoff
  
  data.frame(sample       = sids_sub,
             exposure     = sub("[0-9]+$", "", sids_sub),
             lenti_pos    = lenti_pos,
             target_depth = target_umis,
             cell_type    = cell_type_label) %>%
    group_by(sample, exposure, target_depth, cell_type) %>%
    summarise(n_cells     = n(),
              n_lenti_pos = sum(lenti_pos),
              pct_lenti   = mean(lenti_pos) * 100,
              .groups     = "drop")
}

# ── Run ───────────────────────────────────────────────────────────────────────
depths <- c(1000, 2000, 3000, 4000, 5000, 7500, 10000, 15000)
cutoff <- 1e-4

cat("LT-HSC with cutoff 1e-4...\n")
hsc_results <- bind_rows(
  original_lenti(counts_hsc, "LT-HSC", cutoff),
  bind_rows(Filter(Negate(is.null),
                   lapply(depths, function(d) lenti_at_depth(counts_hsc, d, "LT-HSC", cutoff))))
)

cat("\nClassical Monocytes with cutoff 1e-4...\n")
mono_results <- bind_rows(
  original_lenti(counts_mono, "Classical Monocytes", cutoff),
  bind_rows(Filter(Negate(is.null),
                   lapply(depths, function(d) lenti_at_depth(counts_mono, d, "Classical Monocytes", cutoff))))
)

all_results <- bind_rows(hsc_results, mono_results) %>%
  mutate(depth_label = ifelse(is.infinite(target_depth), "Original",
                              as.character(target_depth)),
         depth_label = factor(depth_label,
                              levels = c(as.character(depths), "Original")))

# ── Plot functions ────────────────────────────────────────────────────────────
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
    labs(title    = paste0("Lenti+ % heatmap across UMI depths - ", cell_type_label),
         subtitle = paste0("Cutoff >= 1e-4 | Excluded where <", min_cells,
                           " cells survive | n = cells used"),
         x = "UMI subsampling depth", y = NULL) +
    theme_classic() +
    theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 9),
          plot.title     = element_text(face = "bold", size = 11),
          plot.subtitle  = element_text(size = 9, color = "grey40"),
          strip.text     = element_text(face = "bold"))
}

plot_line <- function(df, cell_type_label) {
  df_ct <- df %>%
    filter(cell_type == cell_type_label,
           !is.infinite(target_depth))
  
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
    labs(title    = paste0("Impact of UMI subsampling on Lenti+ detection - ", cell_type_label),
         subtitle = "Mean +/- SEM | Cutoff >= 1e-4 | Excluded where <50 cells survive",
         x = "UMI subsampling depth", y = "% Lenti+ cells",
         color = "Exposure") +
    theme_classic() +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 9),
          plot.title      = element_text(face = "bold", size = 11),
          plot.subtitle   = element_text(size = 9, color = "grey40"),
          legend.position = "top")
}

plot_per_sample <- function(df, cell_type_label) {
  df_ct <- df %>%
    filter(cell_type == cell_type_label) %>%
    mutate(is_original = is.infinite(target_depth))
  
  df_lines <- df_ct %>% filter(!is_original)
  df_orig  <- df_ct %>% filter(is_original)
  
  ggplot() +
    geom_line(data = df_lines,
              aes(x = depth_label, y = pct_lenti,
                  color = exposure, group = sample),
              linewidth = 0.8, alpha = 0.8) +
    geom_point(data = df_lines,
               aes(x = depth_label, y = pct_lenti,
                   color = exposure, group = sample),
               size = 2, alpha = 0.9) +
    geom_point(data = df_orig,
               aes(x = depth_label, y = pct_lenti,
                   color = exposure, group = sample),
               size = 4, shape = 18, alpha = 0.9) +
    scale_color_manual(values = c("HB" = "#F44336", "LB" = "#2196F3")) +
    facet_wrap(~ exposure, nrow = 1) +
    labs(title    = paste0("Lenti+ % per sample across UMI depths - ", cell_type_label),
         subtitle = "Lines connect subsampled depths (>=50 cells) | Diamond = original | Cutoff >= 1e-4",
         x = "UMI subsampling depth", y = "% Lenti+ cells",
         color = "Exposure") +
    theme_classic() +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 9),
          plot.title      = element_text(face = "bold", size = 11),
          plot.subtitle   = element_text(size = 9, color = "grey40"),
          legend.position = "none",
          strip.text      = element_text(face = "bold"))
}

# ── Generate and save plots ───────────────────────────────────────────────────
p_heat_hsc  <- plot_heatmap_safe(all_results, "LT-HSC")
p_heat_mono <- plot_heatmap_safe(all_results, "Classical Monocytes")
p_line_hsc  <- plot_line(all_results, "LT-HSC")
p_line_mono <- plot_line(all_results, "Classical Monocytes")
p_samp_hsc  <- plot_per_sample(all_results, "LT-HSC")
p_samp_mono <- plot_per_sample(all_results, "Classical Monocytes")

ggsave("lenti_1e4_heatmap_hsc.pdf",          plot = p_heat_hsc,              width = 13, height = 7)
ggsave("lenti_1e4_heatmap_mono.pdf",         plot = p_heat_mono,             width = 13, height = 7)
ggsave("lenti_1e4_line_combined.pdf",        plot = p_line_hsc / p_line_mono,width = 10, height = 12)
ggsave("lenti_1e4_persample_hsc.pdf",        plot = p_samp_hsc,              width = 12, height = 6)
ggsave("lenti_1e4_persample_mono.pdf",       plot = p_samp_mono,             width = 12, height = 6)

write.csv(all_results, "lenti_1e4_subsampling_results.csv", row.names = FALSE)

cat("All done!\n")

