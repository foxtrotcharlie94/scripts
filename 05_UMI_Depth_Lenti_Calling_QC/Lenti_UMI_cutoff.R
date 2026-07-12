library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

setwd("C:/Users/fc809/Downloads")

# ── Load counts ───────────────────────────────────────────────────────────────
counts_hsc  <- read.table("LT-HSCs_gene_nUMI_per_cell_for_Foivos.txt",
                          header = TRUE, row.names = 1, sep = "\t")
counts_mono <- read.table("classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt",
                          header = TRUE, row.names = 1, sep = "\t")

threshold_pos <- 4e-5
threshold_neg <- 1e-5
plot_floor    <- 5e-8

get_lenti_df <- function(counts, cell_type_label) {
  total_umis  <- colSums(counts)
  lenti_umis  <- as.numeric(counts["Lenti", ])
  lenti_ratio <- lenti_umis / total_umis
  sample_ids  <- sub("_.*", "", colnames(counts))
  data.frame(
    cell         = colnames(counts),
    sample       = sample_ids,
    exposure     = sub("[0-9]+$", "", sample_ids),
    total_umis   = total_umis,
    lenti_umis   = lenti_umis,
    lenti_ratio  = lenti_ratio,
    lenti_status = ifelse(lenti_ratio >= threshold_pos, "Lenti+",
                          ifelse(lenti_ratio <= threshold_neg, "Lenti-", "Unclassified")),
    cell_type    = cell_type_label
  )
}

df_hsc  <- get_lenti_df(counts_hsc,  "LT-HSC")
df_mono <- get_lenti_df(counts_mono, "Classical Monocytes")
df_all  <- bind_rows(df_hsc, df_mono) %>%
  mutate(cell_type = factor(cell_type, levels = c("LT-HSC", "Classical Monocytes")))

# UMI bin definitions (shared across approaches)
umi_breaks <- c(0, seq(1000, 20000, 1000), Inf)
umi_labels <- c("<1k", paste0(seq(1, 19), "-", seq(2, 20), "k"), ">20k")
bins_mid   <- c(500, seq(1500, 19500, 1000), 21000)

# ══════════════════════════════════════════════════════════════════════════════
# APPROACH 1: Binomial detection probability
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Approach 1: Binomial ──\n")

umi_range    <- c(100, 200, 500, seq(1000, 50000, 500))
true_ratios  <- c(4e-5, 8e-5, 1.5e-4, 3e-4, 1e-3)
ratio_labels <- c("4e-5 (threshold)", "8e-5", "1.5e-4", "3e-4", "1e-3")

binom_df <- expand.grid(umi_depth = umi_range, true_ratio = true_ratios) %>%
  mutate(prob_detect = (1 - (1 - true_ratio)^umi_depth) * 100,
         ratio_label = factor(true_ratio, levels = true_ratios,
                              labels = ratio_labels))

min_umi_df <- expand.grid(
  true_ratio = exp(seq(log(1e-5), log(1e-2), length.out = 300)),
  power      = c(0.50, 0.80, 0.90)) %>%
  mutate(min_umi     = ceiling(log(1 - power) / log(1 - true_ratio)),
         power_label = paste0(power * 100, "% power"))

min_umi_80 <- ceiling(log(0.20) / log(1 - threshold_pos))
cat("80% power at 4e-5:", min_umi_80, "UMIs\n")

p1a <- ggplot(binom_df, aes(x = umi_depth, y = prob_detect,
                            color = ratio_label, group = ratio_label)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 80, linetype = "dashed",
             color = "grey50", linewidth = 0.7) +
  geom_vline(xintercept = min_umi_80, linetype = "dashed",
             color = "#8B0000", linewidth = 0.9) +
  annotate("text", x = min_umi_80 * 1.6, y = 8,
           label = paste0("80%: ", scales::comma(min_umi_80), " UMIs"),
           size = 3.2, color = "#8B0000", fontface = "bold", hjust = 0) +
  annotate("text", x = 47000, y = 81.5, label = "80%",
           size = 3, color = "grey40") +
  scale_x_log10(labels = scales::comma) +
  scale_color_brewer(palette = "Set1") +
  labs(title    = "Binomial detection probability",
       subtitle = "P(detecting >=1 Lenti UMI) | Vertical line = 80% at 4e-5 threshold",
       x = "Total UMIs per cell", y = "Detection probability (%)",
       color = "True Lenti ratio") +
  theme_classic() +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40"))

p1b <- ggplot(min_umi_df, aes(x = true_ratio, y = min_umi,
                              color = power_label, group = power_label)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = threshold_pos, linetype = "dashed",
             color = "#2CA02C", linewidth = 0.8) +
  geom_point(data = data.frame(true_ratio  = threshold_pos,
                               min_umi     = min_umi_80,
                               power_label = "80% power"),
             size = 4, shape = 21, fill = "white", stroke = 1.8) +
  annotate("text", x = threshold_pos * 1.6, y = min_umi_80 * 0.5,
           label = paste0("80%: ", scales::comma(min_umi_80), " UMIs"),
           size = 3.2, color = "#8B0000", fontface = "bold") +
  annotate("text", x = threshold_pos * 0.35, y = 45000,
           label = "4e-5\nthreshold", size = 3, color = "#2CA02C") +
  scale_x_log10(labels = scales::scientific) +
  scale_y_log10(labels = scales::comma) +
  scale_color_manual(values = c("50% power" = "#FF7F0E",
                                "80% power" = "#8B0000",
                                "90% power" = "#1F77B4")) +
  labs(title    = "Min UMIs for detection power",
       subtitle = "Circle marks 80% crossing at 4e-5",
       x = "True Lenti ratio", y = "Min UMIs per cell required",
       color = "Power level") +
  theme_classic() +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40"))

ggsave("detection_approach1_binomial.pdf",
       plot = p1a | p1b, width = 14, height = 6)
cat("Approach 1 saved\n")

# ══════════════════════════════════════════════════════════════════════════════
# APPROACH 2: Empirical consistency from subsampling
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Approach 2: Empirical consistency ──\n")

downsample_cell <- function(cell_counts, target) {
  total <- sum(cell_counts)
  if (total <= target) return(cell_counts)
  sampled <- sample(rep(seq_along(cell_counts), cell_counts), target)
  tabulate(sampled, nbins = length(cell_counts))
}

depths <- c(1000, 2000, 3000, 4000, 5000, 7500, 10000, 15000)

empirical_consistency <- function(counts, df_full, cell_type_label) {
  lenti_pos_cells <- df_full %>%
    filter(cell_type == cell_type_label, lenti_status == "Lenti+") %>%
    pull(cell)
  sample_ids <- sub("_.*", "", colnames(counts))
  total_umis <- colSums(counts)
  
  bind_rows(lapply(depths, function(d) {
    set.seed(42)
    keep <- colnames(counts) %in% lenti_pos_cells & total_umis >= d
    valid_samples <- names(which(table(sample_ids[keep]) >= 10))
    keep <- keep & sample_ids %in% valid_samples
    if (sum(keep) < 10) return(NULL)
    
    counts_sub    <- counts[, keep]
    sids_sub      <- sample_ids[keep]
    counts_ds     <- apply(counts_sub, 2, downsample_cell, target = d)
    rownames(counts_ds) <- rownames(counts_sub)
    
    lenti_ratio_ds <- as.numeric(counts_ds["Lenti", ]) / colSums(counts_ds)
    
    data.frame(depth           = d,
               sample          = sids_sub,
               exposure        = sub("[0-9]+$", "", sids_sub),
               still_lenti_pos = lenti_ratio_ds >= threshold_pos,
               cell_type       = cell_type_label) %>%
      group_by(depth, exposure, cell_type) %>%
      summarise(n_cells      = n(),
                pct_retained = mean(still_lenti_pos) * 100,
                .groups      = "drop")
  }))
}

cat("  LT-HSC empirical consistency...\n")
consist_hsc  <- empirical_consistency(counts_hsc,  df_all, "LT-HSC")
cat("  Classical Monocytes empirical consistency...\n")
consist_mono <- empirical_consistency(counts_mono, df_all, "Classical Monocytes")
consist_all  <- bind_rows(consist_hsc, consist_mono) %>%
  mutate(cell_type = factor(cell_type, levels = c("LT-HSC", "Classical Monocytes")))

summary_consist <- consist_all %>%
  group_by(cell_type, depth, exposure) %>%
  summarise(mean_pct = mean(pct_retained),
            sem_pct  = sd(pct_retained) / sqrt(n()),
            .groups  = "drop")

depth_80 <- summary_consist %>%
  filter(mean_pct >= 80) %>%
  group_by(cell_type, exposure) %>%
  summarise(depth_80 = min(depth), .groups = "drop")
cat("Depth at 80% consistency:\n"); print(depth_80)

p2 <- ggplot(summary_consist, aes(x = factor(depth), y = mean_pct,
                                  color = exposure, group = exposure)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_pct - sem_pct,
                    ymax = mean_pct + sem_pct),
                width = 0.2, linewidth = 0.7) +
  geom_hline(yintercept = 80, linetype = "dashed",
             color = "#8B0000", linewidth = 0.9) +
  geom_text(data = depth_80,
            aes(x = factor(depth_80),
                y = ifelse(exposure == "HB", 75, 68),
                label = paste0(scales::comma(depth_80), "\nUMIs"),
                color = exposure),
            size = 2.8, fontface = "bold", show.legend = FALSE) +
  annotate("text", x = 8.2, y = 81, label = "80%",
           size = 3, color = "#8B0000", fontface = "bold") +
  scale_color_manual(values = c("HB" = "#F44336", "LB" = "#2196F3")) +
  facet_wrap(~ cell_type) +
  labs(title    = "Empirical Lenti+ classification consistency",
       subtitle = "% of originally-Lenti+ cells still called Lenti+ after downsampling\nLabels = UMI depth where 80% consistency first reached",
       x = "UMI subsampling depth", y = "% Lenti+ cells retained",
       color = "Exposure") +
  theme_classic() +
  theme(plot.title      = element_text(face = "bold"),
        plot.subtitle   = element_text(size = 9, color = "grey40"),
        strip.text      = element_text(face = "bold"),
        legend.position = "top",
        axis.text.x     = element_text(angle = 45, hjust = 1))

ggsave("detection_approach2_empirical.pdf", plot = p2, width = 12, height = 6)
cat("Approach 2 saved\n")

# ══════════════════════════════════════════════════════════════════════════════
# APPROACH 3: Wilson CI per cell
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Approach 3: Wilson CI ──\n")

wilson_ci <- function(k, n, conf = 0.95) {
  z      <- qnorm(1 - (1 - conf) / 2)
  denom  <- 1 + z^2 / n
  center <- (k/n + z^2/(2*n)) / denom
  margin <- z * sqrt((k/n) * (1 - k/n) / n + z^2/(4*n^2)) / denom
  list(lower = pmax(0, center - margin), upper = center + margin)
}

add_ci <- function(df) {
  ci <- wilson_ci(df$lenti_umis, df$total_umis)
  df %>% mutate(
    ci_lower   = ci$lower,
    ci_upper   = ci$upper,
    conf_class = case_when(
      lenti_umis == 0          ~ "Confident Lenti-",  # exact zero, no ambiguity
      ci_lower > threshold_pos ~ "Confident Lenti+",
      TRUE                     ~ "Uncertain"
    ))
}

df_hsc_ci  <- add_ci(df_hsc)
df_mono_ci <- add_ci(df_mono)
df_all_ci  <- bind_rows(df_hsc_ci, df_mono_ci) %>%
  mutate(cell_type = factor(cell_type, levels = c("LT-HSC", "Classical Monocytes")))

# Scatter — zero-ratio cells plotted at floor
set.seed(42)
df_plot_ci <- df_all_ci %>%
  mutate(lenti_ratio_plot = ifelse(lenti_ratio == 0, plot_floor, lenti_ratio)) %>%
  group_by(cell_type) %>%
  group_modify(~ slice_sample(.x, n = min(5000, nrow(.x)))) %>%
  ungroup()

p3a <- ggplot(df_plot_ci, aes(x = total_umis, y = lenti_ratio_plot,
                              color = conf_class)) +
  geom_point(size = 0.4, alpha = 0.4) +
  geom_hline(yintercept = threshold_pos, linetype = "dashed",
             color = "#2CA02C", linewidth = 0.8) +
  geom_hline(yintercept = plot_floor, linetype = "dotted",
             color = "grey60", linewidth = 0.5) +
  annotate("text", x = 600, y = plot_floor * 1.6,
           label = "0 Lenti UMIs (floor)", size = 2.5,
           color = "grey50", hjust = 0) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::scientific,
                limits = c(plot_floor * 0.5, NA)) +
  scale_color_manual(values = c("Confident Lenti+" = "#2CA02C",
                                "Confident Lenti-" = "#9467BD",
                                "Uncertain"        = "#CCCCCC")) +
  facet_wrap(~ cell_type) +
  labs(title    = "Per-cell Lenti classification confidence (Wilson CI)",
       subtitle = "Zero Lenti UMIs plotted at floor (5e-8) | 5000 cells sampled per cell type",
       x = "Total UMIs per cell", y = "Observed Lenti ratio",
       color = "Classification") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

# Stacked bar — uncertainty computed among nonzero-Lenti cells
df_binned <- df_all_ci %>%
  mutate(umi_bin = cut(total_umis, breaks = umi_breaks, labels = umi_labels)) %>%
  group_by(cell_type, umi_bin) %>%
  summarise(pct_conf_pos = mean(conf_class == "Confident Lenti+") * 100,
            pct_uncertain = mean(conf_class == "Uncertain") * 100,
            pct_conf_neg  = mean(conf_class == "Confident Lenti-") * 100,
            n_cells       = n(),
            .groups       = "drop") %>%
  pivot_longer(cols = c(pct_conf_pos, pct_uncertain, pct_conf_neg),
               names_to = "category", values_to = "pct") %>%
  mutate(category = recode(category,
                           "pct_conf_pos"  = "Confident Lenti+",
                           "pct_uncertain" = "Uncertain",
                           "pct_conf_neg"  = "Confident Lenti-"),
         category = factor(category,
                           levels = c("Confident Lenti+",
                                      "Uncertain",
                                      "Confident Lenti-")))

# Annotation: first bin where <20% of nonzero-Lenti cells are uncertain
cutoff_3 <- df_all_ci %>%
  mutate(umi_bin = cut(total_umis, breaks = umi_breaks, labels = umi_labels)) %>%
  filter(lenti_umis > 0) %>%
  group_by(cell_type, umi_bin) %>%
  summarise(pct_unc = mean(conf_class == "Uncertain") * 100,
            n       = n(),
            .groups = "drop") %>%
  filter(pct_unc < 20, n >= 10) %>%
  group_by(cell_type) %>%
  summarise(first_bin = first(umi_bin), .groups = "drop")
cat("Approach 3 - <20% uncertain (nonzero-Lenti cells) from:\n"); print(cutoff_3)

p3b <- ggplot(df_binned, aes(x = umi_bin, y = pct, fill = category)) +
  geom_bar(stat = "identity", width = 0.85) +
  geom_segment(data = cutoff_3,
               aes(x = first_bin, xend = first_bin, y = 101, yend = 107),
               inherit.aes = FALSE, color = "#8B0000", linewidth = 0.8,
               arrow = arrow(length = unit(0.15, "cm"), ends = "first")) +
  geom_label(data = cutoff_3,
             aes(x = first_bin, y = 109,
                 label = paste0("<20% uncertain\n(nonzero cells)\nfrom ", first_bin)),
             inherit.aes = FALSE, size = 2.5, fontface = "bold",
             color = "#8B0000", fill = "white", label.size = 0) +
  scale_fill_manual(values = c("Confident Lenti+" = "#2CA02C",
                               "Uncertain"        = "#CCCCCC",
                               "Confident Lenti-" = "#9467BD")) +
  scale_y_continuous(limits = c(0, 125)) +
  facet_wrap(~ cell_type) +
  labs(title    = "Classification confidence by UMI depth bin",
       subtitle = "Wilson 95% CI | Zeros = Confident Lenti- | Arrow = first bin <20% uncertain among nonzero-Lenti cells",
       x = "Total UMIs per cell (binned)", y = "% of cells",
       fill = "Classification") +
  theme_classic() +
  theme(plot.title      = element_text(face = "bold"),
        plot.subtitle   = element_text(size = 9, color = "grey40"),
        strip.text      = element_text(face = "bold"),
        axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "top")

ggsave("detection_approach3_CI.pdf",
       plot = p3a / p3b, width = 14, height = 12)
cat("Approach 3 saved\n")

# ══════════════════════════════════════════════════════════════════════════════
# APPROACH 4: Beta mixture model (fixed pi = 0.5)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Approach 4: Beta mixture ──\n")

fit_beta_mixture <- function(df, cell_type_label) {
  ratios <- df %>%
    filter(cell_type == cell_type_label,
           lenti_ratio > 0, lenti_ratio < 1) %>%
    pull(lenti_ratio)
  cat("  Fitting Beta mixture for", cell_type_label,
      "- n =", length(ratios), "\n")
  
  set.seed(42)
  log_r <- log10(ratios)
  km    <- kmeans(log_r, centers = 2, nstart = 20)
  
  beta_params <- function(m, v) {
    v     <- min(v, m*(1-m)*0.99)
    alpha <- m * (m*(1-m)/v - 1)
    beta  <- (1-m) * (m*(1-m)/v - 1)
    c(max(0.01, alpha), max(0.01, beta))
  }
  
  grp1 <- ratios[km$cluster == which.min(km$centers)]
  grp2 <- ratios[km$cluster == which.max(km$centers)]
  p1   <- beta_params(mean(grp1), var(grp1))
  p2   <- beta_params(mean(grp2), var(grp2))
  
  for (i in 1:50) {
    d1  <- tryCatch(0.5 * dbeta(ratios, p1[1], p1[2]),
                    error = function(e) rep(1e-300, length(ratios)))
    d2  <- tryCatch(0.5 * dbeta(ratios, p2[1], p2[2]),
                    error = function(e) rep(1e-300, length(ratios)))
    tot <- d1 + d2 + 1e-300
    r1  <- d1 / tot;  r2 <- d2 / tot
    m1  <- sum(r1 * ratios) / sum(r1)
    m2  <- sum(r2 * ratios) / sum(r2)
    v1  <- sum(r1 * (ratios - m1)^2) / sum(r1)
    v2  <- sum(r2 * (ratios - m2)^2) / sum(r2)
    p1  <- beta_params(m1, v1)
    p2  <- beta_params(m2, v2)
  }
  cat("  Comp1 (Lenti-): alpha=", round(p1[1],3), "beta=", round(p1[2],3),
      "| mean=", round(p1[1]/(p1[1]+p1[2]),8), "\n")
  cat("  Comp2 (Lenti+): alpha=", round(p2[1],3), "beta=", round(p2[2],3),
      "| mean=", round(p2[1]/(p2[1]+p2[2]),6), "\n")
  list(ratios = ratios, p1 = p1, p2 = p2)
}

fit_hsc  <- fit_beta_mixture(df_all, "LT-HSC")
fit_mono <- fit_beta_mixture(df_all, "Classical Monocytes")

plot_mixture <- function(fit, cell_type_label) {
  ratios  <- fit$ratios
  r_min   <- max(min(ratios), 1e-8)
  r_max   <- quantile(ratios, 0.999)
  r_range <- exp(seq(log(r_min), log(r_max), length.out = 500))
  
  mix_df <- data.frame(
    ratio   = r_range,
    comp1   = 0.5 * dbeta(r_range, fit$p1[1], fit$p1[2]),
    comp2   = 0.5 * dbeta(r_range, fit$p2[1], fit$p2[2])
  ) %>% mutate(mixture = comp1 + comp2)
  
  ggplot() +
    geom_histogram(aes(x = ratios, y = after_stat(density)),
                   bins = 80, fill = "grey80", color = "white", alpha = 0.7) +
    geom_line(data = mix_df, aes(x = ratio, y = comp1,
                                 color = "Lenti- component"),
              linewidth = 1.2) +
    geom_line(data = mix_df, aes(x = ratio, y = comp2,
                                 color = "Lenti+ component"),
              linewidth = 1.2) +
    geom_line(data = mix_df, aes(x = ratio, y = mixture, color = "Mixture"),
              linewidth = 1.5, linetype = "dashed") +
    geom_vline(xintercept = threshold_pos, color = "#2CA02C",
               linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = threshold_neg, color = "#9467BD",
               linetype = "dashed", linewidth = 0.8) +
    scale_x_log10(labels = scales::scientific, limits = c(1e-8, NA)) +
    scale_color_manual(values = c("Lenti- component" = "#9467BD",
                                  "Lenti+ component"  = "#2CA02C",
                                  "Mixture"           = "black")) +
    labs(title    = paste0("Beta mixture - ", cell_type_label),
         subtitle = "pi fixed = 0.5 | Purple dashed = threshold_neg | Green dashed = threshold_pos",
         x = "Lenti UMIs / Total UMIs", y = "Density", color = NULL) +
    theme_classic() +
    theme(plot.title      = element_text(face = "bold"),
          plot.subtitle   = element_text(size = 9, color = "grey40"),
          legend.position = "top")
}

add_posterior <- function(df, fit, cell_type_label) {
  df_pos <- df %>%
    filter(cell_type == cell_type_label,
           lenti_ratio > 0, lenti_ratio < 1)
  d1 <- tryCatch(0.5 * dbeta(df_pos$lenti_ratio, fit$p1[1], fit$p1[2]),
                 error = function(e) rep(1e-300, nrow(df_pos)))
  d2 <- tryCatch(0.5 * dbeta(df_pos$lenti_ratio, fit$p2[1], fit$p2[2]),
                 error = function(e) rep(1e-300, nrow(df_pos)))
  df_pos$prob_lenti_pos <- d2 / (d1 + d2 + 1e-300)
  df_pos
}

df_hsc_post  <- add_posterior(df_hsc,  fit_hsc,  "LT-HSC")
df_mono_post <- add_posterior(df_mono, fit_mono, "Classical Monocytes")
df_posterior <- bind_rows(df_hsc_post, df_mono_post) %>%
  mutate(cell_type = factor(cell_type,
                            levels = c("LT-HSC", "Classical Monocytes")))

# Stacked bar — uncertainty among nonzero-Lenti cells only
df_post_binned <- df_posterior %>%
  mutate(umi_bin = cut(total_umis, breaks = umi_breaks, labels = umi_labels)) %>%
  group_by(cell_type, umi_bin) %>%
  summarise(pct_high_conf = mean(prob_lenti_pos > 0.8 |
                                   prob_lenti_pos < 0.2) * 100,
            pct_uncertain = mean(prob_lenti_pos >= 0.2 &
                                   prob_lenti_pos <= 0.8) * 100,
            n_cells       = n(),
            .groups       = "drop")

cutoff_4 <- df_post_binned %>%
  filter(pct_high_conf >= 80) %>%
  group_by(cell_type) %>%
  summarise(first_bin = first(umi_bin), .groups = "drop")
cat("Approach 4 - >80% confident (nonzero-Lenti) from:\n"); print(cutoff_4)

p4_mix <- plot_mixture(fit_hsc,  "LT-HSC") |
  plot_mixture(fit_mono, "Classical Monocytes")

p4_bar <- df_post_binned %>%
  pivot_longer(c(pct_high_conf, pct_uncertain),
               names_to = "cat", values_to = "pct") %>%
  mutate(cat = recode(cat,
                      "pct_high_conf" = "High confidence",
                      "pct_uncertain" = "Uncertain"),
         cat = factor(cat, levels = c("High confidence", "Uncertain"))) %>%
  ggplot(aes(x = umi_bin, y = pct, fill = cat)) +
  geom_bar(stat = "identity", width = 0.85) +
  geom_hline(yintercept = 80, linetype = "dashed",
             color = "#8B0000", linewidth = 0.9) +
  geom_segment(data = cutoff_4,
               aes(x = first_bin, xend = first_bin, y = 101, yend = 107),
               inherit.aes = FALSE, color = "#8B0000", linewidth = 0.8,
               arrow = arrow(length = unit(0.15, "cm"), ends = "first")) +
  geom_label(data = cutoff_4,
             aes(x = first_bin, y = 109,
                 label = paste0(">80% confident\nfrom ", first_bin)),
             inherit.aes = FALSE, size = 2.8, fontface = "bold",
             color = "#8B0000", fill = "white", label.size = 0) +
  annotate("text", x = length(umi_labels), y = 81.5, label = "80%",
           size = 3, color = "#8B0000", fontface = "bold", hjust = 1) +
  scale_fill_manual(values = c("High confidence" = "#2166AC",
                               "Uncertain"       = "#AAAAAA")) +
  scale_y_continuous(limits = c(0, 125)) +
  facet_wrap(~ cell_type) +
  labs(title    = "Posterior confidence by UMI bin",
       subtitle = "Among nonzero-Lenti cells | High confidence = P(Lenti+) > 0.8 or < 0.2 | pi fixed = 0.5",
       x = "Total UMIs per cell (binned)", y = "% cells", fill = NULL) +
  theme_classic() +
  theme(plot.title      = element_text(face = "bold"),
        plot.subtitle   = element_text(size = 9, color = "grey40"),
        strip.text      = element_text(face = "bold"),
        axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "top")

ggsave("detection_approach4_mixture.pdf",
       plot = p4_mix / p4_bar, width = 14, height = 12)
cat("Approach 4 saved\n")

# ══════════════════════════════════════════════════════════════════════════════
# APPROACH 5: Summary overlay
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Approach 5: Summary ──\n")

umi_x <- seq(1000, 20000, 1000)

# Approach 1: binomial at 1e-4
prob_detect <- (1 - (1 - 1e-4)^umi_x) * 100

# Approach 2: average across cell types and exposures
consist_avg    <- summary_consist %>%
  group_by(depth) %>%
  summarise(mean_all = mean(mean_pct), .groups = "drop") %>%
  arrange(depth)
consist_interp <- approx(consist_avg$depth, consist_avg$mean_all,
                         xout = umi_x, rule = 2)$y

# Approach 3: % confident among nonzero-Lenti cells
conf3_by_bin <- df_all_ci %>%
  mutate(umi_bin = cut(total_umis, breaks = umi_breaks, labels = umi_labels)) %>%
  filter(lenti_umis > 0) %>%
  group_by(umi_bin) %>%
  summarise(pct_conf = mean(conf_class != "Uncertain") * 100,
            n        = n(),
            .groups  = "drop") %>%
  filter(n >= 10) %>%
  mutate(bin_mid = bins_mid[match(umi_bin, umi_labels)])
conf3_interp <- approx(conf3_by_bin$bin_mid, conf3_by_bin$pct_conf,
                       xout = umi_x, rule = 2)$y

# Approach 4: % high posterior among nonzero-Lenti cells
conf4_by_bin <- df_posterior %>%
  mutate(umi_bin = cut(total_umis, breaks = umi_breaks, labels = umi_labels)) %>%
  filter(lenti_umis > 0) %>%
  group_by(umi_bin) %>%
  summarise(pct_conf = mean(prob_lenti_pos > 0.8 |
                              prob_lenti_pos < 0.2) * 100,
            n        = n(),
            .groups  = "drop") %>%
  filter(n >= 10) %>%
  mutate(bin_mid = bins_mid[match(umi_bin, umi_labels)])
conf4_interp <- approx(conf4_by_bin$bin_mid, conf4_by_bin$pct_conf,
                       xout = umi_x, rule = 2)$y

summary_df <- data.frame(
  umi_depth = rep(umi_x, 4),
  value     = c(prob_detect, consist_interp, conf3_interp, conf4_interp),
  approach  = rep(c("Approach 1: Binomial (ratio=1e-4)",
                    "Approach 2: Empirical consistency",
                    "Approach 3: Confident CI",
                    "Approach 4: High posterior (Beta mixture)"),
                  each = length(umi_x))
)

cols_app <- c("Approach 1: Binomial (ratio=1e-4)"        = "#1F77B4",
              "Approach 2: Empirical consistency"         = "#FF7F0E",
              "Approach 3: Confident CI"                  = "#2CA02C",
              "Approach 4: High posterior (Beta mixture)" = "#9467BD")

crossings <- summary_df %>%
  group_by(approach) %>%
  arrange(umi_depth) %>%
  summarise(
    cross_umi = {
      v <- value; u <- umi_depth
      if (v[1] >= 80) NA_real_   # already above 80% at 1k — exclude from crossing
      else {
        idx <- which(v >= 80)[1]
        if (is.na(idx)) NA_real_
        else u[idx-1] + (80 - v[idx-1]) * (u[idx] - u[idx-1]) /
          (v[idx] - v[idx-1])
      }
    },
    starts_above = value[1] >= 80,
    .groups      = "drop"
  )

crossings_finite <- crossings %>% filter(!is.na(cross_umi))
conservative_cut <- if (nrow(crossings_finite) > 0)
  max(crossings_finite$cross_umi) else NA
cat("80% crossings:\n"); print(crossings)
cat("Conservative cutoff:", round(conservative_cut), "UMIs\n")

# Stagger label y positions
y_pos_vec <- c(72, 62, 52, 42)

p5 <- ggplot(summary_df, aes(x = umi_depth, y = value,
                             color = approach, group = approach)) +
  { if (!is.na(conservative_cut))
    geom_rect(aes(xmin = conservative_cut, xmax = 20000,
                  ymin = 80, ymax = 100),
              fill = "#8B0000", alpha = 0.06, inherit.aes = FALSE) } +
  geom_line(linewidth = 2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 80, linetype = "dashed",
             color = "#8B0000", linewidth = 1.2) +
  annotate("text", x = 20500, y = 81, label = "80%",
           size = 3.5, color = "#8B0000", fontface = "bold", hjust = 0) +
  # Vertical crossing lines only for approaches that cross within the range
  geom_vline(data = crossings_finite,
             aes(xintercept = cross_umi, color = approach),
             linetype = "dotted", linewidth = 0.9, show.legend = FALSE) +
  # Crossing labels — staggered vertically
  geom_label(data = crossings_finite %>%
               mutate(y_pos = y_pos_vec[seq_len(nrow(.))]),
             aes(x = cross_umi, y = y_pos,
                 label = paste0(round(cross_umi/1000, 1), "k"),
                 color = approach),
             size = 3, fontface = "bold", fill = "white",
             label.size = 0.3, show.legend = FALSE) +
  # Footnote for approaches already ≥80% at 1k
  { if (any(crossings$starts_above, na.rm = TRUE)) {
    above_labels <- crossings %>%
      filter(starts_above) %>%
      pull(approach) %>%
      gsub("Approach [0-9]+: ", "", .) %>%
      paste(collapse = ", ")
    annotate("text", x = 1000, y = 5,
             label = paste0("≥80% at 1k UMIs: ", above_labels),
             size = 2.8, color = "grey40", hjust = 0, fontface = "italic")
  }
  } +
  scale_color_manual(values = cols_app) +
  scale_x_continuous(breaks = umi_x,
                     labels = paste0(umi_x/1000, "k"),
                     limits = c(500, 21500)) +
  scale_y_continuous(limits = c(0, 108)) +
  labs(title    = "Summary: Recommended UMI cutoff across all approaches",
       subtitle = paste0(
         "Approaches 3 & 4 computed among cells with >0 Lenti UMIs only\n",
         "Dashed = 80% target | Dotted verticals = 80% crossings | ",
         if (!is.na(conservative_cut))
           paste0("Conservative cutoff: ~", round(conservative_cut/1000, 1), "k UMIs")
         else "No single conservative cutoff within range"),
       x = "Total UMIs per cell (minimum cutoff)",
       y = "% / Probability", color = NULL) +
  theme_classic() +
  theme(plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(size = 9, color = "grey40"),
        legend.position = "bottom",
        legend.text     = element_text(size = 9),
        axis.text.x     = element_text(angle = 45, hjust = 1, size = 9))

ggsave("detection_approach5_summary.pdf", plot = p5, width = 14, height = 8)
cat("Approach 5 saved\n")

# Save per-cell scores
df_scores <- df_all_ci %>%
  left_join(bind_rows(df_hsc_post  %>% select(cell, prob_lenti_pos),
                      df_mono_post %>% select(cell, prob_lenti_pos)),
            by = "cell")
write.csv(df_scores, "lenti_detection_per_cell_scores.csv", row.names = FALSE)

cat("\nAll detection analyses complete!\n")

