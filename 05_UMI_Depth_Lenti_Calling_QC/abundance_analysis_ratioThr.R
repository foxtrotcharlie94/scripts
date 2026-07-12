# ============================================================
# Abundance analysis — LentiPos vs LentiNeg defined by ratio:
#   LentiPos: LentiAll / total UMIs >= 5e-5
#   LentiNeg: LentiAll / total UMIs <= 1e-5
#   (cells in between are ambiguous and dropped)
# Cluster labels from harmony_clusters_res0.5 UMAP; clusters
# 9, 10, 16 not annotated -> shown as their numbers.
# ============================================================

setwd("C:/Users/fc809/Downloads/")
tryCatch(detach("package:ensembldb", unload = TRUE), error = function(e) NULL)

library(tidyverse)
library(data.table)

# ── Thresholds ─────────────────────────────────────────────────────────────────
THR_POS <- 5e-5
THR_NEG <- 1e-5

# ── Output folder ──────────────────────────────────────────────────────────────
out_dir <- "C:/Users/fc809/Downloads/abundance_analysis_ratioThr/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Cluster labels from UMAP (res 0.5) ────────────────────────────────────────
cluster_labels <- c(
  "0"  = "Classical Ly6Chigh\nmonocytes",
  "1"  = "MoP",
  "2"  = "LT-HSCs + MPP1",
  "3"  = "MHCII+\nmonocytes",
  "4"  = "GMPs",
  "5"  = "Non-classical\nLy6Clow monocytes",
  "6"  = "MPP2\n(Mk/Ery primed)",
  "7"  = "MoPs\nS",
  "8"  = "MoPs\nG2/M",
  "11" = "MDPs (S)",
  "12" = "Promyelocyte",
  "13" = "Mk-biased HSCs",
  "14" = "Myelocyte",
  "15" = "Metamyelocyte",
  "17" = "IFNg-stimulated\nmonocytes",
  "18" = "PreDC",
  "19" = "Cycling PreDC",
  "20" = "MEP",
  "21" = "Erythroblast",
  "22" = "Mast cells"
)
# Clusters 9, 10, 16 not annotated -> fall through to their number
cluster_labels_flat <- gsub("\n", " ", cluster_labels)

# Clusters excluded from all plots
exclude_clusters <- c("9", "10", "20", "21", "22")

# ── Load per-cell ratios and classify ──────────────────────────────────────────
ratio_csv <- "C:/Users/fc809/Downloads/lenti_ratio_per_cluster/per_cell_lenti_ratio.csv"

cells <- fread(ratio_csv) |>
  as_tibble() |>
  mutate(
    sample    = str_extract(cell, "^[HL]B[0-9]+"),
    condition = str_extract(sample, "^[HL]B"),
    cluster   = as.character(cluster),
    genotype  = case_when(
      ratio >= THR_POS ~ "LentiPos",
      ratio <= THR_NEG ~ "LentiNeg",
      TRUE             ~ NA_character_
    )
  ) |>
  filter(!is.na(sample))

cat("Total cells:", nrow(cells), "\n")
cat("Genotype classification (incl. ambiguous as NA):\n")
print(table(cells$genotype, useNA = "ifany"))
cat("\nCells per sample x genotype (classified only):\n")
print(table(cells$sample, cells$genotype, useNA = "no"))

cells_class <- cells |> filter(!is.na(genotype))

# ── Per (sample, genotype) cluster proportions ─────────────────────────────────
suppressWarnings({
  uniq <- unique(cells_class$cluster); num <- as.numeric(uniq)
  all_clusters <- if (any(is.na(num))) sort(uniq) else uniq[order(num)]
})

# Label lookup: full label, fallback to cluster number
label_lookup <- setNames(
  ifelse(all_clusters %in% names(cluster_labels),
         cluster_labels[all_clusters], all_clusters),
  all_clusters)
label_lookup_flat <- setNames(
  ifelse(all_clusters %in% names(cluster_labels_flat),
         cluster_labels_flat[all_clusters], all_clusters),
  all_clusters)

props_df <- cells_class |>
  count(sample, condition, genotype, cluster) |>
  complete(nesting(sample, condition, genotype),
           cluster = all_clusters, fill = list(n = 0)) |>
  group_by(sample, genotype) |>
  mutate(proportion = n / sum(n)) |>
  ungroup() |>
  mutate(
    donor   = sample,
    cluster = factor(cluster, levels = all_clusters)
  ) |>
  filter(!as.character(cluster) %in% exclude_clusters) |>
  mutate(cluster = droplevels(cluster))

fwrite(props_df, file.path(out_dir, "props_per_sample_genotype_cluster.csv"))

cat("\nSamples x genotypes present:\n")
print(props_df |> distinct(sample, genotype) |> arrange(sample, genotype))

# ── Shared plot builder ────────────────────────────────────────────────────────
build_plot <- function(mean_df, indiv_df, pval_df, fill_var, colors,
                       title, subtitle) {
  
  plot_df <- mean_df %>% left_join(pval_df, by = "cluster")
  
  star_df <- plot_df %>%
    group_by(cluster, sig_label) %>%
    summarise(y_pos = max(mean_prop), .groups = "drop") %>%
    dplyr::filter(sig_label != "")
  
  y_nudge <- max(indiv_df$proportion) * 0.12
  
  ggplot(plot_df, aes(x = cluster, y = mean_prop, fill = .data[[fill_var]])) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7, color = NA) +
    geom_point(
      data = indiv_df,
      aes(x = cluster, y = proportion,
          fill = .data[[fill_var]], group = .data[[fill_var]]),
      position = position_dodge(width = 0.75),
      shape = 21, size = 2.8, color = "white", stroke = 0.4, alpha = 0.85,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = star_df,
      aes(x = cluster, y = y_pos + y_nudge, label = sig_label),
      inherit.aes = FALSE, size = 10, fontface = "bold", color = "grey20"
    ) +
    scale_x_discrete(labels = function(x) label_lookup[x]) +
    scale_fill_manual(values = colors, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                       labels = scales::label_number(accuracy = 0.01)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Mean Proportion") +
    theme_bw(base_size = 18) +
    theme(
      plot.title         = element_text(face = "bold", size = 20, hjust = 0.5),
      plot.subtitle      = element_text(size = 15, hjust = 0.5, color = "grey40",
                                        margin = margin(b = 8)),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 13, lineheight = 0.85),
      axis.text.y        = element_text(size = 16),
      axis.title.y       = element_text(size = 17),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = c(0.95, 0.88),
      legend.background  = element_rect(color = "grey80", linewidth = 0.4),
      legend.text        = element_text(size = 16),
      legend.key.size    = unit(1.1, "lines")
    )
}

# ── FDR helpers (NA-pair safe for paired) ──────────────────────────────────────
compute_fdr <- function(data, group_var, ref_level) {
  clusters <- levels(droplevels(data$cluster))
  pvals <- map_dbl(clusters, function(cl) {
    d <- data %>% dplyr::filter(cluster == cl)
    g1 <- d %>% dplyr::filter(.data[[group_var]] == ref_level) %>% pull(proportion)
    g2 <- d %>% dplyr::filter(.data[[group_var]] != ref_level) %>% pull(proportion)
    if (length(g1) < 2 | length(g2) < 2) return(NA_real_)
    tryCatch(t.test(g1, g2, paired = FALSE)$p.value,
             error = function(e) NA_real_)
  })
  tibble(cluster = factor(clusters, levels = clusters), pval = pvals) %>%
    mutate(
      fdr       = p.adjust(pval, method = "BH"),
      sig_label = case_when(
        is.na(fdr)  ~ "",
        fdr < 0.001 ~ "***",
        fdr < 0.01  ~ "**",
        fdr < 0.05  ~ "*",
        TRUE        ~ ""
      )
    )
}

compute_fdr_paired <- function(data, group_var, level1, level2) {
  clusters <- levels(droplevels(data$cluster))
  pvals <- map_dbl(clusters, function(cl) {
    d <- data %>% dplyr::filter(cluster == cl) %>%
      dplyr::select(donor, all_of(group_var), proportion) %>%
      pivot_wider(names_from = all_of(group_var), values_from = proportion) %>%
      drop_na(all_of(c(level1, level2)))
    if (nrow(d) < 3) return(NA_real_)
    tryCatch(t.test(d[[level1]], d[[level2]], paired = TRUE)$p.value,
             error = function(e) NA_real_)
  })
  tibble(cluster = factor(clusters, levels = clusters), pval = pvals) %>%
    mutate(
      fdr       = p.adjust(pval, method = "BH"),
      sig_label = case_when(
        is.na(fdr)  ~ "",
        fdr < 0.001 ~ "***",
        fdr < 0.01  ~ "**",
        fdr < 0.05  ~ "*",
        TRUE        ~ ""
      )
    )
}

sub_geno     <- sprintf("LentiPos (ratio>=%g) vs LentiNeg (ratio<=%g) | Paired t-test, BH-adj  (* p<0.05, ** p<0.01, *** p<0.001)",
                        THR_POS, THR_NEG)
sub_cond_pos <- sprintf("HB vs LB (LentiPos: ratio>=%g) | Unpaired t-test, BH-adj  (* p<0.05, ** p<0.01, *** p<0.001)", THR_POS)
sub_cond_neg <- sprintf("HB vs LB (LentiNeg: ratio<=%g) | Unpaired t-test, BH-adj  (* p<0.05, ** p<0.01, *** p<0.001)", THR_NEG)

# ── Plot 1: HB — LentiPos vs LentiNeg (paired) ────────────────────────────────
hb_data <- props_df %>% dplyr::filter(condition == "HB")
pval_hb <- compute_fdr_paired(hb_data, "genotype", "LentiPos", "LentiNeg")
mean_hb <- hb_data %>% group_by(cluster, genotype) %>%
  summarise(mean_prop = mean(proportion), .groups = "drop")

p_hb <- build_plot(
  mean_df = mean_hb, indiv_df = hb_data, pval_df = pval_hb,
  fill_var = "genotype",
  colors   = c("LentiNeg" = "#C5C6C7", "LentiPos" = "#4472C4"),
  title    = "HB Samples: Mean Cell Proportions per Cluster",
  subtitle = sub_geno
)

# ── Plot 2: LB — LentiPos vs LentiNeg (paired) ────────────────────────────────
lb_data <- props_df %>% dplyr::filter(condition == "LB")
pval_lb <- compute_fdr_paired(lb_data, "genotype", "LentiPos", "LentiNeg")
mean_lb <- lb_data %>% group_by(cluster, genotype) %>%
  summarise(mean_prop = mean(proportion), .groups = "drop")

p_lb <- build_plot(
  mean_df = mean_lb, indiv_df = lb_data, pval_df = pval_lb,
  fill_var = "genotype",
  colors   = c("LentiNeg" = "#C5C6C7", "LentiPos" = "#4472C4"),
  title    = "LB Samples: Mean Cell Proportions per Cluster",
  subtitle = sub_geno
)

# ── Plot 3: LentiPos — HB vs LB (unpaired) ────────────────────────────────────
lpos_data <- props_df %>% dplyr::filter(genotype == "LentiPos")
pval_lpos <- compute_fdr(lpos_data, "condition", "HB")
mean_lpos <- lpos_data %>% group_by(cluster, condition) %>%
  summarise(mean_prop = mean(proportion), .groups = "drop")

p_lpos <- build_plot(
  mean_df = mean_lpos, indiv_df = lpos_data, pval_df = pval_lpos,
  fill_var = "condition",
  colors   = c("HB" = "#2ECC71", "LB" = "#9B59B6"),
  title    = "LentiPos: Mean Cell Proportions per Cluster",
  subtitle = sub_cond_pos
)

# ── Plot 4: LentiNeg — HB vs LB (unpaired) ────────────────────────────────────
lneg_data <- props_df %>% dplyr::filter(genotype == "LentiNeg")
pval_lneg <- compute_fdr(lneg_data, "condition", "HB")
mean_lneg <- lneg_data %>% group_by(cluster, condition) %>%
  summarise(mean_prop = mean(proportion), .groups = "drop")

p_lneg <- build_plot(
  mean_df = mean_lneg, indiv_df = lneg_data, pval_df = pval_lneg,
  fill_var = "condition",
  colors   = c("HB" = "#2ECC71", "LB" = "#9B59B6"),
  title    = "LentiNeg: Mean Cell Proportions per Cluster",
  subtitle = sub_cond_neg
)

# ── Save plots ─────────────────────────────────────────────────────────────────
ggsave(file.path(out_dir, "HB_cell_proportions_LentiPos_vs_LentiNeg.pdf"), p_hb,   width = 14, height = 7)
ggsave(file.path(out_dir, "HB_cell_proportions_LentiPos_vs_LentiNeg.png"), p_hb,   width = 14, height = 7, dpi = 160)
ggsave(file.path(out_dir, "LB_cell_proportions_LentiPos_vs_LentiNeg.pdf"), p_lb,   width = 14, height = 7)
ggsave(file.path(out_dir, "LB_cell_proportions_LentiPos_vs_LentiNeg.png"), p_lb,   width = 14, height = 7, dpi = 160)
ggsave(file.path(out_dir, "LentiPos_cell_proportions_HB_vs_LB.pdf"),       p_lpos, width = 14, height = 7)
ggsave(file.path(out_dir, "LentiPos_cell_proportions_HB_vs_LB.png"),       p_lpos, width = 14, height = 7, dpi = 160)
ggsave(file.path(out_dir, "LentiNeg_cell_proportions_HB_vs_LB.pdf"),       p_lneg, width = 14, height = 7)
ggsave(file.path(out_dir, "LentiNeg_cell_proportions_HB_vs_LB.png"),       p_lneg, width = 14, height = 7, dpi = 160)

message("Saved: all four cell proportion plots")

# ══════════════════════════════════════════════════════════════════════════════
# Lenti coverage in clusters (recomputed from per-cell CSV, ALL cells)
# ══════════════════════════════════════════════════════════════════════════════
df_cov <- cells |>
  filter(!cluster %in% exclude_clusters) |>
  group_by(cluster) |>
  summarise(
    total_cells = n(),
    lenti_umis  = sum(lenti),
    .groups = "drop"
  ) |>
  arrange(suppressWarnings(as.numeric(cluster))) |>
  mutate(label = label_lookup_flat[cluster])

fwrite(df_cov, file.path(out_dir, "lenti_coverage_per_cluster.csv"))

cor_test <- cor.test(df_cov$total_cells, df_cov$lenti_umis, method = "pearson")
r_val   <- round(cor_test$estimate, 3)
p_val   <- signif(cor_test$p.value, 2)
p_label <- ifelse(p_val < 0.001,
                  paste0("p = ", formatC(p_val, format = "e", digits = 1)),
                  paste0("p = ", p_val))
annot   <- paste0("r = ", r_val, "\n", p_label)

fit      <- lm(lenti_umis ~ 0 + total_cells, data = df_cov)
pred_df  <- tibble(total_cells = seq(0, max(df_cov$total_cells) * 1.05, length.out = 200))
pred_out <- predict(fit, newdata = pred_df, interval = "confidence")
pred_df  <- bind_cols(pred_df, as_tibble(pred_out))

p_corr <- ggplot(df_cov, aes(x = total_cells, y = lenti_umis)) +
  geom_ribbon(data = pred_df, aes(x = total_cells, ymin = lwr, ymax = upr),
              inherit.aes = FALSE, fill = "grey90") +
  geom_line(data = pred_df, aes(x = total_cells, y = fit),
            inherit.aes = FALSE, color = "grey60", linewidth = 0.8) +
  geom_point(size = 4.5, color = "#2C3E50") +
  ggrepel::geom_text_repel(aes(label = label), size = 5, color = "grey30",
                           box.padding = 0.4, max.overlaps = 30,
                           segment.color = "grey70", segment.size = 0.3) +
  annotate("text", x = Inf, y = -Inf, label = annot,
           hjust = 1.1, vjust = -0.8, size = 7,
           fontface = "bold", color = "#2C3E50") +
  scale_x_continuous(labels = scales::comma, limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(labels = scales::comma, limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Total LentiAll UMIs per cluster",
    x     = "Total cells per cluster",
    y     = "Total LentiAll UMIs per cluster"
  ) +
  theme_bw(base_size = 18) +
  theme(
    plot.title  = element_text(face = "bold", size = 20),
    axis.text   = element_text(size = 16),
    axis.title  = element_text(size = 17),
    panel.grid  = element_line(color = "grey93")
  )

ggsave(file.path(out_dir, "lenti_umi_correlation.pdf"), p_corr, width = 11, height = 7.5)
ggsave(file.path(out_dir, "lenti_umi_correlation.png"), p_corr, width = 11, height = 7.5, dpi = 160)

message("Saved: lenti_umi_correlation.pdf / .png")
message("\nAll outputs in: ", out_dir)

