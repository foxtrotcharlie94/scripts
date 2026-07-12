setwd("C:/Users/fc809/Downloads/")
tryCatch(detach("package:ensembldb", unload = TRUE), error = function(e) NULL)

library(tidyverse)

# ── Cluster labels ─────────────────────────────────────────────────────────────
cluster_labels <- c(
  "0"  = "Classical Monocytes\nLy6Chigh",
  "1"  = "Promyelocytes (GrP)\nS phase",
  "2"  = "Transitioning Monocytes\nLy6Cmed",
  "3"  = "G1 LT-HSCs",
  "4"  = "MPP",
  "5"  = "NC Monocytes\nLy6Clow",
  "6"  = "cMoPs\nS phase",
  "7"  = "moDCs",
  "8"  = "cMoPs\nG2/M",
  "9"  = "MDPs",
  "10" = "GMPs\nS phase",
  "11" = "GMPs\nG2/M phase",
  "12" = "MDPs\n(cycling)",
  "13" = "Promyelocytes (GrP)\nG2/M",
  "14" = "MPP2",
  "15" = "Late Mitotic\nPromyelocytes (GrP)",
  "16" = "Pre-Neu",
  "17" = "Inflammatory\nMonocytes",
  "18" = "Neutrophils",
  "19" = "Gr-GMPs",
  "20" = "Pro-pDCs",
  "21" = "G0 LT-HSCs",
  "22" = "Pro-pDCs\ncycling",
  "23" = "CFU-E",
  "24" = "Erythroblasts /\nBasophils",
  "25" = "Mast Cell P"
)

exclude_clusters <- c("3", "4", "14", "21", "23", "16", "18", "24", "25")
kept_labels      <- cluster_labels[!names(cluster_labels) %in% exclude_clusters]

# ── Load and reshape ───────────────────────────────────────────────────────────
raw <- read.delim(
  "number_of_cells_per_cluster_for_each_sampleLenti_res0.6.txt",
  row.names = 1, check.names = FALSE
)

props <- raw / rowSums(raw)

props_df <- props %>%
  rownames_to_column("sample_id") %>%
  mutate(
    condition = ifelse(grepl("HB", sample_id), "HB", "LB"),
    genotype  = ifelse(grepl("LentiPos", sample_id), "LentiPos", "LentiNeg"),
    donor     = str_extract(sample_id, "^[HL]B[0-9]+")
  ) %>%
  pivot_longer(cols = -c(sample_id, condition, genotype, donor),
               names_to = "cluster", values_to = "proportion") %>%
  mutate(cluster = factor(cluster, levels = as.character(0:25))) %>%
  dplyr::filter(!cluster %in% exclude_clusters)

# ── Shared plot builder ────────────────────────────────────────────────────────
build_plot <- function(mean_df, indiv_df, pval_df, fill_var, colors,
                       title, subtitle) {
  
  plot_df <- mean_df %>%
    left_join(pval_df, by = "cluster") %>%
    mutate(label = factor(cluster_labels[as.character(cluster)], levels = kept_labels))
  
  indiv_df <- indiv_df %>%
    mutate(label = factor(cluster_labels[as.character(cluster)], levels = kept_labels))
  
  star_df <- plot_df %>%
    group_by(cluster, label, sig_label) %>%
    summarise(y_pos = max(mean_prop), .groups = "drop") %>%
    dplyr::filter(sig_label != "")
  
  y_nudge <- max(indiv_df$proportion) * 0.12
  
  ggplot(plot_df, aes(x = label, y = mean_prop, fill = .data[[fill_var]])) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7, color = NA) +
    geom_point(
      data = indiv_df,
      aes(x = label, y = proportion,
          fill = .data[[fill_var]], group = .data[[fill_var]]),
      position = position_dodge(width = 0.75),
      shape = 21, size = 2.8, color = "white", stroke = 0.4, alpha = 0.85,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = star_df,
      aes(x = label, y = y_pos + y_nudge, label = sig_label),
      inherit.aes = FALSE, size = 10, fontface = "bold", color = "grey20"
    ) +
    scale_fill_manual(values = colors, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                       labels = scales::label_number(accuracy = 0.01)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Mean Proportion") +
    theme_bw(base_size = 18) +
    theme(
      plot.title         = element_text(face = "bold", size = 20, hjust = 0.5),
      plot.subtitle      = element_text(size = 15, hjust = 0.5, color = "grey40",
                                        margin = margin(b = 8)),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 16, lineheight = 0.85),
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

# ── Helper: compute FDR ────────────────────────────────────────────────────────
compute_fdr <- function(data, group_var, ref_level) {
  clusters <- levels(droplevels(data$cluster))
  pvals <- map_dbl(clusters, function(cl) {
    d <- data %>% dplyr::filter(cluster == cl)
    g1 <- d %>% dplyr::filter(.data[[group_var]] == ref_level) %>% pull(proportion)
    g2 <- d %>% dplyr::filter(.data[[group_var]] != ref_level) %>% pull(proportion)
    if (length(g1) < 2 | length(g2) < 2) return(NA_real_)
    t.test(g1, g2, paired = FALSE)$p.value
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
      pivot_wider(names_from = all_of(group_var), values_from = proportion)
    if (nrow(d) < 3) return(NA_real_)
    t.test(d[[level1]], d[[level2]], paired = TRUE)$p.value
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

# ── Plot 1: HB — LentiPos vs LentiNeg (paired) ────────────────────────────────
hb_data <- props_df %>% dplyr::filter(condition == "HB")
pval_hb  <- compute_fdr_paired(hb_data, "genotype", "LentiPos", "LentiNeg")
mean_hb  <- hb_data %>% group_by(cluster, genotype) %>%
  summarise(mean_prop = mean(proportion), .groups = "drop")

p_hb <- build_plot(
  mean_df   = mean_hb,
  indiv_df  = hb_data,
  pval_df   = pval_hb,
  fill_var  = "genotype",
  colors    = c("LentiNeg" = "#C5C6C7", "LentiPos" = "#4472C4"),
  title     = "HB Samples: Mean Cell Proportions per Cluster",
  subtitle  = "LentiPos vs LentiNeg | Paired t-test, BH-adjusted (* FDR<0.05, ** FDR<0.01, *** FDR<0.001)"
)

# ── Plot 2: LB — LentiPos vs LentiNeg (paired) ────────────────────────────────
lb_data <- props_df %>% dplyr::filter(condition == "LB")
pval_lb  <- compute_fdr_paired(lb_data, "genotype", "LentiPos", "LentiNeg")
mean_lb  <- lb_data %>% group_by(cluster, genotype) %>%
  summarise(mean_prop = mean(proportion), .groups = "drop")

p_lb <- build_plot(
  mean_df   = mean_lb,
  indiv_df  = lb_data,
  pval_df   = pval_lb,
  fill_var  = "genotype",
  colors    = c("LentiNeg" = "#C5C6C7", "LentiPos" = "#4472C4"),
  title     = "LB Samples: Mean Cell Proportions per Cluster",
  subtitle  = "LentiPos vs LentiNeg | Paired t-test, BH-adjusted (* FDR<0.05, ** FDR<0.01, *** FDR<0.001)"
)

# ── Plot 3: LentiPos — HB vs LB (unpaired) ────────────────────────────────────
lpos_data <- props_df %>% dplyr::filter(genotype == "LentiPos")
pval_lpos  <- compute_fdr(lpos_data, "condition", "HB")
mean_lpos  <- lpos_data %>% group_by(cluster, condition) %>%
  summarise(mean_prop = mean(proportion), .groups = "drop")

p_lpos <- build_plot(
  mean_df   = mean_lpos,
  indiv_df  = lpos_data,
  pval_df   = pval_lpos,
  fill_var  = "condition",
  colors    = c("HB" = "#2ECC71", "LB" = "#9B59B6"),
  title     = "LentiPos: Mean Cell Proportions per Cluster",
  subtitle  = "HB vs LB | Unpaired t-test, BH-adjusted (* FDR<0.05, ** FDR<0.01, *** FDR<0.001)"
)

# ── Plot 4: LentiNeg — HB vs LB (unpaired) ────────────────────────────────────
lneg_data <- props_df %>% dplyr::filter(genotype == "LentiNeg")
pval_lneg  <- compute_fdr(lneg_data, "condition", "HB")
mean_lneg  <- lneg_data %>% group_by(cluster, condition) %>%
  summarise(mean_prop = mean(proportion), .groups = "drop")

p_lneg <- build_plot(
  mean_df   = mean_lneg,
  indiv_df  = lneg_data,
  pval_df   = pval_lneg,
  fill_var  = "condition",
  colors    = c("HB" = "#2ECC71", "LB" = "#9B59B6"),
  title     = "LentiNeg: Mean Cell Proportions per Cluster",
  subtitle  = "HB vs LB | Unpaired t-test, BH-adjusted (* FDR<0.05, ** FDR<0.01, *** FDR<0.001)"
)

# ── Save all four plots ────────────────────────────────────────────────────────
ggsave("HB_cell_proportions_LentiPos_vs_LentiNeg.pdf", p_hb,   width = 14, height = 6.5)
ggsave("HB_cell_proportions_LentiPos_vs_LentiNeg.png", p_hb,   width = 14, height = 6.5, dpi = 160)
ggsave("LB_cell_proportions_LentiPos_vs_LentiNeg.pdf", p_lb,   width = 14, height = 6.5)
ggsave("LB_cell_proportions_LentiPos_vs_LentiNeg.png", p_lb,   width = 14, height = 6.5, dpi = 160)
ggsave("LentiPos_cell_proportions_HB_vs_LB.pdf",       p_lpos, width = 14, height = 6.5)
ggsave("LentiPos_cell_proportions_HB_vs_LB.png",       p_lpos, width = 14, height = 6.5, dpi = 160)
ggsave("LentiNeg_cell_proportions_HB_vs_LB.pdf",       p_lneg, width = 14, height = 6.5)
ggsave("LentiNeg_cell_proportions_HB_vs_LB.png",       p_lneg, width = 14, height = 6.5, dpi = 160)

message("Saved: all four cell proportion plots")

# ══════════════════════════════════════════════════════════════════════════════
# Lenti coverage in clusters
# ══════════════════════════════════════════════════════════════════════════════

# ── Cluster labels (single-line for scatter) ──────────────────────────────────
cluster_labels_flat <- c(
  "0"  = "Classical Monocytes Ly6Chigh",
  "1"  = "Promyelocytes (GrP) S phase",
  "2"  = "Transitioning Monocytes Ly6Cmed",
  "3"  = "G1 LT-HSCs",
  "4"  = "MPP",
  "5"  = "NC Monocytes Ly6Clow",
  "6"  = "cMoPs S phase",
  "7"  = "moDCs",
  "8"  = "cMoPs G2/M",
  "9"  = "MDPs",
  "10" = "GMPs S phase",
  "11" = "GMPs G2/M phase",
  "12" = "MDPs (cycling)",
  "13" = "Promyelocytes (GrP) G2/M",
  "14" = "MPP2",
  "15" = "Late Mitotic Promyelocytes (GrP)",
  "16" = "Pre-Neu",
  "17" = "Inflammatory Monocytes",
  "18" = "Neutrophils",
  "19" = "Gr-GMPs",
  "20" = "Pro-pDCs",
  "21" = "G0 LT-HSCs",
  "22" = "Pro-pDCs cycling",
  "23" = "CFU-E",
  "24" = "Erythroblasts / Basophils",
  "25" = "Mast Cell P"
)

total_cells <- colSums(raw)

lenti_umis <- c(
  "0"  = 61925, "1"  = 49186, "2"  = 30256, "3"  = 38534,
  "4"  = 36043, "5"  = 6350,  "6"  = 22744, "7"  = 14909,
  "8"  = 9879,  "9"  = 989,   "10" = 44377, "11" = 50867,
  "12" = 976,   "13" = 14832, "14" = 1714,  "15" = 15452,
  "16" = 4533,  "17" = 7296,  "18" = 1480,  "19" = 39094,
  "20" = 9193,  "21" = 9572,  "22" = 7564,  "23" = 8437,
  "24" = 532,   "25" = 723
)

clusters_ord <- as.character(0:25)
df <- tibble(
  cluster     = clusters_ord,
  total_cells = as.numeric(total_cells[clusters_ord]),
  lenti_umis  = as.numeric(lenti_umis[clusters_ord]),
  label       = cluster_labels_flat[clusters_ord]
)

cor_test <- cor.test(df$total_cells, df$lenti_umis, method = "pearson")
r_val    <- round(cor_test$estimate, 3)
p_val    <- signif(cor_test$p.value, 2)
p_label  <- ifelse(p_val < 0.001,
                   paste0("p = ", formatC(p_val, format = "e", digits = 1)),
                   paste0("p = ", p_val))
annot    <- paste0("r = ", r_val, "\n", p_label)

fit      <- lm(lenti_umis ~ 0 + total_cells, data = df)
pred_df  <- tibble(total_cells = seq(0, max(df$total_cells) * 1.05, length.out = 200))
pred_out <- predict(fit, newdata = pred_df, interval = "confidence")
pred_df  <- bind_cols(pred_df, as_tibble(pred_out))

p_corr <- ggplot(df, aes(x = total_cells, y = lenti_umis)) +
  geom_ribbon(data = pred_df, aes(x = total_cells, ymin = lwr, ymax = upr),
              inherit.aes = FALSE, fill = "grey90") +
  geom_line(data = pred_df, aes(x = total_cells, y = fit),
            inherit.aes = FALSE, color = "grey60", linewidth = 0.8) +
  geom_point(size = 4.5, color = "#2C3E50") +
  ggrepel::geom_text_repel(aes(label = label), size = 5, color = "grey30",
                           box.padding = 0.4, max.overlaps = 20,
                           segment.color = "grey70", segment.size = 0.3) +
  annotate("text", x = Inf, y = -Inf, label = annot,
           hjust = 1.1, vjust = -0.8, size = 7,
           fontface = "bold", color = "#2C3E50") +
  scale_x_continuous(labels = scales::comma, limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(labels = scales::comma, limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Total Lenti UMIs per cluster",
    x     = "Total cells per cluster",
    y     = "Total Lenti UMIs per cluster"
  ) +
  theme_bw(base_size = 18) +
  theme(
    plot.title  = element_text(face = "bold", size = 20),
    axis.text   = element_text(size = 16),
    axis.title  = element_text(size = 17),
    panel.grid  = element_line(color = "grey93")
  )

ggsave("lenti_umi_correlation.pdf", p_corr, width = 10, height = 7)
ggsave("lenti_umi_correlation.png", p_corr, width = 10, height = 7, dpi = 160)

message("Saved: lenti_umi_correlation.pdf / .png")

