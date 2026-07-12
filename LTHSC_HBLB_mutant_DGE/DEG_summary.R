# ══════════════════════════════════════════════════════════════════════════════
#  DEG Comparison: All samples vs excluding HB3/HB4
#  Both inputs are limma log-CPM results (CSV)
#  Plots per comparison:
#    1. Bar: DEG counts at FDR < 0.2 side by side
#    2. Scatter: logFC concordance, coloured by significance
#    3. Scatter: FDR concordance (-log10)
#  One combined PDF per dataset
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(patchwork)
library(ggrepel)

setwd("C:/Users/fc809/Downloads/LT-HSCs (1)")

datasets <- list(
  list(
    name     = "LT-HSCs",
    all_dir  = "allSamples_analysis",
    excl_dir = "excluding_H3_H4",
    out      = "summary_plots"
  ),
  list(
    name     = "Classical Monocytes",
    all_dir  = "../Classical_Monocytes (1)/allSamples_analysis",
    excl_dir = "../Classical_Monocytes (1)/excluding_H3_H4",
    out      = "../Classical_Monocytes (1)/summary_plots"
  )
)

comparisons <- list(
  list(folder = "LentiPos_vs_LentiNeg_in_HB", label = "LentiPos vs LentiNeg (HB)"),
  list(folder = "LentiPos_vs_LentiNeg_in_LB", label = "LentiPos vs LentiNeg (LB)"),
  list(folder = "LentiPosHB_vs_LentiPosLB",   label = "LentiPos HB vs LB"),
  list(folder = "LentiNegHB_vs_LentiNegLB",   label = "LentiNeg HB vs LB")
)

FDR_THRESH <- 0.2

# ── Helper: find limma CSV in a given base/folder ─────────────────────────────
find_limma <- function(base, folder) {
  p     <- file.path(base, folder)
  files <- list.files(p, pattern = "_limma_.*results\\.csv$", full.names = TRUE)
  if (length(files) > 0) return(files[1])
  NULL
}

# ── Helper: load and standardise a limma CSV ──────────────────────────────────
load_limma <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)
  if (!"adj.P.Val" %in% names(df) && "FDR"  %in% names(df)) df <- dplyr::rename(df, adj.P.Val = "FDR")
  if (!"adj.P.Val" %in% names(df) && "padj" %in% names(df)) df <- dplyr::rename(df, adj.P.Val = "padj")
  if (!"gene"      %in% names(df))                           df <- dplyr::rename(df, gene = 1)
  df %>% dplyr::select(all_of(c("gene", "logFC", "adj.P.Val")))
}

# ══════════════════════════════════════════════════════════════════════════════
for (ds in datasets) {
  dir.create(ds$out, showWarnings = FALSE, recursive = TRUE)
  cat("\n══", ds$name, "══\n")
  
  all_plots  <- list()
  count_rows <- list()
  
  for (comp in comparisons) {
    
    f_all  <- find_limma(ds$all_dir,  comp$folder)
    f_excl <- find_limma(ds$excl_dir, comp$folder)
    
    if (is.null(f_all))  { cat("  [skip] all-samples result not found:", comp$folder, "\n"); next }
    if (is.null(f_excl)) { cat("  [skip] excl-HB3/HB4 result not found:", comp$folder, "\n"); next }
    
    cat("  Reading:", comp$label, "\n")
    
    all_res  <- load_limma(f_all)  %>% dplyr::rename(logFC_all  = "logFC", fdr_all  = "adj.P.Val")
    excl_res <- load_limma(f_excl) %>% dplyr::rename(logFC_excl = "logFC", fdr_excl = "adj.P.Val")
    
    # ── Merge ─────────────────────────────────────────────────────────────────
    merged <- inner_join(all_res, excl_res, by = "gene") %>%
      mutate(
        sig = case_when(
          fdr_all  < FDR_THRESH & fdr_excl < FDR_THRESH ~ "Sig in both",
          fdr_excl < FDR_THRESH                          ~ "Sig excl. HB3/HB4 only",
          fdr_all  < FDR_THRESH                          ~ "Sig all samples only",
          TRUE                                           ~ "NS"
        ),
        sig = factor(sig, levels = c("Sig in both", "Sig excl. HB3/HB4 only",
                                     "Sig all samples only", "NS"))
      )
    
    # ── Count table ───────────────────────────────────────────────────────────
    count_rows[[comp$folder]] <- tibble(
      comparison = comp$label,
      source     = c("All samples", "Excl. HB3/HB4"),
      n_up       = c(sum(all_res$fdr_all   < FDR_THRESH & all_res$logFC_all   > 0, na.rm = TRUE),
                     sum(excl_res$fdr_excl < FDR_THRESH & excl_res$logFC_excl > 0, na.rm = TRUE)),
      n_down     = c(sum(all_res$fdr_all   < FDR_THRESH & all_res$logFC_all   < 0, na.rm = TRUE),
                     sum(excl_res$fdr_excl < FDR_THRESH & excl_res$logFC_excl < 0, na.rm = TRUE))
    )
    
    r_lfc <- round(cor(merged$logFC_all, merged$logFC_excl, use = "complete.obs"), 3)
    
    sig_colors <- c(
      "Sig in both"             = "#E74C3C",
      "Sig excl. HB3/HB4 only" = "#9B59B6",
      "Sig all samples only"    = "#F39C12",
      "NS"                      = "grey80"
    )
    
    top_discord <- merged %>%
      filter(sig != "NS") %>%
      mutate(diff = abs(logFC_all - logFC_excl)) %>%
      slice_max(diff, n = 8)
    
    # Plot A: logFC scatter
    pA <- ggplot(merged, aes(logFC_all, logFC_excl, color = sig)) +
      geom_point(alpha = 0.3, size = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "grey30", linewidth = 0.6) +
      geom_text_repel(data = top_discord, aes(label = gene),
                      size = 3.5, max.overlaps = 8,
                      segment.size = 0.3, segment.color = "grey60") +
      annotate("text", x = -Inf, y = Inf,
               label = paste0("r = ", r_lfc),
               hjust = -0.15, vjust = 1.4, size = 4.5,
               fontface = "bold", color = "black") +
      scale_color_manual(values = sig_colors, name = NULL) +
      labs(title = comp$label,
           x = "logFC — all samples",
           y = "logFC — excl. HB3/HB4") +
      theme_bw(base_size = 12) +
      theme(plot.title      = element_text(face = "bold", size = 12, color = "black"),
            axis.text       = element_text(size = 11, color = "black"),
            axis.title      = element_text(size = 11, color = "black"),
            legend.position = "bottom",
            legend.text     = element_text(size = 9))
    
    # Plot B: -log10 FDR scatter
    pB <- ggplot(merged, aes(-log10(fdr_all + 1e-10), -log10(fdr_excl + 1e-10), color = sig)) +
      geom_point(alpha = 0.3, size = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "grey30", linewidth = 0.6) +
      geom_hline(yintercept = -log10(FDR_THRESH), linetype = "dotted",
                 color = "grey50", linewidth = 0.5) +
      geom_vline(xintercept = -log10(FDR_THRESH), linetype = "dotted",
                 color = "grey50", linewidth = 0.5) +
      scale_color_manual(values = sig_colors, name = NULL) +
      labs(title = comp$label,
           x = "-log10(FDR) — all samples",
           y = "-log10(FDR) — excl. HB3/HB4") +
      theme_bw(base_size = 12) +
      theme(plot.title      = element_text(face = "bold", size = 12, color = "black"),
            axis.text       = element_text(size = 11, color = "black"),
            axis.title      = element_text(size = 11, color = "black"),
            legend.position = "bottom",
            legend.text     = element_text(size = 9))
    
    all_plots[[comp$folder]] <- list(lfc = pA, fdr = pB)
  }
  
  # ── Bar plot: DEG counts side by side ─────────────────────────────────────
  if (length(count_rows) > 0) {
    count_df <- bind_rows(count_rows) %>%
      pivot_longer(c(n_up, n_down), names_to = "direction", values_to = "n") %>%
      mutate(
        direction  = recode(direction, n_up = "Up", n_down = "Down"),
        direction  = factor(direction, levels = c("Up", "Down")),
        comparison = factor(comparison, levels = map_chr(comparisons, "label")),
        source     = factor(source, levels = c("All samples", "Excl. HB3/HB4"))
      )
    
    p_bar <- ggplot(count_df, aes(x = source, y = n, fill = direction)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
      geom_text(aes(label = n),
                position = position_dodge(width = 0.7),
                vjust = -0.4, size = 3.5, fontface = "bold", color = "black") +
      facet_wrap(~ comparison, nrow = 1) +
      scale_fill_manual(values = c(Up = "#E74C3C", Down = "#3498DB"), name = "Direction") +
      labs(
        title    = paste0(ds$name, " — DEG Counts: All samples vs Excl. HB3/HB4 (FDR < ", FDR_THRESH, ")"),
        subtitle = "Left = all samples  |  Right = excluding HB3/HB4",
        x = NULL, y = "Number of DEGs"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title         = element_text(face = "bold", size = 13, color = "black"),
        plot.subtitle      = element_text(size = 10, color = "grey30"),
        axis.text.x        = element_text(size = 9,  color = "black", lineheight = 1.2),
        axis.text.y        = element_text(size = 10, color = "black"),
        axis.title.y       = element_text(size = 11, color = "black"),
        strip.text         = element_text(size = 10, face = "bold"),
        strip.background   = element_rect(fill = "grey92", color = NA),
        panel.grid.major.x = element_blank(),
        legend.text        = element_text(size = 10)
      )
    
    ggsave(file.path(ds$out, "DEG_comparison_counts.pdf"), p_bar, width = 14, height = 6)
    ggsave(file.path(ds$out, "DEG_comparison_counts.png"), p_bar, width = 14, height = 6, dpi = 160)
    cat("  Saved: DEG_comparison_counts\n")
  }
  
  # ── logFC and FDR scatter grids ───────────────────────────────────────────
  if (length(all_plots) > 0) {
    lfc_panels <- map(all_plots, "lfc")
    fdr_panels <- map(all_plots, "fdr")
    n <- length(lfc_panels)
    
    p_lfc_grid <- wrap_plots(lfc_panels, nrow = 1) +
      plot_annotation(
        title = paste0(ds$name, " — logFC Concordance: All samples vs Excl. HB3/HB4"),
        theme = theme(plot.title = element_text(face = "bold", size = 14, color = "black"))
      )
    
    p_fdr_grid <- wrap_plots(fdr_panels, nrow = 1) +
      plot_annotation(
        title = paste0(ds$name, " — FDR Concordance: All samples vs Excl. HB3/HB4"),
        theme = theme(plot.title = element_text(face = "bold", size = 14, color = "black"))
      )
    
    ggsave(file.path(ds$out, "DEG_comparison_logFC.pdf"), p_lfc_grid, width = 5.5 * n, height = 6)
    ggsave(file.path(ds$out, "DEG_comparison_logFC.png"), p_lfc_grid, width = 5.5 * n, height = 6, dpi = 150)
    ggsave(file.path(ds$out, "DEG_comparison_FDR.pdf"),   p_fdr_grid, width = 5.5 * n, height = 6)
    ggsave(file.path(ds$out, "DEG_comparison_FDR.png"),   p_fdr_grid, width = 5.5 * n, height = 6, dpi = 150)
    cat("  Saved: DEG_comparison_logFC + DEG_comparison_FDR\n")
  }
}

cat("\n✓ Done.\n")

