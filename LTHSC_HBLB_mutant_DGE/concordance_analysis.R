# ══════════════════════════════════════════════════════════════════════════════
#  Concordance Analysis: All samples vs Excl. HB3/HB4
#  Plots:
#    1. logFC scatter per comparison (gene-level concordance)
#    2. GSEA NES scatter per comparison x database (pathway-level concordance)
#  For both LT-HSCs and Classical Monocytes
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(ggrepel)
library(patchwork)

setwd("C:/Users/fc809/Downloads")

DEG_FDR     <- 0.2
PATHWAY_FDR <- 0.05

datasets <- list(
  list(
    name     = "LT-HSCs",
    all_dir  = "LT-HSCs (1)/allSamples_analysis",
    excl_dir = "LT-HSCs (1)/excluding_H3_H4",
    out      = "LT-HSCs (1)/concordance_analysis"
  ),
  list(
    name     = "Classical Monocytes",
    all_dir  = "Classical_Monocytes (1)/allSamples_analysis",
    excl_dir = "Classical_Monocytes (1)/excluding_H3_H4",
    out      = "Classical_Monocytes (1)/concordance_analysis"
  )
)

comparisons <- list(
  list(folder = "LentiPos_vs_LentiNeg_in_HB", title = "LentiPos vs LentiNeg (HB)"),
  list(folder = "LentiPos_vs_LentiNeg_in_LB", title = "LentiPos vs LentiNeg (LB)"),
  list(folder = "LentiPosHB_vs_LentiPosLB",   title = "LentiPos HB vs LB"),
  list(folder = "LentiNegHB_vs_LentiNegLB",   title = "LentiNeg HB vs LB")
)

db_tags <- c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Hallmark")

# ── Helpers ───────────────────────────────────────────────────────────────────
find_limma <- function(base, folder) {
  p     <- file.path(base, folder)
  files <- list.files(p, pattern = "_limma_.*results\\.csv$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  df <- read_csv(files[1], show_col_types = FALSE)
  if (!"adj.P.Val" %in% names(df) && "FDR"  %in% names(df)) df <- dplyr::rename(df, adj.P.Val = "FDR")
  if (!"adj.P.Val" %in% names(df) && "padj" %in% names(df)) df <- dplyr::rename(df, adj.P.Val = "padj")
  if (!"gene"      %in% names(df))                           df <- dplyr::rename(df, gene = 1)
  df
}

find_gsea <- function(base, folder, db_tag) {
  p     <- file.path(base, folder)
  files <- list.files(p, pattern = paste0("_GSEA_", db_tag, "\\.csv$"), full.names = TRUE)
  if (length(files) == 0) return(NULL)
  read_csv(files[1], show_col_types = FALSE)
}

sig_colors <- c(
  "Sig in both"   = "#E74C3C",
  "Sig excl only" = "#9B59B6",
  "Sig all only"  = "#F39C12",
  "NS"            = "grey75"
)

# ══════════════════════════════════════════════════════════════════════════════
for (ds in datasets) {
  dir.create(ds$out, showWarnings = FALSE, recursive = TRUE)
  cat("\n══", ds$name, "══\n")
  
  # ── Plot 1: logFC concordance ───────────────────────────────────────────────
  lfc_plots <- list()
  
  for (comp in comparisons) {
    res_all  <- find_limma(ds$all_dir,  comp$folder)
    res_excl <- find_limma(ds$excl_dir, comp$folder)
    
    if (is.null(res_all) || is.null(res_excl)) {
      cat("  [skip logFC]", comp$title, "\n"); next
    }
    
    merged <- inner_join(
      res_all  %>% dplyr::select(gene, logFC, adj.P.Val) %>%
        dplyr::rename(logFC_all  = "logFC", fdr_all  = "adj.P.Val"),
      res_excl %>% dplyr::select(gene, logFC, adj.P.Val) %>%
        dplyr::rename(logFC_excl = "logFC", fdr_excl = "adj.P.Val"),
      by = "gene"
    ) %>%
      mutate(sig = case_when(
        fdr_excl < DEG_FDR & fdr_all < DEG_FDR ~ "Sig in both",
        fdr_excl < DEG_FDR                      ~ "Sig excl only",
        fdr_all  < DEG_FDR                      ~ "Sig all only",
        TRUE                                    ~ "NS"
      ))
    
    r  <- round(cor(merged$logFC_excl, merged$logFC_all, method = "pearson"),  3)
    rs <- round(cor(merged$logFC_excl, merged$logFC_all, method = "spearman"), 3)
    
    top_discord <- merged %>%
      filter(sig != "NS") %>%
      mutate(diff = abs(logFC_excl - logFC_all)) %>%
      slice_max(diff, n = 8)
    
    p <- ggplot(merged, aes(logFC_all, logFC_excl, color = sig)) +
      geom_point(alpha = 0.35, size = 0.9) +
      geom_abline(slope = 1, intercept = 0, color = "grey30",
                  linetype = "dashed", linewidth = 0.6) +
      geom_text_repel(data = top_discord, aes(label = gene),
                      size = 4, max.overlaps = 10,
                      segment.size = 0.3, segment.color = "grey60") +
      annotate("text", x = -Inf, y = Inf,
               label = paste0("r = ", r, "\nρ = ", rs),
               hjust = -0.1, vjust = 1.3, size = 4.5,
               fontface = "bold", color = "black") +
      scale_color_manual(values = sig_colors, name = NULL) +
      labs(
        title    = comp$title,
        subtitle = paste0(
          sum(merged$sig == "Sig in both"),   " sig in both  |  ",
          sum(merged$sig == "Sig excl only"), " excl-only  |  ",
          sum(merged$sig == "Sig all only"),  " all-only",
          "  (FDR < ", DEG_FDR, ")"
        ),
        x = "logFC — all samples",
        y = "logFC — excl. HB3/HB4"
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title      = element_text(face = "bold", size = 13, color = "black"),
        plot.subtitle   = element_text(size = 10, color = "black"),
        axis.text       = element_text(size = 12, color = "black"),
        axis.title      = element_text(size = 13, color = "black"),
        legend.position = "bottom",
        legend.text     = element_text(size = 11)
      )
    
    lfc_plots[[comp$folder]] <- p
    cat("  logFC:", comp$title, "— r =", r, "\n")
  }
  
  if (length(lfc_plots) > 0) {
    n    <- length(lfc_plots)
    ncol <- min(4, n)
    nrow <- ceiling(n / ncol)
    combined_lfc <- wrap_plots(lfc_plots, ncol = ncol) +
      plot_annotation(
        title = paste0(ds$name, " — Gene-level logFC Concordance: All samples vs Excl. HB3/HB4"),
        theme = theme(plot.title = element_text(face = "bold", size = 15, color = "black"))
      )
    ggsave(file.path(ds$out, "concordance_logFC.pdf"),
           combined_lfc, width = 6.5 * ncol, height = 5.5 * nrow)
    ggsave(file.path(ds$out, "concordance_logFC.png"),
           combined_lfc, width = 6.5 * ncol, height = 5.5 * nrow, dpi = 160)
    cat("  Saved: concordance_logFC\n")
  }
  
  # ── Plot 2: GSEA NES concordance ───────────────────────────────────────────
  for (db_tag in db_tags) {
    nes_plots <- list()
    
    for (comp in comparisons) {
      gsea_all  <- find_gsea(ds$all_dir,  comp$folder, db_tag)
      gsea_excl <- find_gsea(ds$excl_dir, comp$folder, db_tag)
      
      if (is.null(gsea_all) || is.null(gsea_excl)) next
      
      all_tbl  <- gsea_all  %>% dplyr::select(ID, Description, NES, p.adjust) %>%
        dplyr::rename(NES_all  = "NES", fdr_all  = "p.adjust")
      excl_tbl <- gsea_excl %>% dplyr::select(ID, Description, NES, p.adjust) %>%
        dplyr::rename(NES_excl = "NES", fdr_excl = "p.adjust")
      
      merged <- full_join(all_tbl, excl_tbl, by = "ID") %>%
        mutate(
          Description = coalesce(Description.x, Description.y),
          NES_all  = replace_na(NES_all,  0),
          NES_excl = replace_na(NES_excl, 0),
          fdr_all  = replace_na(fdr_all,  1),
          fdr_excl = replace_na(fdr_excl, 1),
          sig = case_when(
            fdr_excl < PATHWAY_FDR & fdr_all < PATHWAY_FDR ~ "Sig in both",
            fdr_excl < PATHWAY_FDR                         ~ "Sig excl only",
            fdr_all  < PATHWAY_FDR                         ~ "Sig all only",
            TRUE                                           ~ "NS"
          )
        ) %>%
        filter(sig != "NS")
      
      if (nrow(merged) < 1) next
      
      r <- round(cor(merged$NES_excl, merged$NES_all, method = "pearson"), 3)
      
      top_discord <- merged %>%
        mutate(diff = abs(NES_excl - NES_all)) %>%
        slice_max(diff, n = 6)
      
      p <- ggplot(merged, aes(NES_all, NES_excl, color = sig)) +
        geom_point(alpha = 0.6, size = 2) +
        geom_abline(slope = 1, intercept = 0, color = "grey30",
                    linetype = "dashed", linewidth = 0.6) +
        geom_text_repel(data = top_discord,
                        aes(label = str_wrap(Description, 30)),
                        size = 3.5, max.overlaps = 8,
                        segment.size = 0.3, segment.color = "grey60") +
        annotate("text", x = -Inf, y = Inf,
                 label = paste0("r = ", r),
                 hjust = -0.15, vjust = 1.4, size = 4.5,
                 fontface = "bold", color = "black") +
        scale_color_manual(values = sig_colors, name = NULL) +
        labs(
          title    = comp$title,
          subtitle = paste0(
            sum(merged$sig == "Sig in both"),   " sig in both  |  ",
            sum(merged$sig == "Sig excl only"), " excl-only  |  ",
            sum(merged$sig == "Sig all only"),  " all-only",
            "  (FDR < ", PATHWAY_FDR, ")"
          ),
          x = "NES — all samples",
          y = "NES — excl. HB3/HB4"
        ) +
        theme_bw(base_size = 13) +
        theme(
          plot.title      = element_text(face = "bold", size = 13, color = "black"),
          plot.subtitle   = element_text(size = 10, color = "black"),
          axis.text       = element_text(size = 12, color = "black"),
          axis.title      = element_text(size = 13, color = "black"),
          legend.position = "bottom",
          legend.text     = element_text(size = 11)
        )
      
      nes_plots[[comp$folder]] <- p
      cat("  NES:", db_tag, "|", comp$title, "— r =", r, "\n")
    }
    
    if (length(nes_plots) > 0) {
      n    <- length(nes_plots)
      ncol <- min(4, n)
      nrow <- ceiling(n / ncol)
      combined_nes <- wrap_plots(nes_plots, ncol = ncol) +
        plot_annotation(
          title = paste0(ds$name, " — GSEA NES Concordance (", db_tag, "): All samples vs Excl. HB3/HB4"),
          theme = theme(plot.title = element_text(face = "bold", size = 15, color = "black"))
        )
      ggsave(file.path(ds$out, paste0("concordance_NES_", db_tag, ".pdf")),
             combined_nes, width = 6.5 * ncol, height = 5.5 * nrow)
      ggsave(file.path(ds$out, paste0("concordance_NES_", db_tag, ".png")),
             combined_nes, width = 6.5 * ncol, height = 5.5 * nrow, dpi = 160)
      cat("  Saved: concordance_NES_", db_tag, "\n")
    }
  }
}

cat("\n✓ Done.\n")
