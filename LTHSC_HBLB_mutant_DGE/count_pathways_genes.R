# ══════════════════════════════════════════════════════════════════════════════
#  Summary Bar Plots
#    1. DEGs (FDR < 0.1): excl. HB3/HB4 vs all samples
#    2. Significant pathways (FDR < 0.05) per comparison x database
#  For both LT-HSCs and Classical Monocytes
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(patchwork)

setwd("C:/Users/fc809/Downloads/LT-HSCs (1)")

# ── Comparisons ───────────────────────────────────────────────────────────────
comparisons <- list(
  list(folder = "LentiPos_vs_LentiNeg_in_HB", label = "LentiPos vs LentiNeg\n(HB)"),
  list(folder = "LentiPos_vs_LentiNeg_in_LB", label = "LentiPos vs LentiNeg\n(LB)"),
  list(folder = "LentiPosHB_vs_LentiPosLB",   label = "LentiPos\nHB vs LB"),
  list(folder = "LentiNegHB_vs_LentiNegLB",   label = "LentiNeg\nHB vs LB")
)
comp_labels <- map_chr(comparisons, "label")

# ── Databases ─────────────────────────────────────────────────────────────────
db_tags <- c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Hallmark")
db_colors <- c(
  GO_BP    = "#3498DB",
  GO_MF    = "#2ECC71",
  GO_CC    = "#1ABC9C",
  KEGG     = "#E74C3C",
  Hallmark = "#9B59B6"
)

# ── Datasets ──────────────────────────────────────────────────────────────────
datasets <- list(
  list(name = "LT-HSCs",            base = ".",                          out = "summary_plots"),
  list(name = "Classical Monocytes", base = "../Classical_Monocytes (1)", out = "../Classical_Monocytes (1)/summary_plots")
)

# ── Helpers ───────────────────────────────────────────────────────────────────
find_file <- function(base, folder, pattern) {
  for (sub in c("results", "")) {
    p     <- if (sub == "") file.path(base, folder) else file.path(base, folder, sub)
    files <- list.files(p, pattern = pattern, full.names = TRUE)
    if (length(files) > 0) return(files[1])
  }
  NULL
}

fix_deg_cols <- function(df) {
  if (!"adj.P.Val" %in% names(df) && "padj"           %in% names(df)) df <- rename(df, adj.P.Val = padj)
  if (!"adj.P.Val" %in% names(df) && "FDR"            %in% names(df)) df <- rename(df, adj.P.Val = FDR)
  if (!"logFC"     %in% names(df) && "log2FoldChange" %in% names(df)) df <- rename(df, logFC = log2FoldChange)
  df
}

# ══════════════════════════════════════════════════════════════════════════════
for (ds in datasets) {
  dir.create(ds$out, showWarnings = FALSE, recursive = TRUE)
  cat("\n══", ds$name, "══\n")
  
  # ── 1. DEG counts: excl HB3/HB4 and all samples ───────────────────────────
  analyses <- list(
    list(label = "Excl. HB3/HB4", base = ds$base),
    list(label = "All samples",   base = file.path(ds$base, "allSamples_analysis"))
  )
  
  deg_rows <- list()
  for (an in analyses) {
    for (comp in comparisons) {
      f <- find_file(an$base, comp$folder, "_limma_voom_results\\.csv$")
      if (is.null(f)) { cat("  [DEG] not found:", an$label, "|", comp$folder, "\n"); next }
      res    <- fix_deg_cols(read_csv(f, show_col_types = FALSE))
      n_up   <- sum(res$adj.P.Val < 0.1 & res$logFC >  0, na.rm = TRUE)
      n_down <- sum(res$adj.P.Val < 0.1 & res$logFC <  0, na.rm = TRUE)
      deg_rows[[paste(an$label, comp$folder)]] <- tibble(
        analysis   = an$label,
        comparison = comp$label,
        direction  = c("Up", "Down"),
        n_genes    = c(n_up, n_down)
      )
      cat("  [DEG]", an$label, "|", comp$folder, "— up:", n_up, "down:", n_down, "\n")
    }
  }
  
  if (length(deg_rows) == 0) {
    cat("  No DEG results found — skipping DEG plot\n")
  } else {
    deg_df <- bind_rows(deg_rows) %>%
      mutate(
        analysis   = factor(analysis,   levels = c("Excl. HB3/HB4", "All samples")),
        comparison = factor(comparison, levels = comp_labels),
        direction  = factor(direction,  levels = c("Up", "Down"))
      )
    
    p_deg <- ggplot(deg_df, aes(x = comparison, y = n_genes, fill = direction)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.65, alpha = 0.85) +
      geom_text(aes(label = n_genes),
                position = position_dodge(width = 0.7),
                vjust = -0.4, size = 3.5, fontface = "bold", color = "black") +
      facet_wrap(~ analysis, nrow = 2) +
      scale_fill_manual(values = c(Up = "#E74C3C", Down = "#3498DB"), name = "Direction") +
      labs(
        title    = paste0(ds$name, " — Differentially Expressed Genes (FDR < 0.1)"),
        subtitle = "Top: excluding HB3/HB4  |  Bottom: all samples",
        x = NULL, y = "Number of DEGs"
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title         = element_text(face = "bold", size = 14, color = "black"),
        plot.subtitle      = element_text(size = 11, color = "grey30"),
        axis.text.x        = element_text(size = 11, color = "black", lineheight = 1.2),
        axis.text.y        = element_text(size = 11, color = "black"),
        axis.title.y       = element_text(size = 12, color = "black"),
        strip.text         = element_text(size = 11, face = "bold"),
        strip.background   = element_rect(fill = "grey92", color = NA),
        legend.text        = element_text(size = 11),
        legend.title       = element_text(size = 11, face = "bold"),
        panel.grid.major.x = element_blank()
      )
    
    ggsave(file.path(ds$out, "DEG_counts_FDR0.1.pdf"), p_deg, width = 11, height = 9)
    ggsave(file.path(ds$out, "DEG_counts_FDR0.1.png"), p_deg, width = 11, height = 9, dpi = 160)
    cat("  Saved: DEG_counts_FDR0.1\n")
  }
  
  # ── 2. Pathway counts (excl HB3/HB4 only) ────────────────────────────────
  pathway_rows <- list()
  for (comp in comparisons) {
    for (db in db_tags) {
      f <- find_file(ds$base, comp$folder, paste0("_GSEA_", db, "\\.csv$"))
      if (is.null(f)) {
        pathway_rows[[paste(comp$folder, db)]] <- tibble(
          comparison = comp$label, database = db,
          n_sig = 0L, n_up = 0L, n_down = 0L
        )
        next
      }
      res <- read_csv(f, show_col_types = FALSE)
      if (!"p.adjust" %in% names(res) && "qvalue" %in% names(res))
        res <- rename(res, p.adjust = qvalue)
      pathway_rows[[paste(comp$folder, db)]] <- tibble(
        comparison = comp$label,
        database   = db,
        n_sig      = sum(res$p.adjust < 0.05, na.rm = TRUE),
        n_up       = sum(res$p.adjust < 0.05 & res$NES > 0, na.rm = TRUE),
        n_down     = sum(res$p.adjust < 0.05 & res$NES < 0, na.rm = TRUE)
      )
      cat("  [GSEA]", comp$folder, db, "— sig:", tail(pathway_rows,1)[[1]]$n_sig, "\n")
    }
  }
  
  pathway_df <- bind_rows(pathway_rows) %>%
    mutate(
      comparison = factor(comparison, levels = comp_labels),
      database   = factor(database,   levels = db_tags)
    )
  
  p_stack <- ggplot(pathway_df, aes(x = comparison, y = n_sig, fill = database)) +
    geom_col(width = 0.65, alpha = 0.85) +
    geom_text(aes(label = ifelse(n_sig > 0, n_sig, "")),
              position = position_stack(vjust = 0.5),
              size = 3.2, fontface = "bold", color = "white") +
    scale_fill_manual(values = db_colors, name = "Database") +
    labs(
      title    = paste0(ds$name, " — Significant Pathways (FDR < 0.05)"),
      subtitle = "Excluding HB3/HB4  |  Stacked by database",
      x = NULL, y = "Number of significant pathways"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title         = element_text(face = "bold", size = 14, color = "black"),
      plot.subtitle      = element_text(size = 11, color = "grey30"),
      axis.text.x        = element_text(size = 11, color = "black", lineheight = 1.2),
      axis.text.y        = element_text(size = 11, color = "black"),
      axis.title.y       = element_text(size = 12, color = "black"),
      legend.text        = element_text(size = 11),
      legend.title       = element_text(size = 11, face = "bold"),
      panel.grid.major.x = element_blank()
    )
  
  pathway_long <- pathway_df %>%
    dplyr::select(comparison, database, n_up, n_down) %>%
    pivot_longer(c(n_up, n_down), names_to = "direction", values_to = "n") %>%
    mutate(
      direction = recode(direction, n_up = "Up (NES > 0)", n_down = "Down (NES < 0)"),
      direction = factor(direction, levels = c("Up (NES > 0)", "Down (NES < 0)"))
    )
  
  p_facet <- ggplot(pathway_long, aes(x = comparison, y = n, fill = direction)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65, alpha = 0.85) +
    geom_text(aes(label = ifelse(n > 0, n, "")),
              position = position_dodge(width = 0.7),
              vjust = -0.4, size = 3.2, fontface = "bold", color = "black") +
    facet_wrap(~ database, scales = "free_y", nrow = 1) +
    scale_fill_manual(
      values = c("Up (NES > 0)" = "#E74C3C", "Down (NES < 0)" = "#3498DB"),
      name = "Direction"
    ) +
    labs(
      title    = paste0(ds$name, " — Significant Pathways by Database (FDR < 0.05)"),
      subtitle = "Excluding HB3/HB4  |  Split by enrichment direction",
      x = NULL, y = "Number of significant pathways"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 14, color = "black"),
      plot.subtitle      = element_text(size = 11, color = "grey30"),
      axis.text.x        = element_text(size = 9, color = "black", angle = 35,
                                        hjust = 1, lineheight = 1.1),
      axis.text.y        = element_text(size = 10, color = "black"),
      axis.title.y       = element_text(size = 11, color = "black"),
      strip.text         = element_text(size = 11, face = "bold"),
      strip.background   = element_rect(fill = "grey92", color = NA),
      legend.text        = element_text(size = 10),
      legend.title       = element_text(size = 10, face = "bold"),
      panel.grid.major.x = element_blank()
    )
  
  ggsave(file.path(ds$out, "pathway_counts_stacked.pdf"),     p_stack, width = 11, height = 6)
  ggsave(file.path(ds$out, "pathway_counts_stacked.png"),     p_stack, width = 11, height = 6, dpi = 160)
  ggsave(file.path(ds$out, "pathway_counts_by_database.pdf"), p_facet, width = 16, height = 6)
  ggsave(file.path(ds$out, "pathway_counts_by_database.png"), p_facet, width = 16, height = 6, dpi = 160)
  cat("  Saved: pathway_counts_stacked + pathway_counts_by_database\n")
}

cat("\n✓ Done.\n")

