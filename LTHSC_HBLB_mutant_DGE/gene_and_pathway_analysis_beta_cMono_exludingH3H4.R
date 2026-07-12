# ══════════════════════════════════════════════════════════════════════════════
#  limma DGE (log-CPM input) + Pathway Analysis — Five comparisons — Classical Monocytes
# HB3 and HB4 EXCLUDED
#  Pathways (GSEA): GO (BP/MF/CC), KEGG, Hallmark (MSigDB)
#  Each pathway plot: top 10 up + top 10 down combined, with counts annotated
# ══════════════════════════════════════════════════════════════════════════════
#
#  Required packages:
#  install.packages(c("tidyverse", "BiocManager", "msigdbr"))
#  BiocManager::install(c("limma", "clusterProfiler", "org.Mm.eg.db"))
#
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(limma)
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(ggrepel)

base_dir     <- "C:/Users/fc809/Downloads/Classical_Monocytes (1)"
results_base <- "C:/Users/fc809/Downloads/Classical_Monocytes (1)/excluding_H3_H4"

# ── Define comparisons ────────────────────────────────────────────────────────
comparisons <- list(
  list(
    folder          = "LentiPos_vs_LentiNeg_in_HB",
    sample_patterns = list(group1 = c("LentiPos_HB", "Pos_HB", "PosHB"),
                           group2 = c("LentiNeg_HB", "Neg_HB", "NegHB")),
    label1 = "LentiPos", label2 = "LentiNeg",
    title  = "LentiPos vs LentiNeg (HB)"
  ),
  list(
    folder          = "LentiPos_vs_LentiNeg_in_LB",
    sample_patterns = list(group1 = c("LentiPos_LB", "Pos_LB", "PosLB"),
                           group2 = c("LentiNeg_LB", "Neg_LB", "NegLB")),
    label1 = "LentiPos", label2 = "LentiNeg",
    title  = "LentiPos vs LentiNeg (LB)"
  ),
  list(
    folder          = "LentiPosHB_vs_LentiPosLB",
    sample_patterns = list(group1 = c("LentiPos_HB", "Pos_HB", "PosHB"),
                           group2 = c("LentiPos_LB", "Pos_LB", "PosLB")),
    label1 = "HB", label2 = "LB",
    title  = "LentiPos HB vs LentiPos LB"
  ),
  list(
    folder          = "LentiPosLB_vs_LentiNegHB",
    sample_patterns = list(group1 = c("LentiPos_LB", "Pos_LB", "PosLB"),
                           group2 = c("LentiNeg_HB", "Neg_HB", "NegHB")),
    label1 = "LentiPos_LB", label2 = "LentiNeg_HB",
    title  = "LentiPos LB vs LentiNeg HB"
  ),
  list(
    folder          = "LentiNegHB_vs_LentiNegLB",
    sample_patterns = list(group1 = c("LentiNeg_HB", "Neg_HB", "NegHB"),
                           group2 = c("LentiNeg_LB", "Neg_LB", "NegLB")),
    label1 = "LentiNeg_HB", label2 = "LentiNeg_LB",
    title  = "LentiNeg HB vs LentiNeg LB"
  )
)

# ── Helper: find DEG file ─────────────────────────────────────────────────────
find_deg_file <- function(folder_path) {
  files <- list.files(folder_path, pattern = "DEG_", full.names = TRUE)
  if (length(files) == 0) stop("No DEG_ file found in: ", folder_path)
  cat("  DEG file:", basename(files[1]), "\n")
  files[1]
}

# ── Helper: assign group labels ───────────────────────────────────────────────
assign_group <- function(col_names, patterns, label1, label2) {
  group <- rep(NA_character_, length(col_names))
  for (pat in patterns$group1)
    group[grepl(pat, col_names, fixed = TRUE)] <- label1
  for (pat in patterns$group2)
    group[grepl(pat, col_names, fixed = TRUE)] <- label2
  group
}

# ── Helper: run limma on log-CPM (for CPM input — no re-normalization) ────────
run_limma_cpm <- function(cpm, group, label1, label2) {
  group  <- factor(group, levels = c(label1, label2))
  logcpm <- log2(cpm + 0.5)                    # log2 transform; +0.5 avoids log(0)
  keep   <- rowMeans(logcpm) > 0.5             # keep genes with mean log2CPM > 0.5
  logcpm <- logcpm[keep, ]
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  fit    <- lmFit(logcpm, design)
  con    <- makeContrasts(contrasts = paste0(label1, "-", label2),
                          levels = design)
  fit    <- contrasts.fit(fit, con)
  fit    <- eBayes(fit, robust = TRUE)
  topTable(fit, n = Inf, sort.by = "P") %>%
    rownames_to_column("gene") %>%
    as_tibble()
}

# ── Helper: volcano plot ───────────────────────────────────────────────────────
make_volcano <- function(res, title, label1, label2,
                         fdr_cut = 0.05, lfc_cut = 1) {
  df <- res %>%
    mutate(
      sig = case_when(
        adj.P.Val < fdr_cut & logFC >  lfc_cut ~ "Up",
        adj.P.Val < fdr_cut & logFC < -lfc_cut ~ "Down",
        TRUE ~ "NS"
      ),
      neg_log10_fdr = -log10(adj.P.Val + 1e-300)
    )
  
  top_labels <- df %>% filter(sig != "NS") %>%
    arrange(adj.P.Val) %>% slice_head(n = 20)
  
  ggplot(df, aes(logFC, neg_log10_fdr, color = sig)) +
    geom_point(alpha = 0.5, size = 1.2) +
    geom_text_repel(data = top_labels, aes(label = gene),
                    size = 2.8, max.overlaps = 20,
                    segment.size = 0.3, segment.color = "grey60") +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut),
               linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_hline(yintercept = -log10(fdr_cut),
               linetype = "dashed", color = "grey50", linewidth = 0.4) +
    scale_color_manual(
      values = c("Up" = "#E74C3C", "Down" = "#3498DB", "NS" = "grey70"),
      name = NULL
    ) +
    labs(
      title    = title,
      subtitle = paste0(
        sum(df$sig == "Up"),   " up in ", label1, "  |  ",
        sum(df$sig == "Down"), " up in ", label2,
        "  (FDR<", fdr_cut, ", |logFC|>", lfc_cut, ")"
      ),
      x = paste0("log2 FC  (", label1, " / ", label2, ")"),
      y = "-log10(FDR)"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title      = element_text(face = "bold", size = 14),
      plot.subtitle   = element_text(size = 12, color = "black"),
      axis.text       = element_text(size = 11),
      axis.title      = element_text(size = 12),
      legend.position = "top",
      legend.text     = element_text(size = 11),
      plot.margin     = margin(t = 10, r = 20, b = 10, l = 10)
    )
}

# ── Helper: GSEA dot plot (top 10 positive + top 10 negative NES) ────────────
# Keywords to exclude from pathway plots (non-haematopoietic tissue terms)
tissue_exclusion_pattern <- paste(
  "brain", "neuron", "neural", "cerebr", "cortex", "hippocamp",
  "liver", "hepat",
  "kidney", "renal",
  "ovar", "uterus", "uterine", "placenta",
  "retina", "optic",
  "cardiac", "cardio", "heart",
  "pulmonary", "lung", "bronch", "alveol",
  "osteoclast",
  "measles",
  "alcoholism",
  "prion",
  "parkinson",
  "sclerosis",
  "neurodegeneration", "neurodegenerative",
  "huntington",
  "alzheimer",
  sep = "|"
)

make_gsea_dotplot <- function(gsea_res, db_name, comparison_title,
                              label1, label2, n = 10) {
  
  df <- gsea_res %>% filter(p.adjust < 0.05)
  if (nrow(df) == 0) return(NULL)
  
  n_up   <- sum(df$NES > 0)
  n_down <- sum(df$NES < 0)
  
  # Filter tissue-irrelevant terms only when > n significant pathways per direction
  filter_tissue <- function(sub_df) {
    if (nrow(sub_df) > n) {
      filtered <- sub_df %>%
        filter(!grepl(tissue_exclusion_pattern, Description, ignore.case = TRUE))
      # Fall back to unfiltered if filtering removes everything
      if (nrow(filtered) > 0) return(filtered)
    }
    sub_df
  }
  
  top_up <- df %>% filter(NES > 0) %>%
    arrange(p.adjust) %>%
    filter_tissue() %>%
    slice_head(n = n) %>%
    mutate(direction = paste0("Enriched in ", label1))
  
  top_down <- df %>% filter(NES < 0) %>%
    arrange(p.adjust) %>%
    filter_tissue() %>%
    slice_head(n = n) %>%
    mutate(direction = paste0("Enriched in ", label2))
  
  combined <- bind_rows(top_up, top_down)
  if (nrow(combined) == 0) return(NULL)
  
  combined <- combined %>%
    mutate(Description = str_wrap(Description, 45)) %>%
    arrange(desc(NES)) %>%
    mutate(Description = factor(Description, levels = unique(Description)))
  
  dir_colors <- setNames(c("#E74C3C", "#3498DB"),
                         c(paste0("Enriched in ", label1),
                           paste0("Enriched in ", label2)))
  
  subtitle <- paste0(
    comparison_title, "\n",
    n_up,   " pathways enriched in ", label1, "  |  ",
    n_down, " pathways enriched in ", label2,
    "  (FDR<0.05, showing top ", n, " each)"
  )
  
  ggplot(combined, aes(x = NES, y = Description,
                       size = setSize, color = direction)) +
    geom_point() +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    scale_color_manual(values = dir_colors, name = NULL) +
    scale_size_continuous(name = "Gene set size", range = c(2, 8)) +
    labs(
      title    = paste0(db_name, " — GSEA"),
      subtitle = subtitle,
      x = "Normalised Enrichment Score (NES)", y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", size = 13),
      plot.subtitle   = element_text(size = 13, color = "black", lineheight = 1.3, face = "bold"),
      axis.text       = element_text(size = 13, color = "black"),
      axis.text.y     = element_text(size = 13, color = "black"),
      axis.title      = element_text(size = 12),
      legend.position = "top",
      legend.text     = element_text(size = 11),
      plot.margin     = margin(t = 10, r = 20, b = 10, l = 10)
    )
}

# ── Helper: symbol to Entrez ──────────────────────────────────────────────────
to_entrez_df <- function(res) {
  # Returns a named numeric vector of logFC, keyed by Entrez ID
  # Excludes Lenti construct and unmapped genes
  res <- res %>% filter(!gene %in% c("Lenti", "lenti"))
  mapped <- suppressMessages(
    bitr(res$gene, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Mm.eg.db, drop = TRUE)
  )
  res %>%
    inner_join(mapped, by = c("gene" = "SYMBOL")) %>%
    arrange(desc(logFC)) %>%
    { setNames(.$logFC, .$ENTREZID) }
}

# ── Helper: save plot ─────────────────────────────────────────────────────────
save_plot <- function(p, path, w = 10, h = 7) {
  if (is.null(p)) return(invisible(NULL))
  ggsave(paste0(path, ".pdf"), p, width = w, height = h)
  ggsave(paste0(path, ".png"), p, width = w, height = h, dpi = 160)
}

# ── Helper: run GSEA for one database ────────────────────────────────────────
run_gsea_db <- function(gene_list, type, ont = NULL, hallmark_df = NULL) {
  if (type == "GO") {
    suppressMessages(
      gseGO(geneList     = gene_list,
            OrgDb        = org.Mm.eg.db,
            ont          = ont,
            keyType      = "ENTREZID",
            minGSSize    = 15,
            maxGSSize    = 500,
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            eps          = 0,
            nPermSimple  = 1000,
            verbose      = FALSE)
    )
  } else if (type == "KEGG") {
    suppressMessages(
      gseKEGG(geneList      = gene_list,
              organism      = "mmu",
              keyType       = "ncbi-geneid",
              minGSSize     = 15,
              maxGSSize     = 500,
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              eps           = 0,
              nPermSimple   = 1000,
              verbose       = FALSE)
    )
  } else if (type == "Hallmark") {
    # TERM2GENE must have cols: (term_name, entrez_id) matching geneList names
    term2gene <- hallmark_df %>%
      dplyr::rename(term = gs_name, gene = ncbi_gene) %>%
      dplyr::mutate(gene = as.character(gene))
    gene_list_chr <- gene_list
    names(gene_list_chr) <- as.character(names(gene_list))
    suppressMessages(
      GSEA(geneList      = gene_list_chr,
           TERM2GENE     = term2gene,
           minGSSize     = 15,
           maxGSSize     = 500,
           pvalueCutoff  = 0.05,
           pAdjustMethod = "BH",
           eps           = 0,
           nPermSimple   = 1000,
           verbose       = FALSE)
    )
  }
}

# ── Helper: full GSEA pathway analysis ───────────────────────────────────────
run_pathway <- function(res, safe_label, out_dir, label1, label2,
                        comparison_title) {
  
  # Ranked gene list: logFC, named by Entrez ID, sorted descending
  gene_list <- to_entrez_df(res)
  cat("  Ranked gene list length:", length(gene_list), "\n")
  
  hallmark_df <- msigdbr(species = "Mus musculus") %>%
    filter(gs_collection == "H") %>%
    dplyr::select(gs_name, ncbi_gene)
  
  dbs <- list(
    list(name = "GO BP",    type = "GO",      ont = "BP"),
    list(name = "GO MF",    type = "GO",      ont = "MF"),
    list(name = "GO CC",    type = "GO",      ont = "CC"),
    list(name = "KEGG",     type = "KEGG",    ont = NULL),
    list(name = "Hallmark", type = "Hallmark",ont = NULL)
  )
  
  for (db in dbs) {
    cat("  Running GSEA:", db$name, "\n")
    gsea_res <- tryCatch(
      run_gsea_db(gene_list, db$type, db$ont, hallmark_df),
      error = function(e) { cat("    ERROR:", conditionMessage(e), "\n"); NULL }
    )
    
    if (is.null(gsea_res)) next
    
    df <- gsea_res@result
    n_sig <- sum(df$p.adjust < 0.05)
    cat("    Significant terms:", n_sig, "\n")
    if (n_sig == 0) next
    
    db_tag <- gsub(" ", "_", db$name)
    
    # Save full results CSV
    write_csv(df, file.path(out_dir, paste0(safe_label, "_GSEA_", db_tag, ".csv")))
    
    # Combined dot plot
    p <- make_gsea_dotplot(df, db$name, comparison_title, label1, label2, n = 10)
    n_up   <- sum(df$p.adjust < 0.05 & df$NES > 0)
    n_down <- sum(df$p.adjust < 0.05 & df$NES < 0)
    n_shown  <- min(n_up, 10) + min(n_down, 10)
    h <- max(8, n_shown * 0.55 + 4)
    save_plot(p, file.path(out_dir, paste0(safe_label, "_GSEA_", db_tag)),
              w = 14, h = h)
  }
}

# ══════════════════════════════════════════════════════════════════════════════
#  Main loop
# ══════════════════════════════════════════════════════════════════════════════
for (comp in comparisons) {
  
  cat("\n══════════════════════════════════════════\n")
  cat("Comparison:", comp$title, "\n")
  
  folder_path <- file.path(base_dir, comp$folder)
  out_dir     <- file.path(results_base, comp$folder)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  deg_file <- find_deg_file(folder_path)
  raw      <- read.delim(deg_file, row.names = 1, check.names = FALSE)
  
  result_cols <- c("logFC","logCPM","F","LR","PValue","FDR",
                   "logCPMrank","unshrunk.logFC","AveExpr","t","B","P.Value")
  sample_cols <- colnames(raw)[
    !colnames(raw) %in% result_cols & sapply(raw, is.numeric)
  ]
  
  group <- assign_group(sample_cols, comp$sample_patterns,
                        comp$label1, comp$label2)
  
  keep        <- !is.na(group) & !grepl("HB3|HB4", sample_cols)  # exclude HB3 and HB4
  sample_cols <- sample_cols[keep]
  group       <- group[keep]
  
  if (length(sample_cols) < 4) {
    warning("Too few samples (", length(sample_cols), ") for: ",
            comp$title, " — skipping")
    next
  }
  
  cat("Samples:\n")
  walk2(sample_cols, group, ~ cat("  ", .x, "→", .y, "\n"))
  
  cpm        <- raw[!rownames(raw) %in% c("Lenti", "lenti"), sample_cols, drop = FALSE]
  safe_label <- gsub("[^A-Za-z0-9_]", "_", comp$title)
  
  res <- run_limma_cpm(cpm, group, comp$label1, comp$label2)
  
  res <- res %>% filter(!gene %in% c("Lenti", "lenti"))  # exclude construct
  n_sig <- sum(res$adj.P.Val < 0.05 & abs(res$logFC) > 1)
  cat("Significant DEGs (adj.P.Val<0.05, |logFC|>1):", n_sig, "\n")
  
  write_csv(res, file.path(out_dir, paste0(safe_label, "_limma_voom_results.csv")))
  
  p_vol <- make_volcano(res, comp$title, comp$label1, comp$label2)
  save_plot(p_vol, file.path(out_dir, paste0(safe_label, "_volcano")), w = 8, h = 6)
  
  cat("Running pathway analysis...\n")
  tryCatch(
    run_pathway(res, safe_label, out_dir,
                comp$label1, comp$label2, comp$title),
    error = function(e) cat("  ERROR in pathway analysis:", conditionMessage(e), "\n")
  )
  
  cat("✓ Done. Outputs in:", out_dir, "\n")
}

cat("\n✓ All comparisons complete.\n")

