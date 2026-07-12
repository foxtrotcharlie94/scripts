## ============================================================
## Jakobsen HSPC datasets × LT-HSC comparisons — cross-dataset GSEA
##
## Six Jakobsen DE tables compared against your LT-HSC comparisons:
##   1. HSPC_TET2mut_vs_WT       (TotalA/B schema)
##   2. HSPC2_TET2mut_vs_WT      (TotalA/B schema)
##   3. CHWT_vs_nonCH            (TotalA/B schema)
##   4. Aksoz_Aged_vs_Young      (limma-style, has a title row)
##   5. Ainciburu_Aged_vs_Young  (DESeq2 schema)
##   6. Inhouse_Aged_vs_Young    (DESeq2 schema)
##
## Significance filter : padj/FDR < 0.05
## Ranking & direction : log2FC / log2FoldChange / logFC
## Many-to-one ortholog dedup : keep mouse symbol with largest |log2FC|
##
## OUTPUT FOLDER PER DATASET
##   <OUTDIR>/<dataset_name>/
##     GSEA_<comparison>.png          (5x — individual GSEA curves)
##     Summary_dotplot.png            (4 core comparisons)
##     LeadingEdge_heatmap_*.png
##     LeadingEdge_Venn_*.png
##     RankedList_scatter.png
##     orthologs.csv
##     gene_sets.csv
##     pathways/                      (optional — see RUN_PATHWAYS flag)
##
## OUTPUT GLOBAL
##   <OUTDIR>/Combined_dotplot.png    (all 6 datasets faceted)
##   <OUTDIR>/Combined_summary.csv    (long-format NES + FDR table)
##   <OUTDIR>/Jakobsen_HSPC_multi_analysis.R   (script copy)
## ============================================================

library(clusterProfiler)
library(enrichplot)
library(babelgene)
library(pheatmap)
library(ggVennDiagram)
library(org.Hs.eg.db)
library(msigdbr)
library(ggplot2)
# dplyr loaded LAST so its verbs (rename, select, filter, ...) take precedence
# over masks from S4Vectors / AnnotationDbi / clusterProfiler.
library(tidyr)
library(dplyr)

set.seed(42)

OUTDIR        <- "C:/Users/fc809/Downloads/Jakobsen_HSPC_multi_GSEA"
SIG_CUTOFF    <- 0.05         # padj/FDR cutoff
RUN_PATHWAYS  <- TRUE         # GO BP + KEGG + Hallmarks on the original HUMAN data

# Path to this script — update if you rename or move it.
SCRIPT_PATH   <- "C:/Users/fc809/Downloads/Jakobsen_HSPC_multi_analysis.R"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

## ── DATASET CONFIG ────────────────────────────────────────────────────────────
## set_a = genes with log2FC > 0 (higher in numerator of comparison)
## set_b = genes with log2FC < 0 (higher in denominator)

datasets <- list(
  HSPC_TET2mut_vs_WT = list(
    file     = "C:/Users/fc809/Downloads/Jakobsen_HSPC_TET2mut_vs_TET2WT.csv",
    skip     = 0,
    gene_col = "Gene", fc_col = "log2FC",         fdr_col = "FDR",
    set_a    = "MUT_up", set_b = "WT_up",
    label    = "TET2-MUT vs TET2-WT (HSPC1)"
  ),
  HSPC2_TET2mut_vs_WT = list(
    file     = "C:/Users/fc809/Downloads/Jakobsen_HSPC2_TET2mut_vs_TET2WT.csv",
    skip     = 0,
    gene_col = "Gene", fc_col = "log2FC",         fdr_col = "FDR",
    set_a    = "MUT_up", set_b = "WT_up",
    label    = "TET2-MUT vs TET2-WT (HSPC2)"
  ),
  CHWT_vs_nonCH = list(
    file     = "C:/Users/fc809/Downloads/Jakobsen_HSPC_CH-WT_vs_non-CH.csv",
    skip     = 0,
    gene_col = "Gene", fc_col = "log2FC",         fdr_col = "FDR",
    set_a    = "CHwt_up", set_b = "nonCH_up",
    label    = "CH-WT vs non-CH (HSPC)"
  ),
  Aksoz_Aged_vs_Young = list(
    file     = "C:/Users/fc809/Downloads/Jakobsen_HSPC_Aksöz_Aged_vs_Young_donors.csv",
    skip     = 1,                                 # title row above header
    gene_col = "Gene", fc_col = "logFC",          fdr_col = "padj",
    set_a    = "Aged_up", set_b = "Young_up",
    label    = "Aged vs Young (Aksöz cohort)",
    selection = "topN", topN = 200                # top/bottom 200 by logFC (no padj filter)
  ),
  Ainciburu_Aged_vs_Young = list(
    file     = "C:/Users/fc809/Downloads/Jakobsen_HSPC_Ainciburu_Aged_vs_Young_donors.csv",
    skip     = 0,
    gene_col = "Gene", fc_col = "log2FoldChange", fdr_col = "padj",
    set_a    = "Aged_up", set_b = "Young_up",
    label    = "Aged vs Young (Ainciburu cohort)",
    selection = "topN", topN = 200                # top 200 by |logFC| per direction (no-op if set < 200)
  ),
  Inhouse_Aged_vs_Young = list(
    file     = "C:/Users/fc809/Downloads/Jakobsen_HSPC_In-house_Aged_vs_Young_donors.csv",
    skip     = 0,
    gene_col = "Gene", fc_col = "log2FoldChange", fdr_col = "padj",
    set_a    = "Aged_up", set_b = "Young_up",
    label    = "Aged vs Young (In-house cohort)",
    selection = "topN", topN = 200                # top 200 by |logFC| per direction
  )
)

## ── LT-HSC COMPARISON FILES ───────────────────────────────────────────────────

HSC_DIR <- "C:/Users/fc809/Downloads/LT-HSCs (1)"

hsc_files <- list(
  LentiPosLB_vs_LentiNegHB   = file.path(HSC_DIR, "LentiPosLB_vs_LentiNegHB",
                                         "DEG_LentiPosLB_vs_LentiNegHB_pseudoBulk_res0.3.sce_0.txt"),
  LentiNegHB_vs_LentiNegLB   = file.path(HSC_DIR, "LentiNegHB_vs_LentiNegLB",
                                         "DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.3.sce_0 - Copy.txt"),
  LentiPosHB_vs_LentiPosLB   = file.path(HSC_DIR, "LentiPosHB_vs_LentiPosLB",
                                         "DEG_LentiPosHB_vs_LentiPosLB_pseudoBulk_res0.3.sce_0.txt"),
  LentiPos_vs_LentiNeg_in_LB = file.path(HSC_DIR, "LentiPos_vs_LentiNeg_in_LB",
                                         "DEG_LentiPos_vs_LentiNeg_in_LB_pseudoBulk_res0.3.sce_0.txt"),
  LentiPos_vs_LentiNeg_in_HB = file.path(HSC_DIR, "LentiPos_vs_LentiNeg_in_HB",
                                         "DEG_LentiPos_vs_LentiNeg_in_HB_pseudoBulk_res0.3.sce_0.txt")
)

core_comparisons <- c(
  "LentiNegHB_vs_LentiNegLB",
  "LentiPosHB_vs_LentiPosLB",
  "LentiPos_vs_LentiNeg_in_LB",
  "LentiPos_vs_LentiNeg_in_HB"
)

clean_labels <- c(
  LentiNegHB_vs_LentiNegLB   = "LentiNeg HB vs LB",
  LentiPosHB_vs_LentiPosLB   = "LentiPos HB vs LB",
  LentiPos_vs_LentiNeg_in_LB = "Lenti+/- in LB",
  LentiPos_vs_LentiNeg_in_HB = "Lenti+/- in HB",
  LentiPosLB_vs_LentiNegHB   = "LentiPos-LB vs LentiNeg-HB"
)

## ── SHARED HELPERS ────────────────────────────────────────────────────────────

detect_cols <- function(df) {
  gene_col <- colnames(df)[1]
  fc_col   <- intersect(c("logFC", "log2FoldChange", "log2FC", "LogFC"),
                        colnames(df))[1]
  if (is.na(fc_col))
    stop("Cannot find logFC column. Columns present: ",
         paste(colnames(df), collapse = ", "))
  list(gene = gene_col, fc = fc_col)
}

make_ranked_hsc <- function(filepath) {
  df   <- read.delim(filepath, stringsAsFactors = FALSE)
  cols <- detect_cols(df)
  df %>%
    filter(.data[[cols$gene]] != "Lenti") %>%
    filter(!is.na(.data[[cols$fc]])) %>%
    arrange(desc(.data[[cols$fc]])) %>%
    distinct(.data[[cols$gene]], .keep_all = TRUE) %>%
    { setNames(.[[cols$fc]], .[[cols$gene]]) }
}

extract_leading_edge <- function(gsea_obj, set_id) {
  df  <- as.data.frame(gsea_obj)
  row <- df[df$ID == set_id, , drop = FALSE]
  if (nrow(row) == 0) return(character(0))
  le <- row$core_enrichment[1]
  if (is.na(le) || le == "") return(character(0))
  trimws(unlist(strsplit(as.character(le), "/")))
}

load_jakobsen <- function(cfg) {
  df <- read.csv(cfg$file, skip = cfg$skip, stringsAsFactors = FALSE,
                 check.names = FALSE)
  # Strip BOM from first column if present
  colnames(df)[1] <- sub("^\ufeff", "", colnames(df)[1])
  
  required <- c(cfg$gene_col, cfg$fc_col, cfg$fdr_col)
  missing  <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("Missing columns in ", basename(cfg$file), ": ",
         paste(missing, collapse = ", "),
         "\nFound: ", paste(colnames(df), collapse = ", "))
  }
  
  # Filter sex-specific genes when the column is present (Ainciburu + In-house)
  if ("SexSpecific" %in% colnames(df)) {
    n_before <- nrow(df)
    is_sex   <- toupper(trimws(as.character(df$SexSpecific))) %in% c("YES", "TRUE", "Y")
    df       <- df[!is_sex, , drop = FALSE]
    cat("  [", basename(cfg$file), "] SexSpecific=Yes filter: removed ",
        n_before - nrow(df), " of ", n_before, " genes\n", sep = "")
  }
  
  out <- data.frame(
    gene   = df[[cfg$gene_col]],
    log2FC = as.numeric(df[[cfg$fc_col]]),
    padj   = as.numeric(df[[cfg$fdr_col]]),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$log2FC) & !is.na(out$padj) & nzchar(out$gene), ]
  out
}

## ── LOAD HSC DEG TABLES + RANKED LISTS — ONCE ─────────────────────────────────

BiocParallel::register(BiocParallel::SerialParam())

cat("=== Column detection check (first HSC file) ===\n")
peek <- read.delim(hsc_files[[1]], nrow = 3, stringsAsFactors = FALSE)
cat("Columns:", paste(colnames(peek), collapse = ", "), "\n\n")

ranked_lists <- list()
deg_tables   <- list()
for (nm in names(hsc_files)) {
  cat("Loading HSC:", nm, "\n")
  deg_tables[[nm]]   <- read.delim(hsc_files[[nm]], stringsAsFactors = FALSE)
  ranked_lists[[nm]] <- make_ranked_hsc(hsc_files[[nm]])
  cat("  ranked list:", length(ranked_lists[[nm]]), "genes\n")
}

## ── GLOBAL ORTHOLOG MAP — ONE LOOKUP FOR ALL DATASETS ─────────────────────────

cat("\n=== Loading all six Jakobsen DE tables to build union gene list ===\n")
jak_tables <- lapply(datasets, load_jakobsen)
for (nm in names(jak_tables)) {
  d <- jak_tables[[nm]]
  cat(sprintf("%-26s n=%d  sig(padj<%.2g)=%d\n",
              nm, nrow(d), SIG_CUTOFF, sum(d$padj < SIG_CUTOFF)))
}

all_human_genes <- unique(unlist(lapply(jak_tables, function(d) d$gene)))
cat("\nUnion of unique human genes across all datasets:",
    length(all_human_genes), "\n")

cat("=== Human → mouse ortholog mapping (babelgene) ===\n")
ortho <- babelgene::orthologs(genes   = all_human_genes,
                              species = "mouse",
                              human   = TRUE)
cat("Mappings returned        :", nrow(ortho), "\n")
cat("Unique human w/ ortholog :", length(unique(ortho$human_symbol)), "\n")
cat("Unique mouse symbols     :", length(unique(ortho$symbol)), "\n\n")
ortho_map <- ortho %>% dplyr::select(human_symbol, mouse_symbol = symbol)

## ── PER-DATASET ANALYSIS FUNCTION ─────────────────────────────────────────────

run_analysis <- function(jak, cfg, outdir) {
  
  set_a <- cfg$set_a
  set_b <- cfg$set_b
  
  message("\n############################################################")
  message("###  Dataset : ", cfg$label)
  message("###  Sets    : ", set_a, "  (log2FC > 0)  /  ",
          set_b, "  (log2FC < 0)")
  message("###  Output  : ", outdir)
  message("############################################################\n")
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  ## Map orthologs (many-to-one dedup uses |log2FC|)
  jak_mouse <- jak %>%
    inner_join(ortho_map, by = c("gene" = "human_symbol"),
               relationship = "many-to-many") %>%
    group_by(mouse_symbol) %>%
    slice_max(abs(log2FC), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::rename(human_gene = gene, gene = mouse_symbol)
  
  cat("After mouse mapping & dedup :", nrow(jak_mouse), "mouse genes\n")
  cat("Significant (padj <", SIG_CUTOFF, "):",
      sum(jak_mouse$padj < SIG_CUTOFF), "mouse genes\n")
  
  ## Build gene sets — two selection strategies
  ##   "padj_split" (default) : padj < SIG_CUTOFF, then split by sign of logFC
  ##   "topN"                 : padj < SIG_CUTOFF, then within each direction
  ##                            take the top N by |logFC| (so set_a = top N most
  ##                            positive sig, set_b = top N most negative sig).
  ##                            Useful when sig set is so large the directional
  ##                            signal gets diluted (e.g. Aksöz cohort).
  selection <- if (is.null(cfg$selection)) "padj_split" else cfg$selection
  topN      <- if (is.null(cfg$topN))      200          else cfg$topN
  
  jak_sig <- jak_mouse %>%
    filter(padj < SIG_CUTOFF) %>%
    distinct(gene, .keep_all = TRUE)
  
  set_a_full <- jak_sig %>% filter(log2FC > 0) %>% arrange(desc(log2FC)) %>% pull(gene)
  set_b_full <- jak_sig %>% filter(log2FC < 0) %>% arrange(log2FC)       %>% pull(gene)
  
  if (selection == "topN") {
    set_a_genes <- head(set_a_full, topN)
    set_b_genes <- head(set_b_full, topN)
    method_label <- paste0("padj<", SIG_CUTOFF,
                           ", top ", topN, " by |log2FC| per direction",
                           " (full sig pool: ", length(set_a_full), "+/",
                           length(set_b_full), "-)")
  } else {
    set_a_genes <- set_a_full
    set_b_genes <- set_b_full
    method_label <- paste0("padj<", SIG_CUTOFF, ", split by sign of log2FC")
  }
  
  N_A <- length(set_a_genes); N_B <- length(set_b_genes)
  cat("Selection method:", method_label, "\n")
  cat("Set sizes:", set_a, "=", N_A, " / ", set_b, "=", N_B, "\n\n")
  
  if (N_A < 5 && N_B < 5) {
    warning("Both gene sets too small for ", cfg$label,
            " — skipping this dataset")
    return(NULL)
  }
  
  term2gene <- data.frame(
    term = c(rep(set_a, N_A), rep(set_b, N_B)),
    gene = c(set_a_genes, set_b_genes)
  )
  
  ## GSEA per HSC comparison
  run_gsea <- function(rl) {
    GSEA(geneList      = rl,
         TERM2GENE     = term2gene,
         minGSSize     = 10,
         maxGSSize     = max(5000, N_A, N_B) + 1,
         pvalueCutoff  = 1,
         pAdjustMethod = "BH",
         eps           = 0,
         verbose       = TRUE)
  }
  
  plot_gsea <- function(gsea_obj, title, outpath) {
    stats <- as.data.frame(gsea_obj)[, c("ID", "NES", "p.adjust")]
    cat("\n---", title, "---\n"); print(stats)
    
    ids_present <- intersect(c(set_a, set_b), stats$ID)
    if (length(ids_present) == 0) {
      message("No gene sets returned for ", title, " — skipping plot")
      return(invisible(NULL))
    }
    
    nes_a <- if (set_a %in% stats$ID) round(stats$NES[stats$ID == set_a], 2) else NA
    fdr_a <- if (set_a %in% stats$ID) signif(stats$p.adjust[stats$ID == set_a], 2) else NA
    nes_b <- if (set_b %in% stats$ID) round(stats$NES[stats$ID == set_b], 2) else NA
    fdr_b <- if (set_b %in% stats$ID) signif(stats$p.adjust[stats$ID == set_b], 2) else NA
    
    p <- gseaplot2(gsea_obj,
                   geneSetID = ids_present,
                   title     = paste0(title,
                                      "\n(", cfg$label,
                                      "  |  ", method_label, "  |  ",
                                      set_a, " n=", N_A, ", ",
                                      set_b, " n=", N_B, ")"),
                   subplots  = 1:3,
                   base_size = 11)
    
    pal <- setNames(c("#B71C1C", "#1A237E"), c(set_a, set_b))
    labels_vec <- setNames(
      c(paste0(set_a, " (n=", N_A, ")  NES=", nes_a, "  FDR=", fdr_a),
        paste0(set_b, " (n=", N_B, ")  NES=", nes_b, "  FDR=", fdr_b)),
      c(set_a, set_b)
    )
    
    p[[1]] <- p[[1]] +
      scale_color_manual(values = pal[ids_present],
                         labels = labels_vec[ids_present]) +
      theme(legend.position = "top", legend.text = element_text(size = 9))
    p[[2]] <- p[[2]] + scale_color_manual(values = pal[ids_present])
    
    png(outpath, width = 1400, height = 1200, res = 150)
    print(p); dev.off()
    message("Saved: ", basename(outpath))
  }
  
  gsea_results <- list()
  for (nm in names(hsc_files)) {
    cat("\n=== GSEA:", nm, "===\n")
    gsea_results[[nm]] <- run_gsea(ranked_lists[[nm]])
    plot_gsea(
      gsea_obj = gsea_results[[nm]],
      title    = paste0(cfg$label, ": ", clean_labels[nm], " (LT-HSCs)"),
      outpath  = file.path(outdir, paste0("GSEA_on_", nm, ".png"))
    )
  }
  message("\nIndividual GSEA plots done.")
  
  ## Per-dataset summary dotplot
  summary_rows <- lapply(core_comparisons, function(nm) {
    df <- as.data.frame(gsea_results[[nm]])
    df <- df[df$ID %in% c(set_a, set_b), c("ID", "NES", "p.adjust")]
    if (nrow(df) == 0) return(NULL)
    df$comparison <- clean_labels[nm]
    df
  })
  gsea_summary <- do.call(rbind, summary_rows)
  
  if (!is.null(gsea_summary) && nrow(gsea_summary) > 0) {
    gsea_summary$comparison <- factor(gsea_summary$comparison,
                                      levels = unname(clean_labels[core_comparisons]))
    gsea_summary$ID <- factor(gsea_summary$ID, levels = c(set_a, set_b))
    id_labels <- c(setNames(paste0(set_a, "\n(n=", N_A, ")"), set_a),
                   setNames(paste0(set_b, "\n(n=", N_B, ")"), set_b))
    
    p_dot <- ggplot(gsea_summary,
                    aes(x = comparison, y = ID, fill = NES,
                        size = -log10(p.adjust))) +
      geom_point(shape = 21, color = "grey30") +
      scale_fill_gradient2(low = "#1A237E", mid = "white", high = "#B71C1C",
                           midpoint = 0, limits = c(-2.5, 2.5),
                           oob = scales::squish, name = "NES") +
      scale_size_continuous(range = c(4, 14), name = "-log10(FDR)") +
      scale_y_discrete(labels = id_labels, expand = expansion(add = 0.8)) +
      theme_classic(base_size = 13) +
      theme(axis.text.x   = element_text(size = 10, angle = 30, hjust = 1,
                                         margin = margin(t = 4)),
            axis.text.y   = element_text(size = 11),
            plot.title    = element_text(size = 12, hjust = 0.5),
            plot.subtitle = element_text(size = 10, hjust = 0.5),
            plot.margin   = margin(t = 15, r = 10, b = 25, l = 15)) +
      labs(x = NULL, y = "Gene set",
           title    = paste0(cfg$label, " across LT-HSC comparisons"),
           subtitle = paste0(method_label, "  |  ",
                             set_a, " n=", N_A, ", ", set_b, " n=", N_B))
    
    ggsave(file.path(outdir, "Summary_dotplot.png"),
           p_dot, width = 8, height = 4.5, dpi = 150)
    message("Saved: Summary_dotplot.png")
  }
  
  ## Leading edge heatmaps + Venn
  cols_info <- detect_cols(deg_tables[[core_comparisons[1]]])
  gene_col  <- cols_info$gene
  fc_col    <- cols_info$fc
  
  cat("\n=== Leading edge sizes ===\n")
  for (nm in core_comparisons) {
    for (sid in c(set_a, set_b)) {
      le <- extract_leading_edge(gsea_results[[nm]], sid)
      cat(nm, "|", sid, ":", length(le), "genes\n")
    }
  }
  
  for (set_id in c(set_a, set_b)) {
    le_genes <- unique(unlist(lapply(core_comparisons, function(nm)
      extract_leading_edge(gsea_results[[nm]], set_id))))
    cat("\nLeading edge union (", set_id, "):", length(le_genes), "genes\n")
    if (length(le_genes) == 0) {
      message("No leading edge genes for ", set_id, " — skipping heatmap/Venn")
      next
    }
    
    fc_mat <- sapply(core_comparisons, function(nm) {
      df <- deg_tables[[nm]]
      fc <- setNames(df[[fc_col]], df[[gene_col]])
      fc[le_genes]
    })
    rownames(fc_mat) <- le_genes
    colnames(fc_mat) <- unname(clean_labels[core_comparisons])
    fc_mat <- fc_mat[rowSums(is.na(fc_mat)) <= 2, , drop = FALSE]
    fc_mat[is.na(fc_mat)] <- 0
    fc_mat <- pmax(pmin(fc_mat, 2), -2)
    
    set_total <- if (set_id == set_a) N_A else N_B
    png(file.path(outdir, paste0("LeadingEdge_heatmap_", set_id, ".png")),
        width = 1000, height = max(700, nrow(fc_mat) * 14 + 250), res = 130)
    pheatmap(
      fc_mat,
      color             = colorRampPalette(c("#1A237E", "white", "#B71C1C"))(100),
      breaks            = seq(-2, 2, length.out = 101),
      cluster_cols      = FALSE,
      clustering_method = "ward.D2",
      fontsize_row      = 7, fontsize_col = 12,
      border_color      = NA,
      main              = paste0("Leading edge — ", set_id,
                                 "  (", nrow(fc_mat), " / ", set_total,
                                 ")  |  ", cfg$label))
    dev.off()
    message("Saved: LeadingEdge_heatmap_", set_id, ".png")
    
    le_list <- setNames(
      lapply(core_comparisons, function(nm)
        extract_leading_edge(gsea_results[[nm]], set_id)),
      unname(clean_labels[core_comparisons])
    )
    if (any(lengths(le_list) > 0)) {
      fill_high <- if (set_id == set_a) "#B71C1C" else "#1A237E"
      p_venn <- ggVennDiagram(le_list, label_alpha = 0, set_size = 3.5) +
        scale_fill_gradient(low = "white", high = fill_high) +
        labs(title = paste0("Leading edge — ", set_id,
                            "  (of ", set_total, ")  |  ", cfg$label)) +
        theme(legend.position = "none",
              plot.title      = element_text(hjust = 0.5, size = 11))
      ggsave(file.path(outdir, paste0("LeadingEdge_Venn_", set_id, ".png")),
             p_venn, width = 7, height = 6, dpi = 150)
      message("Saved: LeadingEdge_Venn_", set_id, ".png")
    }
  }
  
  ## Ranked-list scatter (highlighting any leading-edge gene from this dataset)
  scatter_pairs <- list(
    list(x = "LentiNegHB_vs_LentiNegLB", y = "LentiPosHB_vs_LentiPosLB",
         xlab = "logFC (LentiNeg HB vs LB)",
         ylab = "logFC (LentiPos HB vs LB)",
         title = "HB vs LB: LentiNeg vs LentiPos"),
    list(x = "LentiPos_vs_LentiNeg_in_LB", y = "LentiPos_vs_LentiNeg_in_HB",
         xlab = "logFC (Lenti+/- in LB)",
         ylab = "logFC (Lenti+/- in HB)",
         title = "Lenti+/- effect: LB vs HB context")
  )
  all_le <- unique(unlist(lapply(core_comparisons, function(nm)
    c(extract_leading_edge(gsea_results[[nm]], set_a),
      extract_leading_edge(gsea_results[[nm]], set_b)))))
  
  plot_list <- lapply(scatter_pairs, function(p) {
    fc_x   <- ranked_lists[[p$x]]
    fc_y   <- ranked_lists[[p$y]]
    shared <- intersect(names(fc_x), names(fc_y))
    df <- data.frame(gene  = shared,
                     x     = as.numeric(fc_x[shared]),
                     y     = as.numeric(fc_y[shared]),
                     is_le = shared %in% all_le)
    df <- df[complete.cases(df), ]
    r <- sprintf("%.3f", cor(df$x, df$y))
    
    ggplot(df, aes(x = x, y = y)) +
      geom_point(data = subset(df, !is_le),
                 color = "grey80", size = 0.4, alpha = 0.5) +
      geom_point(data = subset(df,  is_le),
                 aes(color = x + y), size = 1.5, alpha = 0.9) +
      scale_color_gradient2(low = "#1A237E", mid = "grey50", high = "#B71C1C",
                            midpoint = 0,
                            name = "Leading edge\n(sum logFC)",
                            breaks = c(-2, 0, 2),
                            labels = c("Down both", "Mixed", "Up both")) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
      geom_smooth(method = "lm", se = FALSE, color = "black",
                  linewidth = 0.6, linetype = "dotted") +
      annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
               label = paste0("r = ", r), size = 4) +
      labs(x = p$xlab, y = p$ylab, title = p$title) +
      theme_classic(base_size = 11) +
      theme(legend.title    = element_text(size = 8),
            legend.text     = element_text(size = 7),
            legend.position = c(0.13, 0.82),
            legend.background = element_rect(fill = alpha("white", 0.7),
                                             color = NA))
  })
  p_scatter <- cowplot::plot_grid(plotlist = plot_list, nrow = 1)
  ggsave(file.path(outdir, "RankedList_scatter.png"),
         p_scatter, width = 11, height = 5, dpi = 150)
  message("Saved: RankedList_scatter.png")
  
  ## Export ortholog table + gene sets
  write.csv(jak_mouse %>% arrange(desc(log2FC)),
            file.path(outdir, "orthologs.csv"),
            row.names = FALSE)
  write.csv(data.frame(set  = c(rep(set_a, N_A), rep(set_b, N_B)),
                       gene = c(set_a_genes, set_b_genes)),
            file.path(outdir, "gene_sets.csv"),
            row.names = FALSE)
  message("Saved: orthologs.csv, gene_sets.csv")
  
  ## ---- Pathway analysis on the original HUMAN data --------------------------
  
  if (RUN_PATHWAYS) {
    pw_dir <- file.path(outdir, "pathways")
    dir.create(pw_dir, showWarnings = FALSE)
    
    ranked_jak_full <- jak %>%
      arrange(desc(log2FC)) %>%
      distinct(gene, .keep_all = TRUE) %>%
      { setNames(.$log2FC, .$gene) }
    
    sym2entrez <- bitr(names(ranked_jak_full),
                       fromType = "SYMBOL",
                       toType   = "ENTREZID",
                       OrgDb    = org.Hs.eg.db)
    ranked_jak_entrez <- ranked_jak_full[sym2entrez$SYMBOL]
    names(ranked_jak_entrez) <- sym2entrez$ENTREZID
    ranked_jak_entrez <- ranked_jak_entrez[!duplicated(names(ranked_jak_entrez))]
    ranked_jak_entrez <- sort(ranked_jak_entrez, decreasing = TRUE)
    
    save_pathway_plots <- function(gsea_obj, prefix, top_n = 20) {
      res_df <- as.data.frame(gsea_obj)
      write.csv(res_df, file.path(pw_dir, paste0(prefix, ".csv")),
                row.names = FALSE)
      cat("\n---", prefix, "--- significant (FDR<0.05):",
          sum(res_df$p.adjust < 0.05), "\n")
      
      sig <- res_df[res_df$p.adjust < 0.05, ]
      if (nrow(sig) == 0) {
        message("No significant pathways for ", prefix, " — using top ", top_n)
        sig <- res_df[order(res_df$pvalue), ][seq_len(min(top_n, nrow(res_df))), ]
      }
      top_pos <- sig[sig$NES > 0, ][order(-sig$NES[sig$NES > 0]), ][seq_len(min(15, sum(sig$NES > 0))), ]
      top_neg <- sig[sig$NES < 0, ][order( sig$NES[sig$NES < 0]), ][seq_len(min(15, sum(sig$NES < 0))), ]
      bar_df  <- rbind(top_pos, top_neg)
      
      if (nrow(bar_df) > 0) {
        bar_df$Description <- ifelse(nchar(bar_df$Description) > 55,
                                     paste0(substr(bar_df$Description, 1, 52), "..."),
                                     bar_df$Description)
        # Truncation can collide → make unique before assigning as factor levels
        bar_df$Description <- make.unique(bar_df$Description, sep = " ")
        bar_df$Description <- factor(bar_df$Description,
                                     levels = bar_df$Description[order(bar_df$NES)])
        p_bar <- ggplot(bar_df, aes(x = NES, y = Description, fill = p.adjust)) +
          geom_col() +
          scale_fill_gradient(low = "#B71C1C", high = "grey80", name = "FDR",
                              limits = c(0, 0.05)) +
          geom_vline(xintercept = 0, color = "grey30") +
          labs(x = "Normalised Enrichment Score", y = NULL,
               title = paste0(prefix, "  |  positive NES = ",
                              set_a, " (", cfg$label, ")")) +
          theme_classic(base_size = 10) +
          theme(plot.title  = element_text(size = 9, hjust = 0.5),
                axis.text.y = element_text(size = 7))
        ggsave(file.path(pw_dir, paste0(prefix, "_bar.png")),
               p_bar, width = 10, height = max(4, nrow(bar_df) * 0.3 + 1.5),
               dpi = 150)
      }
      
      tryCatch({
        p_dot <- dotplot(gsea_obj, showCategory = top_n, split = ".sign") +
          facet_grid(. ~ .sign) +
          theme(strip.text  = element_text(size = 9),
                axis.text.y = element_text(size = 7))
        ggsave(file.path(pw_dir, paste0(prefix, "_dot.png")),
               p_dot, width = 10,
               height = max(5, min(top_n, nrow(res_df)) * 0.28 + 2), dpi = 150)
      }, error = function(e) message("dotplot failed for ", prefix, ": ", e$message))
    }
    
    cat("\n=== GO BP (human) ===\n")
    gsea_go <- gseGO(geneList = ranked_jak_full, OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL", ont = "BP",
                     minGSSize = 15, maxGSSize = 500,
                     pvalueCutoff = 1, pAdjustMethod = "BH",
                     eps = 0, verbose = TRUE)
    gsea_go <- setReadable(gsea_go, OrgDb = org.Hs.eg.db, keyType = "SYMBOL")
    save_pathway_plots(gsea_go, "GSEA_GOBP")
    
    cat("\n=== KEGG (human) ===\n")
    gsea_kegg <- gseKEGG(geneList = ranked_jak_entrez, organism = "hsa",
                         minGSSize = 15, maxGSSize = 500,
                         pvalueCutoff = 1, pAdjustMethod = "BH",
                         eps = 0, verbose = TRUE)
    gsea_kegg <- setReadable(gsea_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    save_pathway_plots(gsea_kegg, "GSEA_KEGG")
    
    cat("\n=== Hallmarks (human) ===\n")
    hallmarks <- msigdbr(species = "Homo sapiens", category = "H") %>%
      dplyr::select(gs_name, gene_symbol)
    gsea_h <- GSEA(geneList = ranked_jak_full, TERM2GENE = hallmarks,
                   minGSSize = 15, maxGSSize = 500,
                   pvalueCutoff = 1, pAdjustMethod = "BH",
                   eps = 0, verbose = TRUE)
    save_pathway_plots(gsea_h, "GSEA_Hallmarks")
  }
  
  ## Return summary for combined plot
  if (!is.null(gsea_summary) && nrow(gsea_summary) > 0) {
    out <- gsea_summary
    out$dataset <- cfg$label
    out$N_set   <- ifelse(out$ID == set_a, N_A, N_B)
    out$ID      <- as.character(out$ID)
    return(out)
  }
  return(NULL)
}

## ── EXECUTE OVER ALL DATASETS ─────────────────────────────────────────────────

combined <- list()
for (nm in names(datasets)) {
  cfg <- datasets[[nm]]
  jak <- jak_tables[[nm]]
  res <- tryCatch(
    run_analysis(jak, cfg, file.path(OUTDIR, nm)),
    error = function(e) {
      message("ERROR in ", nm, ": ", e$message)
      NULL
    }
  )
  combined[[nm]] <- res
}

## ── COMBINED FACETED DOTPLOT ──────────────────────────────────────────────────

combined_df <- bind_rows(combined, .id = "dataset_key")
write.csv(combined_df, file.path(OUTDIR, "Combined_summary.csv"),
          row.names = FALSE)
message("Saved: Combined_summary.csv")

if (nrow(combined_df) > 0) {
  combined_df$dataset    <- factor(combined_df$dataset,
                                   levels = unname(sapply(datasets, `[[`, "label")))
  combined_df$comparison <- factor(combined_df$comparison,
                                   levels = unname(clean_labels[core_comparisons]))
  # Order direction: set_a (positive direction) above set_b (negative)
  # Build a per-dataset ordered ID factor
  id_levels <- unique(do.call(c, lapply(datasets, function(d) c(d$set_a, d$set_b))))
  combined_df$ID <- factor(combined_df$ID, levels = id_levels)
  
  # Annotate set sizes in y-axis labels per facet — easiest via ID-with-N relabel
  combined_df$ID_label <- paste0(combined_df$ID, "\n(n=", combined_df$N_set, ")")
  # Keep ordering consistent within facet
  combined_df$ID_label <- factor(combined_df$ID_label,
                                 levels = unique(combined_df$ID_label[
                                   order(combined_df$dataset, combined_df$ID)]))
  
  p_combined <- ggplot(combined_df,
                       aes(x = comparison, y = ID_label, fill = NES,
                           size = -log10(p.adjust))) +
    geom_point(shape = 21, color = "grey30") +
    facet_grid(dataset ~ ., scales = "free_y", space = "free_y",
               switch = "y") +
    scale_fill_gradient2(low = "#1A237E", mid = "white", high = "#B71C1C",
                         midpoint = 0, limits = c(-2.5, 2.5),
                         oob = scales::squish, name = "NES") +
    scale_size_continuous(range = c(3, 11), name = "-log10(FDR)") +
    theme_classic(base_size = 11) +
    theme(axis.text.x      = element_text(size = 9, angle = 30, hjust = 1),
          axis.text.y      = element_text(size = 8),
          strip.text.y.left = element_text(angle = 0, size = 8, hjust = 1),
          strip.background = element_rect(fill = "grey95", color = NA),
          strip.placement  = "outside",
          panel.spacing.y  = unit(0.4, "lines"),
          plot.title       = element_text(size = 12, hjust = 0.5),
          plot.subtitle    = element_text(size = 9, hjust = 0.5)) +
    labs(x = NULL, y = NULL,
         title    = "Jakobsen HSPC signatures across LT-HSC comparisons",
         subtitle = paste0("padj<", SIG_CUTOFF,
                           " split-by-sign for TET2 + CH datasets; ",
                           "top 200 by |logFC| per direction for aging cohorts",
                           "  |  SexSpecific genes filtered when annotated"))
  
  n_facets <- length(unique(combined_df$dataset))
  ggsave(file.path(OUTDIR, "Combined_dotplot.png"),
         p_combined, width = 9, height = 1.6 * n_facets + 2, dpi = 150)
  message("Saved: Combined_dotplot.png")
}

## ── COPY SCRIPT FOR PROVENANCE ────────────────────────────────────────────────

if (file.exists(SCRIPT_PATH)) {
  file.copy(SCRIPT_PATH, file.path(OUTDIR, basename(SCRIPT_PATH)),
            overwrite = TRUE)
  message("Script copied to: ", file.path(OUTDIR, basename(SCRIPT_PATH)))
} else {
  warning("SCRIPT_PATH not found — script not copied: ", SCRIPT_PATH)
}

message("\n############################################################")
message("###  ALL DONE")
message("###  Output root: ", OUTDIR)
message("############################################################")