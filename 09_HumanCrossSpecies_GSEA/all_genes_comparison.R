#!/usr/bin/env Rscript
# =============================================================================
# GSEA for DGE output (human genes) — Hallmark + GO:BP
# -----------------------------------------------------------------------------
# Handles two input formats automatically:
#   * DESeq2 (.tsv): gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
#   * MAST   (.csv): "", primerid, Pr(>Chisq), coef, ci.hi, ci.lo, fdr, label
#
# Runs each file twice, with two ranking metrics:
#   metric 1  "logFC"        : effect size              (DESeq2 log2FoldChange | MAST coef)
#   metric 2  "logFC_signedP": effect * -log10(padj)    (DESeq2 padj           | MAST fdr)
#
# Gene sets: MSigDB Hallmark (H) and GO Biological Process (C5, GO:BP),
#            human gene symbols, pulled via the msigdbr package.
#
# Output: one results table (.csv) per file × metric × collection, written to
#         <OUTPUT_DIR>, plus a combined "GSEA_all_results.csv".
# =============================================================================

# ----------------------------- CONFIG ----------------------------------------
# Folder to scan for DGE tables (searched recursively). Edit if needed.
INPUT_DIR  <- "C:/Users/fc809/Downloads/DGE_tables/DGE_tables"   # holds the DESeq2_* / MAST_* subfolders
OUTPUT_DIR <- "C:/Users/fc809/Downloads/DGE_tables/GSEA_results" # results written here (created if missing)

# fgsea parameters
MIN_SIZE   <- 15
MAX_SIZE   <- 500
SPECIES    <- "Homo sapiens"

# padj floor: padj of 0 -> -log10(0) = Inf. Replace 0 with the smallest
# non-zero padj in that file before taking -log10 (avoids Inf weights).
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  # ---- install missing packages -------------------------------------------
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  if (!requireNamespace("fgsea", quietly = TRUE))        BiocManager::install("fgsea", update = FALSE, ask = FALSE)
  if (!requireNamespace("msigdbr", quietly = TRUE))      install.packages("msigdbr", repos = "https://cloud.r-project.org")
  if (!requireNamespace("data.table", quietly = TRUE))   install.packages("data.table", repos = "https://cloud.r-project.org")
  library(fgsea)
  library(msigdbr)
  library(data.table)
})

set.seed(42)

# ----------------------- build gene set collections ---------------------------
# msigdbr changed its argument names across versions; try the new API first,
# then fall back to the legacy one.
get_sets <- function(collection, subcollection = NULL) {
  df <- tryCatch(
    {
      if (is.null(subcollection)) {
        msigdbr(species = SPECIES, collection = collection)
      } else {
        msigdbr(species = SPECIES, collection = collection, subcollection = subcollection)
      }
    },
    error = function(e) {
      if (is.null(subcollection)) {
        msigdbr(species = SPECIES, category = collection)
      } else {
        msigdbr(species = SPECIES, category = collection, subcategory = subcollection)
      }
    }
  )
  # column holding the gene set name / gene symbol also varies by version
  name_col <- intersect(c("gs_name"), names(df))[1]
  sym_col  <- intersect(c("gene_symbol", "human_gene_symbol"), names(df))[1]
  split(df[[sym_col]], df[[name_col]])
}

message("Loading gene sets (Hallmark + GO:BP) via msigdbr ...")
gene_sets <- list(
  Hallmark = get_sets("H"),
  GOBP     = tryCatch(get_sets("C5", "GO:BP"), error = function(e) get_sets("C5", "BP"))
)
for (nm in names(gene_sets))
  message(sprintf("  %-8s: %d gene sets", nm, length(gene_sets[[nm]])))

# --------------------------- file discovery -----------------------------------
files <- list.files(INPUT_DIR, pattern = "\\.(tsv|csv)$", recursive = TRUE, full.names = TRUE)
if (length(files) == 0) stop(sprintf("No .tsv/.csv files found under '%s'. Check INPUT_DIR.", INPUT_DIR))
message(sprintf("Found %d input file(s).", length(files)))

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------- parse a file into a tidy data.frame ----------------------
# Returns data.frame with columns: gene, logFC, padj  (or NULL if unrecognized)
parse_dge <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "tsv") {
    d <- fread(path, sep = "\t", data.table = FALSE)
  } else {
    d <- fread(path, sep = ",", data.table = FALSE)
  }
  cn <- names(d)
  
  # ---- DESeq2 format ----
  if (all(c("gene", "log2FoldChange") %in% cn)) {
    padj <- if ("padj" %in% cn) d$padj else d$pvalue
    return(data.frame(gene = as.character(d$gene),
                      logFC = as.numeric(d$log2FoldChange),
                      padj  = as.numeric(padj),
                      stringsAsFactors = FALSE))
  }
  # ---- MAST format ----
  if ("primerid" %in% cn && "coef" %in% cn) {
    padj <- if ("fdr" %in% cn) d$fdr else d[["Pr(>Chisq)"]]
    return(data.frame(gene = as.character(d$primerid),
                      logFC = as.numeric(d$coef),
                      padj  = as.numeric(padj),
                      stringsAsFactors = FALSE))
  }
  warning(sprintf("Unrecognized format, skipping: %s", path))
  NULL
}

# ----------------- build a named, sorted ranking vector -----------------------
make_ranks <- function(df, metric) {
  df <- df[!is.na(df$gene) & !is.na(df$logFC), , drop = FALSE]
  df <- df[!duplicated(df$gene), , drop = FALSE]          # keep first occurrence per gene
  
  if (metric == "logFC") {
    v <- df$logFC
  } else { # logFC_signedP : effect * -log10(padj)
    p <- df$padj
    p[is.na(p)] <- 1                                      # missing padj -> weight 0
    nz <- min(p[p > 0], na.rm = TRUE)                     # smallest non-zero padj
    p[p == 0] <- nz                                       # floor zeros to avoid Inf
    v <- df$logFC * -log10(p)
  }
  names(v) <- df$gene
  v <- v[is.finite(v)]
  sort(v, decreasing = TRUE)
}

# --------------------------------- run ----------------------------------------
metrics <- c("logFC", "logFC_signedP")
all_results <- list()

for (path in files) {
  df <- parse_dge(path)
  if (is.null(df)) next
  rel  <- sub(paste0("^", INPUT_DIR, "/?"), "", path)
  tag  <- gsub("[/ ]+", "_", tools::file_path_sans_ext(rel))   # safe label per file
  
  for (metric in metrics) {
    ranks <- make_ranks(df, metric)
    if (length(ranks) < MIN_SIZE) {
      warning(sprintf("Too few ranked genes (%d) for %s [%s]; skipping.", length(ranks), rel, metric))
      next
    }
    for (coll in names(gene_sets)) {
      message(sprintf("GSEA: %-55s | %-13s | %s", rel, metric, coll))
      res <- fgsea(pathways = gene_sets[[coll]],
                   stats    = ranks,
                   minSize  = MIN_SIZE,
                   maxSize  = MAX_SIZE,
                   eps      = 0)
      res <- res[order(res$padj), ]
      # leadingEdge is a list column -> collapse to a string for CSV export
      res$leadingEdge <- vapply(res$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
      res$file        <- rel
      res$metric      <- metric
      res$collection  <- coll
      
      out <- file.path(OUTPUT_DIR, sprintf("%s__%s__%s.csv", tag, metric, coll))
      fwrite(res, out)
      all_results[[length(all_results) + 1]] <- res
    }
  }
}

# --------------------------- combined output ----------------------------------
if (length(all_results)) {
  combined <- rbindlist(all_results, use.names = TRUE, fill = TRUE)
  setcolorder(combined, c("file", "metric", "collection", "pathway",
                          "pval", "padj", "NES", "ES", "size"))
  fwrite(combined, file.path(OUTPUT_DIR, "GSEA_all_results.csv"))
  message(sprintf("\nDone. %d result tables + combined file written to '%s/'.",
                  length(all_results), OUTPUT_DIR))
} else {
  message("No results produced — check inputs.")
}

###################

#!/usr/bin/env Rscript
# =============================================================================
# Visuals for the human CD14 monocyte GSEA (Hallmark + GO:BP)
# -----------------------------------------------------------------------------
# Scoped to the FOUR files of interest only. Reads the result CSVs produced by
# run_gsea.R (GSEA_results/) for dotplots / barplots / cross-contrast heatmap;
# recomputes the ranked lists + gene sets only for the running-ES enrichment
# plots (the curves aren't stored in the CSVs).
#
# Produces, in GSEA_results/figures/ :
#   dotplot_<collection>.pdf          top pathways, x=NES, size=-log10FDR
#   barplot_<collection>.pdf          top pathways as horizontal NES bars
#   heatmap_<collection>.pdf          pathways x 4 contrasts, NES + * (FDR<0.05)
#   enrichment_<label>.pdf            running-ES for top pathways (per file)
# =============================================================================

# ============================== CONFIG =======================================
GSEA_DIR <- "C:/Users/fc809/Downloads/DGE_tables/GSEA_results"   # run_gsea.R output
FIG_DIR  <- file.path(GSEA_DIR, "figures")

# label -> (relative key as stored in the CSV 'file' column)  +  raw DGE path
TARGETS <- list(
  "TET2 CHIP vs Control"    = list(
    key = "DESeq2_Cohort_CCUS_CHIP_vs_Control/CD14+_Monocytes_TET2_CHIP_vs_Control",
    raw = "C:/Users/fc809/Downloads/DGE_tables/DGE_tables/DESeq2_Cohort_CCUS_CHIP_vs_Control/CD14+_Monocytes_TET2_CHIP_vs_Control.tsv"),
  "TET2 CH vs none"         = list(
    key = "DESeq2_Cohort_CH_vs_Control/CD14+_Monocytes_TET2_vs_none",
    raw = "C:/Users/fc809/Downloads/DGE_tables/DGE_tables/DESeq2_Cohort_CH_vs_Control/CD14+_Monocytes_TET2_vs_none.tsv"),
  "TET2 Mutant vs Control"  = list(
    key = "MAST_DE_Mutant_vs_Control/TET2 CD14 Mutant vs Control",
    raw = "C:/Users/fc809/Downloads/DGE_tables/DGE_tables/MAST_DE_Mutant_vs_Control/TET2 CD14 Mutant vs Control.csv"),
  "TET2 Mutant vs Wildtype" = list(
    key = "MAST_DE_Mutant_vs_Wildtype/TET2 CD14 Mutant vs Wildtype",
    raw = "C:/Users/fc809/Downloads/DGE_tables/DGE_tables/MAST_DE_Mutant_vs_Wildtype/TET2 CD14 Mutant vs Wildtype.csv")
)

COLLECTIONS  <- c("Hallmark", "GOBP")
METRICS      <- c("logFC", "logFC_signedP")
TOP_DOT      <- 20     # pathways per panel in dot/bar plots
TOP_ENRICH   <- 10      # top significant pathways to draw running-ES per file
FDR_CUTOFF   <- 0.05
MIN_SIZE     <- 15
MAX_SIZE     <- 500
# =============================================================================

suppressPackageStartupMessages({
  need <- function(p, bioc = FALSE) if (!requireNamespace(p, quietly = TRUE)) {
    if (bioc) { if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); BiocManager::install(p, update=FALSE, ask=FALSE) }
    else install.packages(p, repos = "https://cloud.r-project.org")
  }
  need("fgsea", bioc = TRUE); need("msigdbr"); need("data.table"); need("ggplot2"); need("patchwork")
  library(fgsea); library(msigdbr); library(data.table); library(ggplot2); library(patchwork)
})
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
COL_POS <- "#C0392B"; COL_NEG <- "#1F3A93"

prettify <- function(x) {
  x <- sub("^HALLMARK_", "", x); x <- sub("^GOBP_", "", x)
  x <- gsub("_", " ", x); tolower(x)
}

# --------- load the per-file result CSVs for the 4 targets -------------------
tag_of <- function(key) gsub("[/ ]+", "_", key)   # mirrors run_gsea.R naming
load_res <- function() {
  rows <- list()
  for (lab in names(TARGETS)) for (mt in METRICS) for (co in COLLECTIONS) {
    f <- file.path(GSEA_DIR, sprintf("%s__%s__%s.csv", tag_of(TARGETS[[lab]]$key), mt, co))
    if (!file.exists(f)) { warning(sprintf("missing result CSV: %s", basename(f))); next }
    d <- fread(f)
    d$label <- lab; d$metric <- mt; d$collection <- co
    rows[[length(rows)+1]] <- d
  }
  if (!length(rows)) stop("No result CSVs found - has run_gsea.R finished these 4 files?")
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}
res <- load_res()
res$plab <- prettify(res$pathway)

# ============================ (1) DOTPLOTS ===================================
for (co in COLLECTIONS) {
  pdf(file.path(FIG_DIR, sprintf("dotplot_%s.pdf", co)), width = 9, height = 7, onefile = TRUE)
  for (lab in names(TARGETS)) for (mt in METRICS) {
    d <- res[res$collection==co & res$label==lab & res$metric==mt & is.finite(res$NES), ]
    if (!nrow(d)) next
    d <- d[order(d$padj), ][seq_len(min(TOP_DOT, nrow(d))), ]
    d$plab <- factor(d$plab, levels = d$plab[order(d$NES)])
    p <- ggplot(d, aes(NES, plab, size = -log10(padj), fill = NES)) +
      geom_vline(xintercept = 0, color = "grey80") +
      geom_point(shape = 21, color = "grey30") +
      scale_fill_gradient2(low = COL_NEG, mid = "white", high = COL_POS, midpoint = 0) +
      scale_size_continuous(range = c(2, 7)) +
      labs(title = sprintf("%s  -  %s  [%s]", lab, co, mt),
           y = NULL, x = "NES", size = "-log10 FDR") +
      theme_bw(base_size = 10)
    print(p)
  }
  dev.off()
}

# ============================ (2) BARPLOTS ===================================
for (co in COLLECTIONS) {
  pdf(file.path(FIG_DIR, sprintf("barplot_%s.pdf", co)), width = 9, height = 7, onefile = TRUE)
  for (lab in names(TARGETS)) for (mt in METRICS) {
    d <- res[res$collection==co & res$label==lab & res$metric==mt & is.finite(res$NES), ]
    if (!nrow(d)) next
    d <- d[order(d$padj), ][seq_len(min(TOP_DOT, nrow(d))), ]
    d$plab <- factor(d$plab, levels = d$plab[order(d$NES)])
    d$sig  <- ifelse(d$padj < FDR_CUTOFF, "FDR<0.05", "ns")
    p <- ggplot(d, aes(NES, plab, fill = NES, alpha = sig)) +
      geom_col(color = "grey40", linewidth = .2) +
      scale_fill_gradient2(low = COL_NEG, mid = "white", high = COL_POS, midpoint = 0) +
      scale_alpha_manual(values = c("FDR<0.05" = 1, "ns" = 0.35)) +
      labs(title = sprintf("%s  -  %s  [%s]", lab, co, mt), y = NULL, x = "NES") +
      theme_bw(base_size = 10)
    print(p)
  }
  dev.off()
}

# ==================== (3) CROSS-CONTRAST NES HEATMAP =========================
for (co in COLLECTIONS) for (mt in METRICS) {
  sub <- res[res$collection==co & res$metric==mt & is.finite(res$NES), ]
  if (!nrow(sub)) next
  top <- unique(unlist(lapply(split(sub, sub$label), function(x)
    x$pathway[order(x$padj)][seq_len(min(15, nrow(x)))])))
  h <- sub[sub$pathway %in% top, ]
  h$plab <- prettify(h$pathway)
  h$star <- ifelse(h$padj < FDR_CUTOFF, "*", "")
  ord <- tapply(h$NES, h$plab, mean); h$plab <- factor(h$plab, levels = names(sort(ord)))
  g <- ggplot(h, aes(label, plab, fill = NES)) +
    geom_tile(color = "white") + geom_text(aes(label = star), vjust = .78) +
    scale_fill_gradient2(low = COL_NEG, mid = "white", high = COL_POS, midpoint = 0) +
    labs(title = sprintf("NES across contrasts - %s [%s]", co, mt),
         subtitle = "* FDR<0.05", x = NULL, y = NULL) +
    theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(FIG_DIR, sprintf("heatmap_%s_%s.pdf", co, mt)), g,
         width = 8.5, height = max(5, 0.22 * nlevels(h$plab) + 2))
}

# ================== (4) ENRICHMENT RUNNING-SCORE PLOTS =======================
# these need the ranked lists + gene sets, so recompute (fast, just 4 files)
pick <- function(cn, cand) { h <- cand[tolower(cand) %in% tolower(cn)]; if (length(h)) cn[match(tolower(h[1]), tolower(cn))] else NA_character_ }
read_dge <- function(path) {
  ext <- tolower(tools::file_ext(path)); d <- fread(path, sep = if (ext=="tsv") "\t" else ",", data.table = FALSE)
  cn <- names(d)
  g <- pick(cn, c("gene","primerid","symbol")); l <- pick(cn, c("log2FoldChange","coef","logFC"))
  a <- pick(cn, c("padj","fdr","FDR")); pv <- pick(cn, c("pvalue","Pr(>Chisq)","PValue"))
  pj <- if (!is.na(a)) d[[a]] else d[[pv]]
  data.frame(gene = as.character(d[[g]]), logFC = as.numeric(d[[l]]), padj = as.numeric(pj))
}
make_ranks <- function(df, metric) {
  df <- df[!is.na(df$gene) & !is.na(df$logFC), ]; df <- df[!duplicated(df$gene), ]
  if (metric == "logFC") v <- df$logFC else {
    p <- df$padj; p[is.na(p)] <- 1; nz <- suppressWarnings(min(p[p>0])); if (is.finite(nz)) p[p==0] <- nz
    v <- df$logFC * -log10(p)
  }
  names(v) <- df$gene; v <- v[is.finite(v)]; sort(v, decreasing = TRUE)
}
get_sets <- function(coll) {
  if (coll == "Hallmark") {
    df <- tryCatch(msigdbr(species="Homo sapiens", collection="H"),
                   error=function(e) msigdbr(species="Homo sapiens", category="H"))
  } else {
    df <- tryCatch(msigdbr(species="Homo sapiens", collection="C5", subcollection="GO:BP"),
                   error=function(e) msigdbr(species="Homo sapiens", category="C5", subcategory="BP"))
  }
  sym <- intersect(c("gene_symbol","human_gene_symbol"), names(df))[1]
  split(df[[sym]], df$gs_name)
}
sets <- list(Hallmark = get_sets("Hallmark"), GOBP = get_sets("GOBP"))

running_es <- function(stats, pathway) {
  stats <- sort(stats, decreasing = TRUE); hits <- names(stats) %in% pathway
  Nh <- sum(hits); n <- length(stats); if (Nh == 0) return(NULL)
  r <- abs(stats)
  list(curve = data.frame(rank = seq_len(n), ES = cumsum(ifelse(hits, r, 0))/sum(r[hits]) -
                            cumsum(ifelse(!hits, 1, 0))/(n - Nh)),
       ticks = which(hits), metric = as.numeric(stats))
}
one_enrichment <- function(ranks, geneset, ttl, nes, fdr) {
  e <- running_es(ranks, geneset); if (is.null(e)) return(NULL); n <- length(ranks)
  col <- if (!is.na(nes) && nes < 0) COL_NEG else COL_POS
  p1 <- ggplot(e$curve, aes(rank, ES)) + geom_hline(yintercept = 0, color = "grey70", linewidth = .3) +
    geom_line(color = col, linewidth = .7) +
    labs(y = "Running ES", title = ttl,
         subtitle = sprintf("NES=%.2f   FDR=%.2g", nes, fdr)) +
    theme_classic(base_size = 11) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p2 <- ggplot(data.frame(rank = e$ticks)) +
    geom_segment(aes(x = rank, xend = rank, y = 0, yend = 1), color = col, linewidth = .25) +
    scale_x_continuous(limits = c(1, n)) + theme_void()
  p3 <- ggplot(data.frame(rank = seq_len(n), m = e$metric), aes(rank, m)) +
    geom_area(fill = "grey75") + labs(x = "Rank in Ordered Dataset", y = "Ranked metric") +
    theme_classic(base_size = 11)
  p1 / p2 / p3 + plot_layout(heights = c(3, 0.7, 1.4))
}

for (lab in names(TARGETS)) {
  df <- read_dge(TARGETS[[lab]]$raw)
  pdf(file.path(FIG_DIR, sprintf("enrichment_%s.pdf", gsub("[^A-Za-z0-9]+","_",lab))),
      width = 8, height = 7, onefile = TRUE)
  for (mt in METRICS) {
    ranks <- make_ranks(df, mt)
    for (co in COLLECTIONS) {
      d <- res[res$collection==co & res$label==lab & res$metric==mt & is.finite(res$NES), ]
      if (!nrow(d)) next
      d <- d[order(d$padj), ][seq_len(min(TOP_ENRICH, nrow(d))), ]
      for (i in seq_len(nrow(d))) {
        gs <- sets[[co]][[d$pathway[i]]]; if (is.null(gs)) next
        pl <- one_enrichment(ranks, gs,
                             sprintf("%s [%s]: %s", lab, mt, prettify(d$pathway[i])), d$NES[i], d$padj[i])
        if (!is.null(pl)) print(pl)
      }
    }
  }
  dev.off()
}

message(sprintf("Done. Figures in %s/", FIG_DIR))

