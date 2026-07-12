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
