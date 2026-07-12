#!/usr/bin/env Rscript
# =============================================================================
# Cross-species signature GSEA
# -----------------------------------------------------------------------------
# Idea: turn each human DGE file into UP / DOWN gene "signatures" (padj<0.05),
#       then ask whether those signatures are enriched in the ranked gene lists
#       from mouse monocyte-subset DGE results.
#
#   pathways (gene sets) = human significant UP / DOWN genes  (human symbols)
#   stats    (ranked list) = mouse DGE per subset, mapped to human orthologs
#
# For every mouse subset the GSEA is run twice:
#       metric 1  "logFC"          : effect size
#       metric 2  "logFC_signedP"  : effect * -log10(padj)
#
# Mouse -> human orthologs via the {babelgene} package (offline, MSigDB/HGNC
# ortholog tables; no internet/biomaRt needed).
# =============================================================================

# ============================== CONFIG =======================================
# --- Human DGE files -> become UP/DOWN signatures (padj or fdr < SIG_CUTOFF).
#     Format is auto-detected (DESeq2 .tsv or MAST .csv).
HUMAN_FILES <- c(
  "C:/Users/fc809/Downloads/DGE_tables/DGE_tables/DESeq2_Cohort_CCUS_CHIP_vs_Control/CD14+_Monocytes_TET2_CHIP_vs_Control.tsv",
  "C:/Users/fc809/Downloads/DGE_tables/DGE_tables/DESeq2_Cohort_CH_vs_Control/CD14+_Monocytes_TET2_vs_none.tsv",
  "C:/Users/fc809/Downloads/DGE_tables/DGE_tables/MAST_DE_Mutant_vs_Control/TET2 CD14 Mutant vs Control.csv",
  "C:/Users/fc809/Downloads/DGE_tables/DGE_tables/MAST_DE_Mutant_vs_Wildtype/TET2 CD14 Mutant vs Wildtype.csv"
)

# --- Mouse DGE files -> ranked lists (one per monocyte subset). EDIT THESE.
#     names() are used as subset labels in the output.
MOUSE_FILES <- c(
  subset1 = "C:/Users/fc809/Downloads/DGE_tables/mouse/SUBSET1.csv",
  subset2 = "C:/Users/fc809/Downloads/DGE_tables/mouse/SUBSET2.csv",
  subset3 = "C:/Users/fc809/Downloads/DGE_tables/mouse/SUBSET3.csv"
)

OUTPUT_DIR <- "C:/Users/fc809/Downloads/DGE_tables/GSEA_signature_results"

SIG_CUTOFF <- 0.05      # padj/fdr threshold to call a human gene "significant"
MIN_SIZE   <- 5         # min overlap between signature and ranked list
MAX_SIZE   <- 5000      # effectively uncapped for custom signatures
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  if (!requireNamespace("fgsea", quietly = TRUE))        BiocManager::install("fgsea", update = FALSE, ask = FALSE)
  if (!requireNamespace("babelgene", quietly = TRUE))    install.packages("babelgene", repos = "https://cloud.r-project.org")
  if (!requireNamespace("data.table", quietly = TRUE))   install.packages("data.table", repos = "https://cloud.r-project.org")
  library(fgsea); library(babelgene); library(data.table)
})
set.seed(42)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------- flexible column detection ---------------------------
pick <- function(cn, candidates) {
  hit <- candidates[tolower(candidates) %in% tolower(cn)]
  if (length(hit)) cn[match(tolower(hit[1]), tolower(cn))] else NA_character_
}

# Read any DGE table -> data.frame(gene, logFC, padj). Works for DESeq2, MAST,
# edgeR, limma, and Seurat::FindMarkers style outputs.
read_dge <- function(path) {
  ext <- tolower(tools::file_ext(path))
  d <- fread(path, sep = if (ext == "tsv") "\t" else ",", data.table = FALSE)
  cn <- names(d)

  gcol <- pick(cn, c("gene","Gene","symbol","SYMBOL","primerid","gene_symbol","names","X","V1"))
  # Seurat/edgeR often store genes as rownames -> first unnamed/`V1`/`X` column
  if (is.na(gcol) && (cn[1] %in% c("V1","X",""))) gcol <- cn[1]
  lcol <- pick(cn, c("log2FoldChange","avg_log2FC","logFC","coef","avg_logFC","log2FC","logfc"))
  acol <- pick(cn, c("padj","p_val_adj","FDR","fdr","adj.P.Val","qvalue","q_value","padjust"))
  pcol <- pick(cn, c("pvalue","p_val","PValue","P.Value","Pr(>Chisq)","pval"))

  if (is.na(gcol) || is.na(lcol))
    stop(sprintf("Could not find gene/logFC columns in %s\n  columns: %s",
                 basename(path), paste(cn, collapse = ", ")))
  padj <- if (!is.na(acol)) d[[acol]] else if (!is.na(pcol)) d[[pcol]] else NA
  data.frame(gene  = as.character(d[[gcol]]),
             logFC = as.numeric(d[[lcol]]),
             padj  = as.numeric(padj),
             stringsAsFactors = FALSE)
}

# ------------------ build UP / DOWN human signatures -------------------------
message("Building human UP/DOWN signatures (", SIG_CUTOFF, " cutoff) ...")
signatures <- list()
for (path in HUMAN_FILES) {
  d <- read_dge(path)
  tag <- gsub("[/ ]+", "_", tools::file_path_sans_ext(basename(path)))
  sig <- d[!is.na(d$padj) & d$padj < SIG_CUTOFF & !is.na(d$logFC), ]
  up   <- unique(sig$gene[sig$logFC > 0])
  down <- unique(sig$gene[sig$logFC < 0])
  signatures[[paste0(tag, "__UP")]]   <- up
  signatures[[paste0(tag, "__DOWN")]] <- down
  message(sprintf("  %-55s UP=%d DOWN=%d", tag, length(up), length(down)))
}
signatures <- signatures[lengths(signatures) > 0]

# ----------- map a mouse ranked list to human-ortholog symbols ---------------
mouse_to_human_ranks <- function(d, metric) {
  d <- d[!is.na(d$gene) & !is.na(d$logFC), , drop = FALSE]

  if (metric == "logFC") {
    val <- d$logFC
  } else {                       # logFC * -log10(padj)
    p <- d$padj; p[is.na(p)] <- 1
    nz <- suppressWarnings(min(p[p > 0]))
    if (is.finite(nz)) p[p == 0] <- nz
    val <- d$logFC * -log10(p)
  }
  d$val <- val
  d <- d[is.finite(d$val), ]

  # mouse symbol -> human ortholog (human = FALSE: input genes are the model species)
  orth <- babelgene::orthologs(genes = unique(d$gene), species = "mouse", human = FALSE)
  map  <- orth[, c("symbol", "human_symbol")]          # symbol = mouse, human_symbol = human
  m <- merge(d[, c("gene","val")], map, by.x = "gene", by.y = "symbol")
  if (!nrow(m)) return(setNames(numeric(0), character(0)))

  # collapse many mouse genes -> one human gene: keep the strongest (max |val|)
  m <- m[order(-abs(m$val)), ]
  m <- m[!duplicated(m$human_symbol), ]
  v <- setNames(m$val, m$human_symbol)
  sort(v, decreasing = TRUE)
}

# --------------------------------- run ---------------------------------------
metrics <- c("logFC", "logFC_signedP")
all_res <- list()

for (i in seq_along(MOUSE_FILES)) {
  label <- names(MOUSE_FILES)[i]
  path  <- MOUSE_FILES[i]
  if (!file.exists(path)) { warning(sprintf("Mouse file not found, skipping: %s", path)); next }
  d <- read_dge(path)

  for (metric in metrics) {
    ranks <- mouse_to_human_ranks(d, metric)
    message(sprintf("Subset %-10s | %-13s | %d human-mapped genes", label, metric, length(ranks)))
    if (length(ranks) < MIN_SIZE) { warning("  too few mapped genes; skipping."); next }

    res <- fgsea(pathways = signatures, stats = ranks,
                 minSize = MIN_SIZE, maxSize = MAX_SIZE, eps = 0)
    res <- res[order(res$padj), ]
    res$leadingEdge <- vapply(res$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
    res$mouse_subset <- label
    res$metric       <- metric
    fwrite(res, file.path(OUTPUT_DIR, sprintf("%s__%s.csv", label, metric)))
    all_res[[length(all_res) + 1]] <- res
  }
}

if (length(all_res)) {
  combined <- rbindlist(all_res, use.names = TRUE, fill = TRUE)
  setcolorder(combined, c("mouse_subset","metric","pathway","pval","padj","NES","ES","size"))
  fwrite(combined, file.path(OUTPUT_DIR, "signature_GSEA_all_results.csv"))
  message(sprintf("\nDone. Tested %d signatures across %d subset x metric runs. Output: %s/",
                  length(signatures), length(all_res), OUTPUT_DIR))
} else {
  message("No results — check that MOUSE_FILES paths are set correctly.")
}
