## ============================================================
## GSE207063 — BM Ly6Chi Classical Monocyte Aging
## Barman et al. 2022, Aging Cell
## MALE MICE ONLY: Young (2-6mo) vs Old (24-30mo)
##
## Cross-dataset gene sets are built using ALL genes with FDR < 0.05,
## split by direction (instead of top/bottom 200 by effect size).
##
## Steps:
##   1.  Download from GEO
##   2.  Load + filter male samples
##   3.  Seurat QC + clustering + UMAP
##   4.  Pseudobulk aggregation
##   5.  edgeR quasi-likelihood DE (Old vs Young)
##   6.  Volcano plot
##   7.  GSEA — GO BP (gseGO)
##   8.  GSEA — KEGG (gseKEGG)
##   9.  GSEA — Hallmarks (GSEA + msigdbr)
##   10. Cross-dataset GSEA (aging <-> HB/LB / Lenti monocyte contrasts)
##   11. H2B-style pathway analysis on Old-vs-Young ranked list
##   12. Plots + console summary
##
## All outputs, including a copy of this script, are written to OUTDIR.
## ============================================================

## ── PACKAGES ──────────────────────────────────────────────────────────────────
## Install if needed:
## install.packages(c("BiocManager","ggrepel","stringr"))
## BiocManager::install(c("GEOquery","Seurat","edgeR",
##                        "clusterProfiler","org.Mm.eg.db","enrichplot","msigdbr"))

# NOTE on package load order:
# dplyr is loaded LAST so its verbs (select, filter, mutate, arrange, rename)
# win masking conflicts with Bioconductor packages (AnnotationDbi::select,
# plyr::rename/mutate, etc.). Don't reorder without reason.

library(GEOquery)
library(Seurat)
library(edgeR)
library(tibble)
library(ggplot2)
library(ggrepel)
library(stringr)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(msigdbr)
library(dplyr)  # must be last

set.seed(42)

OUTDIR           <- "C:/Users/fc809/Downloads/AgingCell_analysis_wFDR"
PADJ_CUTOFF      <- 0.05   # stringent cutoff (used by pathway-analysis displays)
RANK_PADJ_CUTOFF <- 0.25   # looser cutoff for building cross-dataset gene sets
# (Old_up/Old_down and HB_up/HB_down)

# Path to this script — update if you rename/move the file.
# At the end of the run, the script is copied into OUTDIR for provenance.
SCRIPT_PATH <- "C:/Users/fc809/Downloads/AgingCell_analysis.R"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
setwd(OUTDIR)

## ── LOAD MONOCYTE DEG TABLES (from separate edgeR pseudobulk pipeline) ───────
## These four tables are used by the cross-dataset GSEA blocks further down.
## The .txt files are tab-delimited edgeR QLF outputs; expected columns are
## X (gene), logFC, and some FDR column (FDR / padj / adj.P.Val).

DEG_ROOT <- "C:/Users/fc809/Downloads/Classical_Monocytes (1)"

deg_files <- c(
  your_mono_lentipos = file.path(DEG_ROOT, "LentiPosHB_vs_LentiPosLB",
                                 "DEG_LentiPosHB_vs_LentiPosLB_pseudoBulk_res0.6.sce_0.txt"),
  your_mono_lentineg = file.path(DEG_ROOT, "LentiNegHB_vs_LentiNegLB",
                                 "DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.6.sce_0.txt"),
  lenti_hb           = file.path(DEG_ROOT, "LentiPos_vs_LentiNeg_in_HB",
                                 "DEG_LentiPos_vs_LentiNeg_in_HB_pseudoBulk_res0.6.sce_0.txt"),
  lenti_lb           = file.path(DEG_ROOT, "LentiPos_vs_LentiNeg_in_LB",
                                 "DEG_LentiPos_vs_LentiNeg_in_LB_pseudoBulk_res0.6.sce_0.txt")
)

load_deg <- function(path) {
  if (!file.exists(path)) stop("DEG file not found: ", path)
  # Try tab-delimited first, fall back to CSV
  df <- tryCatch(
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e)
      read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  )
  # Repair empty / NA column names. R's write.table with row.names=TRUE
  # produces a header line missing the first field, which read.delim keeps
  # as "". First unnamed column is almost always the gene symbol -> "X".
  bad <- is.na(names(df)) | names(df) == ""
  if (any(bad)) {
    if (bad[1]) {
      names(df)[1] <- "X"
      bad[1] <- FALSE
    }
    # Any remaining empties get generic names
    names(df)[bad] <- paste0("V", seq_len(sum(bad)))
  }
  # Normalise gene column to "X" (that's what the cross-dataset code expects)
  if (!"X" %in% colnames(df)) {
    gc <- intersect(c("gene", "Gene", "symbol", "Symbol", "gene_name"),
                    colnames(df))
    if (length(gc) > 0) {
      names(df)[names(df) == gc[1]] <- "X"
    } else if (!identical(rownames(df), as.character(seq_len(nrow(df))))) {
      df$X <- rownames(df)
    }
  }
  df
}

for (nm in names(deg_files)) assign(nm, load_deg(deg_files[[nm]]))

cat("\nLoaded monocyte DEG tables:\n")
for (nm in names(deg_files)) {
  d <- get(nm)
  cat("  ", nm, ":", nrow(d), "rows  |  cols:",
      paste(head(colnames(d), 8), collapse = ", "), "\n")
}
cat("\n")

## ── 1. DOWNLOAD FROM GEO ─────────────────────────────────────────────────────
message("=== Step 1: Downloading GSE207063 ===")

gse   <- getGEO("GSE207063", GSEMatrix = TRUE, AnnotGPL = FALSE)
pheno <- pData(gse[[1]])
message("Samples found:")
print(pheno[, c("title", "geo_accession")])

getGEOSuppFiles("GSE207063", baseDir = OUTDIR)

# Extract any tarballs
tar_files <- list.files(file.path(OUTDIR, "GSE207063"),
                        pattern = "\\.tar$", full.names = TRUE)
for (f in tar_files) {
  message("Extracting: ", basename(f))
  untar(f, exdir = file.path(OUTDIR, "GSE207063"))
}

message("Downloaded files:")
print(basename(list.files(file.path(OUTDIR, "GSE207063"),
                          recursive = TRUE, full.names = FALSE)))

## ── 2. LOCATE + LOAD MALE SAMPLES ─────────────────────────────────────────────
message("=== Step 2: Loading male samples ===")

# Find valid 10x directories (matrix + barcodes + features)
find_10x_dirs <- function(base_dir) {
  Filter(function(d) {
    f <- list.files(d)
    any(grepl("matrix\\.mtx", f)) &&
      any(grepl("barcodes\\.tsv", f)) &&
      any(grepl("features\\.tsv|genes\\.tsv", f))
  }, list.dirs(base_dir, recursive = TRUE))
}

is_male <- function(nm) {
  grepl("male|_m_|_m$|^m_", tolower(nm)) &
    !grepl("female|_f_|_f$|^f_", tolower(nm))
}

sample_dirs <- find_10x_dirs(file.path(OUTDIR, "GSE207063"))
h5_files    <- list.files(file.path(OUTDIR, "GSE207063"),
                          pattern = "\\.h5$", recursive = TRUE, full.names = TRUE)

if (length(sample_dirs) > 0) {
  paths     <- sample_dirs[is_male(basename(sample_dirs))]
  names_vec <- basename(paths)
} else {
  paths     <- h5_files[is_male(basename(h5_files))]
  names_vec <- sub("\\.h5$", "", basename(paths))
}

message("Male sample paths found: ", length(paths))
print(names_vec)


## ── LOAD MALE DATA ────────────────────────────────────────────────────────────
male_dir <- file.path(OUTDIR, "male_10x")
dir.create(male_dir, showWarnings = FALSE)

GEO_DIR <- file.path(OUTDIR, "GSE207063")

# Copy (or overwrite) the gz files with the names Read10X expects
copy_map <- c(
  "GSE207063_Male_Young_Old_barcodes.tsv.gz" = "barcodes.tsv.gz",
  "GSE207063_Male_Young_Old_features.tsv.gz" = "features.tsv.gz",
  "GSE207063_Male_Young_Old_matrix.mtx.gz"   = "matrix.mtx.gz"
)

for (src_name in names(copy_map)) {
  src <- file.path(GEO_DIR,  src_name)
  dst <- file.path(male_dir, copy_map[[src_name]])
  if (!file.exists(src)) stop("Missing GEO file: ", src)
  ok <- file.copy(src, dst, overwrite = TRUE)
  message("Copied ", src_name, " -> ", copy_map[[src_name]],
          "  (success=", ok, ", size=", file.info(dst)$size, " bytes)")
}

# Sanity check all three files are in place and non-empty before Read10X
for (f in copy_map) {
  p <- file.path(male_dir, f)
  if (!file.exists(p) || file.info(p)$size == 0)
    stop("male_10x file missing/empty: ", p)
}

# Load with Read10X — returns a LIST with "Gene Expression" + "Antibody Capture" (HTOs)
counts <- Read10X(male_dir)
stopifnot(is.list(counts),
          "Gene Expression"  %in% names(counts),
          "Antibody Capture" %in% names(counts))

gex_mat <- counts[["Gene Expression"]]
hto_mat <- counts[["Antibody Capture"]]

cat("Gene Expression: ", nrow(gex_mat), "features x", ncol(gex_mat), "cells\n")
cat("Antibody Capture:", nrow(hto_mat), "features x", ncol(hto_mat), "cells\n")

# ── Build Seurat object with LIGHT pre-filter so HTODemux sees the max cell pool ─
# min.features is set to 50 (not 200) — real monocytes with low RNA yield can
# still have perfectly good HTO signal. We apply the strict RNA QC later,
# after demux, so as not to throw away cells pre-emptively.
so <- CreateSeuratObject(counts = gex_mat, min.cells = 3, min.features = 50)
cat("After CreateSeuratObject (min.features=50):", ncol(so), "cells\n")

barcode_suffixes <- sub(".*-", "-", colnames(so))
print(table(barcode_suffixes))

# Attach HTO assay — keep only cells that passed the light RNA filter
hto_shared <- hto_mat[, colnames(so), drop = FALSE]
so[["HTO"]] <- CreateAssayObject(counts = hto_shared)

# Normalize HTO + demultiplex — guarded so re-running the block doesn't
# double-CLR-normalize and break HTODemux thresholds.
# positive.quantile = 0.95 is looser than Seurat's default 0.99 and gives
# far fewer false-negative singlets; bump to 0.99 only if you see spurious
# doublets or want to be conservative.
if (!"HTO_classification.global" %in% colnames(so[[]])) {
  so <- NormalizeData(so, assay = "HTO", normalization.method = "CLR")
  so <- HTODemux(so, assay = "HTO", positive.quantile = 0.95)
} else {
  message("HTO already demuxed on this object — skipping NormalizeData/HTODemux")
}

cat("\nHTODemux classification:\n")
print(table(so$HTO_classification.global))

# Keep only singlets
if (!identical(unique(as.character(so$HTO_classification.global)), "Singlet")) {
  so <- subset(so, subset = HTO_classification.global == "Singlet")
}
message("Singlets kept: ", ncol(so))

# Assign age group and sample ID from hashtag
hashtag_to_age <- c(
  Hashtag1 = "Young", Hashtag2 = "Young", Hashtag3 = "Young",
  Hashtag4 = "Young", Hashtag5 = "Young",
  Hashtag6 = "Old",   Hashtag7 = "Old",   Hashtag8 = "Old",
  Hashtag9 = "Old",   Hashtag10 = "Old"
)

hashtag_to_sample <- c(
  Hashtag1 = "M_Young_1", Hashtag2 = "M_Young_2", Hashtag3 = "M_Young_3",
  Hashtag4 = "M_Young_4", Hashtag5 = "M_Young_5",
  Hashtag6 = "M_Old_1",   Hashtag7 = "M_Old_2",   Hashtag8 = "M_Old_3",
  Hashtag9 = "M_Old_4",   Hashtag10 = "M_Old_5"
)

head(so$hash.ID)
class(so$hash.ID)

meta <- data.frame(
  age_group = hashtag_to_age[as.character(so$hash.ID)],
  sample    = hashtag_to_sample[as.character(so$hash.ID)],
  row.names = colnames(so)
)

so <- AddMetaData(so, meta)
print(table(so$sample, so$age_group))

# ── QC + clustering ──────────────────────────────────────────────────────────
so[["pct.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")

cat("\nBefore RNA QC — cells per sample:\n")
print(table(so$sample))
cat("nFeature_RNA summary:\n"); print(summary(so$nFeature_RNA))
cat("pct.mt summary:\n");       print(summary(so$pct.mt))

so <- subset(so, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & pct.mt < 20)
message("Cells after QC: ", ncol(so))

cat("\nAfter RNA QC — cells per sample:\n")
print(table(so$sample, so$age_group))

so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
so <- ScaleData(so, features = rownames(so), verbose = FALSE)
so <- RunPCA(so, npcs = 30, verbose = FALSE)
so <- RunUMAP(so, dims = 1:20, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:20, verbose = FALSE)
so <- FindClusters(so, resolution = 0.4, verbose = FALSE)

saveRDS(so, file.path(OUTDIR, "GSE207063_seurat_male.rds"))

# UMAP
p1 <- DimPlot(so, group.by = "age_group",
              cols = c("Old" = "#C0392B", "Young" = "#2980B9"), pt.size = 0.6) +
  ggtitle("Age Group") + theme_classic()
p2 <- DimPlot(so, group.by = "seurat_clusters", label = TRUE,
              label.size = 4, pt.size = 0.6) +
  ggtitle("Clusters") + theme_classic()

png(file.path(OUTDIR, "GSE207063_umap_male.png"), width = 1400, height = 600, res = 120)
print(p1 | p2)
dev.off()
message("UMAP saved.")

# ── Pseudobulk aggregation ────────────────────────────────────────────────────
counts_mat <- GetAssayData(so, layer = "counts")
sample_ids <- unique(so$sample)

agg_mat <- sapply(sample_ids, function(s) {
  Matrix::rowSums(counts_mat[, so$sample == s, drop = FALSE])
})
colnames(agg_mat) <- sample_ids

coldata <- data.frame(
  sample    = sample_ids,
  age_group = factor(ifelse(grepl("Old", sample_ids), "Old", "Young"),
                     levels = c("Young", "Old")),
  row.names = sample_ids
)

keep <- rowSums(agg_mat >= 10) >= 4   # at least 4 of 10 samples
agg_mat <- agg_mat[keep, ]
message("Genes after filtering: ", nrow(agg_mat))

# ── edgeR QL ─────────────────────────────────────────────────────────────────
y      <- DGEList(counts = round(agg_mat), group = coldata$age_group)
y      <- calcNormFactors(y, method = "TMM")
design <- model.matrix(~ age_group, data = coldata)
y      <- estimateDisp(y, design, robust = TRUE)
fit    <- glmQLFit(y, design, robust = TRUE)
qlf    <- glmQLFTest(fit, coef = "age_groupOld")

res_df <- topTags(qlf, n = Inf, sort.by = "PValue")$table %>%
  rownames_to_column("gene") %>%
  dplyr::rename(log2FC = logFC, pval = PValue) %>%
  dplyr::mutate(direction = case_when(
    FDR < 0.05 & log2FC >  0.5 ~ "Up_Old",
    FDR < 0.05 & log2FC < -0.5 ~ "Up_Young",
    TRUE ~ "NS"
  ))

message("DE summary:")
print(table(res_df$direction))

write.csv(res_df, file.path(OUTDIR, "GSE207063_DEGs_male_OldVsYoung.csv"), row.names = FALSE)

###PATHWAY_ANALYSIS####

# (Packages already loaded at top of script — do NOT re-library here,
# it re-masks dplyr verbs and breaks downstream select/filter/mutate.)

# ── Gene ID mapping ───────────────────────────────────────────────────────────
entrez_map <- bitr(res_df$gene,
                   fromType = "SYMBOL",
                   toType   = "ENTREZID",
                   OrgDb    = org.Mm.eg.db)

res_annot <- res_df %>%
  left_join(entrez_map, by = c("gene" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID), !is.na(pval)) %>%
  distinct(gene, .keep_all = TRUE)

message("Mapped ", nrow(res_annot), " / ", nrow(res_df), " genes")

# ── Ranked gene list ──────────────────────────────────────────────────────────
ranked_vec <- res_annot %>%
  arrange(desc(abs(log2FC))) %>%
  distinct(ENTREZID, .keep_all = TRUE) %>%
  arrange(desc(log2FC)) %>%
  { setNames(.$log2FC, .$ENTREZID) }

message("Ranked list: ", length(ranked_vec), " genes")

# ── GO BP ─────────────────────────────────────────────────────────────────────
gsea_go <- gseGO(
  geneList      = ranked_vec,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  eps           = 0
)
gsea_go_simp <- simplify(gsea_go, cutoff = 0.7, by = "p.adjust", select_fun = min)
gsea_go_df   <- as.data.frame(gsea_go_simp) %>%
  mutate(direction = ifelse(NES > 0, "Enriched_Old", "Enriched_Young")) %>%
  arrange(p.adjust)
message("GO BP significant: ", nrow(gsea_go_df))
write.csv(gsea_go_df, file.path(OUTDIR, "GSE207063_GSEA_GO_BP.csv"), row.names = FALSE)

# ── KEGG ──────────────────────────────────────────────────────────────────────
gsea_kegg <- gseKEGG(
  geneList      = ranked_vec,
  organism      = "mmu",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  eps           = 0
)
gsea_kegg    <- setReadable(gsea_kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
gsea_kegg_df <- as.data.frame(gsea_kegg) %>%
  mutate(direction = ifelse(NES > 0, "Enriched_Old", "Enriched_Young")) %>%
  arrange(p.adjust)
message("KEGG significant: ", nrow(gsea_kegg_df))
write.csv(gsea_kegg_df, file.path(OUTDIR, "GSE207063_GSEA_KEGG.csv"), row.names = FALSE)

# ── Hallmarks ─────────────────────────────────────────────────────────────────

h_raw <- msigdbr(species = "Mus musculus", category = "H")

h_sets <- h_raw %>%
  dplyr::select(gs_name, ncbi_gene) %>%
  dplyr::mutate(ncbi_gene = as.character(ncbi_gene))

gsea_hallmarks <- GSEA(
  geneList      = ranked_vec,
  TERM2GENE     = h_sets,
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0
)
gsea_hallmarks <- setReadable(gsea_hallmarks, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
gsea_h_df <- as.data.frame(gsea_hallmarks) %>%
  mutate(direction   = ifelse(NES > 0, "Enriched_Old", "Enriched_Young"),
         Description = gsub("HALLMARK_|_", " ", ID)) %>%
  arrange(p.adjust)
message("Hallmarks significant (padj<0.05): ", sum(gsea_h_df$p.adjust < 0.05))
write.csv(gsea_h_df, file.path(OUTDIR, "GSE207063_GSEA_Hallmarks.csv"), row.names = FALSE)

message("\n--- Hallmarks significant ---")
gsea_h_df %>% filter(p.adjust < 0.05) %>% arrange(NES) %>%
  dplyr::select(Description, NES, p.adjust, setSize) %>% print()

# ── Quick console summary ─────────────────────────────────────────────────────
message("\n--- GO BP top 10 ---")
gsea_go_df %>% dplyr::select(Description, NES, p.adjust, setSize) %>% head(10) %>% print()

message("\n--- KEGG all significant ---")
gsea_kegg_df %>% arrange(NES) %>% dplyr::select(Description, NES, p.adjust, setSize) %>% print()

message("\n--- Hallmarks significant ---")
gsea_h_df %>% filter(p.adjust < 0.05) %>% arrange(NES) %>%
  dplyr::select(Description, NES, p.adjust, setSize) %>% print()

# ── Check what we have ────────────────────────────────────────────────────────
message("GO BP significant: ", nrow(gsea_go_df))
message("KEGG significant:  ", nrow(gsea_kegg_df))

# ── GO BP bar plot ────────────────────────────────────────────────────────────
if (nrow(gsea_go_df) > 0) {
  plot_df <- gsea_go_df %>%
    slice_max(abs(NES), n = 30) %>%
    arrange(NES) %>%
    mutate(
      label = str_wrap(Description, width = 45),
      label = factor(label, levels = label)
    )
  
  p_go_bar <- ggplot(plot_df, aes(x = NES, y = label, fill = direction)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = c("Enriched_Old" = "#C0392B", "Enriched_Young" = "#2980B9"),
                      labels = c("Enriched_Old" = "Old", "Enriched_Young" = "Young")) +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey40") +
    labs(title    = "GO BP GSEA — Old vs Young BM Monocytes",
         subtitle = paste0("Top 30 of ", nrow(gsea_go_df), " significant terms"),
         x = "NES", y = NULL, fill = NULL) +
    theme_classic(base_size = 11) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = 8))
  
  ggsave(file.path(OUTDIR, "GSE207063_GSEA_bar_GOBP.png"),
         p_go_bar, width = 11,
         height = max(6, 0.38 * nrow(plot_df) + 2),
         dpi = 180, limitsize = FALSE)
  message("Saved: GSE207063_GSEA_bar_GOBP.png")
  
  # Dotplot
  p_go_dot <- dotplot(gsea_go_simp, showCategory = 20, split = ".sign", font.size = 9) +
    facet_grid(~.sign) +
    ggtitle("GO BP — Old vs Young BM Monocytes") +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(OUTDIR, "GSE207063_GSEA_dotplot_GOBP.png"),
         p_go_dot, width = 12, height = 9, dpi = 180)
  message("Saved: GSE207063_GSEA_dotplot_GOBP.png")
  
  # Enrichment map
  tryCatch({
    gsea_go_net <- pairwise_termsim(gsea_go_simp)
    p_emap <- emapplot(gsea_go_net, showCategory = 30, color = "NES") +
      ggtitle("GO BP enrichment map — Old vs Young")
    ggsave(file.path(OUTDIR, "GSE207063_GSEA_emap_GOBP.png"),
           p_emap, width = 12, height = 10, dpi = 180)
    message("Saved: GSE207063_GSEA_emap_GOBP.png")
  }, error = function(e) message("emapplot skipped: ", conditionMessage(e)))
  
  # Enrichment curves — top 6
  top_go_ids <- gsea_go_df %>% slice_max(abs(NES), n = 6) %>% pull(ID)
  png(file.path(OUTDIR, "GSE207063_GSEA_curves_GOBP.png"),
      width = 1800, height = 1200, res = 130)
  print(gseaplot2(gsea_go_simp, geneSetID = top_go_ids,
                  title = "Top GO BP gene sets", base_size = 10))
  dev.off()
  message("Saved: GSE207063_GSEA_curves_GOBP.png")
  
} else message("No significant GO BP terms — skipping plots")

# ── KEGG bar plot ─────────────────────────────────────────────────────────────
if (nrow(gsea_kegg_df) > 0) {
  plot_df_k <- gsea_kegg_df %>%
    arrange(NES) %>%
    mutate(
      label = str_wrap(Description, width = 45),
      label = factor(label, levels = label)
    )
  
  p_kegg_bar <- ggplot(plot_df_k, aes(x = NES, y = label, fill = direction)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = c("Enriched_Old" = "#C0392B", "Enriched_Young" = "#2980B9"),
                      labels = c("Enriched_Old" = "Old", "Enriched_Young" = "Young")) +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey40") +
    labs(title = "KEGG GSEA — Old vs Young BM Monocytes",
         x = "NES", y = NULL, fill = NULL) +
    theme_classic(base_size = 11) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = 8))
  
  ggsave(file.path(OUTDIR, "GSE207063_GSEA_bar_KEGG.png"),
         p_kegg_bar, width = 11,
         height = max(6, 0.38 * nrow(plot_df_k) + 2),
         dpi = 180, limitsize = FALSE)
  message("Saved: GSE207063_GSEA_bar_KEGG.png")
  
  # Dotplot
  p_kegg_dot <- dotplot(gsea_kegg, showCategory = 20, split = ".sign", font.size = 9) +
    facet_grid(~.sign) +
    ggtitle("KEGG — Old vs Young BM Monocytes") +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(OUTDIR, "GSE207063_GSEA_dotplot_KEGG.png"),
         p_kegg_dot, width = 12, height = 9, dpi = 180)
  message("Saved: GSE207063_GSEA_dotplot_KEGG.png")
  
  # Enrichment curves — top 6
  top_kegg_ids <- gsea_kegg_df %>% slice_max(abs(NES), n = 6) %>% pull(ID)
  png(file.path(OUTDIR, "GSE207063_GSEA_curves_KEGG.png"),
      width = 1800, height = 1200, res = 130)
  print(gseaplot2(gsea_kegg, geneSetID = top_kegg_ids,
                  title = "Top KEGG pathways", base_size = 10))
  dev.off()
  message("Saved: GSE207063_GSEA_curves_KEGG.png")
  
} else message("No significant KEGG pathways — skipping plots")

message("\nAll done. Outputs in: ", OUTDIR)

# ── Hallmarks bar plot ────────────────────────────────────────────────────────
sig_h <- gsea_h_df %>% filter(p.adjust < 0.05) %>% arrange(NES) %>%
  mutate(
    label = gsub("HALLMARK_", "", ID),
    label = gsub("_", " ", label),
    label = str_wrap(label, width = 40),
    label = factor(label, levels = label)
  )

p_bar <- ggplot(sig_h, aes(x = NES, y = label, fill = direction)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = c("Enriched_Old" = "#C0392B", "Enriched_Young" = "#2980B9"),
                    labels = c("Enriched_Old" = "Old", "Enriched_Young" = "Young")) +
  geom_text(aes(label = formatC(p.adjust, format = "e", digits = 1)),
            hjust = ifelse(sig_h$NES > 0, -0.1, 1.1), size = 2.8) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey40") +
  xlim(min(sig_h$NES) * 1.35, max(sig_h$NES) * 1.35) +
  labs(title    = "Hallmarks GSEA — Old vs Young BM Monocytes",
       subtitle = "Pseudobulk edgeR | ranked by log2FC",
       x = "Normalized Enrichment Score (NES)", y = NULL, fill = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 9))

ggsave(file.path(OUTDIR, "GSE207063_GSEA_bar_Hallmarks.png"),
       p_bar, width = 10, height = 7, dpi = 200)
message("Saved: GSE207063_GSEA_bar_Hallmarks.png")

# ── Dotplot ───────────────────────────────────────────────────────────────────
p_dot <- dotplot(gsea_hallmarks, showCategory = 20, split = ".sign", font.size = 9) +
  facet_grid(~.sign) +
  ggtitle("Hallmarks — Old vs Young BM Monocytes") +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(OUTDIR, "GSE207063_GSEA_dotplot_Hallmarks.png"),
       p_dot, width = 12, height = 8, dpi = 180)
message("Saved: GSE207063_GSEA_dotplot_Hallmarks.png")

# ── Enrichment curves for top 6 by |NES| ─────────────────────────────────────
top_ids <- gsea_h_df %>% filter(p.adjust < 0.05) %>%
  slice_max(abs(NES), n = 6) %>% pull(ID)

png(file.path(OUTDIR, "GSE207063_GSEA_curves_Hallmarks.png"),
    width = 1800, height = 1200, res = 130)
print(gseaplot2(gsea_hallmarks, geneSetID = top_ids,
                title = "Top Hallmark gene sets — Old vs Young", base_size = 10))
dev.off()
message("Saved: GSE207063_GSEA_curves_Hallmarks.png")

# (No library() reloads needed — all packages are loaded at the top.)

# ── Build gene sets from Old vs Young (GSE207063) ─────────────────────────────
# Use ALL genes with FDR < RANK_PADJ_CUTOFF, split by direction:
#   Old_up   = FDR < RANK_PADJ_CUTOFF AND higher in Old   (log2FC > 0)
#   Old_down = FDR < RANK_PADJ_CUTOFF AND higher in Young (log2FC < 0)

aging_sig <- res_df %>%
  dplyr::filter(!is.na(FDR), !is.na(log2FC), FDR < RANK_PADJ_CUTOFF) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

old_up <- aging_sig %>%
  dplyr::filter(log2FC > 0) %>%
  dplyr::arrange(desc(log2FC)) %>%
  dplyr::pull(gene)

old_down <- aging_sig %>%
  dplyr::filter(log2FC < 0) %>%
  dplyr::arrange(log2FC) %>%
  dplyr::pull(gene)

N_OLD_UP   <- length(old_up)
N_OLD_DOWN <- length(old_down)

cat("Aging signature (GSE207063) — FDR <", RANK_PADJ_CUTOFF, ":\n")
cat("  Total significant genes :", nrow(aging_sig), "\n")
cat("  Old_up   (higher in Old)   :", N_OLD_UP,   "genes\n")
cat("  Old_down (higher in Young) :", N_OLD_DOWN, "genes\n\n")

aging_sets <- data.frame(
  term = c(rep("Old_up",   length(old_up)),
           rep("Old_down", length(old_down))),
  gene = c(old_up, old_down)
)

# ── Cross-dataset GSEA: aging signature on HB-vs-LB ranked lists ─────────────
# Run for BOTH LentiNeg and LentiPos HB-vs-LB contrasts (mirrors H2B pipeline).

run_aging_on_hblb <- function(deg_tbl, tag) {
  ranked_hb <- deg_tbl %>%
    dplyr::arrange(desc(logFC)) %>%
    dplyr::distinct(X, .keep_all = TRUE) %>%
    { setNames(.$logFC, .$X) }
  
  cat("\n[", tag, "] Ranked HB vs LB list:", length(ranked_hb), "genes\n")
  
  gsea_res <- GSEA(
    geneList      = ranked_hb,
    TERM2GENE     = aging_sets,
    minGSSize     = 10,
    maxGSSize     = max(5000, N_OLD_UP, N_OLD_DOWN) + 1,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    eps           = 0,
    verbose       = TRUE
  )
  
  res_stats <- as.data.frame(gsea_res)[, c("ID", "NES", "pvalue", "p.adjust")]
  cat("\n--- [", tag, "] aging sets on HB-vs-LB ---\n")
  print(res_stats)
  
  p <- gseaplot2(gsea_res,
                 geneSetID = c("Old_up", "Old_down"),
                 title     = paste0("Aging signature enrichment in ", tag,
                                    " HB vs LB monocytes\n",
                                    "(Aging sets: FDR<", RANK_PADJ_CUTOFF,
                                    "  |  Old_up n=", N_OLD_UP,
                                    ", Old_down n=", N_OLD_DOWN, ")"),
                 subplots  = 1:3,
                 base_size = 11)
  
  p[[1]] <- p[[1]] +
    scale_color_manual(
      values = c("Old_up" = "#1A237E", "Old_down" = "#B71C1C"),
      labels = c(
        "Old_up"   = paste0("Old-up (n=",   N_OLD_UP,
                            ")   NES=", round(res_stats$NES[res_stats$ID=="Old_up"],   2),
                            "  FDR=",   signif(res_stats$p.adjust[res_stats$ID=="Old_up"],   2)),
        "Old_down" = paste0("Old-down (n=", N_OLD_DOWN,
                            ") NES=",   round(res_stats$NES[res_stats$ID=="Old_down"], 2),
                            "  FDR=",   signif(res_stats$p.adjust[res_stats$ID=="Old_down"], 2))
      )
    ) +
    theme(legend.position = "top", legend.text = element_text(size = 9))
  
  p[[2]] <- p[[2]] +
    scale_color_manual(values = c("Old_up" = "#1A237E", "Old_down" = "#B71C1C"))
  
  outfile <- file.path(OUTDIR, paste0("GSEA_OldvsYoung_on_HBvsLB_", tag, ".png"))
  png(outfile, width = 1400, height = 1200, res = 150)
  print(p)
  dev.off()
  message("Saved: ", basename(outfile))
  
  invisible(gsea_res)
}

gsea_res_lentineg <- run_aging_on_hblb(your_mono_lentineg, "LentiNeg")
gsea_res_lentipos <- run_aging_on_hblb(your_mono_lentipos, "LentiPos")

# ── Ranked lists ──────────────────────────────────────────────────────────────
ranked_lenti_hb <- lenti_hb %>%
  filter(X != "Lenti") %>%          # remove the lenti transgene itself
  arrange(desc(logFC)) %>%
  distinct(X, .keep_all = TRUE) %>%
  { setNames(.$logFC, .$X) }

ranked_lenti_lb <- lenti_lb %>%
  filter(X != "Lenti") %>%
  arrange(desc(logFC)) %>%
  distinct(X, .keep_all = TRUE) %>%
  { setNames(.$logFC, .$X) }

cat("Ranked Lenti HB:", length(ranked_lenti_hb), "genes\n")
cat("Ranked Lenti LB:", length(ranked_lenti_lb), "genes\n")

# ── GSEA — Lenti+ vs Lenti- in HB ────────────────────────────────────────────
gsea_lenti_hb <- GSEA(
  geneList      = ranked_lenti_hb,
  TERM2GENE     = aging_sets,
  minGSSize     = 10,
  maxGSSize     = max(5000, N_OLD_UP, N_OLD_DOWN) + 1,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = TRUE
)

# ── GSEA — Lenti+ vs Lenti- in LB ────────────────────────────────────────────
gsea_lenti_lb <- GSEA(
  geneList      = ranked_lenti_lb,
  TERM2GENE     = aging_sets,
  minGSSize     = 10,
  maxGSSize     = max(5000, N_OLD_UP, N_OLD_DOWN) + 1,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = TRUE
)

# ── Print stats ───────────────────────────────────────────────────────────────
cat("\n--- Lenti HB ---\n")
print(as.data.frame(gsea_lenti_hb)[, c("ID", "NES", "pvalue", "p.adjust")])

cat("\n--- Lenti LB ---\n")
print(as.data.frame(gsea_lenti_lb)[, c("ID", "NES", "pvalue", "p.adjust")])

# ── Plot helper ───────────────────────────────────────────────────────────────
plot_aging_gsea <- function(gsea_obj, title, outpath) {
  stats <- as.data.frame(gsea_obj)[, c("ID", "NES", "p.adjust")]
  
  p <- gseaplot2(gsea_obj,
                 geneSetID = c("Old_up", "Old_down"),
                 title     = paste0(title,
                                    "\n(Aging sets: FDR<", PADJ_CUTOFF,
                                    "  |  Old_up n=", N_OLD_UP,
                                    ", Old_down n=", N_OLD_DOWN, ")"),
                 subplots  = 1:3,
                 base_size = 11)
  
  p[[1]] <- p[[1]] +
    scale_color_manual(
      values = c("Old_up" = "#1A237E", "Old_down" = "#B71C1C"),
      labels = c(
        "Old_up"   = paste0("Old-up (n=",   N_OLD_UP,
                            ")   NES=", round(stats$NES[stats$ID == "Old_up"],   2),
                            "  FDR=",   signif(stats$p.adjust[stats$ID == "Old_up"],   2)),
        "Old_down" = paste0("Old-down (n=", N_OLD_DOWN,
                            ") NES=",   round(stats$NES[stats$ID == "Old_down"], 2),
                            "  FDR=",   signif(stats$p.adjust[stats$ID == "Old_down"], 2))
      )
    ) +
    theme(legend.position = "top", legend.text = element_text(size = 9))
  
  p[[2]] <- p[[2]] +
    scale_color_manual(values = c("Old_up" = "#1A237E", "Old_down" = "#B71C1C"))
  
  png(outpath, width = 1400, height = 1200, res = 150)
  print(p)
  dev.off()
  message("Saved: ", basename(outpath))
}

plot_aging_gsea(gsea_lenti_hb,
                "Aging signature enrichment in Lenti+ vs Lenti- (HB monocytes)",
                file.path(OUTDIR, "GSEA_OldvsYoung_on_LentiHB.png"))

plot_aging_gsea(gsea_lenti_lb,
                "Aging signature enrichment in Lenti+ vs Lenti- (LB monocytes)",
                file.path(OUTDIR, "GSEA_OldvsYoung_on_LentiLB.png"))
###switch####

# ── Cross-dataset GSEA: HB-vs-LB signature on aging ranked list ──────────────
# Runs for BOTH LentiNeg and LentiPos HB-vs-LB contrasts.
# Uses ALL genes with FDR < RANK_PADJ_CUTOFF from the HB-vs-LB DEG table,
# split by direction of logFC.

run_hblb_on_aging <- function(deg_tbl, tag) {
  
  # Ranked list from Old vs Young — built inside so the function is
  # self-contained and can be called without needing ranked_aging in scope.
  ranked_aging <- res_df %>%
    dplyr::arrange(desc(log2FC)) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    { setNames(.$log2FC, .$gene) }
  
  fdr_col <- intersect(c("FDR", "padj", "adj.P.Val", "adj_p"), colnames(deg_tbl))[1]
  if (is.na(fdr_col)) {
    stop("[", tag, "] can't find FDR column. Columns present: ",
         paste(colnames(deg_tbl), collapse = ", "))
  }
  cat("\n[", tag, "] FDR column in DEG table:", fdr_col, "\n")
  
  hb_sig <- deg_tbl %>%
    dplyr::filter(!is.na(.data[[fdr_col]]), !is.na(logFC),
                  .data[[fdr_col]] < RANK_PADJ_CUTOFF) %>%
    dplyr::distinct(X, .keep_all = TRUE)
  
  hb_up   <- hb_sig %>% dplyr::filter(logFC > 0) %>%
    dplyr::arrange(desc(logFC)) %>% dplyr::pull(X)
  hb_down <- hb_sig %>% dplyr::filter(logFC < 0) %>%
    dplyr::arrange(logFC)        %>% dplyr::pull(X)
  
  n_up   <- length(hb_up)
  n_down <- length(hb_down)
  
  cat("[", tag, "] HB vs LB signature — FDR <", RANK_PADJ_CUTOFF, ":\n")
  cat("  Total significant:", nrow(hb_sig), "\n")
  cat("  HB_up   (higher in HB):", n_up,   "\n")
  cat("  HB_down (higher in LB):", n_down, "\n")
  
  # Guard: GSEA needs at least one set >= minGSSize (10)
  if (n_up < 10 && n_down < 10) {
    warning("[", tag, "] Not enough significant genes to run GSEA ",
            "(need >=10 in at least one direction). Skipping.")
    return(invisible(NULL))
  }
  
  hb_sets <- data.frame(
    term = c(rep("HB_up", n_up), rep("HB_down", n_down)),
    gene = c(hb_up, hb_down)
  )
  
  gsea_obj <- GSEA(
    geneList      = ranked_aging,
    TERM2GENE     = hb_sets,
    minGSSize     = 10,
    maxGSSize     = max(5000, n_up, n_down) + 1,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    eps           = 0,
    verbose       = TRUE
  )
  
  stats <- as.data.frame(gsea_obj)[, c("ID", "NES", "p.adjust")]
  cat("\n--- [", tag, "] HB-vs-LB sets on aging ranked list ---\n")
  print(stats)
  
  # Only plot directions that actually made it through GSEA (>= minGSSize)
  plot_ids <- intersect(c("HB_up", "HB_down"), stats$ID)
  if (length(plot_ids) == 0) {
    warning("[", tag, "] GSEA returned no usable gene sets. Skipping plot.")
    return(invisible(gsea_obj))
  }
  
  p <- gseaplot2(gsea_obj,
                 geneSetID = plot_ids,
                 title     = paste0(tag, " HB vs LB signature enrichment ",
                                    "in Old vs Young monocytes\n",
                                    "(HB sets: FDR<", RANK_PADJ_CUTOFF,
                                    "  |  HB_up n=", n_up,
                                    ", HB_down n=", n_down, ")"),
                 subplots  = 1:3,
                 base_size = 11)
  
  # Build colour / label maps dynamically so we don't reference missing IDs
  col_map <- c("HB_up" = "#1A237E", "HB_down" = "#B71C1C")[plot_ids]
  lab_map <- setNames(sapply(plot_ids, function(id) {
    n  <- if (id == "HB_up") n_up else n_down
    pretty <- if (id == "HB_up") "HB-up" else "HB-down"
    paste0(pretty, " (n=", n,
           ")  NES=", round(stats$NES[stats$ID == id], 2),
           "  FDR=",  signif(stats$p.adjust[stats$ID == id], 2))
  }), plot_ids)
  
  p[[1]] <- p[[1]] +
    scale_color_manual(values = col_map, labels = lab_map) +
    theme(legend.position = "top", legend.text = element_text(size = 9))
  
  p[[2]] <- p[[2]] + scale_color_manual(values = col_map)
  
  outfile <- file.path(OUTDIR, paste0("GSEA_HBvsLB_on_OldvsYoung_", tag, ".png"))
  png(outfile, width = 1400, height = 1200, res = 150)
  print(p)
  dev.off()
  message("Saved: ", basename(outfile))
  
  invisible(gsea_obj)
}

gsea_hb_on_aging_lentineg <- run_hblb_on_aging(your_mono_lentineg, "LentiNeg")
gsea_hb_on_aging_lentipos <- run_hblb_on_aging(your_mono_lentipos, "LentiPos")


## ============================================================
## H2B-STYLE PATHWAY ANALYSIS ON Old-vs-Young RANKED LIST
##
## Ranked list: res_df log2FC (positive = higher in Old monocytes)
## Three analyses, consistent bar + dotplot + CSV output:
##   A. GSEA — GO Biological Process (gseGO)
##   B. GSEA — KEGG (gseKEGG)
##   C. GSEA — MSigDB Hallmarks
##
## OUTPUTS (written to OUTDIR/AgingCell_pathways/):
##   Aging_GSEA_GOBP_bar.png      + .csv + _dot.png
##   Aging_GSEA_KEGG_bar.png      + .csv + _dot.png
##   Aging_GSEA_Hallmarks_bar.png + .csv + _dot.png
## ============================================================

# (All packages loaded at top — no library() reloads here to keep
# dplyr verbs unmasked.)

PW_DIR <- file.path(OUTDIR, "AgingCell_pathways")
dir.create(PW_DIR, showWarnings = FALSE)

## ── Ranked list: positive = higher in OLD monocytes ──────────────────────────
ranked_aging_full <- res_df %>%
  filter(!is.na(log2FC)) %>%
  arrange(desc(log2FC)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  { setNames(.$log2FC, .$gene) }

cat("Full aging ranked list:", length(ranked_aging_full), "genes\n")

## ── Symbol -> Entrez for KEGG ─────────────────────────────────────────────────
sym2entrez_age <- bitr(names(ranked_aging_full),
                       fromType = "SYMBOL",
                       toType   = "ENTREZID",
                       OrgDb    = org.Mm.eg.db)
cat("Mapped", nrow(sym2entrez_age), "of", length(ranked_aging_full),
    "genes to Entrez\n")

ranked_aging_entrez <- ranked_aging_full[sym2entrez_age$SYMBOL]
names(ranked_aging_entrez) <- sym2entrez_age$ENTREZID
ranked_aging_entrez <- ranked_aging_entrez[!duplicated(names(ranked_aging_entrez))]
ranked_aging_entrez <- sort(ranked_aging_entrez, decreasing = TRUE)

## ── Plot helper (same style as H2B pipeline) ──────────────────────────────────
save_pathway_plots <- function(gsea_obj, prefix, top_n = 20) {
  
  res_df_pw <- as.data.frame(gsea_obj)
  write.csv(res_df_pw, file.path(PW_DIR, paste0(prefix, ".csv")), row.names = FALSE)
  cat("\n---", prefix, "--- significant (FDR<0.05):",
      sum(res_df_pw$p.adjust < 0.05), "\n")
  
  sig <- res_df_pw[res_df_pw$p.adjust < 0.05, ]
  if (nrow(sig) == 0) {
    message("No significant pathways for ", prefix,
            " at FDR<0.05 — using top ", top_n)
    sig <- res_df_pw[order(res_df_pw$pvalue), ][seq_len(min(top_n, nrow(res_df_pw))), ]
  }
  
  # Barplot: top enriched (Old) and depleted (Young)
  top_pos <- sig[sig$NES > 0, ][order(-sig$NES[sig$NES > 0]), ][seq_len(min(15, sum(sig$NES > 0))), ]
  top_neg <- sig[sig$NES < 0, ][order( sig$NES[sig$NES < 0]), ][seq_len(min(15, sum(sig$NES < 0))), ]
  bar_df  <- rbind(top_pos, top_neg)
  
  if (nrow(bar_df) > 0) {
    bar_df$Description <- ifelse(nchar(bar_df$Description) > 55,
                                 paste0(substr(bar_df$Description, 1, 52), "..."),
                                 bar_df$Description)
    bar_df$Description <- factor(bar_df$Description,
                                 levels = bar_df$Description[order(bar_df$NES)])
    
    p_bar <- ggplot(bar_df, aes(x = NES, y = Description, fill = p.adjust)) +
      geom_col() +
      scale_fill_gradient(low = "#B71C1C", high = "grey80", name = "FDR",
                          limits = c(0, 0.05)) +
      geom_vline(xintercept = 0, color = "grey30") +
      labs(x = "Normalised Enrichment Score",
           y = NULL,
           title = paste0(prefix,
                          "  |  positive NES = enriched in Old monocytes")) +
      theme_classic(base_size = 10) +
      theme(plot.title  = element_text(size = 9, hjust = 0.5),
            axis.text.y = element_text(size = 7))
    
    ggsave(file.path(PW_DIR, paste0(prefix, "_bar.png")),
           p_bar, width = 10, height = max(4, nrow(bar_df) * 0.3 + 1.5), dpi = 150)
    message("Saved: ", prefix, "_bar.png")
  }
  
  # Dotplot
  tryCatch({
    p_dot <- dotplot(gsea_obj, showCategory = top_n, split = ".sign") +
      facet_grid(. ~ .sign) +
      theme(strip.text  = element_text(size = 9),
            axis.text.y = element_text(size = 7))
    ggsave(file.path(PW_DIR, paste0(prefix, "_dot.png")),
           p_dot, width = 10,
           height = max(5, min(top_n, nrow(res_df_pw)) * 0.28 + 2), dpi = 150)
    message("Saved: ", prefix, "_dot.png")
  }, error = function(e) message("dotplot failed for ", prefix, ": ", e$message))
}

## ── A. GO Biological Process ─────────────────────────────────────────────────
cat("\n=== GO BP ===\n")
gsea_go_pw <- gseGO(
  geneList      = ranked_aging_full,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = TRUE
)
gsea_go_pw <- setReadable(gsea_go_pw, OrgDb = org.Mm.eg.db, keyType = "SYMBOL")
save_pathway_plots(gsea_go_pw, "Aging_GSEA_GOBP")

## ── B. KEGG ───────────────────────────────────────────────────────────────────
cat("\n=== KEGG ===\n")
gsea_kegg_pw <- gseKEGG(
  geneList      = ranked_aging_entrez,
  organism      = "mmu",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = TRUE
)
gsea_kegg_pw <- setReadable(gsea_kegg_pw, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
save_pathway_plots(gsea_kegg_pw, "Aging_GSEA_KEGG")

## ── C. MSigDB Hallmarks ───────────────────────────────────────────────────────
cat("\n=== Hallmarks ===\n")
hallmarks_pw <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

gsea_hallmarks_pw <- GSEA(
  geneList      = ranked_aging_full,
  TERM2GENE     = hallmarks_pw,
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = TRUE
)
save_pathway_plots(gsea_hallmarks_pw, "Aging_GSEA_Hallmarks")

message("\n=== Pathway analysis complete. Results in: ", PW_DIR, " ===")

## ── COPY THIS SCRIPT INTO OUTDIR FOR PROVENANCE ───────────────────────────────

if (file.exists(SCRIPT_PATH)) {
  script_dest <- file.path(OUTDIR, basename(SCRIPT_PATH))
  file.copy(SCRIPT_PATH, script_dest, overwrite = TRUE)
  message("Script copied to: ", script_dest)
} else {
  warning("SCRIPT_PATH not found — script not copied into OUTDIR: ", SCRIPT_PATH)
}

## ============================================================
## SUMMARY DOTPLOT — all four cross-dataset GSEA results on one figure
##
## Panel A (forward): Aging sets (Old_up/Old_down) tested on HB-vs-LB ranked lists
## Panel B (reverse): HB-vs-LB sets (HB_up/HB_down) tested on the Old-vs-Young ranked list
## ============================================================

# Helper: pull NES + p.adjust from a GSEA object as a tidy frame
extract_gsea <- function(gsea_obj, comparison, fallback_ids) {
  if (is.null(gsea_obj)) {
    return(data.frame(
      ID         = fallback_ids,
      NES        = NA_real_,
      p.adjust   = NA_real_,
      comparison = comparison
    ))
  }
  d <- as.data.frame(gsea_obj)[, c("ID", "NES", "p.adjust")]
  # Fill in any missing IDs (e.g. a direction dropped by minGSSize)
  missing_ids <- setdiff(fallback_ids, d$ID)
  if (length(missing_ids) > 0) {
    d <- rbind(d, data.frame(ID = missing_ids, NES = NA, p.adjust = NA))
  }
  d$comparison <- comparison
  d
}

## Panel A — forward (aging sets on HB-vs-LB)
fwd_df <- rbind(
  extract_gsea(gsea_res_lentineg, "LentiNeg HB vs LB", c("Old_up", "Old_down")),
  extract_gsea(gsea_res_lentipos, "LentiPos HB vs LB", c("Old_up", "Old_down"))
)
fwd_df$comparison <- factor(fwd_df$comparison,
                            levels = c("LentiNeg HB vs LB", "LentiPos HB vs LB"))
fwd_df$ID <- factor(fwd_df$ID, levels = c("Old_up", "Old_down"))

fwd_labels <- c(
  Old_up   = paste0("Old_up\n(n=",   N_OLD_UP,   ")"),
  Old_down = paste0("Old_down\n(n=", N_OLD_DOWN, ")")
)

## Panel B — reverse (HB-vs-LB sets on aging)
# HB set sizes are computed inside the run_hblb_on_aging function;
# recompute here so we can label the axis correctly.
hb_sizes <- function(deg_tbl) {
  fdr_col <- intersect(c("FDR", "padj", "adj.P.Val", "adj_p"),
                       colnames(deg_tbl))[1]
  sig <- deg_tbl %>%
    dplyr::filter(!is.na(.data[[fdr_col]]), !is.na(logFC),
                  .data[[fdr_col]] < RANK_PADJ_CUTOFF) %>%
    dplyr::distinct(X, .keep_all = TRUE)
  c(up   = sum(sig$logFC > 0),
    down = sum(sig$logFC < 0))
}
hb_ln <- hb_sizes(your_mono_lentineg)
hb_lp <- hb_sizes(your_mono_lentipos)

rev_df <- rbind(
  extract_gsea(gsea_hb_on_aging_lentineg, "OldvsYoung — LentiNeg HB sets",
               c("HB_up", "HB_down")),
  extract_gsea(gsea_hb_on_aging_lentipos, "OldvsYoung — LentiPos HB sets",
               c("HB_up", "HB_down"))
)
rev_df$comparison <- factor(rev_df$comparison,
                            levels = c("OldvsYoung — LentiNeg HB sets",
                                       "OldvsYoung — LentiPos HB sets"))
rev_df$ID <- factor(rev_df$ID, levels = c("HB_up", "HB_down"))

# For the reverse plot, gene set sizes depend on which LentiNeg/Pos contrast
# the sets came from — show both sizes per row label
rev_labels <- c(
  HB_up   = paste0("HB_up\n(LentiNeg n=", hb_ln["up"],
                   ", LentiPos n=", hb_lp["up"], ")"),
  HB_down = paste0("HB_down\n(LentiNeg n=", hb_ln["down"],
                   ", LentiPos n=", hb_lp["down"], ")")
)

cat("\nSummary dotplot — forward:\n"); print(fwd_df)
cat("\nSummary dotplot — reverse:\n"); print(rev_df)

## Shared aesthetics (match H2B summary dotplot style)
nes_scale <- scale_fill_gradient2(low = "#1A237E", mid = "white", high = "#B71C1C",
                                  midpoint = 0, limits = c(-2.5, 2.5),
                                  oob = scales::squish, name = "NES")
size_scale <- scale_size_continuous(range = c(4, 14), name = "-log10(FDR)")

base_theme <- theme_classic(base_size = 13) +
  theme(axis.text.x   = element_text(size = 10, angle = 30, hjust = 1,
                                     margin = margin(t = 4)),
        axis.text.y   = element_text(size = 11),
        plot.title    = element_text(size = 12, hjust = 0.5,
                                     margin = margin(b = 8)),
        plot.subtitle = element_text(size = 9, hjust = 0.5,
                                     margin = margin(b = 8)),
        plot.margin   = margin(t = 15, r = 10, b = 15, l = 15))

p_fwd <- ggplot(fwd_df,
                aes(x = comparison, y = ID, fill = NES,
                    size = -log10(p.adjust))) +
  geom_point(shape = 21, color = "grey30") +
  nes_scale + size_scale +
  scale_y_discrete(labels = fwd_labels, expand = expansion(add = 0.8)) +
  base_theme +
  labs(x = NULL, y = "Aging gene set",
       title    = "Aging signature tested on HB vs LB monocyte contrasts",
       subtitle = paste0("Aging sets: GSE207063 DEGs with FDR < ", RANK_PADJ_CUTOFF))

p_rev <- ggplot(rev_df,
                aes(x = comparison, y = ID, fill = NES,
                    size = -log10(p.adjust))) +
  geom_point(shape = 21, color = "grey30") +
  nes_scale + size_scale +
  scale_y_discrete(labels = rev_labels, expand = expansion(add = 0.8)) +
  base_theme +
  labs(x = NULL, y = "HB-vs-LB gene set",
       title    = "HB vs LB signatures tested on Old vs Young monocytes",
       subtitle = paste0("HB sets: pseudobulk DEGs with FDR < ", RANK_PADJ_CUTOFF))

# Stack with patchwork, shared legend on the right
summary_plot <- (p_fwd / p_rev) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(OUTDIR, "Summary_dotplot_crossdataset.png"),
       summary_plot, width = 9, height = 8, dpi = 150)
message("Saved: Summary_dotplot_crossdataset.png")


## ============================================================
## H2B-STYLE PATHWAY ANALYSIS ON Old-vs-Young RANKED LIST
##
## Ranked list: res_df log2FC (positive = higher in Old monocytes)
## Three analyses, consistent bar + dotplot + CSV output:
##   A. GSEA — GO Biological Process (gseGO)
##   B. GSEA — KEGG (gseKEGG)
##   C. GSEA — MSigDB Hallmarks
##
## OUTPUTS (written to OUTDIR/AgingCell_pathways/):
##   Aging_GSEA_GOBP_bar.png      + .csv + _dot.png
##   Aging_GSEA_KEGG_bar.png      + .csv + _dot.png
##   Aging_GSEA_Hallmarks_bar.png + .csv + _dot.png
## ============================================================

# (All packages loaded at top — no library() reloads here to keep
# dplyr verbs unmasked.)

PW_DIR <- file.path(OUTDIR, "AgingCell_pathways")
dir.create(PW_DIR, showWarnings = FALSE)

## ── Ranked list: positive = higher in OLD monocytes ──────────────────────────
ranked_aging_full <- res_df %>%
  filter(!is.na(log2FC)) %>%
  arrange(desc(log2FC)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  { setNames(.$log2FC, .$gene) }

cat("Full aging ranked list:", length(ranked_aging_full), "genes\n")

## ── Symbol -> Entrez for KEGG ─────────────────────────────────────────────────
sym2entrez_age <- bitr(names(ranked_aging_full),
                       fromType = "SYMBOL",
                       toType   = "ENTREZID",
                       OrgDb    = org.Mm.eg.db)
cat("Mapped", nrow(sym2entrez_age), "of", length(ranked_aging_full),
    "genes to Entrez\n")

ranked_aging_entrez <- ranked_aging_full[sym2entrez_age$SYMBOL]
names(ranked_aging_entrez) <- sym2entrez_age$ENTREZID
ranked_aging_entrez <- ranked_aging_entrez[!duplicated(names(ranked_aging_entrez))]
ranked_aging_entrez <- sort(ranked_aging_entrez, decreasing = TRUE)

## ── Plot helper (same style as H2B pipeline) ──────────────────────────────────
save_pathway_plots <- function(gsea_obj, prefix, top_n = 20) {
  
  res_df_pw <- as.data.frame(gsea_obj)
  write.csv(res_df_pw, file.path(PW_DIR, paste0(prefix, ".csv")), row.names = FALSE)
  cat("\n---", prefix, "--- significant (FDR<0.05):",
      sum(res_df_pw$p.adjust < 0.05), "\n")
  
  sig <- res_df_pw[res_df_pw$p.adjust < 0.05, ]
  if (nrow(sig) == 0) {
    message("No significant pathways for ", prefix,
            " at FDR<0.05 — using top ", top_n)
    sig <- res_df_pw[order(res_df_pw$pvalue), ][seq_len(min(top_n, nrow(res_df_pw))), ]
  }
  
  # Barplot: top enriched (Old) and depleted (Young)
  top_pos <- sig[sig$NES > 0, ][order(-sig$NES[sig$NES > 0]), ][seq_len(min(15, sum(sig$NES > 0))), ]
  top_neg <- sig[sig$NES < 0, ][order( sig$NES[sig$NES < 0]), ][seq_len(min(15, sum(sig$NES < 0))), ]
  bar_df  <- rbind(top_pos, top_neg)
  
  if (nrow(bar_df) > 0) {
    bar_df$Description <- ifelse(nchar(bar_df$Description) > 55,
                                 paste0(substr(bar_df$Description, 1, 52), "..."),
                                 bar_df$Description)
    bar_df$Description <- factor(bar_df$Description,
                                 levels = bar_df$Description[order(bar_df$NES)])
    
    p_bar <- ggplot(bar_df, aes(x = NES, y = Description, fill = p.adjust)) +
      geom_col() +
      scale_fill_gradient(low = "#B71C1C", high = "grey80", name = "FDR",
                          limits = c(0, 0.05)) +
      geom_vline(xintercept = 0, color = "grey30") +
      labs(x = "Normalised Enrichment Score",
           y = NULL,
           title = paste0(prefix,
                          "  |  positive NES = enriched in Old monocytes")) +
      theme_classic(base_size = 10) +
      theme(plot.title  = element_text(size = 9, hjust = 0.5),
            axis.text.y = element_text(size = 7))
    
    ggsave(file.path(PW_DIR, paste0(prefix, "_bar.png")),
           p_bar, width = 10, height = max(4, nrow(bar_df) * 0.3 + 1.5), dpi = 150)
    message("Saved: ", prefix, "_bar.png")
  }
  
  # Dotplot
  tryCatch({
    p_dot <- dotplot(gsea_obj, showCategory = top_n, split = ".sign") +
      facet_grid(. ~ .sign) +
      theme(strip.text  = element_text(size = 9),
            axis.text.y = element_text(size = 7))
    ggsave(file.path(PW_DIR, paste0(prefix, "_dot.png")),
           p_dot, width = 10,
           height = max(5, min(top_n, nrow(res_df_pw)) * 0.28 + 2), dpi = 150)
    message("Saved: ", prefix, "_dot.png")
  }, error = function(e) message("dotplot failed for ", prefix, ": ", e$message))
}

## ── A. GO Biological Process ─────────────────────────────────────────────────
cat("\n=== GO BP ===\n")
gsea_go_pw <- gseGO(
  geneList      = ranked_aging_full,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = TRUE
)
gsea_go_pw <- setReadable(gsea_go_pw, OrgDb = org.Mm.eg.db, keyType = "SYMBOL")
save_pathway_plots(gsea_go_pw, "Aging_GSEA_GOBP")

## ── B. KEGG ───────────────────────────────────────────────────────────────────
cat("\n=== KEGG ===\n")
gsea_kegg_pw <- gseKEGG(
  geneList      = ranked_aging_entrez,
  organism      = "mmu",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = TRUE
)
gsea_kegg_pw <- setReadable(gsea_kegg_pw, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
save_pathway_plots(gsea_kegg_pw, "Aging_GSEA_KEGG")

## ── C. MSigDB Hallmarks ───────────────────────────────────────────────────────
cat("\n=== Hallmarks ===\n")
hallmarks_pw <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

gsea_hallmarks_pw <- GSEA(
  geneList      = ranked_aging_full,
  TERM2GENE     = hallmarks_pw,
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = TRUE
)
save_pathway_plots(gsea_hallmarks_pw, "Aging_GSEA_Hallmarks")

message("\n=== Pathway analysis complete. Results in: ", PW_DIR, " ===")

## ── COPY THIS SCRIPT INTO OUTDIR FOR PROVENANCE ───────────────────────────────

if (file.exists(SCRIPT_PATH)) {
  script_dest <- file.path(OUTDIR, basename(SCRIPT_PATH))
  file.copy(SCRIPT_PATH, script_dest, overwrite = TRUE)
  message("Script copied to: ", script_dest)
} else {
  warning("SCRIPT_PATH not found — script not copied into OUTDIR: ", SCRIPT_PATH)
}

## ============================================================
## SUMMARY DOTPLOT — all four cross-dataset GSEA results on one figure
##
## Panel A (forward): Aging sets (Old_up/Old_down) tested on HB-vs-LB ranked lists
## Panel B (reverse): HB-vs-LB sets (HB_up/HB_down) tested on the Old-vs-Young ranked list
## ============================================================

# Helper: pull NES + p.adjust from a GSEA object as a tidy frame
extract_gsea <- function(gsea_obj, comparison, fallback_ids) {
  if (is.null(gsea_obj)) {
    return(data.frame(
      ID         = fallback_ids,
      NES        = NA_real_,
      p.adjust   = NA_real_,
      comparison = comparison
    ))
  }
  d <- as.data.frame(gsea_obj)[, c("ID", "NES", "p.adjust")]
  # Fill in any missing IDs (e.g. a direction dropped by minGSSize)
  missing_ids <- setdiff(fallback_ids, d$ID)
  if (length(missing_ids) > 0) {
    d <- rbind(d, data.frame(ID = missing_ids, NES = NA, p.adjust = NA))
  }
  d$comparison <- comparison
  d
}

## Panel A — forward (aging sets on HB-vs-LB)
fwd_df <- rbind(
  extract_gsea(gsea_res_lentineg, "LentiNeg HB vs LB", c("Old_up", "Old_down")),
  extract_gsea(gsea_res_lentipos, "LentiPos HB vs LB", c("Old_up", "Old_down"))
)
fwd_df$comparison <- factor(fwd_df$comparison,
                            levels = c("LentiNeg HB vs LB", "LentiPos HB vs LB"))
fwd_df$ID <- factor(fwd_df$ID, levels = c("Old_up", "Old_down"))

fwd_labels <- c(
  Old_up   = paste0("Old_up\n(n=",   N_OLD_UP,   ")"),
  Old_down = paste0("Old_down\n(n=", N_OLD_DOWN, ")")
)

## Panel B — reverse (HB-vs-LB sets on aging)
# HB set sizes are computed inside the run_hblb_on_aging function;
# recompute here so we can label the axis correctly.
hb_sizes <- function(deg_tbl) {
  fdr_col <- intersect(c("FDR", "padj", "adj.P.Val", "adj_p"),
                       colnames(deg_tbl))[1]
  sig <- deg_tbl %>%
    dplyr::filter(!is.na(.data[[fdr_col]]), !is.na(logFC),
                  .data[[fdr_col]] < RANK_PADJ_CUTOFF) %>%
    dplyr::distinct(X, .keep_all = TRUE)
  c(up   = sum(sig$logFC > 0),
    down = sum(sig$logFC < 0))
}
hb_ln <- hb_sizes(your_mono_lentineg)
hb_lp <- hb_sizes(your_mono_lentipos)

rev_df <- rbind(
  extract_gsea(gsea_hb_on_aging_lentineg, "OldvsYoung — LentiNeg HB sets",
               c("HB_up", "HB_down")),
  extract_gsea(gsea_hb_on_aging_lentipos, "OldvsYoung — LentiPos HB sets",
               c("HB_up", "HB_down"))
)
rev_df$comparison <- factor(rev_df$comparison,
                            levels = c("OldvsYoung — LentiNeg HB sets",
                                       "OldvsYoung — LentiPos HB sets"))
rev_df$ID <- factor(rev_df$ID, levels = c("HB_up", "HB_down"))

# For the reverse plot, gene set sizes depend on which LentiNeg/Pos contrast
# the sets came from — show both sizes per row label
rev_labels <- c(
  HB_up   = paste0("HB_up\n(LentiNeg n=", hb_ln["up"],
                   ", LentiPos n=", hb_lp["up"], ")"),
  HB_down = paste0("HB_down\n(LentiNeg n=", hb_ln["down"],
                   ", LentiPos n=", hb_lp["down"], ")")
)

cat("\nSummary dotplot — forward:\n"); print(fwd_df)
cat("\nSummary dotplot — reverse:\n"); print(rev_df)

## Shared aesthetics (match H2B summary dotplot style)
nes_scale <- scale_fill_gradient2(low = "#1A237E", mid = "white", high = "#B71C1C",
                                  midpoint = 0, limits = c(-2.5, 2.5),
                                  oob = scales::squish, name = "NES")
size_scale <- scale_size_continuous(range = c(4, 14), name = "-log10(FDR)")

base_theme <- theme_classic(base_size = 13) +
  theme(axis.text.x   = element_text(size = 10, angle = 30, hjust = 1,
                                     margin = margin(t = 4)),
        axis.text.y   = element_text(size = 11),
        plot.title    = element_text(size = 12, hjust = 0.5,
                                     margin = margin(b = 8)),
        plot.subtitle = element_text(size = 9, hjust = 0.5,
                                     margin = margin(b = 8)),
        plot.margin   = margin(t = 15, r = 10, b = 15, l = 15))

p_fwd <- ggplot(fwd_df,
                aes(x = comparison, y = ID, fill = NES,
                    size = -log10(p.adjust))) +
  geom_point(shape = 21, color = "grey30") +
  nes_scale + size_scale +
  scale_y_discrete(labels = fwd_labels, expand = expansion(add = 0.8)) +
  base_theme +
  labs(x = NULL, y = "Aging gene set",
       title    = "Aging signature tested on HB vs LB monocyte contrasts",
       subtitle = paste0("Aging sets: GSE207063 DEGs with FDR < ", RANK_PADJ_CUTOFF))

p_rev <- ggplot(rev_df,
                aes(x = comparison, y = ID, fill = NES,
                    size = -log10(p.adjust))) +
  geom_point(shape = 21, color = "grey30") +
  nes_scale + size_scale +
  scale_y_discrete(labels = rev_labels, expand = expansion(add = 0.8)) +
  base_theme +
  labs(x = NULL, y = "HB-vs-LB gene set",
       title    = "HB vs LB signatures tested on Old vs Young monocytes",
       subtitle = paste0("HB sets: pseudobulk DEGs with FDR < ", RANK_PADJ_CUTOFF))

# Stack with patchwork, shared legend on the right
summary_plot <- (p_fwd / p_rev) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(OUTDIR, "Summary_dotplot_crossdataset.png"),
       summary_plot, width = 9, height = 8, dpi = 150)
message("Saved: Summary_dotplot_crossdataset.png")


## ============================================================
## H2B-STYLE PATHWAY ANALYSIS ON Old-vs-Young RANKED LIST
##
## Ranked list: res_df log2FC (positive = higher in Old monocytes)
## Three analyses, consistent bar + dotplot + CSV output:
##   A. GSEA — GO Biological Process (gseGO)
##   B. GSEA — KEGG (gseKEGG)
##   C. GSEA — MSigDB Hallmarks
##
## OUTPUTS (written to OUTDIR/AgingCell_pathways/):
##   Aging_GSEA_GOBP_bar.png      + .csv + _dot.png
##   Aging_GSEA_KEGG_bar.png      + .csv + _dot.png
##   Aging_GSEA_Hallmarks_bar.png + .csv + _dot.png
## ============================================================

# (All packages loaded at top — no library() reloads here to keep
# dplyr verbs unmasked.)

PW_DIR <- file.path(OUTDIR, "AgingCell_pathways")
dir.create(PW_DIR, showWarnings = FALSE)

## ── Ranked list: positive = higher in OLD monocytes ──────────────────────────
ranked_aging_full <- res_df %>%
  filter(!is.na(log2FC)) %>%
  arrange(desc(log2FC)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  { setNames(.$log2FC, .$gene) }

cat("Full aging ranked list:", length(ranked_aging_full), "genes\n")

## ── Symbol -> Entrez for KEGG ─────────────────────────────────────────────────
sym2entrez_age <- bitr(names(ranked_aging_full),
                       fromType = "SYMBOL",
                       toType   = "ENTREZID",
                       OrgDb    = org.Mm.eg.db)
cat("Mapped", nrow(sym2entrez_age), "of", length(ranked_aging_full),
    "genes to Entrez\n")

ranked_aging_entrez <- ranked_aging_full[sym2entrez_age$SYMBOL]
names(ranked_aging_entrez) <- sym2entrez_age$ENTREZID
ranked_aging_entrez <- ranked_aging_entrez[!duplicated(names(ranked_aging_entrez))]
ranked_aging_entrez <- sort(ranked_aging_entrez, decreasing = TRUE)

## ── Plot helper (same style as H2B pipeline) ──────────────────────────────────
save_pathway_plots <- function(gsea_obj, prefix, top_n = 20) {
  
  res_df_pw <- as.data.frame(gsea_obj)
  write.csv(res_df_pw, file.path(PW_DIR, paste0(prefix, ".csv")), row.names = FALSE)
  cat("\n---", prefix, "--- significant (FDR<0.05):",
      sum(res_df_pw$p.adjust < 0.05), "\n")
  
  sig <- res_df_pw[res_df_pw$p.adjust < 0.05, ]
  if (nrow(sig) == 0) {
    message("No significant pathways for ", prefix,
            " at FDR<0.05 — using top ", top_n)
    sig <- res_df_pw[order(res_df_pw$pvalue), ][seq_len(min(top_n, nrow(res_df_pw))), ]
  }
  
  # Barplot: top enriched (Old) and depleted (Young)
  top_pos <- sig[sig$NES > 0, ][order(-sig$NES[sig$NES > 0]), ][seq_len(min(15, sum(sig$NES > 0))), ]
  top_neg <- sig[sig$NES < 0, ][order( sig$NES[sig$NES < 0]), ][seq_len(min(15, sum(sig$NES < 0))), ]
  bar_df  <- rbind(top_pos, top_neg)
  
  if (nrow(bar_df) > 0) {
    bar_df$Description <- ifelse(nchar(bar_df$Description) > 55,
                                 paste0(substr(bar_df$Description, 1, 52), "..."),
                                 bar_df$Description)
    bar_df$Description <- factor(bar_df$Description,
                                 levels = bar_df$Description[order(bar_df$NES)])
    
    p_bar <- ggplot(bar_df, aes(x = NES, y = Description, fill = p.adjust)) +
      geom_col() +
      scale_fill_gradient(low = "#B71C1C", high = "grey80", name = "FDR",
                          limits = c(0, 0.05)) +
      geom_vline(xintercept = 0, color = "grey30") +
      labs(x = "Normalised Enrichment Score",
           y = NULL,
           title = paste0(prefix,
                          "  |  positive NES = enriched in Old monocytes")) +
      theme_classic(base_size = 10) +
      theme(plot.title  = element_text(size = 9, hjust = 0.5),
            axis.text.y = element_text(size = 7))
    
    ggsave(file.path(PW_DIR, paste0(prefix, "_bar.png")),
           p_bar, width = 10, height = max(4, nrow(bar_df) * 0.3 + 1.5), dpi = 150)
    message("Saved: ", prefix, "_bar.png")
  }
  
  # Dotplot
  tryCatch({
    p_dot <- dotplot(gsea_obj, showCategory = top_n, split = ".sign") +
      facet_grid(. ~ .sign) +
      theme(strip.text  = element_text(size = 9),
            axis.text.y = element_text(size = 7))
    ggsave(file.path(PW_DIR, paste0(prefix, "_dot.png")),
           p_dot, width = 10,
           height = max(5, min(top_n, nrow(res_df_pw)) * 0.28 + 2), dpi = 150)
    message("Saved: ", prefix, "_dot.png")
  }, error = function(e) message("dotplot failed for ", prefix, ": ", e$message))
}

## ── A. GO Biological Process ─────────────────────────────────────────────────
cat("\n=== GO BP ===\n")
gsea_go_pw <- gseGO(
  geneList      = ranked_aging_full,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = TRUE
)
gsea_go_pw <- setReadable(gsea_go_pw, OrgDb = org.Mm.eg.db, keyType = "SYMBOL")
save_pathway_plots(gsea_go_pw, "Aging_GSEA_GOBP")

## ── B. KEGG ───────────────────────────────────────────────────────────────────
cat("\n=== KEGG ===\n")
gsea_kegg_pw <- gseKEGG(
  geneList      = ranked_aging_entrez,
  organism      = "mmu",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = TRUE
)
gsea_kegg_pw <- setReadable(gsea_kegg_pw, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
save_pathway_plots(gsea_kegg_pw, "Aging_GSEA_KEGG")

## ── C. MSigDB Hallmarks ───────────────────────────────────────────────────────
cat("\n=== Hallmarks ===\n")
hallmarks_pw <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

gsea_hallmarks_pw <- GSEA(
  geneList      = ranked_aging_full,
  TERM2GENE     = hallmarks_pw,
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = TRUE
)
save_pathway_plots(gsea_hallmarks_pw, "Aging_GSEA_Hallmarks")

message("\n=== Pathway analysis complete. Results in: ", PW_DIR, " ===")

## ── COPY THIS SCRIPT INTO OUTDIR FOR PROVENANCE ───────────────────────────────

if (file.exists(SCRIPT_PATH)) {
  script_dest <- file.path(OUTDIR, basename(SCRIPT_PATH))
  file.copy(SCRIPT_PATH, script_dest, overwrite = TRUE)
  message("Script copied to: ", script_dest)
} else {
  warning("SCRIPT_PATH not found — script not copied into OUTDIR: ", SCRIPT_PATH)
}

## ============================================================
## SUMMARY DOTPLOT — all four cross-dataset GSEA results on one figure
##
## Panel A (forward): Aging sets (Old_up/Old_down) tested on HB-vs-LB ranked lists
## Panel B (reverse): HB-vs-LB sets (HB_up/HB_down) tested on the Old-vs-Young ranked list
## ============================================================

# Helper: pull NES + p.adjust from a GSEA object as a tidy frame
extract_gsea <- function(gsea_obj, comparison, fallback_ids) {
  if (is.null(gsea_obj)) {
    return(data.frame(
      ID         = fallback_ids,
      NES        = NA_real_,
      p.adjust   = NA_real_,
      comparison = comparison
    ))
  }
  d <- as.data.frame(gsea_obj)[, c("ID", "NES", "p.adjust")]
  # Fill in any missing IDs (e.g. a direction dropped by minGSSize)
  missing_ids <- setdiff(fallback_ids, d$ID)
  if (length(missing_ids) > 0) {
    d <- rbind(d, data.frame(ID = missing_ids, NES = NA, p.adjust = NA))
  }
  d$comparison <- comparison
  d
}

## Panel A — forward (aging sets on HB-vs-LB)
fwd_df <- rbind(
  extract_gsea(gsea_res_lentineg, "LentiNeg HB vs LB", c("Old_up", "Old_down")),
  extract_gsea(gsea_res_lentipos, "LentiPos HB vs LB", c("Old_up", "Old_down"))
)
fwd_df$comparison <- factor(fwd_df$comparison,
                            levels = c("LentiNeg HB vs LB", "LentiPos HB vs LB"))
fwd_df$ID <- factor(fwd_df$ID, levels = c("Old_up", "Old_down"))

fwd_labels <- c(
  Old_up   = paste0("Old_up\n(n=",   N_OLD_UP,   ")"),
  Old_down = paste0("Old_down\n(n=", N_OLD_DOWN, ")")
)

## Panel B — reverse (HB-vs-LB sets on aging)
# HB set sizes are computed inside the run_hblb_on_aging function;
# recompute here so we can label the axis correctly.
hb_sizes <- function(deg_tbl) {
  fdr_col <- intersect(c("FDR", "padj", "adj.P.Val", "adj_p"),
                       colnames(deg_tbl))[1]
  sig <- deg_tbl %>%
    dplyr::filter(!is.na(.data[[fdr_col]]), !is.na(logFC),
                  .data[[fdr_col]] < RANK_PADJ_CUTOFF) %>%
    dplyr::distinct(X, .keep_all = TRUE)
  c(up   = sum(sig$logFC > 0),
    down = sum(sig$logFC < 0))
}
hb_ln <- hb_sizes(your_mono_lentineg)
hb_lp <- hb_sizes(your_mono_lentipos)

rev_df <- rbind(
  extract_gsea(gsea_hb_on_aging_lentineg, "OldvsYoung — LentiNeg HB sets",
               c("HB_up", "HB_down")),
  extract_gsea(gsea_hb_on_aging_lentipos, "OldvsYoung — LentiPos HB sets",
               c("HB_up", "HB_down"))
)
rev_df$comparison <- factor(rev_df$comparison,
                            levels = c("OldvsYoung — LentiNeg HB sets",
                                       "OldvsYoung — LentiPos HB sets"))
rev_df$ID <- factor(rev_df$ID, levels = c("HB_up", "HB_down"))

# For the reverse plot, gene set sizes depend on which LentiNeg/Pos contrast
# the sets came from — show both sizes per row label
rev_labels <- c(
  HB_up   = paste0("HB_up\n(LentiNeg n=", hb_ln["up"],
                   ", LentiPos n=", hb_lp["up"], ")"),
  HB_down = paste0("HB_down\n(LentiNeg n=", hb_ln["down"],
                   ", LentiPos n=", hb_lp["down"], ")")
)

cat("\nSummary dotplot — forward:\n"); print(fwd_df)
cat("\nSummary dotplot — reverse:\n"); print(rev_df)

## Shared aesthetics (match H2B summary dotplot style)
nes_scale <- scale_fill_gradient2(low = "#1A237E", mid = "white", high = "#B71C1C",
                                  midpoint = 0, limits = c(-2.5, 2.5),
                                  oob = scales::squish, name = "NES")
size_scale <- scale_size_continuous(range = c(4, 14), name = "-log10(FDR)")

base_theme <- theme_classic(base_size = 13) +
  theme(axis.text.x   = element_text(size = 10, angle = 30, hjust = 1,
                                     margin = margin(t = 4)),
        axis.text.y   = element_text(size = 11),
        plot.title    = element_text(size = 12, hjust = 0.5,
                                     margin = margin(b = 8)),
        plot.subtitle = element_text(size = 9, hjust = 0.5,
                                     margin = margin(b = 8)),
        plot.margin   = margin(t = 15, r = 10, b = 15, l = 15))

p_fwd <- ggplot(fwd_df,
                aes(x = comparison, y = ID, fill = NES,
                    size = -log10(p.adjust))) +
  geom_point(shape = 21, color = "grey30") +
  nes_scale + size_scale +
  scale_y_discrete(labels = fwd_labels, expand = expansion(add = 0.8)) +
  base_theme +
  labs(x = NULL, y = "Aging gene set",
       title    = "Aging signature tested on HB vs LB monocyte contrasts",
       subtitle = paste0("Aging sets: GSE207063 DEGs with FDR < ", RANK_PADJ_CUTOFF))

p_rev <- ggplot(rev_df,
                aes(x = comparison, y = ID, fill = NES,
                    size = -log10(p.adjust))) +
  geom_point(shape = 21, color = "grey30") +
  nes_scale + size_scale +
  scale_y_discrete(labels = rev_labels, expand = expansion(add = 0.8)) +
  base_theme +
  labs(x = NULL, y = "HB-vs-LB gene set",
       title    = "HB vs LB signatures tested on Old vs Young monocytes",
       subtitle = paste0("HB sets: pseudobulk DEGs with FDR < ", RANK_PADJ_CUTOFF))

# Stack with patchwork, shared legend on the right
summary_plot <- (p_fwd / p_rev) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(OUTDIR, "Summary_dotplot_crossdataset.png"),
       summary_plot, width = 9, height = 8, dpi = 150)
message("Saved: Summary_dotplot_crossdataset.png")


## ============================================================
## H2B-STYLE PATHWAY ANALYSIS ON Old-vs-Young RANKED LIST
##
## Ranked list: res_df log2FC (positive = higher in Old monocytes)
## Three analyses, consistent bar + dotplot + CSV output:
##   A. GSEA — GO Biological Process (gseGO)
##   B. GSEA — KEGG (gseKEGG)
##   C. GSEA — MSigDB Hallmarks
##
## OUTPUTS (written to OUTDIR/AgingCell_pathways/):
##   Aging_GSEA_GOBP_bar.png      + .csv + _dot.png
##   Aging_GSEA_KEGG_bar.png      + .csv + _dot.png
##   Aging_GSEA_Hallmarks_bar.png + .csv + _dot.png
## ============================================================

# (All packages loaded at top — no library() reloads here to keep
# dplyr verbs unmasked.)

PW_DIR <- file.path(OUTDIR, "AgingCell_pathways")
dir.create(PW_DIR, showWarnings = FALSE)

## ── Ranked list: positive = higher in OLD monocytes ──────────────────────────
ranked_aging_full <- res_df %>%
  filter(!is.na(log2FC)) %>%
  arrange(desc(log2FC)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  { setNames(.$log2FC, .$gene) }

cat("Full aging ranked list:", length(ranked_aging_full), "genes\n")

## ── Symbol -> Entrez for KEGG ─────────────────────────────────────────────────
sym2entrez_age <- bitr(names(ranked_aging_full),
                       fromType = "SYMBOL",
                       toType   = "ENTREZID",
                       OrgDb    = org.Mm.eg.db)
cat("Mapped", nrow(sym2entrez_age), "of", length(ranked_aging_full),
    "genes to Entrez\n")

ranked_aging_entrez <- ranked_aging_full[sym2entrez_age$SYMBOL]
names(ranked_aging_entrez) <- sym2entrez_age$ENTREZID
ranked_aging_entrez <- ranked_aging_entrez[!duplicated(names(ranked_aging_entrez))]
ranked_aging_entrez <- sort(ranked_aging_entrez, decreasing = TRUE)

## ── Plot helper (same style as H2B pipeline) ──────────────────────────────────
save_pathway_plots <- function(gsea_obj, prefix, top_n = 20) {
  
  res_df_pw <- as.data.frame(gsea_obj)
  write.csv(res_df_pw, file.path(PW_DIR, paste0(prefix, ".csv")), row.names = FALSE)
  cat("\n---", prefix, "--- significant (FDR<0.05):",
      sum(res_df_pw$p.adjust < 0.05), "\n")
  
  sig <- res_df_pw[res_df_pw$p.adjust < 0.05, ]
  if (nrow(sig) == 0) {
    message("No significant pathways for ", prefix,
            " at FDR<0.05 — using top ", top_n)
    sig <- res_df_pw[order(res_df_pw$pvalue), ][seq_len(min(top_n, nrow(res_df_pw))), ]
  }
  
  # Barplot: top enriched (Old) and depleted (Young)
  top_pos <- sig[sig$NES > 0, ][order(-sig$NES[sig$NES > 0]), ][seq_len(min(15, sum(sig$NES > 0))), ]
  top_neg <- sig[sig$NES < 0, ][order( sig$NES[sig$NES < 0]), ][seq_len(min(15, sum(sig$NES < 0))), ]
  bar_df  <- rbind(top_pos, top_neg)
  
  if (nrow(bar_df) > 0) {
    bar_df$Description <- ifelse(nchar(bar_df$Description) > 55,
                                 paste0(substr(bar_df$Description, 1, 52), "..."),
                                 bar_df$Description)
    bar_df$Description <- factor(bar_df$Description,
                                 levels = bar_df$Description[order(bar_df$NES)])
    
    p_bar <- ggplot(bar_df, aes(x = NES, y = Description, fill = p.adjust)) +
      geom_col() +
      scale_fill_gradient(low = "#B71C1C", high = "grey80", name = "FDR",
                          limits = c(0, 0.05)) +
      geom_vline(xintercept = 0, color = "grey30") +
      labs(x = "Normalised Enrichment Score",
           y = NULL,
           title = paste0(prefix,
                          "  |  positive NES = enriched in Old monocytes")) +
      theme_classic(base_size = 10) +
      theme(plot.title  = element_text(size = 9, hjust = 0.5),
            axis.text.y = element_text(size = 7))
    
    ggsave(file.path(PW_DIR, paste0(prefix, "_bar.png")),
           p_bar, width = 10, height = max(4, nrow(bar_df) * 0.3 + 1.5), dpi = 150)
    message("Saved: ", prefix, "_bar.png")
  }
  
  # Dotplot
  tryCatch({
    p_dot <- dotplot(gsea_obj, showCategory = top_n, split = ".sign") +
      facet_grid(. ~ .sign) +
      theme(strip.text  = element_text(size = 9),
            axis.text.y = element_text(size = 7))
    ggsave(file.path(PW_DIR, paste0(prefix, "_dot.png")),
           p_dot, width = 10,
           height = max(5, min(top_n, nrow(res_df_pw)) * 0.28 + 2), dpi = 150)
    message("Saved: ", prefix, "_dot.png")
  }, error = function(e) message("dotplot failed for ", prefix, ": ", e$message))
}

## ── A. GO Biological Process ─────────────────────────────────────────────────
cat("\n=== GO BP ===\n")
gsea_go_pw <- gseGO(
  geneList      = ranked_aging_full,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = TRUE
)
gsea_go_pw <- setReadable(gsea_go_pw, OrgDb = org.Mm.eg.db, keyType = "SYMBOL")
save_pathway_plots(gsea_go_pw, "Aging_GSEA_GOBP")

## ── B. KEGG ───────────────────────────────────────────────────────────────────
cat("\n=== KEGG ===\n")
gsea_kegg_pw <- gseKEGG(
  geneList      = ranked_aging_entrez,
  organism      = "mmu",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = TRUE
)
gsea_kegg_pw <- setReadable(gsea_kegg_pw, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
save_pathway_plots(gsea_kegg_pw, "Aging_GSEA_KEGG")

## ── C. MSigDB Hallmarks ───────────────────────────────────────────────────────
cat("\n=== Hallmarks ===\n")
hallmarks_pw <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

gsea_hallmarks_pw <- GSEA(
  geneList      = ranked_aging_full,
  TERM2GENE     = hallmarks_pw,
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = TRUE
)
save_pathway_plots(gsea_hallmarks_pw, "Aging_GSEA_Hallmarks")

message("\n=== Pathway analysis complete. Results in: ", PW_DIR, " ===")

## ── COPY THIS SCRIPT INTO OUTDIR FOR PROVENANCE ───────────────────────────────

if (file.exists(SCRIPT_PATH)) {
  script_dest <- file.path(OUTDIR, basename(SCRIPT_PATH))
  file.copy(SCRIPT_PATH, script_dest, overwrite = TRUE)
  message("Script copied to: ", script_dest)
} else {
  warning("SCRIPT_PATH not found — script not copied into OUTDIR: ", SCRIPT_PATH)
}
## ============================================================
## SUMMARY DOTPLOT — all four cross-dataset GSEA results on one figure
##
## Panel A (forward): Aging sets (Old_up/Old_down) tested on HB-vs-LB ranked lists
## Panel B (reverse): HB-vs-LB sets (HB_up/HB_down) tested on the Old-vs-Young ranked list
## ============================================================

# Helper: pull NES + p.adjust from a GSEA object as a tidy frame
extract_gsea <- function(gsea_obj, comparison, fallback_ids) {
  if (is.null(gsea_obj)) {
    return(data.frame(
      ID         = fallback_ids,
      NES        = NA_real_,
      p.adjust   = NA_real_,
      comparison = comparison
    ))
  }
  d <- as.data.frame(gsea_obj)[, c("ID", "NES", "p.adjust")]
  # Fill in any missing IDs (e.g. a direction dropped by minGSSize)
  missing_ids <- setdiff(fallback_ids, d$ID)
  if (length(missing_ids) > 0) {
    d <- rbind(d, data.frame(ID = missing_ids, NES = NA, p.adjust = NA))
  }
  d$comparison <- comparison
  d
}

## Panel A — forward (aging sets on HB-vs-LB)
fwd_df <- rbind(
  extract_gsea(gsea_res_lentineg, "LentiNeg HB vs LB", c("Old_up", "Old_down")),
  extract_gsea(gsea_res_lentipos, "LentiPos HB vs LB", c("Old_up", "Old_down"))
)
fwd_df$comparison <- factor(fwd_df$comparison,
                            levels = c("LentiNeg HB vs LB", "LentiPos HB vs LB"))
fwd_df$ID <- factor(fwd_df$ID, levels = c("Old_up", "Old_down"))

fwd_labels <- c(
  Old_up   = paste0("Old_up\n(n=",   N_OLD_UP,   ")"),
  Old_down = paste0("Old_down\n(n=", N_OLD_DOWN, ")")
)

## Panel B — reverse (HB-vs-LB sets on aging)
# HB set sizes are computed inside the run_hblb_on_aging function;
# recompute here so we can label the axis correctly.
hb_sizes <- function(deg_tbl) {
  fdr_col <- intersect(c("FDR", "padj", "adj.P.Val", "adj_p"),
                       colnames(deg_tbl))[1]
  sig <- deg_tbl %>%
    dplyr::filter(!is.na(.data[[fdr_col]]), !is.na(logFC),
                  .data[[fdr_col]] < RANK_PADJ_CUTOFF) %>%
    dplyr::distinct(X, .keep_all = TRUE)
  c(up   = sum(sig$logFC > 0),
    down = sum(sig$logFC < 0))
}
hb_ln <- hb_sizes(your_mono_lentineg)
hb_lp <- hb_sizes(your_mono_lentipos)

rev_df <- rbind(
  extract_gsea(gsea_hb_on_aging_lentineg, "OldvsYoung — LentiNeg HB sets",
               c("HB_up", "HB_down")),
  extract_gsea(gsea_hb_on_aging_lentipos, "OldvsYoung — LentiPos HB sets",
               c("HB_up", "HB_down"))
)
rev_df$comparison <- factor(rev_df$comparison,
                            levels = c("OldvsYoung — LentiNeg HB sets",
                                       "OldvsYoung — LentiPos HB sets"))
rev_df$ID <- factor(rev_df$ID, levels = c("HB_up", "HB_down"))

# For the reverse plot, gene set sizes depend on which LentiNeg/Pos contrast
# the sets came from — show both sizes per row label
rev_labels <- c(
  HB_up   = paste0("HB_up\n(LentiNeg n=", hb_ln["up"],
                   ", LentiPos n=", hb_lp["up"], ")"),
  HB_down = paste0("HB_down\n(LentiNeg n=", hb_ln["down"],
                   ", LentiPos n=", hb_lp["down"], ")")
)

cat("\nSummary dotplot — forward:\n"); print(fwd_df)
cat("\nSummary dotplot — reverse:\n"); print(rev_df)

## Shared aesthetics (match H2B summary dotplot style)
nes_scale <- scale_fill_gradient2(low = "#1A237E", mid = "white", high = "#B71C1C",
                                  midpoint = 0, limits = c(-2.5, 2.5),
                                  oob = scales::squish, name = "NES")
size_scale <- scale_size_continuous(range = c(4, 14), name = "-log10(FDR)")

base_theme <- theme_classic(base_size = 13) +
  theme(axis.text.x   = element_text(size = 10, angle = 30, hjust = 1,
                                     margin = margin(t = 4)),
        axis.text.y   = element_text(size = 11),
        plot.title    = element_text(size = 12, hjust = 0.5,
                                     margin = margin(b = 8)),
        plot.subtitle = element_text(size = 9, hjust = 0.5,
                                     margin = margin(b = 8)),
        plot.margin   = margin(t = 15, r = 10, b = 15, l = 15))

p_fwd <- ggplot(fwd_df,
                aes(x = comparison, y = ID, fill = NES,
                    size = -log10(p.adjust))) +
  geom_point(shape = 21, color = "grey30") +
  nes_scale + size_scale +
  scale_y_discrete(labels = fwd_labels, expand = expansion(add = 0.8)) +
  base_theme +
  labs(x = NULL, y = "Aging gene set",
       title    = "Aging signature tested on HB vs LB monocyte contrasts",
       subtitle = paste0("Aging sets: GSE207063 DEGs with FDR < ", RANK_PADJ_CUTOFF))

p_rev <- ggplot(rev_df,
                aes(x = comparison, y = ID, fill = NES,
                    size = -log10(p.adjust))) +
  geom_point(shape = 21, color = "grey30") +
  nes_scale + size_scale +
  scale_y_discrete(labels = rev_labels, expand = expansion(add = 0.8)) +
  base_theme +
  labs(x = NULL, y = "HB-vs-LB gene set",
       title    = "HB vs LB signatures tested on Old vs Young monocytes",
       subtitle = paste0("HB sets: pseudobulk DEGs with FDR < ", RANK_PADJ_CUTOFF))

# Stack with patchwork, shared legend on the right
summary_plot <- (p_fwd / p_rev) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(OUTDIR, "Summary_dotplot_crossdataset.png"),
       summary_plot, width = 9, height = 8, dpi = 150)
message("Saved: Summary_dotplot_crossdataset.png")


## ============================================================
## H2B-STYLE PATHWAY ANALYSIS ON Old-vs-Young RANKED LIST
##
## Ranked list: res_df log2FC (positive = higher in Old monocytes)
## Three analyses, consistent bar + dotplot + CSV output:
##   A. GSEA — GO Biological Process (gseGO)
##   B. GSEA — KEGG (gseKEGG)
##   C. GSEA — MSigDB Hallmarks
##
## OUTPUTS (written to OUTDIR/AgingCell_pathways/):
##   Aging_GSEA_GOBP_bar.png      + .csv + _dot.png
##   Aging_GSEA_KEGG_bar.png      + .csv + _dot.png
##   Aging_GSEA_Hallmarks_bar.png + .csv + _dot.png
## ============================================================

# (All packages loaded at top — no library() reloads here to keep
# dplyr verbs unmasked.)

PW_DIR <- file.path(OUTDIR, "AgingCell_pathways")
dir.create(PW_DIR, showWarnings = FALSE)

## ── Ranked list: positive = higher in OLD monocytes ──────────────────────────
ranked_aging_full <- res_df %>%
  filter(!is.na(log2FC)) %>%
  arrange(desc(log2FC)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  { setNames(.$log2FC, .$gene) }

cat("Full aging ranked list:", length(ranked_aging_full), "genes\n")

## ── Symbol -> Entrez for KEGG ─────────────────────────────────────────────────
sym2entrez_age <- bitr(names(ranked_aging_full),
                       fromType = "SYMBOL",
                       toType   = "ENTREZID",
                       OrgDb    = org.Mm.eg.db)
cat("Mapped", nrow(sym2entrez_age), "of", length(ranked_aging_full),
    "genes to Entrez\n")

ranked_aging_entrez <- ranked_aging_full[sym2entrez_age$SYMBOL]
names(ranked_aging_entrez) <- sym2entrez_age$ENTREZID
ranked_aging_entrez <- ranked_aging_entrez[!duplicated(names(ranked_aging_entrez))]
ranked_aging_entrez <- sort(ranked_aging_entrez, decreasing = TRUE)

## ── Plot helper (same style as H2B pipeline) ──────────────────────────────────
save_pathway_plots <- function(gsea_obj, prefix, top_n = 20) {
  
  res_df_pw <- as.data.frame(gsea_obj)
  write.csv(res_df_pw, file.path(PW_DIR, paste0(prefix, ".csv")), row.names = FALSE)
  cat("\n---", prefix, "--- significant (FDR<0.05):",
      sum(res_df_pw$p.adjust < 0.05), "\n")
  
  sig <- res_df_pw[res_df_pw$p.adjust < 0.05, ]
  if (nrow(sig) == 0) {
    message("No significant pathways for ", prefix,
            " at FDR<0.05 — using top ", top_n)
    sig <- res_df_pw[order(res_df_pw$pvalue), ][seq_len(min(top_n, nrow(res_df_pw))), ]
  }
  
  # Barplot: top enriched (Old) and depleted (Young)
  top_pos <- sig[sig$NES > 0, ][order(-sig$NES[sig$NES > 0]), ][seq_len(min(15, sum(sig$NES > 0))), ]
  top_neg <- sig[sig$NES < 0, ][order( sig$NES[sig$NES < 0]), ][seq_len(min(15, sum(sig$NES < 0))), ]
  bar_df  <- rbind(top_pos, top_neg)
  
  if (nrow(bar_df) > 0) {
    bar_df$Description <- ifelse(nchar(bar_df$Description) > 55,
                                 paste0(substr(bar_df$Description, 1, 52), "..."),
                                 bar_df$Description)
    bar_df$Description <- factor(bar_df$Description,
                                 levels = bar_df$Description[order(bar_df$NES)])
    
    p_bar <- ggplot(bar_df, aes(x = NES, y = Description, fill = p.adjust)) +
      geom_col() +
      scale_fill_gradient(low = "#B71C1C", high = "grey80", name = "FDR",
                          limits = c(0, 0.05)) +
      geom_vline(xintercept = 0, color = "grey30") +
      labs(x = "Normalised Enrichment Score",
           y = NULL,
           title = paste0(prefix,
                          "  |  positive NES = enriched in Old monocytes")) +
      theme_classic(base_size = 10) +
      theme(plot.title  = element_text(size = 9, hjust = 0.5),
            axis.text.y = element_text(size = 7))
    
    ggsave(file.path(PW_DIR, paste0(prefix, "_bar.png")),
           p_bar, width = 10, height = max(4, nrow(bar_df) * 0.3 + 1.5), dpi = 150)
    message("Saved: ", prefix, "_bar.png")
  }
  
  # Dotplot
  tryCatch({
    p_dot <- dotplot(gsea_obj, showCategory = top_n, split = ".sign") +
      facet_grid(. ~ .sign) +
      theme(strip.text  = element_text(size = 9),
            axis.text.y = element_text(size = 7))
    ggsave(file.path(PW_DIR, paste0(prefix, "_dot.png")),
           p_dot, width = 10,
           height = max(5, min(top_n, nrow(res_df_pw)) * 0.28 + 2), dpi = 150)
    message("Saved: ", prefix, "_dot.png")
  }, error = function(e) message("dotplot failed for ", prefix, ": ", e$message))
}

## ── A. GO Biological Process ─────────────────────────────────────────────────
cat("\n=== GO BP ===\n")
gsea_go_pw <- gseGO(
  geneList      = ranked_aging_full,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = TRUE
)
gsea_go_pw <- setReadable(gsea_go_pw, OrgDb = org.Mm.eg.db, keyType = "SYMBOL")
save_pathway_plots(gsea_go_pw, "Aging_GSEA_GOBP")

## ── B. KEGG ───────────────────────────────────────────────────────────────────
cat("\n=== KEGG ===\n")
gsea_kegg_pw <- gseKEGG(
  geneList      = ranked_aging_entrez,
  organism      = "mmu",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = TRUE
)
gsea_kegg_pw <- setReadable(gsea_kegg_pw, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
save_pathway_plots(gsea_kegg_pw, "Aging_GSEA_KEGG")

## ── C. MSigDB Hallmarks ───────────────────────────────────────────────────────
cat("\n=== Hallmarks ===\n")
hallmarks_pw <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

gsea_hallmarks_pw <- GSEA(
  geneList      = ranked_aging_full,
  TERM2GENE     = hallmarks_pw,
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = TRUE
)
save_pathway_plots(gsea_hallmarks_pw, "Aging_GSEA_Hallmarks")

message("\n=== Pathway analysis complete. Results in: ", PW_DIR, " ===")

## ── COPY THIS SCRIPT INTO OUTDIR FOR PROVENANCE ───────────────────────────────

if (file.exists(SCRIPT_PATH)) {
  script_dest <- file.path(OUTDIR, basename(SCRIPT_PATH))
  file.copy(SCRIPT_PATH, script_dest, overwrite = TRUE)
  message("Script copied to: ", script_dest)
} else {
  warning("SCRIPT_PATH not found — script not copied into OUTDIR: ", SCRIPT_PATH)
}


#

#
#
#
#
#
#

# ── SUMMARY DOTPLOT — aging sets across four contrasts ──────────────────────
extract_gsea <- function(gsea_obj, comparison, fallback_ids) {
  if (is.null(gsea_obj)) {
    return(data.frame(ID = fallback_ids, NES = NA_real_,
                      p.adjust = NA_real_, comparison = comparison))
  }
  d <- as.data.frame(gsea_obj)[, c("ID", "NES", "p.adjust")]
  missing_ids <- setdiff(fallback_ids, d$ID)
  if (length(missing_ids) > 0) {
    d <- rbind(d, data.frame(ID = missing_ids, NES = NA, p.adjust = NA))
  }
  d$comparison <- comparison
  d
}

summary_df <- rbind(
  extract_gsea(gsea_res_lentineg, "LentiNeg HB vs LB",  c("Old_up", "Old_down")),
  extract_gsea(gsea_res_lentipos, "LentiPos HB vs LB",  c("Old_up", "Old_down")),
  extract_gsea(gsea_lenti_lb,     "Lenti+/- in LB",     c("Old_up", "Old_down")),
  extract_gsea(gsea_lenti_hb,     "Lenti+/- in HB",     c("Old_up", "Old_down"))
)

summary_df$comparison <- factor(
  summary_df$comparison,
  levels = c("LentiNeg HB vs LB", "LentiPos HB vs LB",
             "Lenti+/- in LB",    "Lenti+/- in HB")
)
summary_df$ID <- factor(summary_df$ID, levels = c("Old_up", "Old_down"))

y_labels <- c(
  Old_up   = paste0("Old_up\n(n=",   N_OLD_UP,   ")"),
  Old_down = paste0("Old_down\n(n=", N_OLD_DOWN, ")")
)

cat("\nSummary dotplot data:\n"); print(summary_df)

p_summary <- ggplot(summary_df,
                    aes(x = comparison, y = ID, fill = NES,
                        size = -log10(p.adjust))) +
  geom_point(shape = 21, color = "grey30") +
  scale_fill_gradient2(low = "#1A237E", mid = "white", high = "#B71C1C",
                       midpoint = 0, limits = c(-2.5, 2.5),
                       oob = scales::squish, name = "NES") +
  scale_size_continuous(
    range  = c(4, 14),
    name   = "FDR",
    breaks = -log10(c(0.5, 0.1, 0.05, 0.01)),
    labels = c("0.5", "0.1", "0.05", "0.01")
  ) +
  scale_y_discrete(labels = y_labels, expand = expansion(add = 0.8)) +
  theme_classic(base_size = 13) +
  theme(axis.text.x    = element_text(size = 10, angle = 30, hjust = 1,
                                      margin = margin(t = 4)),
        axis.text.y    = element_text(size = 11),
        plot.title     = element_text(size = 12, hjust = 0.5,
                                      margin = margin(b = 10)),
        plot.subtitle  = element_text(size = 10, hjust = 0.5,
                                      margin = margin(b = 8)),
        plot.margin    = margin(t = 15, r = 10, b = 25, l = 15)) +
  labs(x = NULL, y = "Aging gene set",
       title    = "Aging signature across classical monocyte comparisons",
       subtitle = paste0("Gene sets: GSE207063 DEGs with FDR < ", RANK_PADJ_CUTOFF,
                         "  |  Old_up n=", N_OLD_UP,
                         ", Old_down n=", N_OLD_DOWN))

ggsave(file.path(OUTDIR, "Summary_dotplot_crossdataset.png"),
       p_summary, width = 9, height = 5, dpi = 150)
message("Saved: Summary_dotplot_crossdataset.png")

#######################

# ── SETUP ───────────────────────────────────────────────────────────────────
OUTDIR2 <- "C:/Users/fc809/Downloads/HBvLBrank_vs_Lentipos_Lentineg"
dir.create(OUTDIR2, showWarnings = FALSE, recursive = TRUE)

PVAL_CUTOFF <- 0.05

# Identify p-value column
pval_col_ln <- intersect(c("PValue", "pvalue", "p.value", "P.Value", "pval"),
                         colnames(your_mono_lentineg))[1]
if (is.na(pval_col_ln)) stop("can't find p-value column. Present: ",
                             paste(colnames(your_mono_lentineg), collapse = ", "))
cat("p-value column in DEG table:", pval_col_ln, "\n")

# ── Build gene sets from LentiNeg HB-vs-LB at p < 0.05 ──────────────────────
hbneg_sig <- your_mono_lentineg %>%
  dplyr::filter(!is.na(.data[[pval_col_ln]]), !is.na(logFC),
                .data[[pval_col_ln]] < PVAL_CUTOFF) %>%
  dplyr::distinct(X, .keep_all = TRUE)

hbneg_up   <- hbneg_sig %>% dplyr::filter(logFC > 0) %>%
  dplyr::arrange(desc(logFC)) %>% dplyr::pull(X)
hbneg_down <- hbneg_sig %>% dplyr::filter(logFC < 0) %>%
  dplyr::arrange(logFC) %>% dplyr::pull(X)

N_HBNEG_UP   <- length(hbneg_up)
N_HBNEG_DOWN <- length(hbneg_down)

cat("HBneg sets (LentiNeg HB vs LB) — p <", PVAL_CUTOFF, ":\n")
cat("  Total significant:", nrow(hbneg_sig), "\n")
cat("  HBneg_up   (higher in HB):", N_HBNEG_UP,   "\n")
cat("  HBneg_down (higher in LB):", N_HBNEG_DOWN, "\n")

hbneg_sets <- data.frame(
  term = c(rep("HBneg_up", N_HBNEG_UP), rep("HBneg_down", N_HBNEG_DOWN)),
  gene = c(hbneg_up, hbneg_down)
)

# ── GSEA per Lenti+/- contrast ──────────────────────────────────────────────
run_hbneg_on_lenti <- function(deg_tbl, tag) {
  if (N_HBNEG_UP < 10 && N_HBNEG_DOWN < 10) {
    warning("[", tag, "] HBneg sets too small. Skipping."); return(invisible(NULL))
  }
  ranked <- deg_tbl %>%
    dplyr::filter(X != "Lenti") %>%
    dplyr::arrange(desc(logFC)) %>%
    dplyr::distinct(X, .keep_all = TRUE) %>%
    { setNames(.$logFC, .$X) }
  cat("\n[", tag, "] Ranked list:", length(ranked), "genes\n")
  
  gsea_obj <- GSEA(geneList = ranked, TERM2GENE = hbneg_sets,
                   minGSSize = 10,
                   maxGSSize = max(5000, N_HBNEG_UP, N_HBNEG_DOWN) + 1,
                   pvalueCutoff = 1, pAdjustMethod = "BH",
                   eps = 0, verbose = TRUE)
  
  stats <- as.data.frame(gsea_obj)[, c("ID", "NES", "pvalue", "p.adjust")]
  cat("--- [", tag, "] ---\n"); print(stats)
  
  plot_ids <- intersect(c("HBneg_up", "HBneg_down"), stats$ID)
  if (length(plot_ids) == 0) return(invisible(gsea_obj))
  
  p <- gseaplot2(gsea_obj, geneSetID = plot_ids,
                 title = paste0("LentiNeg HB-vs-LB signature in Lenti+/- within ", tag,
                                "\n(HBneg sets: p<", PVAL_CUTOFF,
                                " | HBneg_up n=", N_HBNEG_UP,
                                ", HBneg_down n=", N_HBNEG_DOWN, ")"),
                 subplots = 1:3, base_size = 11)
  
  col_map <- c("HBneg_up" = "#1A237E", "HBneg_down" = "#B71C1C")[plot_ids]
  lab_map <- setNames(sapply(plot_ids, function(id) {
    n <- if (id == "HBneg_up") N_HBNEG_UP else N_HBNEG_DOWN
    pretty <- if (id == "HBneg_up") "HBneg-up" else "HBneg-down"
    paste0(pretty, " (n=", n,
           ") NES=", round(stats$NES[stats$ID == id], 2),
           " p=",    signif(stats$pvalue[stats$ID == id], 2))
  }), plot_ids)
  
  p[[1]] <- p[[1]] +
    scale_color_manual(values = col_map, labels = lab_map) +
    theme(legend.position = "top", legend.text = element_text(size = 9))
  p[[2]] <- p[[2]] + scale_color_manual(values = col_map)
  
  outfile <- file.path(OUTDIR2, paste0("GSEA_HBneg_on_Lenti_", tag, ".png"))
  png(outfile, width = 1400, height = 1200, res = 150); print(p); dev.off()
  message("Saved: ", basename(outfile))
  invisible(gsea_obj)
}

gsea_hbneg_on_lb <- run_hbneg_on_lenti(lenti_lb, "LB")
gsea_hbneg_on_hb <- run_hbneg_on_lenti(lenti_hb, "HB")

# ── Summary dotplot (size = raw p) ──────────────────────────────────────────
extract_gsea <- function(gsea_obj, comparison, fallback_ids) {
  if (is.null(gsea_obj)) {
    return(data.frame(ID = fallback_ids, NES = NA_real_,
                      pvalue = NA_real_, p.adjust = NA_real_,
                      comparison = comparison))
  }
  d <- as.data.frame(gsea_obj)[, c("ID", "NES", "pvalue", "p.adjust")]
  missing_ids <- setdiff(fallback_ids, d$ID)
  if (length(missing_ids) > 0) {
    d <- rbind(d, data.frame(ID = missing_ids, NES = NA,
                             pvalue = NA, p.adjust = NA))
  }
  d$comparison <- comparison
  d
}

hbneg_summary_df <- rbind(
  extract_gsea(gsea_hbneg_on_lb, "Lenti+/- in LB", c("HBneg_up", "HBneg_down")),
  extract_gsea(gsea_hbneg_on_hb, "Lenti+/- in HB", c("HBneg_up", "HBneg_down"))
)
hbneg_summary_df$comparison <- factor(hbneg_summary_df$comparison,
                                      levels = c("Lenti+/- in LB", "Lenti+/- in HB"))
hbneg_summary_df$ID <- factor(hbneg_summary_df$ID,
                              levels = c("HBneg_up", "HBneg_down"))

hbneg_y_labels <- c(
  HBneg_up   = paste0("HBneg_up\n(n=",   N_HBNEG_UP,   ")"),
  HBneg_down = paste0("HBneg_down\n(n=", N_HBNEG_DOWN, ")")
)

cat("\nHBneg summary dotplot data:\n"); print(hbneg_summary_df)

p_hbneg_summary <- ggplot(hbneg_summary_df,
                          aes(x = comparison, y = ID, fill = NES,
                              size = -log10(pvalue))) +
  geom_point(shape = 21, color = "grey30") +
  scale_fill_gradient2(low = "#1A237E", mid = "white", high = "#B71C1C",
                       midpoint = 0, limits = c(-2.5, 2.5),
                       oob = scales::squish, name = "NES") +
  scale_size_continuous(range = c(4, 14), name = "p",
                        breaks = -log10(c(0.5, 0.1, 0.05, 0.01, 0.001)),
                        labels = c("0.5", "0.1", "0.05", "0.01", "0.001")) +
  scale_y_discrete(labels = hbneg_y_labels, expand = expansion(add = 0.8)) +
  theme_classic(base_size = 13) +
  theme(axis.text.x   = element_text(size = 10, angle = 30, hjust = 1,
                                     margin = margin(t = 4)),
        axis.text.y   = element_text(size = 11),
        plot.title    = element_text(size = 12, hjust = 0.5,
                                     margin = margin(b = 10)),
        plot.subtitle = element_text(size = 10, hjust = 0.5,
                                     margin = margin(b = 8)),
        plot.margin   = margin(t = 15, r = 10, b = 25, l = 15)) +
  labs(x = NULL, y = "HBneg gene set",
       title    = "LentiNeg HB-vs-LB signature across Lenti+/- contrasts",
       subtitle = paste0("Gene sets: LentiNeg HB-vs-LB DEGs with p < ", PVAL_CUTOFF,
                         "  |  HBneg_up n=", N_HBNEG_UP,
                         ", HBneg_down n=", N_HBNEG_DOWN))

ggsave(file.path(OUTDIR2, "Summary_dotplot_HBneg_on_Lenti.png"),
       p_hbneg_summary, width = 7, height = 5, dpi = 150)
message("Saved: Summary_dotplot_HBneg_on_Lenti.png")

# Copy script for provenance
if (file.exists(SCRIPT_PATH)) {
  file.copy(SCRIPT_PATH, file.path(OUTDIR2, basename(SCRIPT_PATH)),
            overwrite = TRUE)
  message("Script copied to: ", OUTDIR2)
}

################
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)

# ── SETUP ───────────────────────────────────────────────────────────────────
OUTDIR_HSC <- "C:/Users/fc809/Downloads/HBvLBrank_vs_Lentipos_Lentineg_LTHSC"
dir.create(OUTDIR_HSC, showWarnings = FALSE, recursive = TRUE)

PVAL_CUTOFF <- 0.05

DEG_ROOT_HSC <- "C:/Users/fc809/Downloads/LT-HSCs (1)"

deg_files_hsc <- c(
  lentineg_hblb = file.path(DEG_ROOT_HSC, "LentiNegHB_vs_LentiNegLB",
                            "DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.3.sce_0 - Copy.txt"),
  lentipos_hblb = file.path(DEG_ROOT_HSC, "LentiPosHB_vs_LentiPosLB",
                            "DEG_LentiPosHB_vs_LentiPosLB_pseudoBulk_res0.3.sce_0.txt"),
  lenti_lb      = file.path(DEG_ROOT_HSC, "LentiPos_vs_LentiNeg_in_LB",
                            "DEG_LentiPos_vs_LentiNeg_in_LB_pseudoBulk_res0.3.sce_0.txt"),
  lenti_hb      = file.path(DEG_ROOT_HSC, "LentiPos_vs_LentiNeg_in_HB",
                            "DEG_LentiPos_vs_LentiNeg_in_HB_pseudoBulk_res0.3.sce_0.txt")
)

# ── LOADER ──────────────────────────────────────────────────────────────────
load_deg <- function(path) {
  if (!file.exists(path)) stop("DEG file not found: ", path)
  df <- tryCatch(
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  )
  bad <- is.na(names(df)) | names(df) == ""
  if (any(bad)) {
    if (bad[1]) { names(df)[1] <- "X"; bad[1] <- FALSE }
    names(df)[bad] <- paste0("V", seq_len(sum(bad)))
  }
  if (!"X" %in% colnames(df)) {
    gc <- intersect(c("gene", "Gene", "symbol", "Symbol", "gene_name"),
                    colnames(df))
    if (length(gc) > 0) {
      names(df)[names(df) == gc[1]] <- "X"
    } else if (!identical(rownames(df), as.character(seq_len(nrow(df))))) {
      df$X <- rownames(df)
    }
  }
  df
}

for (nm in names(deg_files_hsc)) assign(nm, load_deg(deg_files_hsc[[nm]]))

cat("\nLoaded LT-HSC DEG tables:\n")
for (nm in names(deg_files_hsc)) {
  d <- get(nm)
  cat("  ", nm, ":", nrow(d), "rows  |  cols:",
      paste(head(colnames(d), 8), collapse = ", "), "\n")
}

# ── Identify p-value column in LentiNeg HB-vs-LB ────────────────────────────
pval_col_ln <- intersect(c("PValue", "pvalue", "p.value", "P.Value", "pval"),
                         colnames(lentineg_hblb))[1]
if (is.na(pval_col_ln)) stop("can't find p-value column. Present: ",
                             paste(colnames(lentineg_hblb), collapse = ", "))
cat("\np-value column in LentiNeg HB-vs-LB:", pval_col_ln, "\n")

# ── Build gene sets from LentiNeg HB-vs-LB at p < 0.05 ──────────────────────
hbneg_sig <- lentineg_hblb %>%
  dplyr::filter(!is.na(.data[[pval_col_ln]]), !is.na(logFC),
                .data[[pval_col_ln]] < PVAL_CUTOFF) %>%
  dplyr::distinct(X, .keep_all = TRUE)

hbneg_up   <- hbneg_sig %>% dplyr::filter(logFC > 0) %>%
  dplyr::arrange(desc(logFC)) %>% dplyr::pull(X)
hbneg_down <- hbneg_sig %>% dplyr::filter(logFC < 0) %>%
  dplyr::arrange(logFC) %>% dplyr::pull(X)

N_HBNEG_UP   <- length(hbneg_up)
N_HBNEG_DOWN <- length(hbneg_down)

cat("\nHBneg sets (LentiNeg HB vs LB, LT-HSCs) — p <", PVAL_CUTOFF, ":\n")
cat("  Total significant:", nrow(hbneg_sig), "\n")
cat("  HBneg_up   (higher in HB):", N_HBNEG_UP,   "\n")
cat("  HBneg_down (higher in LB):", N_HBNEG_DOWN, "\n")

hbneg_sets <- data.frame(
  term = c(rep("HBneg_up", N_HBNEG_UP), rep("HBneg_down", N_HBNEG_DOWN)),
  gene = c(hbneg_up, hbneg_down)
)

# ── GSEA per Lenti+/- contrast ──────────────────────────────────────────────
run_hbneg_on_lenti <- function(deg_tbl, tag) {
  if (N_HBNEG_UP < 10 && N_HBNEG_DOWN < 10) {
    warning("[", tag, "] HBneg sets too small. Skipping."); return(invisible(NULL))
  }
  ranked <- deg_tbl %>%
    dplyr::filter(X != "Lenti") %>%
    dplyr::arrange(desc(logFC)) %>%
    dplyr::distinct(X, .keep_all = TRUE) %>%
    { setNames(.$logFC, .$X) }
  cat("\n[", tag, "] Ranked list:", length(ranked), "genes\n")
  
  gsea_obj <- GSEA(geneList = ranked, TERM2GENE = hbneg_sets,
                   minGSSize = 10,
                   maxGSSize = max(5000, N_HBNEG_UP, N_HBNEG_DOWN) + 1,
                   pvalueCutoff = 1, pAdjustMethod = "BH",
                   eps = 0, verbose = TRUE)
  
  stats <- as.data.frame(gsea_obj)[, c("ID", "NES", "pvalue", "p.adjust")]
  cat("--- [", tag, "] ---\n"); print(stats)
  
  plot_ids <- intersect(c("HBneg_up", "HBneg_down"), stats$ID)
  if (length(plot_ids) == 0) return(invisible(gsea_obj))
  
  p <- gseaplot2(gsea_obj, geneSetID = plot_ids,
                 title = paste0("LentiNeg HB-vs-LB signature (LT-HSCs) in Lenti+/- within ", tag,
                                "\n(HBneg sets: p<", PVAL_CUTOFF,
                                " | HBneg_up n=", N_HBNEG_UP,
                                ", HBneg_down n=", N_HBNEG_DOWN, ")"),
                 subplots = 1:3, base_size = 11)
  
  col_map <- c("HBneg_up" = "#1A237E", "HBneg_down" = "#B71C1C")[plot_ids]
  lab_map <- setNames(sapply(plot_ids, function(id) {
    n <- if (id == "HBneg_up") N_HBNEG_UP else N_HBNEG_DOWN
    pretty <- if (id == "HBneg_up") "HBneg-up" else "HBneg-down"
    paste0(pretty, " (n=", n,
           ") NES=", round(stats$NES[stats$ID == id], 2),
           " FDR=",  signif(stats$p.adjust[stats$ID == id], 2))
  }), plot_ids)
  
  p[[1]] <- p[[1]] +
    scale_color_manual(values = col_map, labels = lab_map) +
    theme(legend.position = "top", legend.text = element_text(size = 9))
  p[[2]] <- p[[2]] + scale_color_manual(values = col_map)
  
  outfile <- file.path(OUTDIR_HSC, paste0("GSEA_HBneg_on_Lenti_", tag, ".png"))
  png(outfile, width = 1400, height = 1200, res = 150); print(p); dev.off()
  message("Saved: ", basename(outfile))
  invisible(gsea_obj)
}

gsea_hbneg_on_lb <- run_hbneg_on_lenti(lenti_lb, "LB")
gsea_hbneg_on_hb <- run_hbneg_on_lenti(lenti_hb, "HB")

# ── Summary dotplot (size = FDR) ────────────────────────────────────────────
extract_gsea <- function(gsea_obj, comparison, fallback_ids) {
  if (is.null(gsea_obj)) {
    return(data.frame(ID = fallback_ids, NES = NA_real_,
                      pvalue = NA_real_, p.adjust = NA_real_,
                      comparison = comparison))
  }
  d <- as.data.frame(gsea_obj)[, c("ID", "NES", "pvalue", "p.adjust")]
  missing_ids <- setdiff(fallback_ids, d$ID)
  if (length(missing_ids) > 0) {
    d <- rbind(d, data.frame(ID = missing_ids, NES = NA,
                             pvalue = NA, p.adjust = NA))
  }
  d$comparison <- comparison
  d
}

hbneg_summary_df <- rbind(
  extract_gsea(gsea_hbneg_on_lb, "Lenti+/- in LB", c("HBneg_up", "HBneg_down")),
  extract_gsea(gsea_hbneg_on_hb, "Lenti+/- in HB", c("HBneg_up", "HBneg_down"))
)
hbneg_summary_df$comparison <- factor(hbneg_summary_df$comparison,
                                      levels = c("Lenti+/- in LB", "Lenti+/- in HB"))
hbneg_summary_df$ID <- factor(hbneg_summary_df$ID,
                              levels = c("HBneg_up", "HBneg_down"))

hbneg_y_labels <- c(
  HBneg_up   = paste0("HBneg_up\n(n=",   N_HBNEG_UP,   ")"),
  HBneg_down = paste0("HBneg_down\n(n=", N_HBNEG_DOWN, ")")
)

cat("\nHBneg summary dotplot data (LT-HSCs):\n"); print(hbneg_summary_df)

p_hbneg_summary <- ggplot(hbneg_summary_df,
                          aes(x = comparison, y = ID, fill = NES,
                              size = -log10(p.adjust))) +
  geom_point(shape = 21, color = "grey30") +
  scale_fill_gradient2(low = "#1A237E", mid = "white", high = "#B71C1C",
                       midpoint = 0, limits = c(-2.5, 2.5),
                       oob = scales::squish, name = "NES") +
  scale_size_continuous(range = c(4, 14), name = "FDR",
                        breaks = -log10(c(0.5, 0.1, 0.05, 0.01, 0.001)),
                        labels = c("0.5", "0.1", "0.05", "0.01", "0.001")) +
  scale_y_discrete(labels = hbneg_y_labels, expand = expansion(add = 0.8)) +
  theme_classic(base_size = 13) +
  theme(axis.text.x   = element_text(size = 10, angle = 30, hjust = 1,
                                     margin = margin(t = 4)),
        axis.text.y   = element_text(size = 11),
        plot.title    = element_text(size = 12, hjust = 0.5,
                                     margin = margin(b = 10)),
        plot.subtitle = element_text(size = 10, hjust = 0.5,
                                     margin = margin(b = 8)),
        plot.margin   = margin(t = 15, r = 10, b = 25, l = 15)) +
  labs(x = NULL, y = "HBneg gene set",
       title    = "LentiNeg HB-vs-LB signature (LT-HSCs) across Lenti+/- contrasts",
       subtitle = paste0("Gene sets: LentiNeg HB-vs-LB DEGs with p < ", PVAL_CUTOFF,
                         "  |  HBneg_up n=", N_HBNEG_UP,
                         ", HBneg_down n=", N_HBNEG_DOWN))

ggsave(file.path(OUTDIR_HSC, "Summary_dotplot_HBneg_on_Lenti.png"),
       p_hbneg_summary, width = 7, height = 5, dpi = 150)
message("Saved: Summary_dotplot_HBneg_on_Lenti.png")

####################################

library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)

# ── SETUP ───────────────────────────────────────────────────────────────────
OUTDIR_HSC2 <- "C:/Users/fc809/Downloads/HBvLBrank_vs_Lentipos_Lentineg_LTHSC_LentiPos"
dir.create(OUTDIR_HSC2, showWarnings = FALSE, recursive = TRUE)

# Path to THIS script — update if you save it with a different name/location.
# The script is copied into OUTDIR_HSC2 at the end for provenance.
SCRIPT_PATH_HSC2 <- "C:/Users/fc809/Downloads/LTHSC_HBpos_on_Lenti.R"

PVAL_CUTOFF <- 0.05

DEG_ROOT_HSC <- "C:/Users/fc809/Downloads/LT-HSCs (1)"

deg_files_hsc <- c(
  lentineg_hblb = file.path(DEG_ROOT_HSC, "LentiNegHB_vs_LentiNegLB",
                            "DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.3.sce_0 - Copy.txt"),
  lentipos_hblb = file.path(DEG_ROOT_HSC, "LentiPosHB_vs_LentiPosLB",
                            "DEG_LentiPosHB_vs_LentiPosLB_pseudoBulk_res0.3.sce_0.txt"),
  lenti_lb      = file.path(DEG_ROOT_HSC, "LentiPos_vs_LentiNeg_in_LB",
                            "DEG_LentiPos_vs_LentiNeg_in_LB_pseudoBulk_res0.3.sce_0.txt"),
  lenti_hb      = file.path(DEG_ROOT_HSC, "LentiPos_vs_LentiNeg_in_HB",
                            "DEG_LentiPos_vs_LentiNeg_in_HB_pseudoBulk_res0.3.sce_0.txt")
)

# ── LOADER ──────────────────────────────────────────────────────────────────
load_deg <- function(path) {
  if (!file.exists(path)) stop("DEG file not found: ", path)
  df <- tryCatch(
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  )
  bad <- is.na(names(df)) | names(df) == ""
  if (any(bad)) {
    if (bad[1]) { names(df)[1] <- "X"; bad[1] <- FALSE }
    names(df)[bad] <- paste0("V", seq_len(sum(bad)))
  }
  if (!"X" %in% colnames(df)) {
    gc <- intersect(c("gene", "Gene", "symbol", "Symbol", "gene_name"),
                    colnames(df))
    if (length(gc) > 0) {
      names(df)[names(df) == gc[1]] <- "X"
    } else if (!identical(rownames(df), as.character(seq_len(nrow(df))))) {
      df$X <- rownames(df)
    }
  }
  df
}

for (nm in names(deg_files_hsc)) assign(nm, load_deg(deg_files_hsc[[nm]]))

cat("\nLoaded LT-HSC DEG tables:\n")
for (nm in names(deg_files_hsc)) {
  d <- get(nm)
  cat("  ", nm, ":", nrow(d), "rows  |  cols:",
      paste(head(colnames(d), 8), collapse = ", "), "\n")
}

# ── Identify p-value column in LentiPos HB-vs-LB ────────────────────────────
pval_col_lp <- intersect(c("PValue", "pvalue", "p.value", "P.Value", "pval"),
                         colnames(lentipos_hblb))[1]
if (is.na(pval_col_lp)) stop("can't find p-value column. Present: ",
                             paste(colnames(lentipos_hblb), collapse = ", "))
cat("\np-value column in LentiPos HB-vs-LB:", pval_col_lp, "\n")

# ── Build gene sets from LentiPos HB-vs-LB at p < 0.05 ──────────────────────
hb_sig <- lentipos_hblb %>%
  dplyr::filter(!is.na(.data[[pval_col_lp]]), !is.na(logFC),
                .data[[pval_col_lp]] < PVAL_CUTOFF) %>%
  dplyr::distinct(X, .keep_all = TRUE)

hb_up   <- hb_sig %>% dplyr::filter(logFC > 0) %>%
  dplyr::arrange(desc(logFC)) %>% dplyr::pull(X)
hb_down <- hb_sig %>% dplyr::filter(logFC < 0) %>%
  dplyr::arrange(logFC) %>% dplyr::pull(X)

N_HB_UP   <- length(hb_up)
N_HB_DOWN <- length(hb_down)

cat("\nHB sets (LentiPos HB vs LB, LT-HSCs) — p <", PVAL_CUTOFF, ":\n")
cat("  Total significant:", nrow(hb_sig), "\n")
cat("  HB_up   (higher in HB):", N_HB_UP,   "\n")
cat("  HB_down (higher in LB):", N_HB_DOWN, "\n")

hb_sets <- data.frame(
  term = c(rep("HB_up", N_HB_UP), rep("HB_down", N_HB_DOWN)),
  gene = c(hb_up, hb_down)
)

# ── GSEA per Lenti+/- contrast ──────────────────────────────────────────────
run_hb_on_lenti <- function(deg_tbl, tag) {
  if (N_HB_UP < 10 && N_HB_DOWN < 10) {
    warning("[", tag, "] HB sets too small. Skipping."); return(invisible(NULL))
  }
  ranked <- deg_tbl %>%
    dplyr::filter(X != "Lenti") %>%
    dplyr::arrange(desc(logFC)) %>%
    dplyr::distinct(X, .keep_all = TRUE) %>%
    { setNames(.$logFC, .$X) }
  cat("\n[", tag, "] Ranked list:", length(ranked), "genes\n")
  
  gsea_obj <- GSEA(geneList = ranked, TERM2GENE = hb_sets,
                   minGSSize = 10,
                   maxGSSize = max(5000, N_HB_UP, N_HB_DOWN) + 1,
                   pvalueCutoff = 1, pAdjustMethod = "BH",
                   eps = 0, verbose = TRUE)
  
  stats <- as.data.frame(gsea_obj)[, c("ID", "NES", "pvalue", "p.adjust")]
  cat("--- [", tag, "] ---\n"); print(stats)
  
  plot_ids <- intersect(c("HB_up", "HB_down"), stats$ID)
  if (length(plot_ids) == 0) return(invisible(gsea_obj))
  
  p <- gseaplot2(gsea_obj, geneSetID = plot_ids,
                 title = paste0("LentiPos HB-vs-LB signature (LT-HSCs) in Lenti+/- within ", tag,
                                "\n(HB sets: p<", PVAL_CUTOFF,
                                " | HB_up n=", N_HB_UP,
                                ", HB_down n=", N_HB_DOWN, ")"),
                 subplots = 1:3, base_size = 11)
  
  col_map <- c("HB_up" = "#1A237E", "HB_down" = "#B71C1C")[plot_ids]
  lab_map <- setNames(sapply(plot_ids, function(id) {
    n <- if (id == "HB_up") N_HB_UP else N_HB_DOWN
    pretty <- if (id == "HB_up") "HB-up" else "HB-down"
    paste0(pretty, " (n=", n,
           ") NES=", round(stats$NES[stats$ID == id], 2),
           " FDR=",  signif(stats$p.adjust[stats$ID == id], 2))
  }), plot_ids)
  
  p[[1]] <- p[[1]] +
    scale_color_manual(values = col_map, labels = lab_map) +
    theme(legend.position = "top", legend.text = element_text(size = 9))
  p[[2]] <- p[[2]] + scale_color_manual(values = col_map)
  
  outfile <- file.path(OUTDIR_HSC2, paste0("GSEA_HBpos_on_Lenti_", tag, ".png"))
  png(outfile, width = 1400, height = 1200, res = 150); print(p); dev.off()
  message("Saved: ", basename(outfile))
  invisible(gsea_obj)
}

gsea_hbpos_on_hb <- run_hb_on_lenti(lenti_hb, "HB")
gsea_hbpos_on_lb <- run_hb_on_lenti(lenti_lb, "LB")

# ── Summary dotplot (size = FDR) ────────────────────────────────────────────
extract_gsea <- function(gsea_obj, comparison, fallback_ids) {
  if (is.null(gsea_obj)) {
    return(data.frame(ID = fallback_ids, NES = NA_real_,
                      pvalue = NA_real_, p.adjust = NA_real_,
                      comparison = comparison))
  }
  d <- as.data.frame(gsea_obj)[, c("ID", "NES", "pvalue", "p.adjust")]
  missing_ids <- setdiff(fallback_ids, d$ID)
  if (length(missing_ids) > 0) {
    d <- rbind(d, data.frame(ID = missing_ids, NES = NA,
                             pvalue = NA, p.adjust = NA))
  }
  d$comparison <- comparison
  d
}

hb_summary_df <- rbind(
  extract_gsea(gsea_hbpos_on_lb, "Lenti+/- in LB", c("HB_up", "HB_down")),
  extract_gsea(gsea_hbpos_on_hb, "Lenti+/- in HB", c("HB_up", "HB_down"))
)
hb_summary_df$comparison <- factor(hb_summary_df$comparison,
                                   levels = c("Lenti+/- in LB", "Lenti+/- in HB"))
hb_summary_df$ID <- factor(hb_summary_df$ID, levels = c("HB_up", "HB_down"))

hb_y_labels <- c(
  HB_up   = paste0("HB_up\n(n=",   N_HB_UP,   ")"),
  HB_down = paste0("HB_down\n(n=", N_HB_DOWN, ")")
)

cat("\nHB summary dotplot data (LT-HSCs, LentiPos signature):\n")
print(hb_summary_df)

p_hb_summary <- ggplot(hb_summary_df,
                       aes(x = comparison, y = ID, fill = NES,
                           size = -log10(p.adjust))) +
  geom_point(shape = 21, color = "grey30") +
  scale_fill_gradient2(low = "#1A237E", mid = "white", high = "#B71C1C",
                       midpoint = 0, limits = c(-2.5, 2.5),
                       oob = scales::squish, name = "NES") +
  scale_size_continuous(range = c(4, 14), name = "FDR",
                        breaks = -log10(c(0.5, 0.1, 0.05, 0.01, 0.001)),
                        labels = c("0.5", "0.1", "0.05", "0.01", "0.001")) +
  scale_y_discrete(labels = hb_y_labels, expand = expansion(add = 0.8)) +
  theme_classic(base_size = 13) +
  theme(axis.text.x   = element_text(size = 10, angle = 30, hjust = 1,
                                     margin = margin(t = 4)),
        axis.text.y   = element_text(size = 11),
        plot.title    = element_text(size = 12, hjust = 0.5,
                                     margin = margin(b = 10)),
        plot.subtitle = element_text(size = 10, hjust = 0.5,
                                     margin = margin(b = 8)),
        plot.margin   = margin(t = 15, r = 10, b = 25, l = 15)) +
  labs(x = NULL, y = "HB gene set",
       title    = "LentiPos HB-vs-LB signature (LT-HSCs) across Lenti+/- contrasts",
       subtitle = paste0("Gene sets: LentiPos HB-vs-LB DEGs with p < ", PVAL_CUTOFF,
                         "  |  HB_up n=", N_HB_UP,
                         ", HB_down n=", N_HB_DOWN))

ggsave(file.path(OUTDIR_HSC2, "Summary_dotplot_HBpos_on_Lenti.png"),
       p_hb_summary, width = 7, height = 5, dpi = 150)
message("Saved: Summary_dotplot_HBpos_on_Lenti.png")

# ── Copy script into OUTDIR_HSC2 for provenance ─────────────────────────────
if (file.exists(SCRIPT_PATH_HSC2)) {
  file.copy(SCRIPT_PATH_HSC2,
            file.path(OUTDIR_HSC2, basename(SCRIPT_PATH_HSC2)),
            overwrite = TRUE)
  message("Script copied to: ", OUTDIR_HSC2)
} else {
  warning("SCRIPT_PATH_HSC2 not found — update the path at top of file: ",
          SCRIPT_PATH_HSC2)
}

####################################################

library(dplyr)
library(ggplot2)
library(ggVennDiagram)

# ── SETUP ───────────────────────────────────────────────────────────────────
OUTDIR_OV <- "C:/Users/fc809/Downloads/HB_signature_overlap_LTHSC"
dir.create(OUTDIR_OV, showWarnings = FALSE, recursive = TRUE)

SCRIPT_PATH_OV <- "C:/Users/fc809/Downloads/LTHSC_HB_signature_overlap.R"

PVAL_CUTOFF <- 0.05

DEG_ROOT_HSC <- "C:/Users/fc809/Downloads/LT-HSCs (1)"

deg_files_hsc <- c(
  lentineg_hblb = file.path(DEG_ROOT_HSC, "LentiNegHB_vs_LentiNegLB",
                            "DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.3.sce_0 - Copy.txt"),
  lentipos_hblb = file.path(DEG_ROOT_HSC, "LentiPosHB_vs_LentiPosLB",
                            "DEG_LentiPosHB_vs_LentiPosLB_pseudoBulk_res0.3.sce_0.txt")
)

# ── LOADER ──────────────────────────────────────────────────────────────────
load_deg <- function(path) {
  if (!file.exists(path)) stop("DEG file not found: ", path)
  df <- tryCatch(
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  )
  bad <- is.na(names(df)) | names(df) == ""
  if (any(bad)) {
    if (bad[1]) { names(df)[1] <- "X"; bad[1] <- FALSE }
    names(df)[bad] <- paste0("V", seq_len(sum(bad)))
  }
  if (!"X" %in% colnames(df)) {
    gc <- intersect(c("gene", "Gene", "symbol", "Symbol", "gene_name"),
                    colnames(df))
    if (length(gc) > 0) {
      names(df)[names(df) == gc[1]] <- "X"
    } else if (!identical(rownames(df), as.character(seq_len(nrow(df))))) {
      df$X <- rownames(df)
    }
  }
  df
}

for (nm in names(deg_files_hsc)) assign(nm, load_deg(deg_files_hsc[[nm]]))

# ── Identify p-value column ─────────────────────────────────────────────────
find_pval_col <- function(df) {
  intersect(c("PValue", "pvalue", "p.value", "P.Value", "pval"), colnames(df))[1]
}
pcol_ln <- find_pval_col(lentineg_hblb)
pcol_lp <- find_pval_col(lentipos_hblb)
if (is.na(pcol_ln) || is.na(pcol_lp)) stop("p-value column not found")
cat("p column in LentiNeg table:", pcol_ln, "\n")
cat("p column in LentiPos table:", pcol_lp, "\n")

# ── Build the four gene sets ────────────────────────────────────────────────
build_sets <- function(df, pcol, tag) {
  sig <- df %>%
    dplyr::filter(!is.na(.data[[pcol]]), !is.na(logFC),
                  .data[[pcol]] < PVAL_CUTOFF) %>%
    dplyr::distinct(X, .keep_all = TRUE)
  up   <- sig %>% dplyr::filter(logFC > 0) %>% dplyr::pull(X)
  down <- sig %>% dplyr::filter(logFC < 0) %>% dplyr::pull(X)
  cat(tag, "— total:", nrow(sig), " up:", length(up), " down:", length(down), "\n")
  list(up = up, down = down)
}

ln_sets <- build_sets(lentineg_hblb, pcol_ln, "LentiNeg HB-vs-LB")
lp_sets <- build_sets(lentipos_hblb, pcol_lp, "LentiPos HB-vs-LB")

# ── Overlap stats ───────────────────────────────────────────────────────────
jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

summarise_overlap <- function(a, b, name_a, name_b) {
  shared <- length(intersect(a, b))
  data.frame(
    comparison  = paste(name_a, "vs", name_b),
    only_a      = length(setdiff(a, b)),
    shared      = shared,
    only_b      = length(setdiff(b, a)),
    jaccard     = round(jaccard(a, b), 3),
    pct_a_in_b  = round(100 * shared / length(a), 1),
    pct_b_in_a  = round(100 * shared / length(b), 1)
  )
}

overlap_tbl <- rbind(
  summarise_overlap(ln_sets$up,   lp_sets$up,   "LentiNeg_up",   "LentiPos_up"),
  summarise_overlap(ln_sets$down, lp_sets$down, "LentiNeg_down", "LentiPos_down"),
  summarise_overlap(ln_sets$up,   lp_sets$down, "LentiNeg_up",   "LentiPos_down"),
  summarise_overlap(ln_sets$down, lp_sets$up,   "LentiNeg_down", "LentiPos_up")
)

cat("\nOverlap summary:\n"); print(overlap_tbl)
write.csv(overlap_tbl,
          file.path(OUTDIR_OV, "Overlap_summary.csv"), row.names = FALSE)

write.csv(data.frame(gene = intersect(ln_sets$up, lp_sets$up)),
          file.path(OUTDIR_OV, "Shared_HB_up.csv"), row.names = FALSE)
write.csv(data.frame(gene = intersect(ln_sets$down, lp_sets$down)),
          file.path(OUTDIR_OV, "Shared_HB_down.csv"), row.names = FALSE)
write.csv(data.frame(gene = setdiff(ln_sets$up, lp_sets$up)),
          file.path(OUTDIR_OV, "LentiNeg_only_HB_up.csv"), row.names = FALSE)
write.csv(data.frame(gene = setdiff(ln_sets$down, lp_sets$down)),
          file.path(OUTDIR_OV, "LentiNeg_only_HB_down.csv"), row.names = FALSE)
write.csv(data.frame(gene = setdiff(lp_sets$up, ln_sets$up)),
          file.path(OUTDIR_OV, "LentiPos_only_HB_up.csv"), row.names = FALSE)
write.csv(data.frame(gene = setdiff(lp_sets$down, ln_sets$down)),
          file.path(OUTDIR_OV, "LentiPos_only_HB_down.csv"), row.names = FALSE)

# ── Venn diagram: HB_up ─────────────────────────────────────────────────────
venn_up <- list(
  `LentiNeg HB_up`  = ln_sets$up,
  `LentiPos HB_up`  = lp_sets$up
)
n_shared_up  <- length(intersect(ln_sets$up, lp_sets$up))
n_only_ln_up <- length(setdiff(ln_sets$up, lp_sets$up))
n_only_lp_up <- length(setdiff(lp_sets$up, ln_sets$up))

p_venn_up <- ggVennDiagram(venn_up, label_alpha = 0, set_size = 3.5) +
  scale_fill_gradient(low = "white", high = "#1A237E") +
  labs(title    = "HB_up gene set overlap (LT-HSCs)",
       subtitle = paste0("LentiNeg-only: ", n_only_ln_up,
                         "  |  Shared: ", n_shared_up,
                         "  |  LentiPos-only: ", n_only_lp_up,
                         "  |  Jaccard: ",
                         round(jaccard(ln_sets$up, lp_sets$up), 3))) +
  theme(legend.position = "none",
        plot.title      = element_text(hjust = 0.5, size = 13),
        plot.subtitle   = element_text(hjust = 0.5, size = 10),
        plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(OUTDIR_OV, "Venn_HB_up.png"),
       p_venn_up, width = 7, height = 6, dpi = 150, bg = "white")
message("Saved: Venn_HB_up.png")

# ── Venn diagram: HB_down ───────────────────────────────────────────────────
venn_down <- list(
  `LentiNeg HB_down` = ln_sets$down,
  `LentiPos HB_down` = lp_sets$down
)
n_shared_dn  <- length(intersect(ln_sets$down, lp_sets$down))
n_only_ln_dn <- length(setdiff(ln_sets$down, lp_sets$down))
n_only_lp_dn <- length(setdiff(lp_sets$down, ln_sets$down))

p_venn_down <- ggVennDiagram(venn_down, label_alpha = 0, set_size = 3.5) +
  scale_fill_gradient(low = "white", high = "#B71C1C") +
  labs(title    = "HB_down gene set overlap (LT-HSCs)",
       subtitle = paste0("LentiNeg-only: ", n_only_ln_dn,
                         "  |  Shared: ", n_shared_dn,
                         "  |  LentiPos-only: ", n_only_lp_dn,
                         "  |  Jaccard: ",
                         round(jaccard(ln_sets$down, lp_sets$down), 3))) +
  theme(legend.position = "none",
        plot.title      = element_text(hjust = 0.5, size = 13),
        plot.subtitle   = element_text(hjust = 0.5, size = 10),
        plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(OUTDIR_OV, "Venn_HB_down.png"),
       p_venn_down, width = 7, height = 6, dpi = 150, bg = "white")
message("Saved: Venn_HB_down.png")

# ── 4-way Venn ──────────────────────────────────────────────────────────────
venn_all <- list(
  `LentiNeg HB_up`   = ln_sets$up,
  `LentiPos HB_up`   = lp_sets$up,
  `LentiNeg HB_down` = ln_sets$down,
  `LentiPos HB_down` = lp_sets$down
)

p_venn_all <- ggVennDiagram(venn_all, label_alpha = 0, set_size = 3) +
  scale_fill_gradient(low = "white", high = "#5E3A8C") +
  labs(title    = "All four HB-signature gene sets (LT-HSCs)",
       subtitle = paste0("Cross-direction overlaps should be near zero ",
                         "(genes shouldn't flip sign)")) +
  theme(legend.position = "none",
        plot.title      = element_text(hjust = 0.5, size = 13),
        plot.subtitle   = element_text(hjust = 0.5, size = 9),
        plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(OUTDIR_OV, "Venn_all4.png"),
       p_venn_all, width = 8, height = 7, dpi = 150, bg = "white")
message("Saved: Venn_all4.png")

# ── Copy script ─────────────────────────────────────────────────────────────
if (file.exists(SCRIPT_PATH_OV)) {
  file.copy(SCRIPT_PATH_OV,
            file.path(OUTDIR_OV, basename(SCRIPT_PATH_OV)),
            overwrite = TRUE)
  message("Script copied to: ", OUTDIR_OV)
} else {
  warning("SCRIPT_PATH_OV not found — update path: ", SCRIPT_PATH_OV)
}

message("\nDone. Outputs in: ", OUTDIR_OV)

########################################################################3

library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)

# ── SETUP ───────────────────────────────────────────────────────────────────
OUTDIR_SH <- "C:/Users/fc809/Downloads/HB_shared_signature_vs_Lenti_LTHSC"
dir.create(OUTDIR_SH, showWarnings = FALSE, recursive = TRUE)

SCRIPT_PATH_SH <- "C:/Users/fc809/Downloads/LTHSC_HB_shared_vs_Lenti.R"

PVAL_CUTOFF <- 0.05

DEG_ROOT_HSC <- "C:/Users/fc809/Downloads/LT-HSCs (1)"

deg_files_hsc <- c(
  lentineg_hblb = file.path(DEG_ROOT_HSC, "LentiNegHB_vs_LentiNegLB",
                            "DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.3.sce_0 - Copy.txt"),
  lentipos_hblb = file.path(DEG_ROOT_HSC, "LentiPosHB_vs_LentiPosLB",
                            "DEG_LentiPosHB_vs_LentiPosLB_pseudoBulk_res0.3.sce_0.txt"),
  lenti_lb      = file.path(DEG_ROOT_HSC, "LentiPos_vs_LentiNeg_in_LB",
                            "DEG_LentiPos_vs_LentiNeg_in_LB_pseudoBulk_res0.3.sce_0.txt"),
  lenti_hb      = file.path(DEG_ROOT_HSC, "LentiPos_vs_LentiNeg_in_HB",
                            "DEG_LentiPos_vs_LentiNeg_in_HB_pseudoBulk_res0.3.sce_0.txt")
)

# ── LOADER ──────────────────────────────────────────────────────────────────
load_deg <- function(path) {
  if (!file.exists(path)) stop("DEG file not found: ", path)
  df <- tryCatch(
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  )
  bad <- is.na(names(df)) | names(df) == ""
  if (any(bad)) {
    if (bad[1]) { names(df)[1] <- "X"; bad[1] <- FALSE }
    names(df)[bad] <- paste0("V", seq_len(sum(bad)))
  }
  if (!"X" %in% colnames(df)) {
    gc <- intersect(c("gene", "Gene", "symbol", "Symbol", "gene_name"),
                    colnames(df))
    if (length(gc) > 0) {
      names(df)[names(df) == gc[1]] <- "X"
    } else if (!identical(rownames(df), as.character(seq_len(nrow(df))))) {
      df$X <- rownames(df)
    }
  }
  df
}

for (nm in names(deg_files_hsc)) assign(nm, load_deg(deg_files_hsc[[nm]]))

find_pval_col <- function(df) {
  intersect(c("PValue", "pvalue", "p.value", "P.Value", "pval"), colnames(df))[1]
}
pcol_ln <- find_pval_col(lentineg_hblb)
pcol_lp <- find_pval_col(lentipos_hblb)
if (is.na(pcol_ln) || is.na(pcol_lp)) stop("p-value column not found")

# ── Build SHARED HB sets (intersection of LentiNeg + LentiPos signatures) ───
build_dir <- function(df, pcol, sign_fn) {
  df %>%
    dplyr::filter(!is.na(.data[[pcol]]), !is.na(logFC),
                  .data[[pcol]] < PVAL_CUTOFF,
                  sign_fn(logFC)) %>%
    dplyr::distinct(X, .keep_all = TRUE) %>%
    dplyr::pull(X)
}

ln_up   <- build_dir(lentineg_hblb, pcol_ln, function(x) x > 0)
ln_down <- build_dir(lentineg_hblb, pcol_ln, function(x) x < 0)
lp_up   <- build_dir(lentipos_hblb, pcol_lp, function(x) x > 0)
lp_down <- build_dir(lentipos_hblb, pcol_lp, function(x) x < 0)

hb_up   <- intersect(ln_up,   lp_up)    # higher in HB in both contrasts
hb_down <- intersect(ln_down, lp_down)  # higher in LB in both contrasts

N_HB_UP   <- length(hb_up)
N_HB_DOWN <- length(hb_down)

cat("\nShared HB signature (intersection of LentiNeg + LentiPos, p<", PVAL_CUTOFF, "):\n", sep = "")
cat("  LentiNeg up:", length(ln_up),   " LentiPos up:", length(lp_up),
    " -> Shared HB_up:",   N_HB_UP,   "\n")
cat("  LentiNeg down:", length(ln_down), " LentiPos down:", length(lp_down),
    " -> Shared HB_down:", N_HB_DOWN, "\n")

hb_sets <- data.frame(
  term = c(rep("HB_up", N_HB_UP), rep("HB_down", N_HB_DOWN)),
  gene = c(hb_up, hb_down)
)

# ── GSEA per Lenti+/- contrast ──────────────────────────────────────────────
run_hb_on_lenti <- function(deg_tbl, tag) {
  if (N_HB_UP < 10 && N_HB_DOWN < 10) {
    warning("[", tag, "] shared HB sets too small. Skipping.")
    return(invisible(NULL))
  }
  ranked <- deg_tbl %>%
    dplyr::filter(X != "Lenti") %>%
    dplyr::arrange(desc(logFC)) %>%
    dplyr::distinct(X, .keep_all = TRUE) %>%
    { setNames(.$logFC, .$X) }
  cat("\n[", tag, "] Ranked list:", length(ranked), "genes\n")
  
  gsea_obj <- GSEA(geneList = ranked, TERM2GENE = hb_sets,
                   minGSSize = 10,
                   maxGSSize = max(5000, N_HB_UP, N_HB_DOWN) + 1,
                   pvalueCutoff = 1, pAdjustMethod = "BH",
                   eps = 0, verbose = TRUE)
  
  stats <- as.data.frame(gsea_obj)[, c("ID", "NES", "pvalue", "p.adjust")]
  cat("--- [", tag, "] ---\n"); print(stats)
  
  plot_ids <- intersect(c("HB_up", "HB_down"), stats$ID)
  if (length(plot_ids) == 0) return(invisible(gsea_obj))
  
  p <- gseaplot2(gsea_obj, geneSetID = plot_ids,
                 title = paste0("Shared HB signature (LT-HSCs) in Lenti+/- within ", tag,
                                "\n(Shared sets: LentiNeg ∩ LentiPos, p<", PVAL_CUTOFF,
                                " | HB_up n=", N_HB_UP,
                                ", HB_down n=", N_HB_DOWN, ")"),
                 subplots = 1:3, base_size = 11)
  
  col_map <- c("HB_up" = "#1A237E", "HB_down" = "#B71C1C")[plot_ids]
  lab_map <- setNames(sapply(plot_ids, function(id) {
    n <- if (id == "HB_up") N_HB_UP else N_HB_DOWN
    pretty <- if (id == "HB_up") "HB-up" else "HB-down"
    paste0(pretty, " (n=", n,
           ") NES=", round(stats$NES[stats$ID == id], 2),
           " FDR=",  signif(stats$p.adjust[stats$ID == id], 2))
  }), plot_ids)
  
  p[[1]] <- p[[1]] +
    scale_color_manual(values = col_map, labels = lab_map) +
    theme(legend.position = "top", legend.text = element_text(size = 9))
  p[[2]] <- p[[2]] + scale_color_manual(values = col_map)
  
  outfile <- file.path(OUTDIR_SH, paste0("GSEA_SharedHB_on_Lenti_", tag, ".png"))
  png(outfile, width = 1400, height = 1200, res = 150, bg = "white")
  print(p); dev.off()
  message("Saved: ", basename(outfile))
  invisible(gsea_obj)
}

gsea_sh_on_lb <- run_hb_on_lenti(lenti_lb, "LB")
gsea_sh_on_hb <- run_hb_on_lenti(lenti_hb, "HB")

# ── Summary dotplot ─────────────────────────────────────────────────────────
extract_gsea <- function(gsea_obj, comparison, fallback_ids) {
  if (is.null(gsea_obj)) {
    return(data.frame(ID = fallback_ids, NES = NA_real_,
                      pvalue = NA_real_, p.adjust = NA_real_,
                      comparison = comparison))
  }
  d <- as.data.frame(gsea_obj)[, c("ID", "NES", "pvalue", "p.adjust")]
  missing_ids <- setdiff(fallback_ids, d$ID)
  if (length(missing_ids) > 0) {
    d <- rbind(d, data.frame(ID = missing_ids, NES = NA,
                             pvalue = NA, p.adjust = NA))
  }
  d$comparison <- comparison
  d
}

hb_summary_df <- rbind(
  extract_gsea(gsea_sh_on_lb, "Lenti+/- in LB", c("HB_up", "HB_down")),
  extract_gsea(gsea_sh_on_hb, "Lenti+/- in HB", c("HB_up", "HB_down"))
)
hb_summary_df$comparison <- factor(hb_summary_df$comparison,
                                   levels = c("Lenti+/- in LB", "Lenti+/- in HB"))
hb_summary_df$ID <- factor(hb_summary_df$ID, levels = c("HB_up", "HB_down"))

hb_y_labels <- c(
  HB_up   = paste0("HB_up\n(n=",   N_HB_UP,   ")"),
  HB_down = paste0("HB_down\n(n=", N_HB_DOWN, ")")
)

cat("\nSummary dotplot data (shared HB signature):\n"); print(hb_summary_df)

p_hb_summary <- ggplot(hb_summary_df,
                       aes(x = comparison, y = ID, fill = NES,
                           size = -log10(p.adjust))) +
  geom_point(shape = 21, color = "grey30") +
  scale_fill_gradient2(low = "#1A237E", mid = "white", high = "#B71C1C",
                       midpoint = 0, limits = c(-2.5, 2.5),
                       oob = scales::squish, name = "NES") +
  scale_size_continuous(range = c(4, 14), name = "FDR",
                        breaks = -log10(c(0.5, 0.1, 0.05, 0.01, 0.001)),
                        labels = c("0.5", "0.1", "0.05", "0.01", "0.001")) +
  scale_y_discrete(labels = hb_y_labels, expand = expansion(add = 0.8)) +
  theme_classic(base_size = 13) +
  theme(axis.text.x   = element_text(size = 10, angle = 30, hjust = 1,
                                     margin = margin(t = 4)),
        axis.text.y   = element_text(size = 11),
        plot.title    = element_text(size = 12, hjust = 0.5,
                                     margin = margin(b = 10)),
        plot.subtitle = element_text(size = 10, hjust = 0.5,
                                     margin = margin(b = 8)),
        plot.margin   = margin(t = 15, r = 10, b = 25, l = 15),
        plot.background = element_rect(fill = "white", color = NA)) +
  labs(x = NULL, y = "HB gene set",
       title    = "Shared HB signature (LT-HSCs) across Lenti+/- contrasts",
       subtitle = paste0("Gene sets: intersection of LentiNeg ∩ LentiPos HB-vs-LB DEGs, p < ",
                         PVAL_CUTOFF,
                         "  |  HB_up n=", N_HB_UP,
                         ", HB_down n=", N_HB_DOWN))

ggsave(file.path(OUTDIR_SH, "Summary_dotplot_SharedHB_on_Lenti.png"),
       p_hb_summary, width = 8, height = 5, dpi = 150, bg = "white")
message("Saved: Summary_dotplot_SharedHB_on_Lenti.png")

# ── Export the shared gene lists for convenience ────────────────────────────
write.csv(data.frame(gene = hb_up),
          file.path(OUTDIR_SH, "Shared_HB_up.csv"), row.names = FALSE)
write.csv(data.frame(gene = hb_down),
          file.path(OUTDIR_SH, "Shared_HB_down.csv"), row.names = FALSE)

# ── Copy script for provenance ──────────────────────────────────────────────
if (file.exists(SCRIPT_PATH_SH)) {
  file.copy(SCRIPT_PATH_SH,
            file.path(OUTDIR_SH, basename(SCRIPT_PATH_SH)),
            overwrite = TRUE)
  message("Script copied to: ", OUTDIR_SH)
} else {
  warning("SCRIPT_PATH_SH not found — update path: ", SCRIPT_PATH_SH)
}

message("\nDone. Outputs in: ", OUTDIR_SH)

##############################

library(dplyr)
library(ggplot2)
library(ggVennDiagram)

# ── SETUP ───────────────────────────────────────────────────────────────────
OUTDIR_MOV <- "C:/Users/fc809/Downloads/HB_signature_overlap_Monocytes"
dir.create(OUTDIR_MOV, showWarnings = FALSE, recursive = TRUE)

SCRIPT_PATH_MOV <- "C:/Users/fc809/Downloads/Mono_HB_signature_overlap.R"

PVAL_CUTOFF <- 0.05

DEG_ROOT_MONO <- "C:/Users/fc809/Downloads/Classical_Monocytes (1)"

deg_files_mono <- c(
  lentineg_hblb = file.path(DEG_ROOT_MONO, "LentiNegHB_vs_LentiNegLB",
                            "DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.6.sce_0.txt"),
  lentipos_hblb = file.path(DEG_ROOT_MONO, "LentiPosHB_vs_LentiPosLB",
                            "DEG_LentiPosHB_vs_LentiPosLB_pseudoBulk_res0.6.sce_0.txt")
)

# ── LOADER ──────────────────────────────────────────────────────────────────
load_deg <- function(path) {
  if (!file.exists(path)) stop("DEG file not found: ", path)
  df <- tryCatch(
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  )
  bad <- is.na(names(df)) | names(df) == ""
  if (any(bad)) {
    if (bad[1]) { names(df)[1] <- "X"; bad[1] <- FALSE }
    names(df)[bad] <- paste0("V", seq_len(sum(bad)))
  }
  if (!"X" %in% colnames(df)) {
    gc <- intersect(c("gene", "Gene", "symbol", "Symbol", "gene_name"),
                    colnames(df))
    if (length(gc) > 0) {
      names(df)[names(df) == gc[1]] <- "X"
    } else if (!identical(rownames(df), as.character(seq_len(nrow(df))))) {
      df$X <- rownames(df)
    }
  }
  df
}

for (nm in names(deg_files_mono)) assign(nm, load_deg(deg_files_mono[[nm]]))

find_pval_col <- function(df) {
  intersect(c("PValue", "pvalue", "p.value", "P.Value", "pval"), colnames(df))[1]
}
pcol_ln <- find_pval_col(lentineg_hblb)
pcol_lp <- find_pval_col(lentipos_hblb)
if (is.na(pcol_ln) || is.na(pcol_lp)) stop("p-value column not found")
cat("p column in LentiNeg table:", pcol_ln, "\n")
cat("p column in LentiPos table:", pcol_lp, "\n")

build_sets <- function(df, pcol, tag) {
  sig <- df %>%
    dplyr::filter(!is.na(.data[[pcol]]), !is.na(logFC),
                  .data[[pcol]] < PVAL_CUTOFF) %>%
    dplyr::distinct(X, .keep_all = TRUE)
  up   <- sig %>% dplyr::filter(logFC > 0) %>% dplyr::pull(X)
  down <- sig %>% dplyr::filter(logFC < 0) %>% dplyr::pull(X)
  cat(tag, "— total:", nrow(sig), " up:", length(up), " down:", length(down), "\n")
  list(up = up, down = down)
}

ln_sets <- build_sets(lentineg_hblb, pcol_ln, "LentiNeg HB-vs-LB")
lp_sets <- build_sets(lentipos_hblb, pcol_lp, "LentiPos HB-vs-LB")

jaccard <- function(a, b) {
  if (length(union(a, b)) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

summarise_overlap <- function(a, b, name_a, name_b) {
  shared <- length(intersect(a, b))
  data.frame(
    comparison  = paste(name_a, "vs", name_b),
    only_a      = length(setdiff(a, b)),
    shared      = shared,
    only_b      = length(setdiff(b, a)),
    jaccard     = round(jaccard(a, b), 3),
    pct_a_in_b  = round(100 * shared / max(1, length(a)), 1),
    pct_b_in_a  = round(100 * shared / max(1, length(b)), 1)
  )
}

overlap_tbl <- rbind(
  summarise_overlap(ln_sets$up,   lp_sets$up,   "LentiNeg_up",   "LentiPos_up"),
  summarise_overlap(ln_sets$down, lp_sets$down, "LentiNeg_down", "LentiPos_down"),
  summarise_overlap(ln_sets$up,   lp_sets$down, "LentiNeg_up",   "LentiPos_down"),
  summarise_overlap(ln_sets$down, lp_sets$up,   "LentiNeg_down", "LentiPos_up")
)

cat("\nOverlap summary:\n"); print(overlap_tbl)
write.csv(overlap_tbl, file.path(OUTDIR_MOV, "Overlap_summary.csv"), row.names = FALSE)

write.csv(data.frame(gene = intersect(ln_sets$up, lp_sets$up)),
          file.path(OUTDIR_MOV, "Shared_HB_up.csv"), row.names = FALSE)
write.csv(data.frame(gene = intersect(ln_sets$down, lp_sets$down)),
          file.path(OUTDIR_MOV, "Shared_HB_down.csv"), row.names = FALSE)
write.csv(data.frame(gene = setdiff(ln_sets$up, lp_sets$up)),
          file.path(OUTDIR_MOV, "LentiNeg_only_HB_up.csv"), row.names = FALSE)
write.csv(data.frame(gene = setdiff(ln_sets$down, lp_sets$down)),
          file.path(OUTDIR_MOV, "LentiNeg_only_HB_down.csv"), row.names = FALSE)
write.csv(data.frame(gene = setdiff(lp_sets$up, ln_sets$up)),
          file.path(OUTDIR_MOV, "LentiPos_only_HB_up.csv"), row.names = FALSE)
write.csv(data.frame(gene = setdiff(lp_sets$down, ln_sets$down)),
          file.path(OUTDIR_MOV, "LentiPos_only_HB_down.csv"), row.names = FALSE)

# ── Venn plots ──────────────────────────────────────────────────────────────
venn_up <- list(
  `LentiNeg HB_up` = ln_sets$up,
  `LentiPos HB_up` = lp_sets$up
)

p_venn_up <- ggVennDiagram(venn_up, label_alpha = 0, set_size = 3.5) +
  scale_fill_gradient(low = "white", high = "#1A237E") +
  labs(title    = "HB_up gene set overlap (classical monocytes)",
       subtitle = paste0("LentiNeg-only: ", length(setdiff(ln_sets$up, lp_sets$up)),
                         "  |  Shared: ", length(intersect(ln_sets$up, lp_sets$up)),
                         "  |  LentiPos-only: ", length(setdiff(lp_sets$up, ln_sets$up)),
                         "  |  Jaccard: ", round(jaccard(ln_sets$up, lp_sets$up), 3))) +
  theme(legend.position = "none",
        plot.title      = element_text(hjust = 0.5, size = 13),
        plot.subtitle   = element_text(hjust = 0.5, size = 10),
        plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(OUTDIR_MOV, "Venn_HB_up.png"),
       p_venn_up, width = 7, height = 6, dpi = 150, bg = "white")
message("Saved: Venn_HB_up.png")

venn_down <- list(
  `LentiNeg HB_down` = ln_sets$down,
  `LentiPos HB_down` = lp_sets$down
)

p_venn_down <- ggVennDiagram(venn_down, label_alpha = 0, set_size = 3.5) +
  scale_fill_gradient(low = "white", high = "#B71C1C") +
  labs(title    = "HB_down gene set overlap (classical monocytes)",
       subtitle = paste0("LentiNeg-only: ", length(setdiff(ln_sets$down, lp_sets$down)),
                         "  |  Shared: ", length(intersect(ln_sets$down, lp_sets$down)),
                         "  |  LentiPos-only: ", length(setdiff(lp_sets$down, ln_sets$down)),
                         "  |  Jaccard: ", round(jaccard(ln_sets$down, lp_sets$down), 3))) +
  theme(legend.position = "none",
        plot.title      = element_text(hjust = 0.5, size = 13),
        plot.subtitle   = element_text(hjust = 0.5, size = 10),
        plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(OUTDIR_MOV, "Venn_HB_down.png"),
       p_venn_down, width = 7, height = 6, dpi = 150, bg = "white")
message("Saved: Venn_HB_down.png")

venn_all <- list(
  `LentiNeg HB_up`   = ln_sets$up,
  `LentiPos HB_up`   = lp_sets$up,
  `LentiNeg HB_down` = ln_sets$down,
  `LentiPos HB_down` = lp_sets$down
)

p_venn_all <- ggVennDiagram(venn_all, label_alpha = 0, set_size = 3) +
  scale_fill_gradient(low = "white", high = "#5E3A8C") +
  labs(title    = "All four HB-signature gene sets (classical monocytes)",
       subtitle = "Cross-direction overlaps should be near zero") +
  theme(legend.position = "none",
        plot.title      = element_text(hjust = 0.5, size = 13),
        plot.subtitle   = element_text(hjust = 0.5, size = 9),
        plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(OUTDIR_MOV, "Venn_all4.png"),
       p_venn_all, width = 8, height = 7, dpi = 150, bg = "white")
message("Saved: Venn_all4.png")

if (file.exists(SCRIPT_PATH_MOV)) {
  file.copy(SCRIPT_PATH_MOV,
            file.path(OUTDIR_MOV, basename(SCRIPT_PATH_MOV)),
            overwrite = TRUE)
  message("Script copied to: ", OUTDIR_MOV)
} else {
  warning("SCRIPT_PATH_MOV not found — update path: ", SCRIPT_PATH_MOV)
}

message("\nDone. Outputs in: ", OUTDIR_MOV)

###########################################


library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)

# ── SETUP ───────────────────────────────────────────────────────────────────
OUTDIR_MSH <- "C:/Users/fc809/Downloads/HB_shared_signature_vs_Lenti_Monocytes"
dir.create(OUTDIR_MSH, showWarnings = FALSE, recursive = TRUE)

SCRIPT_PATH_MSH <- "C:/Users/fc809/Downloads/Mono_HB_shared_vs_Lenti.R"

PVAL_CUTOFF <- 0.05

DEG_ROOT_MONO <- "C:/Users/fc809/Downloads/Classical_Monocytes (1)"

deg_files_mono <- c(
  lentineg_hblb = file.path(DEG_ROOT_MONO, "LentiNegHB_vs_LentiNegLB",
                            "DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.6.sce_0.txt"),
  lentipos_hblb = file.path(DEG_ROOT_MONO, "LentiPosHB_vs_LentiPosLB",
                            "DEG_LentiPosHB_vs_LentiPosLB_pseudoBulk_res0.6.sce_0.txt"),
  lenti_lb      = file.path(DEG_ROOT_MONO, "LentiPos_vs_LentiNeg_in_LB",
                            "DEG_LentiPos_vs_LentiNeg_in_LB_pseudoBulk_res0.6.sce_0.txt"),
  lenti_hb      = file.path(DEG_ROOT_MONO, "LentiPos_vs_LentiNeg_in_HB",
                            "DEG_LentiPos_vs_LentiNeg_in_HB_pseudoBulk_res0.6.sce_0.txt")
)

load_deg <- function(path) {
  if (!file.exists(path)) stop("DEG file not found: ", path)
  df <- tryCatch(
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  )
  bad <- is.na(names(df)) | names(df) == ""
  if (any(bad)) {
    if (bad[1]) { names(df)[1] <- "X"; bad[1] <- FALSE }
    names(df)[bad] <- paste0("V", seq_len(sum(bad)))
  }
  if (!"X" %in% colnames(df)) {
    gc <- intersect(c("gene", "Gene", "symbol", "Symbol", "gene_name"),
                    colnames(df))
    if (length(gc) > 0) {
      names(df)[names(df) == gc[1]] <- "X"
    } else if (!identical(rownames(df), as.character(seq_len(nrow(df))))) {
      df$X <- rownames(df)
    }
  }
  df
}

for (nm in names(deg_files_mono)) assign(nm, load_deg(deg_files_mono[[nm]]))

find_pval_col <- function(df) {
  intersect(c("PValue", "pvalue", "p.value", "P.Value", "pval"), colnames(df))[1]
}
pcol_ln <- find_pval_col(lentineg_hblb)
pcol_lp <- find_pval_col(lentipos_hblb)

build_dir <- function(df, pcol, sign_fn) {
  df %>%
    dplyr::filter(!is.na(.data[[pcol]]), !is.na(logFC),
                  .data[[pcol]] < PVAL_CUTOFF,
                  sign_fn(logFC)) %>%
    dplyr::distinct(X, .keep_all = TRUE) %>%
    dplyr::pull(X)
}

ln_up   <- build_dir(lentineg_hblb, pcol_ln, function(x) x > 0)
ln_down <- build_dir(lentineg_hblb, pcol_ln, function(x) x < 0)
lp_up   <- build_dir(lentipos_hblb, pcol_lp, function(x) x > 0)
lp_down <- build_dir(lentipos_hblb, pcol_lp, function(x) x < 0)

hb_up   <- intersect(ln_up,   lp_up)
hb_down <- intersect(ln_down, lp_down)

N_HB_UP   <- length(hb_up)
N_HB_DOWN <- length(hb_down)

cat("\nShared HB signature (classical monocytes, p<", PVAL_CUTOFF, ")\n", sep = "")
cat("  LentiNeg up:",   length(ln_up),   " LentiPos up:",   length(lp_up),
    " -> Shared HB_up:",   N_HB_UP,   "\n")
cat("  LentiNeg down:", length(ln_down), " LentiPos down:", length(lp_down),
    " -> Shared HB_down:", N_HB_DOWN, "\n")

if (N_HB_UP + N_HB_DOWN == 0)
  stop("Shared sets are empty — can't run GSEA.")

hb_sets <- data.frame(
  term = c(rep("HB_up", N_HB_UP), rep("HB_down", N_HB_DOWN)),
  gene = c(hb_up, hb_down)
)

run_hb_on_lenti <- function(deg_tbl, tag) {
  if (N_HB_UP < 10 && N_HB_DOWN < 10) {
    warning("[", tag, "] shared HB sets too small. Skipping.")
    return(invisible(NULL))
  }
  ranked <- deg_tbl %>%
    dplyr::filter(X != "Lenti") %>%
    dplyr::arrange(desc(logFC)) %>%
    dplyr::distinct(X, .keep_all = TRUE) %>%
    { setNames(.$logFC, .$X) }
  cat("\n[", tag, "] Ranked list:", length(ranked), "genes\n")
  
  gsea_obj <- GSEA(geneList = ranked, TERM2GENE = hb_sets,
                   minGSSize = 10,
                   maxGSSize = max(5000, N_HB_UP, N_HB_DOWN) + 1,
                   pvalueCutoff = 1, pAdjustMethod = "BH",
                   eps = 0, verbose = TRUE)
  
  stats <- as.data.frame(gsea_obj)[, c("ID", "NES", "pvalue", "p.adjust")]
  cat("--- [", tag, "] ---\n"); print(stats)
  
  plot_ids <- intersect(c("HB_up", "HB_down"), stats$ID)
  if (length(plot_ids) == 0) return(invisible(gsea_obj))
  
  p <- gseaplot2(gsea_obj, geneSetID = plot_ids,
                 title = paste0("Shared HB signature (monocytes) in Lenti+/- within ", tag,
                                "\n(Shared sets: LentiNeg ∩ LentiPos, p<", PVAL_CUTOFF,
                                " | HB_up n=", N_HB_UP,
                                ", HB_down n=", N_HB_DOWN, ")"),
                 subplots = 1:3, base_size = 11)
  
  col_map <- c("HB_up" = "#1A237E", "HB_down" = "#B71C1C")[plot_ids]
  lab_map <- setNames(sapply(plot_ids, function(id) {
    n <- if (id == "HB_up") N_HB_UP else N_HB_DOWN
    pretty <- if (id == "HB_up") "HB-up" else "HB-down"
    paste0(pretty, " (n=", n,
           ") NES=", round(stats$NES[stats$ID == id], 2),
           " FDR=",  signif(stats$p.adjust[stats$ID == id], 2))
  }), plot_ids)
  
  p[[1]] <- p[[1]] +
    scale_color_manual(values = col_map, labels = lab_map) +
    theme(legend.position = "top", legend.text = element_text(size = 9))
  p[[2]] <- p[[2]] + scale_color_manual(values = col_map)
  
  outfile <- file.path(OUTDIR_MSH, paste0("GSEA_SharedHB_on_Lenti_", tag, ".png"))
  png(outfile, width = 1400, height = 1200, res = 150, bg = "white")
  print(p); dev.off()
  message("Saved: ", basename(outfile))
  invisible(gsea_obj)
}

gsea_sh_on_lb <- run_hb_on_lenti(lenti_lb, "LB")
gsea_sh_on_hb <- run_hb_on_lenti(lenti_hb, "HB")

extract_gsea <- function(gsea_obj, comparison, fallback_ids) {
  if (is.null(gsea_obj)) {
    return(data.frame(ID = fallback_ids, NES = NA_real_,
                      pvalue = NA_real_, p.adjust = NA_real_,
                      comparison = comparison))
  }
  d <- as.data.frame(gsea_obj)[, c("ID", "NES", "pvalue", "p.adjust")]
  missing_ids <- setdiff(fallback_ids, d$ID)
  if (length(missing_ids) > 0) {
    d <- rbind(d, data.frame(ID = missing_ids, NES = NA,
                             pvalue = NA, p.adjust = NA))
  }
  d$comparison <- comparison
  d
}

hb_summary_df <- rbind(
  extract_gsea(gsea_sh_on_lb, "Lenti+/- in LB", c("HB_up", "HB_down")),
  extract_gsea(gsea_sh_on_hb, "Lenti+/- in HB", c("HB_up", "HB_down"))
)
hb_summary_df$comparison <- factor(hb_summary_df$comparison,
                                   levels = c("Lenti+/- in LB", "Lenti+/- in HB"))
hb_summary_df$ID <- factor(hb_summary_df$ID, levels = c("HB_up", "HB_down"))

hb_y_labels <- c(
  HB_up   = paste0("HB_up\n(n=",   N_HB_UP,   ")"),
  HB_down = paste0("HB_down\n(n=", N_HB_DOWN, ")")
)

cat("\nSummary dotplot data (shared HB signature, monocytes):\n")
print(hb_summary_df)

p_hb_summary <- ggplot(hb_summary_df,
                       aes(x = comparison, y = ID, fill = NES,
                           size = -log10(p.adjust))) +
  geom_point(shape = 21, color = "grey30") +
  scale_fill_gradient2(low = "#1A237E", mid = "white", high = "#B71C1C",
                       midpoint = 0, limits = c(-2.5, 2.5),
                       oob = scales::squish, name = "NES") +
  scale_size_continuous(range = c(4, 14), name = "FDR",
                        breaks = -log10(c(0.5, 0.1, 0.05, 0.01, 0.001)),
                        labels = c("0.5", "0.1", "0.05", "0.01", "0.001")) +
  scale_y_discrete(labels = hb_y_labels, expand = expansion(add = 0.8)) +
  theme_classic(base_size = 13) +
  theme(axis.text.x   = element_text(size = 10, angle = 30, hjust = 1,
                                     margin = margin(t = 4)),
        axis.text.y   = element_text(size = 11),
        plot.title    = element_text(size = 12, hjust = 0.5,
                                     margin = margin(b = 10)),
        plot.subtitle = element_text(size = 10, hjust = 0.5,
                                     margin = margin(b = 8)),
        plot.margin   = margin(t = 15, r = 10, b = 25, l = 15),
        plot.background = element_rect(fill = "white", color = NA)) +
  labs(x = NULL, y = "HB gene set",
       title    = "Shared HB signature (classical monocytes) across Lenti+/- contrasts",
       subtitle = paste0("Gene sets: intersection of LentiNeg ∩ LentiPos HB-vs-LB DEGs, p < ",
                         PVAL_CUTOFF,
                         "  |  HB_up n=", N_HB_UP,
                         ", HB_down n=", N_HB_DOWN))

ggsave(file.path(OUTDIR_MSH, "Summary_dotplot_SharedHB_on_Lenti.png"),
       p_hb_summary, width = 8, height = 5, dpi = 150, bg = "white")
message("Saved: Summary_dotplot_SharedHB_on_Lenti.png")

write.csv(data.frame(gene = hb_up),
          file.path(OUTDIR_MSH, "Shared_HB_up.csv"), row.names = FALSE)
write.csv(data.frame(gene = hb_down),
          file.path(OUTDIR_MSH, "Shared_HB_down.csv"), row.names = FALSE)

if (file.exists(SCRIPT_PATH_MSH)) {
  file.copy(SCRIPT_PATH_MSH,
            file.path(OUTDIR_MSH, basename(SCRIPT_PATH_MSH)),
            overwrite = TRUE)
  message("Script copied to: ", OUTDIR_MSH)
} else {
  warning("SCRIPT_PATH_MSH not found — update path: ", SCRIPT_PATH_MSH)
}

message("\nDone. Outputs in: ", OUTDIR_MSH)
\
################################################


library(dplyr)
library(ggplot2)
library(readxl)
library(clusterProfiler)
library(enrichplot)
library(babelgene)   # install.packages("babelgene") if needed

# ── SETUP ───────────────────────────────────────────────────────────────────
OUTDIR_TET2 <- "C:/Users/fc809/Downloads/TET2_signature_vs_HBvsLB"
dir.create(OUTDIR_TET2, showWarnings = FALSE, recursive = TRUE)

SCRIPT_PATH_TET2 <- "C:/Users/fc809/Downloads/TET2_vs_HBvsLB.R"

PADJ_CUTOFF <- 0.05

# Paths to TET2 DEG (Jakobsen) Excel files
TET2_HSC  <- "C:/Users/fc809/Downloads/Jakobsen_hspc.xlsx"
TET2_MONO <- "C:/Users/fc809/Downloads/Jakobsen_monocytes.xlsx"

DEG_ROOT_HSC  <- "C:/Users/fc809/Downloads/LT-HSCs (1)"
DEG_ROOT_MONO <- "C:/Users/fc809/Downloads/Classical_Monocytes (1)"

# Your DEG files
deg_files_hsc <- c(
  hsc_lentineg_hblb = file.path(DEG_ROOT_HSC, "LentiNegHB_vs_LentiNegLB",
                                "DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.3.sce_0 - Copy.txt"),
  hsc_lentipos_hblb = file.path(DEG_ROOT_HSC, "LentiPosHB_vs_LentiPosLB",
                                "DEG_LentiPosHB_vs_LentiPosLB_pseudoBulk_res0.3.sce_0.txt")
)
deg_files_mono <- c(
  mono_lentineg_hblb = file.path(DEG_ROOT_MONO, "LentiNegHB_vs_LentiNegLB",
                                 "DEG_LentiNegHB_vs_LentiNegLB_pseudoBulk_res0.6.sce_0.txt"),
  mono_lentipos_hblb = file.path(DEG_ROOT_MONO, "LentiPosHB_vs_LentiPosLB",
                                 "DEG_LentiPosHB_vs_LentiPosLB_pseudoBulk_res0.6.sce_0.txt")
)

# ── LOADERS ─────────────────────────────────────────────────────────────────
load_deg <- function(path) {
  if (!file.exists(path)) stop("DEG file not found: ", path)
  df <- tryCatch(
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  )
  bad <- is.na(names(df)) | names(df) == ""
  if (any(bad)) {
    if (bad[1]) { names(df)[1] <- "X"; bad[1] <- FALSE }
    names(df)[bad] <- paste0("V", seq_len(sum(bad)))
  }
  if (!"X" %in% colnames(df)) {
    gc <- intersect(c("gene", "Gene", "symbol", "Symbol", "gene_name"), colnames(df))
    if (length(gc) > 0) names(df)[names(df) == gc[1]] <- "X"
  }
  df
}

load_tet2 <- function(path) {
  if (!file.exists(path)) stop("TET2 file not found: ", path)
  df <- as.data.frame(read_excel(path))
  cat("\nColumns in", basename(path), ":", paste(colnames(df), collapse = ", "), "\n")
  
  gene_col <- intersect(c("Gene", "gene", "symbol", "Symbol", "gene_name",
                          "GeneName", "Gene_Name", "GENE", "X"),
                        colnames(df))[1]
  z_col    <- intersect(c("z_score_ranking", "z", "Z", "z_score", "Z_score",
                          "zscore", "Zscore", "stat", "statistic"),
                        colnames(df))[1]
  fdr_col  <- intersect(c("FDR", "fdr", "padj", "p.adjust", "p_adj",
                          "adj.P.Val", "adj_pval"),
                        colnames(df))[1]
  
  if (is.na(gene_col)) stop("can't find gene column. Cols: ",
                            paste(colnames(df), collapse = ", "))
  if (is.na(z_col))    stop("can't find z-score column. Cols: ",
                            paste(colnames(df), collapse = ", "))
  if (is.na(fdr_col))  stop("can't find FDR column. Cols: ",
                            paste(colnames(df), collapse = ", "))
  
  cat("  -> using gene='", gene_col, "', z='", z_col, "', FDR='", fdr_col, "'\n", sep = "")
  
  data.frame(
    gene_human = as.character(df[[gene_col]]),
    z          = as.numeric(df[[z_col]]),
    FDR        = as.numeric(df[[fdr_col]]),
    stringsAsFactors = FALSE
  )
}

tet2_hsc_raw  <- load_tet2(TET2_HSC)
tet2_mono_raw <- load_tet2(TET2_MONO)

# ── HUMAN -> MOUSE ortholog mapping ─────────────────────────────────────────
human_to_mouse <- function(human_genes) {
  o <- babelgene::orthologs(genes = unique(human_genes),
                            species = "mouse",
                            human   = TRUE)
  o <- o %>%
    dplyr::filter(support_n >= 3) %>%
    dplyr::group_by(human_symbol) %>% dplyr::filter(dplyr::n() == 1) %>%
    dplyr::group_by(symbol)       %>% dplyr::filter(dplyr::n() == 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(human_symbol, mouse_symbol = symbol)
  o
}

build_tet2_sets <- function(tet2_df, tag) {
  sig <- tet2_df %>%
    dplyr::filter(!is.na(FDR), !is.na(z), FDR < PADJ_CUTOFF) %>%
    dplyr::distinct(gene_human, .keep_all = TRUE)
  
  cat("\n[", tag, "] TET2 DEGs at FDR <", PADJ_CUTOFF, ":", nrow(sig), "\n")
  cat("  Up   (z>0):", sum(sig$z > 0), "\n")
  cat("  Down (z<0):", sum(sig$z < 0), "\n")
  
  ortho <- human_to_mouse(sig$gene_human)
  cat("  One-to-one ortholog mapping:", nrow(ortho), "of", nrow(sig),
      "(", round(100 * nrow(ortho) / nrow(sig), 1), "%)\n")
  
  sig_m <- sig %>%
    dplyr::inner_join(ortho, by = c("gene_human" = "human_symbol")) %>%
    dplyr::distinct(mouse_symbol, .keep_all = TRUE)
  
  up   <- sig_m %>% dplyr::filter(z > 0) %>%
    dplyr::arrange(desc(z)) %>% dplyr::pull(mouse_symbol)
  down <- sig_m %>% dplyr::filter(z < 0) %>%
    dplyr::arrange(z)        %>% dplyr::pull(mouse_symbol)
  
  cat("  After ortholog mapping —  TET2_up:", length(up),
      ",  TET2_down:", length(down), "\n")
  
  list(up = up, down = down,
       n_up = length(up), n_down = length(down))
}

tet2_hsc_sets  <- build_tet2_sets(tet2_hsc_raw,  "HSC")
tet2_mono_sets <- build_tet2_sets(tet2_mono_raw, "Monocyte")

# Save mapped gene lists
write.csv(data.frame(gene = tet2_hsc_sets$up),
          file.path(OUTDIR_TET2, "TET2_up_HSC_mouse.csv"), row.names = FALSE)
write.csv(data.frame(gene = tet2_hsc_sets$down),
          file.path(OUTDIR_TET2, "TET2_down_HSC_mouse.csv"), row.names = FALSE)
write.csv(data.frame(gene = tet2_mono_sets$up),
          file.path(OUTDIR_TET2, "TET2_up_Monocyte_mouse.csv"), row.names = FALSE)
write.csv(data.frame(gene = tet2_mono_sets$down),
          file.path(OUTDIR_TET2, "TET2_down_Monocyte_mouse.csv"), row.names = FALSE)

# ── Load YOUR DEG files ─────────────────────────────────────────────────────
for (nm in names(deg_files_hsc))  assign(nm, load_deg(deg_files_hsc[[nm]]))
for (nm in names(deg_files_mono)) assign(nm, load_deg(deg_files_mono[[nm]]))

# ── GSEA wrapper ────────────────────────────────────────────────────────────
make_tet2_term2gene <- function(sets) {
  data.frame(
    term = c(rep("TET2_up",   sets$n_up),
             rep("TET2_down", sets$n_down)),
    gene = c(sets$up, sets$down)
  )
}

run_tet2_gsea <- function(deg_tbl, tet2_sets, tag, cell_label, outfile_tag) {
  if (tet2_sets$n_up < 10 && tet2_sets$n_down < 10) {
    warning("[", tag, "] TET2 sets too small. Skipping.")
    return(invisible(NULL))
  }
  ranked <- deg_tbl %>%
    dplyr::filter(X != "Lenti") %>%
    dplyr::arrange(desc(logFC)) %>%
    dplyr::distinct(X, .keep_all = TRUE) %>%
    { setNames(.$logFC, .$X) }
  cat("\n[", tag, "] Ranked list:", length(ranked), "genes\n")
  
  gsea_obj <- GSEA(
    geneList      = ranked,
    TERM2GENE     = make_tet2_term2gene(tet2_sets),
    minGSSize     = 10,
    maxGSSize     = max(5000, tet2_sets$n_up, tet2_sets$n_down) + 1,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    eps           = 0,
    verbose       = TRUE
  )
  
  stats <- as.data.frame(gsea_obj)[, c("ID", "NES", "pvalue", "p.adjust")]
  cat("--- [", tag, "] ---\n"); print(stats)
  
  plot_ids <- intersect(c("TET2_up", "TET2_down"), stats$ID)
  if (length(plot_ids) == 0) return(invisible(gsea_obj))
  
  p <- gseaplot2(gsea_obj, geneSetID = plot_ids,
                 title = paste0("Jakobsen TET2 ", cell_label, " signature in ", tag,
                                "\n(Sets: human DEGs FDR<", PADJ_CUTOFF, " -> mouse orthologs",
                                " | TET2_up n=", tet2_sets$n_up,
                                ", TET2_down n=", tet2_sets$n_down, ")"),
                 subplots = 1:3, base_size = 11)
  
  col_map <- c("TET2_up" = "#1A237E", "TET2_down" = "#B71C1C")[plot_ids]
  lab_map <- setNames(sapply(plot_ids, function(id) {
    n <- if (id == "TET2_up") tet2_sets$n_up else tet2_sets$n_down
    pretty <- if (id == "TET2_up") "TET2-up" else "TET2-down"
    paste0(pretty, " (n=", n,
           ") NES=", round(stats$NES[stats$ID == id], 2),
           " FDR=",  signif(stats$p.adjust[stats$ID == id], 2))
  }), plot_ids)
  
  p[[1]] <- p[[1]] +
    scale_color_manual(values = col_map, labels = lab_map) +
    theme(legend.position = "top", legend.text = element_text(size = 9))
  p[[2]] <- p[[2]] + scale_color_manual(values = col_map)
  
  outfile <- file.path(OUTDIR_TET2, paste0("GSEA_TET2_", outfile_tag, ".png"))
  png(outfile, width = 1400, height = 1200, res = 150, bg = "white")
  print(p); dev.off()
  message("Saved: ", basename(outfile))
  invisible(gsea_obj)
}

# ── Run GSEA ────────────────────────────────────────────────────────────────
gsea_results <- list(
  hsc_lentineg  = run_tet2_gsea(hsc_lentineg_hblb, tet2_hsc_sets,
                                "LT-HSC LentiNeg HB-vs-LB", "HSC",
                                "HSC_on_HSC_LentiNeg"),
  hsc_lentipos  = run_tet2_gsea(hsc_lentipos_hblb, tet2_hsc_sets,
                                "LT-HSC LentiPos HB-vs-LB", "HSC",
                                "HSC_on_HSC_LentiPos"),
  mono_lentineg = run_tet2_gsea(mono_lentineg_hblb, tet2_mono_sets,
                                "Monocyte LentiNeg HB-vs-LB", "Monocyte",
                                "Mono_on_Mono_LentiNeg"),
  mono_lentipos = run_tet2_gsea(mono_lentipos_hblb, tet2_mono_sets,
                                "Monocyte LentiPos HB-vs-LB", "Monocyte",
                                "Mono_on_Mono_LentiPos")
)

# ── Summary dotplot ─────────────────────────────────────────────────────────
extract_gsea <- function(gsea_obj, comparison, fallback_ids) {
  if (is.null(gsea_obj)) {
    return(data.frame(ID = fallback_ids, NES = NA_real_,
                      pvalue = NA_real_, p.adjust = NA_real_,
                      comparison = comparison))
  }
  d <- as.data.frame(gsea_obj)[, c("ID", "NES", "pvalue", "p.adjust")]
  missing_ids <- setdiff(fallback_ids, d$ID)
  if (length(missing_ids) > 0) {
    d <- rbind(d, data.frame(ID = missing_ids, NES = NA,
                             pvalue = NA, p.adjust = NA))
  }
  d$comparison <- comparison
  d
}

summary_df <- rbind(
  extract_gsea(gsea_results$hsc_lentineg,  "HSC LentiNeg HB-vs-LB",
               c("TET2_up", "TET2_down")),
  extract_gsea(gsea_results$hsc_lentipos,  "HSC LentiPos HB-vs-LB",
               c("TET2_up", "TET2_down")),
  extract_gsea(gsea_results$mono_lentineg, "Mono LentiNeg HB-vs-LB",
               c("TET2_up", "TET2_down")),
  extract_gsea(gsea_results$mono_lentipos, "Mono LentiPos HB-vs-LB",
               c("TET2_up", "TET2_down"))
)
summary_df$comparison <- factor(summary_df$comparison,
                                levels = c("HSC LentiNeg HB-vs-LB", "HSC LentiPos HB-vs-LB",
                                           "Mono LentiNeg HB-vs-LB", "Mono LentiPos HB-vs-LB"))
summary_df$ID <- factor(summary_df$ID, levels = c("TET2_up", "TET2_down"))

y_labels <- c(
  TET2_up   = paste0("TET2_up\n(HSC n=",   tet2_hsc_sets$n_up,
                     ", Mono n=",          tet2_mono_sets$n_up,   ")"),
  TET2_down = paste0("TET2_down\n(HSC n=", tet2_hsc_sets$n_down,
                     ", Mono n=",          tet2_mono_sets$n_down, ")")
)

cat("\nSummary dotplot data:\n"); print(summary_df)
write.csv(summary_df,
          file.path(OUTDIR_TET2, "TET2_GSEA_summary.csv"), row.names = FALSE)

p_summary <- ggplot(summary_df,
                    aes(x = comparison, y = ID, fill = NES,
                        size = -log10(p.adjust))) +
  geom_point(shape = 21, color = "grey30") +
  scale_fill_gradient2(low = "#1A237E", mid = "white", high = "#B71C1C",
                       midpoint = 0, limits = c(-2.5, 2.5),
                       oob = scales::squish, name = "NES") +
  scale_size_continuous(range = c(4, 14), name = "FDR",
                        breaks = -log10(c(0.5, 0.1, 0.05, 0.01, 0.001)),
                        labels = c("0.5", "0.1", "0.05", "0.01", "0.001")) +
  scale_y_discrete(labels = y_labels, expand = expansion(add = 0.8)) +
  theme_classic(base_size = 13) +
  theme(axis.text.x   = element_text(size = 10, angle = 30, hjust = 1,
                                     margin = margin(t = 4)),
        axis.text.y   = element_text(size = 11),
        plot.title    = element_text(size = 12, hjust = 0.5,
                                     margin = margin(b = 10)),
        plot.subtitle = element_text(size = 10, hjust = 0.5,
                                     margin = margin(b = 8)),
        plot.margin   = margin(t = 15, r = 10, b = 25, l = 15),
        plot.background = element_rect(fill = "white", color = NA)) +
  labs(x = NULL, y = "TET2 gene set",
       title    = "Jakobsen TET2-mut signature across HB-vs-LB contrasts",
       subtitle = paste0("TET2 sets from human DEGs (FDR<", PADJ_CUTOFF,
                         ") mapped to mouse orthologs."))

ggsave(file.path(OUTDIR_TET2, "Summary_dotplot_TET2_on_HBvsLB.png"),
       p_summary, width = 9, height = 5, dpi = 150, bg = "white")
message("Saved: Summary_dotplot_TET2_on_HBvsLB.png")

# ── Copy script ─────────────────────────────────────────────────────────────
if (file.exists(SCRIPT_PATH_TET2)) {
  file.copy(SCRIPT_PATH_TET2,
            file.path(OUTDIR_TET2, basename(SCRIPT_PATH_TET2)),
            overwrite = TRUE)
  message("Script copied to: ", OUTDIR_TET2)
} else {
  warning("SCRIPT_PATH_TET2 not found — update path: ", SCRIPT_PATH_TET2)
}

message("\nDone. Outputs in: ", OUTDIR_TET2)
