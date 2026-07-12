## ============================================================
## GSE207063 — BM Ly6Chi Classical Monocyte Aging
## Barman et al. 2022, Aging Cell
## MALE MICE ONLY: Young (2-6mo) vs Old (24-30mo)
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
##   10. Plots + console summary
##
## OUTPUTS:
##   GSE207063_seurat_male.rds
##   GSE207063_DEGs_male_OldVsYoung.csv
##   GSE207063_GSEA_GO_BP.csv
##   GSE207063_GSEA_KEGG.csv
##   GSE207063_GSEA_Hallmarks.csv
##   GSE207063_umap_male.png
##   GSE207063_volcano_male.png
##   GSE207063_GSEA_bar_GOBP.png
##   GSE207063_GSEA_bar_KEGG.png
##   GSE207063_GSEA_bar_Hallmarks.png
##   GSE207063_GSEA_dotplot_GOBP.png
##   GSE207063_GSEA_dotplot_KEGG.png
##   GSE207063_GSEA_dotplot_Hallmarks.png
##   GSE207063_GSEA_curves_GOBP.png
##   GSE207063_GSEA_curves_Hallmarks.png
##   GSE207063_GSEA_emap_GOBP.png
## ============================================================

## ── PACKAGES ──────────────────────────────────────────────────────────────────
## Install if needed:
## install.packages(c("BiocManager","ggrepel","stringr"))
## BiocManager::install(c("GEOquery","Seurat","edgeR",
##                        "clusterProfiler","org.Mm.eg.db","enrichplot","msigdbr"))

library(GEOquery)
library(Seurat)
library(edgeR)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(stringr)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(msigdbr)

set.seed(42)

OUTDIR <- "C:/Users/fc809/Downloads/GSE207063"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
setwd(OUTDIR)

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

# Just copy the gz files with the names Read10X expects
file.copy(file.path(GEO_DIR, "GSE207063_Male_Young_Old_barcodes.tsv.gz"),
          file.path(male_dir, "barcodes.tsv.gz"))
file.copy(file.path(GEO_DIR, "GSE207063_Male_Young_Old_features.tsv.gz"),
          file.path(male_dir, "features.tsv.gz"))
file.copy(file.path(GEO_DIR, "GSE207063_Male_Young_Old_matrix.mtx.gz"),
          file.path(male_dir, "matrix.mtx.gz"))

counts <- Read10X(male_dir)
so <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)
message("Loaded ", ncol(so), " cells")

barcode_suffixes <- sub(".*-", "-", colnames(so))
print(table(barcode_suffixes))
# Load with Read10X
counts <- Read10X(male_dir)
so <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)
message("Loaded ", ncol(so), " cells")

# Add HTO assay
so[["HTO"]] <- CreateAssayObject(counts = counts[["Antibody Capture"]][, colnames(so)])

# Normalize HTO
so <- NormalizeData(so, assay = "HTO", normalization.method = "CLR")

# Demultiplex
so <- HTODemux(so, assay = "HTO", positive.quantile = 0.99)

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
so <- subset(so, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & pct.mt < 20)
message("Cells after QC: ", ncol(so))

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
  rename(log2FC = logFC, pval = PValue) %>%
  mutate(direction = case_when(
    FDR < 0.05 & log2FC >  0.5 ~ "Up_Old",
    FDR < 0.05 & log2FC < -0.5 ~ "Up_Young",
    TRUE ~ "NS"
  ))

message("DE summary:")
print(table(res_df$direction))

write.csv(res_df, file.path(OUTDIR, "GSE207063_DEGs_male_OldVsYoung.csv"), row.names = FALSE)

###PATHWAY_ANALYSIS####

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(msigdbr)

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

h_sets <- h_raw %>%
  dplyr::select(gs_name, ncbi_gene) %>%
  mutate(ncbi_gene = as.character(ncbi_gene))

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
gsea_go_df %>% select(Description, NES, p.adjust, setSize) %>% head(10) %>% print()

message("\n--- KEGG all significant ---")
gsea_kegg_df %>% arrange(NES) %>% select(Description, NES, p.adjust, setSize) %>% print()

message("\n--- Hallmarks significant ---")
gsea_h_df %>% filter(p.adjust < 0.05) %>% arrange(NES) %>%
  select(Description, NES, p.adjust, setSize) %>% print()

# ── Fix select conflict ───────────────────────────────────────────────────────
library(dplyr)

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

library(ggplot2)
library(dplyr)
library(ggrepel)
library(clusterProfiler)

library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(enrichplot)
library(patchwork)

# ── Build gene sets from Old vs Young (GSE207063) ─────────────────────────────
old_up   <- res_df %>% arrange(desc(log2FC)) %>% head(200) %>% pull(gene)
old_down <- res_df %>% arrange(log2FC)       %>% head(200) %>% pull(gene)

aging_sets <- data.frame(
  term = c(rep("Old_up",   length(old_up)),
           rep("Old_down", length(old_down))),
  gene = c(old_up, old_down)
)

# ── Ranked list from your HB vs LB data ──────────────────────────────────────
ranked_hb <- your_mono %>%
  arrange(desc(logFC)) %>%
  distinct(X, .keep_all = TRUE) %>%
  { setNames(.$logFC, .$X) }

cat("Ranked HB vs LB list:", length(ranked_hb), "genes\n")

# ── GSEA ──────────────────────────────────────────────────────────────────────
gsea_res <- GSEA(
  geneList      = ranked_hb,
  TERM2GENE     = aging_sets,
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = TRUE
)

# Print results
print(as.data.frame(gsea_res)[, c("ID", "NES", "pvalue", "p.adjust", "setSize")])

# ── Enrichment curves ─────────────────────────────────────────────────────────
png("C:/Users/fc809/Downloads/GSEA_OldvsYoung_on_HBvsLB.png",
    width = 1400, height = 1200, res = 150)
print(
  gseaplot2(gsea_res,
            geneSetID = c("Old_up", "Old_down"),
            title     = "Aging signature enrichment in HB vs LB monocytes",
            subplots  = 1:3,
            base_size = 11)
)
dev.off()

message("Saved: GSEA_OldvsYoung_on_HBvsLB.png")

library(ggplot2)

# Get NES and FDR for annotation
res_stats <- as.data.frame(gsea_res)[, c("ID", "NES", "pvalue", "p.adjust")]
print(res_stats)

# Build plot and recolor
p <- gseaplot2(gsea_res,
               geneSetID = c("Old_up", "Old_down"),
               title     = "Aging signature enrichment in HB vs LB monocytes",
               subplots  = 1:3,
               base_size = 11)

# Recolor — gseaplot2 returns a patchwork; modify the top panel layers
p[[1]] <- p[[1]] +
  scale_color_manual(
    values = c("Old_up" = "#1A237E", "Old_down" = "#B71C1C"),
    labels = c(
      "Old_up"   = paste0("Old-up   NES=",   round(res_stats$NES[res_stats$ID=="Old_up"],   2),
                          "  FDR=", signif(res_stats$p.adjust[res_stats$ID=="Old_up"],   2)),
      "Old_down" = paste0("Old-down NES=",   round(res_stats$NES[res_stats$ID=="Old_down"], 2),
                          "  FDR=", signif(res_stats$p.adjust[res_stats$ID=="Old_down"], 2))
    )
  ) +
  theme(legend.position = "top",
        legend.text     = element_text(size = 9))

# Recolor tick marks panel
p[[2]] <- p[[2]] +
  scale_color_manual(values = c("Old_up" = "#1A237E", "Old_down" = "#B71C1C"))

png("C:/Users/fc809/Downloads/GSEA_OldvsYoung_on_HBvsLB.png",
    width = 1400, height = 1200, res = 150)
print(p)
dev.off()

message("Saved.")

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
  maxGSSize     = 500,
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
  maxGSSize     = 500,
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
                 title     = title,
                 subplots  = 1:3,
                 base_size = 11)
  
  p[[1]] <- p[[1]] +
    scale_color_manual(
      values = c("Old_up" = "#1A237E", "Old_down" = "#B71C1C"),
      labels = c(
        "Old_up"   = paste0("Old-up   NES=", round(stats$NES[stats$ID == "Old_up"],   2),
                            "  FDR=", signif(stats$p.adjust[stats$ID == "Old_up"],   2)),
        "Old_down" = paste0("Old-down NES=", round(stats$NES[stats$ID == "Old_down"], 2),
                            "  FDR=", signif(stats$p.adjust[stats$ID == "Old_down"], 2))
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
                "C:/Users/fc809/Downloads/GSEA_OldvsYoung_on_LentiHB.png")

plot_aging_gsea(gsea_lenti_lb,
                "Aging signature enrichment in Lenti+ vs Lenti- (LB monocytes)",
                "C:/Users/fc809/Downloads/GSEA_OldvsYoung_on_LentiLB.png")
###switch####

# ── Build gene sets from HB vs LB ────────────────────────────────────────────
hb_up   <- your_mono %>% arrange(desc(logFC)) %>% head(200) %>% pull(X)
hb_down <- your_mono %>% arrange(logFC)       %>% head(200) %>% pull(X)

hb_sets <- data.frame(
  term = c(rep("HB_up",   length(hb_up)),
           rep("HB_down", length(hb_down))),
  gene = c(hb_up, hb_down)
)

# ── Ranked list from Old vs Young (GSE207063) ─────────────────────────────────
ranked_aging <- res_df %>%
  arrange(desc(log2FC)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  { setNames(.$log2FC, .$gene) }

# ── GSEA ──────────────────────────────────────────────────────────────────────
gsea_hb_on_aging <- GSEA(
  geneList      = ranked_aging,
  TERM2GENE     = hb_sets,
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = TRUE
)

# ── Plot ──────────────────────────────────────────────────────────────────────
stats <- as.data.frame(gsea_hb_on_aging)[, c("ID", "NES", "p.adjust")]
print(stats)

p <- gseaplot2(gsea_hb_on_aging,
               geneSetID = c("HB_up", "HB_down"),
               title     = "HB vs LB signature enrichment in Old vs Young monocytes",
               subplots  = 1:3,
               base_size = 11)

p[[1]] <- p[[1]] +
  scale_color_manual(
    values = c("HB_up" = "#1A237E", "HB_down" = "#B71C1C"),
    labels = c(
      "HB_up"   = paste0("HB-up   NES=", round(stats$NES[stats$ID == "HB_up"],   2),
                         "  FDR=", signif(stats$p.adjust[stats$ID == "HB_up"],   2)),
      "HB_down" = paste0("HB-down NES=", round(stats$NES[stats$ID == "HB_down"], 2),
                         "  FDR=", signif(stats$p.adjust[stats$ID == "HB_down"], 2))
    )
  ) +
  theme(legend.position = "top", legend.text = element_text(size = 9))

p[[2]] <- p[[2]] +
  scale_color_manual(values = c("HB_up" = "#1A237E", "HB_down" = "#B71C1C"))

png("C:/Users/fc809/Downloads/GSEA_HBvsLB_on_OldvsYoung.png",
    width = 1400, height = 1200, res = 150)
print(p)
dev.off()

message("Saved: GSEA_HBvsLB_on_OldvsYoung.png")

