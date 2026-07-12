## ============================================================
##  Inflammation gene expression heatmap — Classical Monocytes
##  Same structure as GSEA master heatmap but rows = genes,
##  values = logFC from pseudobulk edgeR
##  4 comparisons × 8 UMI-filtering cases
## ============================================================

library(data.table)
library(tidyverse)
library(edgeR)
library(ComplexHeatmap)
library(circlize)

setwd("C:/Users/fc809/Downloads/")

## ── Parameters ────────────────────────────────────────────────────────────────
MATRIX_FILE  <- "classical_monocytes_gene_nUMI_per_cell_for_Foivos.txt"
THRESH_FILE  <- "Lenti_detection_thresholds_SCTv2.csv"
LENTI_FREQ   <- 4e-5
FIXED_THRESH  <- 13700
FIXED_THRESH2 <- 9227
MIN_CELLS     <- 100
FDR_CUT       <- 0.05    # shade genes non-significant at this FDR

## ── Inflammation gene panel (~100 genes) ─────────────────────────────────────
# Named list: category → genes (order preserved for heatmap row grouping)
gene_categories <- list(
  "Cytokines"             = c("Il1b","Il1a","Il6","Il10","Il12b","Il12a","Il18",
                              "Il23a","Il33","Il17a","Il4","Il13","Il22","Il34",
                              "Il15","Il27","Tnf","Tgfb1","Ifng","Csf1","Csf2","Csf3"),
  "Chemokines & receptors"= c("Ccl2","Ccl3","Ccl4","Ccl5","Ccl7","Ccl8","Ccl12",
                              "Cxcl1","Cxcl2","Cxcl5","Cxcl9","Cxcl10","Cxcl16",
                              "Ccr1","Ccr2","Ccr5","Cx3cr1","Cxcr4"),
  "PRRs & innate"         = c("Tlr2","Tlr4","Tlr7","Tlr9","Myd88","Ticam1",
                              "Nod1","Nod2","Aim2","Sting1"),
  "Transcription factors" = c("Nfkb1","Nfkb2","Rela","Relb","Ikbke",
                              "Irf1","Irf3","Irf5","Irf7","Irf8",
                              "Stat1","Stat3","Stat6","Jun","Fos"),
  "Inflammasome"          = c("Nlrp3","Nlrc4","Casp1","Casp4","Pycard","Gsdmd"),
  "Prostaglandin/ROS"     = c("Ptgs1","Ptgs2","Ptges","Alox5","Alox15","Pla2g4a",
                              "Cybb","Ncf1","Ncf4","Sod1","Sod2","Cat","Gpx1","Hmox1"),
  "DAMPs / alarmins"      = c("Hmgb1","S100a8","S100a9","S100a12","Lcn2","Thbs1"),
  "Macrophage polaris."   = c("Nos2","Arg1","Mrc1","Cd163","Chil3","Retnla",
                              "Cd80","Cd86","Cd68","Adgre1"),
  "Anti-inflammatory"     = c("Il1rn","Socs1","Socs3","Tnfaip3","Nfkbia"),
  "Coagulation/vascular"  = c("F3","Selp","Icam1","Vcam1","Thbd"),
  "Complement"            = c("C1qa","C1qb","C3","C5ar1")
)

# Flat ordered gene list (preserves category order)
inflammation_genes <- unlist(gene_categories, use.names=FALSE)
gene_category_map  <- rep(names(gene_categories), sapply(gene_categories, length))
names(gene_category_map) <- inflammation_genes

## ── Load thresholds & matrix ──────────────────────────────────────────────────
thresholds <- read.csv(THRESH_FILE)
thresh_vec <- setNames(thresholds$UMI_Threshold, thresholds$Sample)

message("Reading matrix...")
mat_raw    <- as.data.frame(fread(MATRIX_FILE))
gene_names <- mat_raw[[1]]
mat        <- as.matrix(mat_raw[, -1])
rownames(mat) <- gene_names
barcodes   <- colnames(mat)
message("  ", nrow(mat), " genes x ", ncol(mat), " cells")

# Check which inflammation genes are present
genes_present <- inflammation_genes[inflammation_genes %in% rownames(mat)]
genes_missing <- setdiff(inflammation_genes, rownames(mat))
if (length(genes_missing)) message("Genes not in matrix: ", paste(genes_missing, collapse=", "))
message("Genes found: ", length(genes_present), "/", length(inflammation_genes))

## ── Per-cell metadata ─────────────────────────────────────────────────────────
lenti_row <- which(grepl("^[Ll]enti$", gene_names))[1]
numi_vec   <- colSums(mat)
lenti_freq <- mat[lenti_row, ] / numi_vec
sample_vec <- sub("_.*", "", barcodes)
cond_vec   <- ifelse(grepl("^HB", sample_vec), "HB", "LB")

meta_all <- data.frame(
  barcode      = barcodes,
  sample       = sample_vec,
  condition    = cond_vec,
  nUMI         = numi_vec,
  lenti_freq   = lenti_freq,
  Lenti_status = ifelse(lenti_freq >= LENTI_FREQ, "LentiPos", "LentiNeg"),
  stringsAsFactors = FALSE
)

## ── Comparisons & cases (same as GSEA script) ─────────────────────────────────
comparisons <- list(
  list(id="LentiPosNeg_LB",  label="LentiPos vs LentiNeg (LB)",
       cell_filter=function(m) m[m$condition=="LB",],
       group_var="Lenti_status", ref_level="LentiNeg", test_level="LentiPos",
       paired=TRUE,  do_match=TRUE),
  list(id="LentiPosNeg_HB",  label="LentiPos vs LentiNeg (HB)",
       cell_filter=function(m) m[m$condition=="HB",],
       group_var="Lenti_status", ref_level="LentiNeg", test_level="LentiPos",
       paired=TRUE,  do_match=TRUE),
  list(id="HBvsLB_LentiNeg", label="HB vs LB (LentiNeg)",
       cell_filter=function(m) m[m$Lenti_status=="LentiNeg",],
       group_var="condition", ref_level="LB", test_level="HB",
       paired=FALSE, do_match=TRUE),
  list(id="HBvsLB_LentiPos", label="HB vs LB (LentiPos)",
       cell_filter=function(m) m[m$Lenti_status=="LentiPos",],
       group_var="condition", ref_level="LB", test_level="HB",
       paired=FALSE, do_match=TRUE)
)

cases <- list(
  "1_all_cells"      = list(label="Case 1: all cells",                           umi_filter="none",   umi_match=FALSE),
  "2_all_matched"    = list(label="Case 2: all cells + UMI matched",              umi_filter="none",   umi_match=TRUE),
  "3_fixed_thresh"   = list(label=paste0("Case 3: nUMI \u2265 ",FIXED_THRESH),    umi_filter="fixed",  umi_match=FALSE),
  "4_fixed_matched"  = list(label=paste0("Case 4: nUMI \u2265 ",FIXED_THRESH," + UMI matched"), umi_filter="fixed",  umi_match=TRUE),
  "5_pooled_thresh"  = list(label=paste0("Case 5: nUMI \u2265 ",FIXED_THRESH2),   umi_filter="fixed2", umi_match=FALSE),
  "6_pooled_matched" = list(label=paste0("Case 6: nUMI \u2265 ",FIXED_THRESH2," + UMI matched"), umi_filter="fixed2", umi_match=TRUE),
  "7_sample_thresh"  = list(label="Case 7: sample-specific SCTv2",                umi_filter="sample", umi_match=FALSE),
  "8_sample_matched" = list(label="Case 8: sample-specific SCTv2 + UMI matched",  umi_filter="sample", umi_match=TRUE)
)

## ── Helpers ───────────────────────────────────────────────────────────────────
filter_min_cells <- function(meta_sub, group_var, paired) {
  counts <- meta_sub %>% group_by(sample, .data[[group_var]]) %>%
    summarise(n=n(), .groups="drop")
  if (paired) {
    wide <- counts %>% pivot_wider(names_from=all_of(group_var), values_from=n, values_fill=0)
    keep <- wide
    for (lv in unique(meta_sub[[group_var]])) keep <- keep[keep[[lv]] >= MIN_CELLS,]
    keep_samples <- keep$sample
  } else {
    keep_samples <- counts %>% dplyr::filter(n >= MIN_CELLS) %>% pull(sample) %>% unique()
  }
  meta_sub[meta_sub$sample %in% keep_samples,]
}

match_umi <- function(meta_sub, group_var, ref_level, test_level) {
  result <- list()
  for (samp in unique(meta_sub$sample)) {
    s    <- meta_sub[meta_sub$sample == samp,]
    ref  <- s[s[[group_var]] == ref_level,]
    test <- s[s[[group_var]] == test_level,]
    if (nrow(ref)==0 || nrow(test)==0) { result[[samp]] <- s; next }
    breaks <- unique(quantile(ref$nUMI, probs=seq(0,1,0.1)))
    if (length(breaks) < 2) {
      keep_test <- test[sample(nrow(test), min(nrow(test), nrow(ref))),]
    } else {
      ref$bin  <- cut(ref$nUMI,  breaks=breaks, include.lowest=TRUE)
      test$bin <- cut(test$nUMI, breaks=breaks, include.lowest=TRUE)
      bin_counts <- table(ref$bin)
      keep_test <- do.call(rbind, lapply(names(bin_counts), function(b) {
        avail <- test[!is.na(test$bin) & test$bin==b,]
        n_want <- bin_counts[b]
        if (nrow(avail)==0 || n_want==0) return(NULL)
        avail[sample(nrow(avail), min(nrow(avail), n_want)),]
      }))
      ref$bin <- NULL; keep_test$bin <- NULL
    }
    result[[samp]] <- rbind(ref, keep_test)
  }
  do.call(rbind, result)
}

run_edger_genes <- function(meta_sub, count_mat, group_var, ref_level,
                            test_level, paired, genes) {
  levels_grp <- c(ref_level, test_level)
  pb_list <- list(); pb_meta <- list()
  for (samp in unique(meta_sub$sample)) {
    for (grp in levels_grp) {
      cells <- meta_sub$barcode[meta_sub$sample==samp & meta_sub[[group_var]]==grp]
      if (length(cells)==0) next
      cn <- paste0(samp,"_",grp)
      pb_list[[cn]] <- rowSums(count_mat[, cells, drop=FALSE])
      pb_meta[[cn]] <- data.frame(sample=samp, group=grp, stringsAsFactors=FALSE)
    }
  }
  if (length(pb_list)==0) return(NULL)
  pb_counts <- do.call(cbind, pb_list)
  pb_df     <- do.call(rbind, pb_meta)
  rownames(pb_df) <- colnames(pb_counts)
  pb_df$group  <- factor(pb_df$group,  levels=levels_grp)
  pb_df$sample <- factor(pb_df$sample)
  n_ref  <- sum(pb_df$group == ref_level)
  n_test <- sum(pb_df$group == test_level)
  if (n_ref < 2 || n_test < 2) return(NULL)
  design <- if (paired) model.matrix(~sample+group, data=pb_df) else
    model.matrix(~group, data=pb_df)
  coef_name <- paste0("group", test_level)
  if (!coef_name %in% colnames(design)) return(NULL)
  dge <- DGEList(counts=pb_counts)
  dge <- dge[filterByExpr(dge, design),]
  dge <- calcNormFactors(dge, method="TMM")
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  res <- glmQLFTest(fit, coef=coef_name)
  tt  <- topTags(res, n=Inf)$table
  tt$gene <- rownames(tt)
  # Return only inflammation genes present after filterByExpr
  tt[tt$gene %in% genes,]
}

## ── Main loop ─────────────────────────────────────────────────────────────────
all_deg_results <- list()

for (cmp in comparisons) {
  message("\n========================================")
  message("Comparison: ", cmp$label)
  
  for (case_id in names(cases)) {
    cfg  <- cases[[case_id]]
    meta <- cmp$cell_filter(meta_all)
    
    if (cfg$umi_filter == "fixed") {
      meta <- meta[meta$nUMI >= FIXED_THRESH,]
    } else if (cfg$umi_filter == "fixed2") {
      meta <- meta[meta$nUMI >= FIXED_THRESH2,]
    } else if (cfg$umi_filter == "sample") {
      meta <- meta[meta$nUMI >= thresh_vec[meta$sample],]
    }
    
    meta <- filter_min_cells(meta, cmp$group_var, cmp$paired)
    if (nrow(meta) == 0) next
    
    if (cfg$umi_match && cmp$do_match) {
      meta <- match_umi(meta, cmp$group_var, cmp$ref_level, cmp$test_level)
      meta <- filter_min_cells(meta, cmp$group_var, cmp$paired)
      if (nrow(meta) == 0) next
    }
    
    deg <- run_edger_genes(meta, mat, cmp$group_var, cmp$ref_level,
                           cmp$test_level, cmp$paired, genes_present)
    if (!is.null(deg) && nrow(deg) > 0) {
      deg$case       <- cfg$label
      deg$comparison <- cmp$label
      deg$col_id     <- paste0(cmp$label, "  |  ", cfg$label)
      all_deg_results[[paste(cmp$id, case_id)]] <- deg
      message("  ", cfg$label, ": ", nrow(deg), " genes tested")
    }
  }
}

all_deg <- bind_rows(all_deg_results)

# Save raw results
write.csv(all_deg, "inflammation_genes_DEG_all_results.csv", row.names=FALSE)
message("\nDEG results saved.")

## ── Build matrix ──────────────────────────────────────────────────────────────
# Column order: comparison blocks × case order
comparison_labels <- sapply(comparisons, `[[`, "label")
case_labels        <- sapply(cases,       `[[`, "label")

col_order <- c()
col_split <- c()
for (cmp in comparisons) {
  for (cs in case_labels) {
    cid <- paste0(cmp$label, "  |  ", cs)
    if (cid %in% unique(all_deg$col_id)) {
      col_order <- c(col_order, cid)
      col_split <- c(col_split, cmp$label)
    }
  }
}

# logFC matrix (all genes × all columns)
lfc_wide <- all_deg %>%
  dplyr::select(gene, col_id, logFC) %>%
  pivot_wider(names_from=col_id, values_from=logFC) %>%
  column_to_rownames("gene")

# FDR matrix for masking
fdr_wide <- all_deg %>%
  dplyr::select(gene, col_id, FDR) %>%
  pivot_wider(names_from=col_id, values_from=FDR) %>%
  column_to_rownames("gene")

# Keep only genes present in matrix
lfc_wide <- lfc_wide[rownames(lfc_wide) %in% genes_present,, drop=FALSE]
fdr_wide <- fdr_wide[rownames(lfc_wide),, drop=FALSE]

# Reorder columns
col_order <- col_order[col_order %in% colnames(lfc_wide)]
lfc_wide  <- lfc_wide[, col_order, drop=FALSE]
fdr_wide  <- fdr_wide[, col_order, drop=FALSE]

# Mask non-significant logFC → NA (shown as grey)
lfc_masked <- lfc_wide
for (r in rownames(lfc_masked))
  for (c in colnames(lfc_masked))
    if (is.na(fdr_wide[r,c]) || fdr_wide[r,c] >= FDR_CUT)
      lfc_masked[r,c] <- NA

# Remove genes with all NA
keep <- rowSums(!is.na(lfc_masked)) > 0
lfc_masked <- lfc_masked[keep,, drop=FALSE]
message("Genes with at least one significant result: ", nrow(lfc_masked))

# Clustering matrix (impute NA → 0)
lfc_clust <- as.matrix(lfc_masked); lfc_clust[is.na(lfc_clust)] <- 0

## ── Draw heatmap ──────────────────────────────────────────────────────────────
max_lfc <- max(abs(lfc_masked), na.rm=TRUE)
if (!is.finite(max_lfc) || max_lfc == 0) max_lfc <- 2

col_fun <- colorRamp2(
  c(-max_lfc, -max_lfc/2, 0, max_lfc/2, max_lfc),
  c("#2980B9", "#7fbbd4", "white", "#e87070", "#C0392B")
)

col_labels <- sub("^.*  \\|  Case \\d+: ", "", colnames(lfc_masked))

# Build row split factor in the order genes appear in lfc_masked
row_cats <- gene_category_map[rownames(lfc_masked)]
row_cats[is.na(row_cats)] <- "Other"
row_split <- factor(row_cats, levels=unique(gene_category_map))

ht <- Heatmap(
  matrix         = as.matrix(lfc_masked),
  col            = col_fun,
  na_col         = "grey92",
  name           = "logFC",
  
  # Row grouping by category — no dendrogram (order fixed by category)
  cluster_rows   = FALSE,
  row_split      = row_split,
  row_gap        = unit(2, "mm"),
  row_title_gp   = gpar(fontsize=11, fontface="bold"),
  row_title_rot  = 0,
  row_names_side = "left",
  row_names_gp   = gpar(fontsize=10, fontface="italic"),
  row_names_max_width = unit(5, "cm"),
  
  cluster_columns  = FALSE,
  column_split     = factor(col_split, levels=comparison_labels),
  column_gap       = unit(6, "mm"),
  column_labels    = col_labels,
  column_names_gp  = gpar(fontsize=10),
  column_names_rot = 45,
  column_names_max_height = unit(7, "cm"),
  column_title_gp  = gpar(fontsize=12, fontface="bold"),
  
  rect_gp = gpar(col="white", lwd=1),
  
  heatmap_legend_param = list(
    title_gp   = gpar(fontsize=12, fontface="bold"),
    labels_gp  = gpar(fontsize=11),
    legend_height = unit(5, "cm"),
    grid_width    = unit(0.8, "cm")
  )
)

## ── Save ──────────────────────────────────────────────────────────────────────
out_stem <- "Inflammation_genes_heatmap_Monocytes"

pdf(paste0(out_stem, ".pdf"),
    width  = max(18, ncol(lfc_masked) * 0.55 + 8),
    height = max(10, nrow(lfc_masked) * 0.35 + 4))

draw(ht,
     column_title    = paste0("Classical Monocytes — Inflammation genes (logFC, FDR < ", FDR_CUT, ")"),
     column_title_gp = gpar(fontsize=14, fontface="bold"),
     padding         = unit(c(2, 2, 2, 5), "cm"),
     merge_legend    = TRUE)

dev.off()
message("Saved: ", out_stem, ".pdf")

