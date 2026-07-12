###############################################################################
## Cluster annotation helpers — runs on the saved <name>_seurat.rds objects.
## For each dataset & each clustering resolution:
##   1) FindAllMarkers            -> top positive markers per cluster (CSV)
##   2) FindConservedMarkers      -> markers conserved across genotype per cluster (CSV)
##   3) DotPlot of canonical HSPC/BM markers -> PNG (fastest way to eyeball identities)
## No re-clustering needed. Edit BASE / DATASETS / RES if desired.
###############################################################################
pkgs <- c("Seurat","dplyr","ggplot2","metap")   # metap=conserved markers (base Seurat otherwise)
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE))
  try(install.packages(p, repos="https://cloud.r-project.org"), silent=TRUE)
suppressPackageStartupMessages({library(Seurat); library(dplyr); library(ggplot2)})

BASE     <- "C:/Users/fc809/Downloads/mutant_proliefration_analysis"
OUTROOT  <- file.path(BASE, "seurat_harmony_output")
DATASETS <- c("Tet2","Asxl1","Jak2","Dnmt3a","Dnmt3a_R878H")   # subset if you like
RES      <- c(0.4, 1.2)          # must match what you clustered at
TOPN     <- 15                    # markers saved per cluster

## canonical mouse markers for the DotPlot (grouped by lineage)
canon <- c(
 "Hlf","Mecom","Procr","Fgd5","Mpl","Meis1","Hoxa9",        # HSC/MPP
 "Cd34","Kit","Flt3","Dntt","Satb1",                        # MPP / LMPP
 "Gata1","Klf1","Car1","Car2","Hba-a1","Gypa",              # Ery/MEP
 "Pf4","Itga2b","Vwf","Gp1bb",                              # MkP
 "Mpo","Elane","Prtn3","Ctsg","Ms4a3",                      # GMP
 "Camp","Ngp","Ltf","S100a8","S100a9",                      # Neutrophil
 "Csf1r","Ly6c2","Ccr2","F13a1",                            # Monocyte
 "Siglech","Bst2","Itgax","Cd74",                           # pDC / cDC
 "Prss34","Mcpt8","Cpa3","Ms4a2",                           # Baso/Mast
 "Prg2","Epx","Il5ra",                                      # Eosinophil
 "Cd79a","Vpreb1","Ebf1","Ms4a1","Jchain","Mzb1",           # B / Plasma
 "Cd3e","Il7r","Ncr1","Nkg7")                               # T / NK

annotate_one <- function(name) {
  od  <- file.path(OUTROOT, name)
  rds <- file.path(od, paste0(name, "_seurat.rds"))
  if (!file.exists(rds)) { message("skip ", name, " (no rds)"); return(invisible()) }
  message("\n==== ", name, " ====")
  obj <- readRDS(rds)
  DefaultAssay(obj) <- "RNA"
  if (inherits(obj[["RNA"]], "Assay5")) obj <- JoinLayers(obj)   # v5: merge split layers
  obj <- NormalizeData(obj, verbose=FALSE)                       # ensure data layer exists
  gg <- canon[canon %in% rownames(obj)]

  for (r in RES) {
    col <- paste0("RNA_snn_res.", r)
    if (!col %in% colnames(obj@meta.data)) { message("  no ", col); next }
    Idents(obj) <- factor(obj@meta.data[[col]])

    ## 1) all markers, robust: low logfc.threshold + per-cluster fallback so the
    ##    central progenitor clusters (dropped at the default 0.25) are all included.
    am <- FindAllMarkers(obj, only.pos=TRUE, min.pct=0.1, logfc.threshold=0.1,
                         return.thresh=1, max.cells.per.ident=2000, verbose=FALSE)
    for (cl in setdiff(levels(Idents(obj)), unique(as.character(am$cluster)))) {
      fm <- tryCatch(FindMarkers(obj, ident.1=cl, only.pos=TRUE, min.pct=0.05,
              logfc.threshold=0, max.cells.per.ident=2000, verbose=FALSE), error=function(e) NULL)
      if (!is.null(fm) && nrow(fm)) { fm$gene<-rownames(fm); fm$cluster<-cl; am<-bind_rows(am, fm) }
    }
    if (nrow(am)) {
      am <- am %>% arrange(cluster, desc(avg_log2FC))
      top <- am %>% group_by(cluster) %>% slice_max(avg_log2FC, n=TOPN, with_ties=FALSE) %>% ungroup()
      write.csv(am,  file.path(od, sprintf("%s_markers_res%s_all.csv",  name, r)), row.names=FALSE)
      write.csv(top, file.path(od, sprintf("%s_markers_res%s_top%d.csv", name, r, TOPN)), row.names=FALSE)
    }

    ## 2) conserved markers across genotype, per cluster (robust to genotype-specific expression)
    cons <- list()
    for (cl in levels(Idents(obj))) {
      cm <- tryCatch(
        FindConservedMarkers(obj, ident.1=cl, grouping.var="genotype",
                             only.pos=TRUE, min.pct=0.25, verbose=FALSE),
        error=function(e) NULL)
      if (!is.null(cm) && nrow(cm)) {
        cm$gene <- rownames(cm); cm$cluster <- cl
        cons[[cl]] <- head(cm, TOPN)
      }
    }
    if (length(cons))
      write.csv(bind_rows(cons), file.path(od, sprintf("%s_conserved_markers_res%s.csv", name, r)),
                row.names=FALSE)

    ## 3) canonical-marker DotPlot
    p <- DotPlot(obj, features=gg, cluster.idents=TRUE) +
         RotatedAxis() + ggtitle(paste0(name, "  (res ", r, ")")) +
         theme(axis.text.x=element_text(size=7))
    ggsave(file.path(od, sprintf("%s_canonical_dotplot_res%s.png", name, r)),
           p, width=16, height=6, dpi=150)
    message("  res ", r, ": wrote markers + conserved + dotplot")
  }
}

for (nm in DATASETS) annotate_one(nm)
message("\nDone. Look in each seurat_harmony_output/<dataset>/ for:")
message("  *_markers_res<r>_top15.csv     (top positive markers per cluster)")
message("  *_conserved_markers_res<r>.csv (conserved across WT & mutant)")
message("  *_canonical_dotplot_res<r>.png (eyeball identities)")
###############################################################################
## Notes
## - presto makes FindAllMarkers ~10-50x faster; if it won't install, Seurat falls
##   back to the (slower) Wilcoxon default automatically.
## - FindConservedMarkers needs 'metap'. If a cluster has cells from only one
##   genotype it is skipped (tryCatch) rather than erroring.
## - grouping.var="genotype" assumes that column exists (it does in these objects).
## - To annotate: open the dotplot, match each cluster's high-marker block to a
##   lineage; cross-check with the top-markers CSV. Conserved markers are the most
##   reliable identity call because they must hold in BOTH genotypes.
###############################################################################
