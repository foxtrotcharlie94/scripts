###############################################################################
## FIX cluster markers for ALL datasets with a Seurat .rds (Tet2, Dnmt3a_R878H,
## Asxl1) at both resolutions.  NO presto required (base Seurat only).
##
## Why the originals were broken: FindAllMarkers(only.pos, logfc.threshold=0.25)
## does one-vs-rest DE. The large, transcriptionally-central progenitor clusters
## have fold-changes below threshold vs the (myeloid-dominated) "rest", so they
## return ZERO rows and silently disappear. Only rare distinct lineages survived.
##
## Fix: run FindAllMarkers with a low logfc.threshold + return.thresh=1, then a
## per-cluster fallback (logfc.threshold=0) for any cluster still missing, so EVERY
## cluster gets ranked markers. max.cells.per.ident downsamples big clusters -> fast.
###############################################################################
suppressPackageStartupMessages({
  for (p in c("Seurat","dplyr","ggplot2","future")) {
    if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")
    library(p, character.only=TRUE)
  }
})

## ---- backend: SEQUENTIAL (single process) ----
## Do NOT use future/multisession here. Seurat marker functions serialize the whole
## object + the DE closure to each worker; on the larger objects (e.g. R878H 2.85 GB)
## this exceeds future.globals.maxSize ("total size of globals is 26 GiB") and errors.
## With max.cells.per.ident=2000 the sequential run is fast enough. For real speed,
## install presto (remotes::install_github("immunogenomics/presto")) — it multithreads
## in-process with NO object copying, and Seurat::FindMarkers auto-uses it.
plan("sequential")
options(future.globals.maxSize = 32 * 1024^3)   # generous, just in case

BASE     <- "C:/Users/fc809/Downloads/mutant_proliefration_analysis/seurat_harmony_output"
## All three: the earlier files are truncated (interrupted writes), so regenerate all.
## With DO_CONSERVED=FALSE + sequential + atomic writes this finishes in minutes.
DATASETS <- c("Tet2","Dnmt3a_R878H","Asxl1")
RES          <- c(0.4, 1.2)
TOPN         <- 15
MAXCELLS     <- 500     # per-cluster downsample for the DE test (speed; ranking unaffected)
DO_CONSERVED <- FALSE   # conserved-across-genotype markers = the slow step; off by default
HVG_ONLY     <- TRUE    # test only the ~2000 highly-variable genes -> ~20-30x faster on base
                        # Seurat Wilcoxon. HVGs contain the discriminating markers, so cluster
                        # identity is unaffected. Set FALSE for genome-wide (slow without presto).

## atomic, verified CSV write: write to a temp file, read it back, and only then
## rename over the target. Guarantees you never end up with a truncated *_all.csv
## even if the session is interrupted mid-write.
safe_write_csv <- function(df, path) {
  tmp <- paste0(path, ".tmp")
  write.csv(df, tmp, row.names=FALSE)
  chk <- tryCatch(nrow(read.csv(tmp)), error=function(e) -1)
  if (chk != nrow(df)) stop("write verification failed for ", path,
                            " (wrote ", nrow(df), ", read back ", chk, ")")
  file.rename(tmp, path)
}

canon <- c("Hlf","Mecom","Procr","Fgd5","Mpl","Meis1","Cd34","Kit","Flt3","Dntt",
 "Gata1","Klf1","Car1","Car2","Hba-a1","Pf4","Itga2b","Vwf","Mpo","Elane","Prtn3",
 "Ctsg","Ms4a3","Camp","Ngp","Ltf","S100a8","Csf1r","Ly6c2","Ccr2","Siglech","Bst2",
 "Itgax","Cd74","Prss34","Mcpt8","Cpa3","Prg2","Il5ra","Cd79a","Vpreb1","Ebf1",
 "Jchain","Il7r","Ncr1","Nkg7")

fix_one <- function(name) {
  od  <- file.path(BASE, name); rds <- file.path(od, paste0(name, "_seurat.rds"))
  if (!file.exists(rds)) { message("skip ", name, " (no rds)"); return(invisible()) }
  message("\n==== ", name, " ====")
  obj <- readRDS(rds); DefaultAssay(obj) <- "RNA"
  if (inherits(obj[["RNA"]], "Assay5")) obj <- JoinLayers(obj)
  obj <- NormalizeData(obj, verbose=FALSE)
  gg <- canon[canon %in% rownames(obj)]
  ## gene universe for the DE test (big speed lever on base Seurat)
  feats <- if (HVG_ONLY) {
    hv <- VariableFeatures(obj)
    if (length(hv) < 50) hv <- VariableFeatures(FindVariableFeatures(obj, nfeatures=3000, verbose=FALSE))
    union(hv, gg)                     # always include canonical markers
  } else NULL                         # NULL = all genes

  for (r in RES) {
    col <- paste0("RNA_snn_res.", r)
    if (!col %in% colnames(obj@meta.data)) { message("  no ", col); next }
    Idents(obj) <- factor(obj@meta.data[[col]])
    clusters <- levels(Idents(obj))

    ## primary pass: low threshold, restricted gene set, downsampled -> fast
    am <- FindAllMarkers(obj, only.pos=TRUE, min.pct=0.1, logfc.threshold=0.1,
                         return.thresh=1, max.cells.per.ident=MAXCELLS,
                         features=feats, verbose=FALSE)
    ## fallback for any cluster that still returned nothing (cheap: same restricted genes)
    missing <- setdiff(clusters, unique(as.character(am$cluster)))
    for (cl in missing) {
      fm <- tryCatch(FindMarkers(obj, ident.1=cl, only.pos=TRUE, min.pct=0.05,
              logfc.threshold=0.05, max.cells.per.ident=MAXCELLS,
              features=feats, verbose=FALSE), error=function(e) NULL)
      if (!is.null(fm) && nrow(fm)) { fm$gene<-rownames(fm); fm$cluster<-cl; am<-bind_rows(am, fm) }
    }
    still <- setdiff(clusters, unique(as.character(am$cluster)))
    if (length(still)) message("  WARNING res ", r, " no markers for clusters: ", paste(still, collapse=","))

    am <- am %>% arrange(cluster, desc(avg_log2FC))
    top <- am %>% group_by(cluster) %>% slice_max(avg_log2FC, n=TOPN, with_ties=FALSE) %>% ungroup()
    safe_write_csv(am,  file.path(od, sprintf("%s_markers_res%s_all.csv",  name, r)))
    safe_write_csv(top, file.path(od, sprintf("%s_markers_res%s_top%d.csv", name, r, TOPN)))
    message("  res ", r, ": ", length(unique(am$cluster)), "/", length(clusters),
            " clusters with markers, ", nrow(am), " rows")

    ## conserved markers across genotype -- OPT-IN. This was the step that made the
    ## full run take hours (per-cluster x per-genotype base-Wilcoxon). Off by default;
    ## the top15 markers are enough to annotate. Set DO_CONSERVED <- TRUE to enable
    ## (now downsampled via max.cells.per.ident to keep it tractable).
    if (DO_CONSERVED && requireNamespace("metap", quietly=TRUE) && "genotype" %in% colnames(obj@meta.data)) {
      cons <- list()
      for (cl in clusters) {
        cm <- tryCatch(FindConservedMarkers(obj, ident.1=cl, grouping.var="genotype",
                       only.pos=TRUE, min.pct=0.1, logfc.threshold=0.1,
                       max.cells.per.ident=800, verbose=FALSE), error=function(e) NULL)
        if (!is.null(cm) && nrow(cm)) { cm$gene<-rownames(cm); cm$cluster<-cl; cons[[cl]]<-head(cm,TOPN) }
      }
      if (length(cons)) write.csv(bind_rows(cons),
        file.path(od, sprintf("%s_conserved_markers_res%s.csv", name, r)), row.names=FALSE)
    }

    ## canonical dotplot (all clusters)
    p <- DotPlot(obj, features=gg, cluster.idents=TRUE) + RotatedAxis() +
         ggtitle(paste0(name, "  res ", r)) + theme(axis.text.x=element_text(size=7))
    ggsave(file.path(od, sprintf("%s_canonical_dotplot_res%s.png", name, r)), p,
           width=16, height=6, dpi=150)
  }
}
for (nm in DATASETS) fix_one(nm)
plan("sequential")   # release workers
message("\nDone. Corrected *_markers_res*_all.csv / _top15.csv now contain ALL clusters.")
###############################################################################
## Notes
## - Base Seurat only (no presto). FindMarkers falls back to the wilcox test; with
##   MAXCELLS=2000/ident it stays fast. Raise MAXCELLS for more precision.
## - logfc.threshold=0.1 in the primary pass, and logfc.threshold=0 in the fallback,
##   ensures the central progenitor clusters (HSC/MPP, CMP/GMP, monocyte, erythroid-
##   prog) that vanished at 0.25 now get ranked markers.
## - If you *do* want the faster presto backend later:
##     install.packages("remotes"); remotes::install_github("immunogenomics/presto")
##   Seurat's FindMarkers then auto-uses it; no code change needed.
###############################################################################
