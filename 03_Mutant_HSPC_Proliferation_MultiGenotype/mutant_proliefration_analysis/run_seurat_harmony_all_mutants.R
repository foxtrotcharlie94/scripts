###############################################################################
## Seurat + Harmony pipeline with cell-cycle regression for all mutant datasets
## Produces, per dataset: UMAP, clusters at TWO resolutions, QC/CC plots, saved object.
##
## Datasets (mouse BM HSPC scRNA-seq):
##   Tet2         GSE209994  Lin- BM,  HTO-multiplexed   (4 WT / 4 KO)
##   Asxl1        GSE158184  Lin- BM,  HTO-multiplexed   (7 WT / 7 KO)
##   Jak2         GSE227026  Lin-cKit+ (LK), CellRanger .h5 (3 WT / 3 mut)
##   Dnmt3a       GSE272266  LSK, unmanipulated arm      (4 WT / 4 Dnmt3a)
##   Dnmt3a-R878H GSE233963  LSK, vehicle arm            (4 WT / 4 R878H)
##
## Batch key for Harmony = individual mouse.  Cell cycle regressed out before PCA.
## Edit BASE if your folder is elsewhere.  Run top-to-bottom (or source()).
###############################################################################

## ------------------------- parallel: use 10 cores -------------------------
## The heavy steps here are C++/BLAS (ScaleData regression solve, RunPCA, Harmony)
## and uwot (UMAP) -- NOT things future/multisession speeds up. So the optimal use
## of 10 cores is a MULTITHREADED BLAS + uwot threads, which is also memory-light
## (no per-worker object copies). Datasets are run sequentially to bound RAM.
NCORES <- 10
Sys.setenv(OMP_NUM_THREADS = NCORES, OPENBLAS_NUM_THREADS = NCORES,
           MKL_NUM_THREADS = NCORES, VECLIB_MAXIMUM_THREADS = NCORES,
           RCPP_PARALLEL_NUM_THREADS = NCORES)

## ------------------------- 0. packages ------------------------------------
pkgs <- c("Seurat","harmony","dplyr","ggplot2","Matrix","patchwork","RhpcBLASctl","metap","hdf5r")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE))
  install.packages(p, repos="https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(dplyr); library(ggplot2)
  library(Matrix); library(patchwork)
})
set.seed(1)
if (!exists("%||%")) `%||%` <- function(a,b) if (is.null(a)) b else a   # null-coalesce
if (requireNamespace("RhpcBLASctl", quietly=TRUE)) {
  RhpcBLASctl::blas_set_num_threads(NCORES)   # PCA / ScaleData / Harmony linear algebra
  RhpcBLASctl::omp_set_num_threads(NCORES)
  message("BLAS threads: ", RhpcBLASctl::blas_get_num_procs(),
          " (needs a threaded BLAS: OpenBLAS/MKL — check with sessionInfo())")
}

## ------------------------- 1. config --------------------------------------
BASE    <- "C:/Users/fc809/Downloads/mutant_proliefration_analysis"   # <- edit if needed
OUTROOT <- file.path(BASE, "seurat_harmony_output")
TMP     <- file.path(BASE, "_extracted")          # tars are unpacked here
dir.create(OUTROOT, showWarnings=FALSE); dir.create(TMP, showWarnings=FALSE)

RES        <- c(0.4, 1.2)     # the TWO clustering resolutions (coarse, fine)
DIMS       <- 1:30            # harmony/UMAP dims
QC_MINGENE <- 500; QC_MAXGENE <- 7000; QC_MAXMT <- 10

## Mouse cell-cycle gene sets (Tirosh/Seurat orthologs) -----------------------
s.genes <- c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6",
  "Cdca7","Dtl","Prim1","Uhrf1","Mlf1ip","Cenpu","Hells","Rfc2","Rpa2","Nasp",
  "Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51",
  "Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn",
  "Pola1","Chaf1b","Brip1","E2f8")
g2m.genes <- c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2",
  "Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Fam64a","Pimreg","Smc4","Ccnb2",
  "Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp",
  "Cdca3","Hn1","Jpt1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5",
  "Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe",
  "Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa")

## ------------------------- 2. helpers -------------------------------------

## read a 10x mtx triplet given a directory + a GSM/file prefix -> sparse matrix
read_triplet <- function(dir, prefix) {
  fm <- list.files(dir, pattern=paste0("^", prefix, ".*matrix\\.mtx\\.gz$"),  full.names=TRUE)[1]
  ff <- list.files(dir, pattern=paste0("^", prefix, ".*features\\.tsv\\.gz$"), full.names=TRUE)[1]
  fb <- list.files(dir, pattern=paste0("^", prefix, ".*barcodes\\.tsv\\.gz$"), full.names=TRUE)[1]
  ReadMtx(mtx=fm, features=ff, cells=fb, feature.column=2)   # rownames = gene symbols
}
## read the feature-type column (3rd col) for a triplet (to split GEX vs HTO/ADT)
read_ftypes <- function(dir, prefix) {
  ff <- list.files(dir, pattern=paste0("^", prefix, ".*features\\.tsv\\.gz$"), full.names=TRUE)[1]
  read.delim(gzfile(ff), header=FALSE)
}

## Core pipeline: takes a Seurat object with meta cols `genotype` and `mouse`,
## regresses cell cycle, runs Harmony over `mouse`, UMAP, clusters at 2 res.
run_pipeline <- function(obj, name) {
  message(">>> ", name, ": ", ncol(obj), " cells before QC")
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern="^mt-")
  obj <- subset(obj, subset = nFeature_RNA >= QC_MINGENE &
                              nFeature_RNA <  QC_MAXGENE &
                              percent.mt   <  QC_MAXMT)
  message("    ", ncol(obj), " cells after QC")

  obj <- NormalizeData(obj, verbose=FALSE)
  obj <- CellCycleScoring(obj, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)
  obj <- FindVariableFeatures(obj, nfeatures=2000, verbose=FALSE)
  ## regress out cell cycle (clustering only); scores kept in metadata
  obj <- ScaleData(obj, vars.to.regress=c("S.Score","G2M.Score"), verbose=FALSE)
  obj <- RunPCA(obj, npcs=max(DIMS), verbose=FALSE)
  obj <- RunHarmony(obj, group.by.vars="mouse", verbose=FALSE)
  obj <- RunUMAP(obj, reduction="harmony", dims=DIMS, n.threads=NCORES, verbose=FALSE)
  obj <- FindNeighbors(obj, reduction="harmony", dims=DIMS, verbose=FALSE)
  for (r in RES) obj <- FindClusters(obj, resolution=r, verbose=FALSE)  # cols: RNA_snn_res.<r>

  ## ---- outputs ----
  od <- file.path(OUTROOT, name); dir.create(od, showWarnings=FALSE, recursive=TRUE)
  rlow  <- paste0("RNA_snn_res.", RES[1])
  rhigh <- paste0("RNA_snn_res.", RES[2])

  p1 <- DimPlot(obj, group.by=rlow,  label=TRUE) + ggtitle(paste0(name, "  res=", RES[1]))
  p2 <- DimPlot(obj, group.by=rhigh, label=TRUE) + ggtitle(paste0(name, "  res=", RES[2]))
  p3 <- DimPlot(obj, group.by="genotype") + ggtitle(paste0(name, "  genotype"))
  p4 <- DimPlot(obj, group.by="Phase")    + ggtitle(paste0(name, "  cell-cycle phase"))
  ggsave(file.path(od, paste0(name,"_UMAP_res",RES[1],".png")), p1, width=7, height=6, dpi=150)
  ggsave(file.path(od, paste0(name,"_UMAP_res",RES[2],".png")), p2, width=7, height=6, dpi=150)
  ggsave(file.path(od, paste0(name,"_UMAP_genotype.png")),      p3, width=7, height=6, dpi=150)
  ggsave(file.path(od, paste0(name,"_UMAP_phase.png")),         p4, width=7, height=6, dpi=150)
  ggsave(file.path(od, paste0(name,"_UMAP_panel.png")),
         (p1|p2)/(p3|p4), width=14, height=12, dpi=150)

  meta <- obj@meta.data[, c("mouse","genotype","Phase","S.Score","G2M.Score",
                            "nFeature_RNA","percent.mt", rlow, rhigh)]
  write.csv(meta, file.path(od, paste0(name,"_metadata_clusters.csv")))

  ## ---- cluster identity: markers + conserved markers + canonical dotplot ----
  om <- obj
  if (inherits(om[["RNA"]], "Assay5")) om <- JoinLayers(om)   # v5: markers need joined layers
  canon <- c("Hlf","Mecom","Procr","Fgd5","Mpl","Meis1","Cd34","Kit","Flt3","Dntt",
             "Gata1","Klf1","Car1","Car2","Hba-a1","Pf4","Itga2b","Vwf","Mpo","Elane",
             "Prtn3","Ctsg","Ms4a3","Camp","Ngp","Ltf","S100a8","Csf1r","Ly6c2","Ccr2",
             "Siglech","Bst2","Itgax","Cd74","Prss34","Mcpt8","Cpa3","Prg2","Il5ra",
             "Cd79a","Vpreb1","Ebf1","Jchain","Cd3e","Il7r","Ncr1","Nkg7")
  canon <- canon[canon %in% rownames(om)]
  for (r in RES) {
    col <- paste0("RNA_snn_res.", r); Idents(om) <- factor(om@meta.data[[col]])
    ## robust markers: low logfc.threshold + per-cluster fallback so the large,
    ## transcriptionally-central progenitor clusters are never dropped (they vanish
    ## at the default logfc.threshold=0.25). Base Seurat, no presto needed.
    am <- FindAllMarkers(om, only.pos=TRUE, min.pct=0.1, logfc.threshold=0.1,
                         return.thresh=1, max.cells.per.ident=2000, verbose=FALSE)
    for (cl in setdiff(levels(Idents(om)), unique(as.character(am$cluster)))) {
      fm <- tryCatch(FindMarkers(om, ident.1=cl, only.pos=TRUE, min.pct=0.05,
              logfc.threshold=0, max.cells.per.ident=2000, verbose=FALSE), error=function(e) NULL)
      if (!is.null(fm) && nrow(fm)) { fm$gene<-rownames(fm); fm$cluster<-cl; am<-bind_rows(am, fm) }
    }
    if (nrow(am)) {
      am <- am %>% arrange(cluster, desc(avg_log2FC))
      write.csv(am, file.path(od, sprintf("%s_markers_res%s_all.csv", name, r)), row.names=FALSE)
      top <- am %>% group_by(cluster) %>% slice_max(avg_log2FC, n=15, with_ties=FALSE) %>% ungroup()
      write.csv(top, file.path(od, sprintf("%s_markers_res%s_top15.csv", name, r)), row.names=FALSE)
    }
    cons <- list()
    for (cl in levels(Idents(om))) {
      cm <- tryCatch(FindConservedMarkers(om, ident.1=cl, grouping.var="genotype",
                     only.pos=TRUE, min.pct=0.1, logfc.threshold=0.1, verbose=FALSE), error=function(e) NULL)
      if (!is.null(cm) && nrow(cm)) { cm$gene<-rownames(cm); cm$cluster<-cl; cons[[cl]]<-head(cm,15) }
    }
    if (length(cons)) write.csv(bind_rows(cons),
        file.path(od, sprintf("%s_conserved_markers_res%s.csv", name, r)), row.names=FALSE)
    ggsave(file.path(od, sprintf("%s_canonical_dotplot_res%s.png", name, r)),
           DotPlot(om, features=canon, cluster.idents=TRUE) + RotatedAxis() +
           ggtitle(paste0(name," (res ",r,")")) + theme(axis.text.x=element_text(size=7)),
           width=16, height=6, dpi=150)
  }
  ## cluster x genotype composition at each resolution
  write.csv(as.matrix(table(obj@meta.data[[rlow]],  obj$genotype)),
            file.path(od, paste0(name,"_composition_res",RES[1],".csv")))
  write.csv(as.matrix(table(obj@meta.data[[rhigh]], obj$genotype)),
            file.path(od, paste0(name,"_composition_res",RES[2],".csv")))
  saveRDS(obj, file.path(od, paste0(name,"_seurat.rds")))
  message("    wrote outputs -> ", od)
  invisible(obj)
}

## HTO demultiplex helper: given full counts (genes+HTO) and the row indices of
## the hashtag features, build a Seurat obj with an HTO assay and call HTODemux.
make_hto_obj <- function(counts, hto_rows, sample_tag, genotype, hto_names=NULL) {
  gex <- counts[-hto_rows, , drop=FALSE]
  hto <- counts[ hto_rows, , drop=FALSE]
  ## give hashtags explicit, unique names (features file often reuses the same name
  ## for both tags, e.g. Tet2's "WT_Unt"). hto_names must match hto_rows order.
  if (!is.null(hto_names)) rownames(hto) <- hto_names
  keep <- Matrix::colSums(gex) > 0 & Matrix::colSums(hto) > 0
  gex <- gex[, keep]; hto <- hto[, keep]
  colnames(gex) <- paste0(sample_tag, "_", colnames(gex))
  colnames(hto) <- colnames(gex)
  o <- CreateSeuratObject(gex, project=sample_tag, min.cells=3, min.features=200)
  hto <- hto[, colnames(o)]
  o[["HTO"]] <- CreateAssayObject(counts=hto)
  o <- NormalizeData(o, assay="HTO", normalization.method="CLR", verbose=FALSE)
  o <- HTODemux(o, assay="HTO", positive.quantile=0.99)
  o <- subset(o, subset = HTO_classification.global == "Singlet")
  o$mouse    <- as.character(o$hash.ID)   # individual mouse = hashtag identity (now uniquely named)
  o$genotype <- genotype
  DefaultAssay(o) <- "RNA"
  o
}

## ------------------------- 3. TET2  (GSE209994, HTO) -----------------------
## 4 loose matrices; last-2 features are hashtag_1/2 -> reps within each file.
load_tet2 <- function() {
  base <- BASE
  ## reps = the two hashtags (in file/hashtag order) -> distinct mouse names
  cfg <- list(
    list(f="GSE209994_WT_Unt_rep1_rep2", g="WT", reps=c("WT_rep1","WT_rep2")),
    list(f="GSE209994_WT_Unt_rep3_rep4", g="WT", reps=c("WT_rep3","WT_rep4")),
    list(f="GSE209994_KO_Unt_rep1_rep2", g="KO", reps=c("KO_rep1","KO_rep2")),
    list(f="GSE209994_KO_Unt_rep3_rep4", g="KO", reps=c("KO_rep3","KO_rep4")))
  objs <- lapply(cfg, function(c) {
    pre <- paste0(c$f, "_filtered_feature_bc_matrix_")
    m  <- ReadMtx(mtx =file.path(base,paste0(pre,"matrix.mtx.gz")),
                  features=file.path(base,paste0(pre,"features.tsv.gz")),
                  cells   =file.path(base,paste0(pre,"barcodes.tsv.gz")), feature.column=2)
    ft <- read.delim(gzfile(file.path(base,paste0(pre,"features.tsv.gz"))), header=FALSE)
    hto_rows <- which(ft$V3 != "Gene Expression")   # 2 hashtag rows, in order 1,2
    make_hto_obj(m, hto_rows, c$f, c$g, hto_names=c$reps)
  })
  merge(objs[[1]], objs[-1])
}

## ------------------------- 4. ASXL1 (GSE158184, HTO) -----------------------
## per-mouse hashtags named e.g. wk04_wt_1 ; keep only this sample's hashtags.
load_asxl1 <- function() {
  base <- BASE
  cfg <- list(
    list(p="GSM4794842_wk04_wt_",  tag="wk04_wt",  g="WT"),
    list(p="GSM4794845_wk36_wt_",  tag="wk36_wt",  g="WT"),
    list(p="GSM4794843_wk04_mut_", tag="wk04_mut", g="KO"),
    list(p="GSM4794844_wk36_mut_", tag="wk36_mut", g="KO"))
  objs <- lapply(cfg, function(c) {
    m  <- ReadMtx(mtx =file.path(base,paste0(c$p,"matrix.mtx.gz")),
                  features=file.path(base,paste0(c$p,"features.tsv.gz")),
                  cells   =file.path(base,paste0(c$p,"barcodes.tsv.gz")), feature.column=2)
    ft <- read.delim(gzfile(file.path(base,paste0(c$p,"features.tsv.gz"))), header=FALSE)
    ## hashtag rows for THIS sample = Antibody Capture whose name starts with the tag
    hto_rows <- which(ft$V3!="Gene Expression" & grepl(paste0("^",c$tag,"_[0-9]+$"), ft$V2))
    make_hto_obj(m, hto_rows, c$tag, c$g)
  })
  merge(objs[[1]], objs[-1])
}

## ------------------------- 5. JAK2 (GSE227026, .h5 in tars) ----------------
load_jak2 <- function() {
  map <- list(GSM7090165="Jak2_WT_1", GSM7090166="Jak2_WT_2", GSM7090167="Jak2_WT_3",
              GSM7090162="Jak2_Mut_1",GSM7090163="Jak2_Mut_2",GSM7090164="Jak2_Mut_3")
  objs <- lapply(names(map), function(gsm) {
    ex <- file.path(TMP, gsm); dir.create(ex, showWarnings=FALSE)
    untar(file.path(BASE, paste0(gsm, ".tar")), exdir=ex)
    h5 <- list.files(ex, pattern="\\.h5$", recursive=TRUE, full.names=TRUE)[1]
    if (is.na(h5)) stop("no .h5 found for ", gsm)
    m  <- Read10X_h5(h5)                       # needs hdf5r; raw (unfiltered)
    if (is.list(m)) m <- m[["Gene Expression"]] %||% m[[1]]   # multi-modal h5 -> take GEX
    o  <- CreateSeuratObject(m, project=map[[gsm]], min.cells=3, min.features=200)
    o$mouse    <- map[[gsm]]
    o$genotype <- ifelse(grepl("WT", map[[gsm]]), "WT", "Jak2mut")
    o
  })
  merge(objs[[1]], objs[-1])
}

## ------------------------- 6. DNMT3A (GSE272266, tar) ----------------------
## unmanipulated arm GSM8397470-477 ; prefixes from the RAW.tar.
load_dnmt3a <- function() {
  ex <- file.path(TMP, "GSE272266"); dir.create(ex, showWarnings=FALSE)
  if (length(list.files(ex))==0) untar(file.path(BASE,"GSE272266_RAW.tar"), exdir=ex)
  map <- list(
    GSM8397470="WT",    GSM8397471="WT",    GSM8397472="WT",    GSM8397473="WT",
    GSM8397474="Dnmt3a",GSM8397475="Dnmt3a",GSM8397476="Dnmt3a",GSM8397477="Dnmt3a")
  objs <- lapply(names(map), function(gsm) {
    m <- read_triplet(ex, gsm)
    o <- CreateSeuratObject(m, project=gsm, min.cells=3, min.features=200)
    o$mouse <- gsm; o$genotype <- map[[gsm]]; o
  })
  merge(objs[[1]], objs[-1])
}

## ------------------------- 7. DNMT3A-R878H (GSE233963, tars) ---------------
load_r878 <- function() {
  map <- list(GSM7439859="WT", GSM7439863="WT", GSM7439867="WT", GSM7439871="WT",
              GSM7439861="R878H",GSM7439865="R878H",GSM7439869="R878H",GSM7439873="R878H")
  objs <- lapply(names(map), function(gsm) {
    ex <- file.path(TMP, gsm); dir.create(ex, showWarnings=FALSE)
    if (length(list.files(ex))==0) untar(file.path(BASE, paste0(gsm,".tar")), exdir=ex)
    m <- read_triplet(ex, gsm)
    o <- CreateSeuratObject(m, project=gsm, min.cells=3, min.features=200)
    o$mouse <- gsm; o$genotype <- map[[gsm]]; o
  })
  merge(objs[[1]], objs[-1])
}

## ------------------------- 8. run everything ------------------------------
datasets <- list(
  Tet2         = load_tet2,
  Asxl1        = load_asxl1,
  Jak2         = load_jak2,
  Dnmt3a       = load_dnmt3a,
  Dnmt3a_R878H = load_r878
)

for (nm in names(datasets)) {
  message("\n=================  ", nm, "  =================")
  obj <- tryCatch(datasets[[nm]](), error=function(e){message("LOAD FAILED: ",conditionMessage(e)); NULL})
  if (!is.null(obj)) tryCatch(run_pipeline(obj, nm),
                              error=function(e) message("PIPELINE FAILED: ", conditionMessage(e)))
}
message("\nAll done. Outputs in: ", OUTROOT)
###############################################################################
## Notes
## - Parallelism (NCORES=10): the real cost here is C++/BLAS (ScaleData regression,
##   RunPCA, Harmony) + uwot (UMAP), so we drive 10 cores via a THREADED BLAS
##   (RhpcBLASctl + OMP/OPENBLAS/MKL env vars) and RunUMAP(n.threads=10). This only
##   helps if R is linked against a threaded BLAS (OpenBLAS or MKL) -- check
##   sessionInfo()$BLAS; stock reference BLAS is single-threaded. On Windows,
##   Microsoft R Open / an MKL build, or radian+OpenBLAS, gives multithreaded BLAS.
##   FindNeighbors (annoy) is largely single-threaded regardless.
## - Deliberately NOT using future::multisession: for these steps it copies the
##   whole object to each worker (RAM blowup) and speeds up almost nothing.
## - The 5 datasets are independent; if you have plenty of RAM you can instead run
##   them concurrently with future.apply::future_lapply(names(datasets), ...) and
##   set BLAS threads to ~2 each. Default here is sequential + threaded BLAS (safe).
## - Two resolutions live in RES (0.4 coarse, 1.2 fine); metadata cols are
##   RNA_snn_res.0.4 and RNA_snn_res.1.2. Change RES and re-run FindClusters only.
## - Harmony batch = individual mouse. For a lighter batch key use `genotype`.
## - Tet2/Asxl1 are hashtag-multiplexed -> HTODemux recovers per-mouse identity;
##   only singlets are kept. If HTODemux is too strict, lower positive.quantile.
## - Jak2 .h5 are RAW (unfiltered) matrices; CreateSeuratObject(min.features=200)
##   plus the QC subset removes empty droplets. For stricter calling use DropletUtils.
## - To also score proliferation on UN-regressed data, use the S.Score/G2M.Score
##   already in meta.data, or AddModuleScore with your proliferation-core gene list
##   BEFORE ScaleData's regression (expression layer is unaffected by regression).
###############################################################################
