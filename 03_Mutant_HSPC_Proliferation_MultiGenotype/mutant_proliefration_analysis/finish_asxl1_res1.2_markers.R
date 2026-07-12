## Finish ONLY the missing piece: Asxl1 res 1.2 cluster markers (no conserved markers,
## which is what made the full run take hours). Fast, single-core, downsampled.
suppressPackageStartupMessages({library(Seurat); library(dplyr); library(ggplot2)})
od <- "C:/Users/fc809/Downloads/mutant_proliefration_analysis/seurat_harmony_output/Asxl1"
obj <- readRDS(file.path(od, "Asxl1_seurat.rds")); DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) obj <- JoinLayers(obj)
obj <- NormalizeData(obj, verbose=FALSE)
Idents(obj) <- factor(obj@meta.data[["RNA_snn_res.1.2"]])
am <- FindAllMarkers(obj, only.pos=TRUE, min.pct=0.1, logfc.threshold=0.1,
                     return.thresh=1, max.cells.per.ident=2000, verbose=FALSE)
for (cl in setdiff(levels(Idents(obj)), unique(as.character(am$cluster)))) {
  fm <- tryCatch(FindMarkers(obj, ident.1=cl, only.pos=TRUE, min.pct=0.05,
          logfc.threshold=0, max.cells.per.ident=2000, verbose=FALSE), error=function(e) NULL)
  if (!is.null(fm) && nrow(fm)) { fm$gene<-rownames(fm); fm$cluster<-cl; am<-bind_rows(am, fm) }
}
am  <- am %>% arrange(cluster, desc(avg_log2FC))
top <- am %>% group_by(cluster) %>% slice_max(avg_log2FC, n=15, with_ties=FALSE) %>% ungroup()
write.csv(am,  file.path(od, "Asxl1_markers_res1.2_all.csv"),   row.names=FALSE)
write.csv(top, file.path(od, "Asxl1_markers_res1.2_top15.csv"), row.names=FALSE)
message("Asxl1 res1.2: ", length(unique(am$cluster)), " clusters, ", nrow(am), " rows -> done")
