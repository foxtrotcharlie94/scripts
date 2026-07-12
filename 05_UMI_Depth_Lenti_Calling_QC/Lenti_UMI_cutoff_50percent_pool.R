## ============================================================
##  CDS size + coordinates for Kmt2d and Pkd1 (mm10)
##  Outputs:
##    cds_sizes_Kmt2d_Pkd1.csv       — total merged CDS bp per gene
##    cds_coordinates_Kmt2d_Pkd1.bed — one row per merged CDS exon block
##
##  NOTE: Run in a fresh R session, or detach ensembldb first:
##    detach("package:ensembldb", unload = TRUE)
## ============================================================

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

## ─── Gene list ────────────────────────────────────────────────────────────────

genes_of_interest <- c("Kmt2d", "Pkd1")

## ─── Map symbols to Entrez IDs ────────────────────────────────────────────────

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

entrez_ids <- AnnotationDbi::mapIds(org.Mm.eg.db,
                                    keys      = genes_of_interest,
                                    column    = "ENTREZID",
                                    keytype   = "SYMBOL",
                                    multiVals = "first")

missing_sym <- genes_of_interest[is.na(entrez_ids)]
if (length(missing_sym))
  message("No Entrez ID found for: ", paste(missing_sym, collapse = ", "))

entrez_ids <- entrez_ids[!is.na(entrez_ids)]

## ─── Get transcript IDs for each gene ────────────────────────────────────────

tx_tbl <- AnnotationDbi::select(txdb,
                                keys    = as.character(entrez_ids),
                                columns = c("GENEID", "TXID"),
                                keytype = "GENEID")
tx_tbl <- unique(tx_tbl)

sym_map       <- setNames(names(entrez_ids), as.character(entrez_ids))
tx_tbl$symbol <- sym_map[tx_tbl$GENEID]

## ─── Get CDS intervals as a plain data frame ─────────────────────────────────

cds_tbl <- AnnotationDbi::select(txdb,
                                 keys    = as.character(tx_tbl$TXID),
                                 columns = c("TXID", "CDSCHROM",
                                             "CDSSTART", "CDSEND", "CDSSTRAND"),
                                 keytype = "TXID")
cds_tbl <- unique(cds_tbl)
cds_tbl <- cds_tbl[!is.na(cds_tbl$CDSSTART), ]

# Join gene symbol using base merge
cds_tbl <- merge(cds_tbl,
                 tx_tbl[, c("TXID", "symbol")],
                 by   = "TXID",
                 all.x = TRUE)
cds_tbl <- cds_tbl[!is.na(cds_tbl$symbol), ]

message("CDS table rows: ", nrow(cds_tbl))
message("Genes with CDS data: ", length(unique(cds_tbl$symbol)))

## ─── Interval merging function ───────────────────────────────────────────────

merge_intervals <- function(starts, ends) {
  if (length(starts) == 0) return(list(starts = integer(0), ends = integer(0)))
  ord    <- order(starts)
  starts <- starts[ord]
  ends   <- ends[ord]
  out_s  <- starts[1]
  out_e  <- ends[1]
  for (i in seq_along(starts)[-1]) {
    if (starts[i] <= out_e[length(out_e)]) {
      out_e[length(out_e)] <- max(out_e[length(out_e)], ends[i])
    } else {
      out_s <- c(out_s, starts[i])
      out_e <- c(out_e, ends[i])
    }
  }
  list(starts = out_s, ends = out_e)
}

## ─── Per-gene loop ───────────────────────────────────────────────────────────

results_bed  <- list()
results_size <- list()

for (sym in names(entrez_ids)) {
  gene_cds <- cds_tbl[cds_tbl$symbol == sym, ]
  
  if (nrow(gene_cds) == 0) {
    message("No CDS found for ", sym)
    next
  }
  
  chrom  <- gene_cds$CDSCHROM[1]
  strand <- gene_cds$CDSSTRAND[1]
  
  merged   <- merge_intervals(gene_cds$CDSSTART, gene_cds$CDSEND)
  total_bp <- sum(merged$ends - merged$starts + 1L)
  
  # BED rows (0-based half-open)
  for (i in seq_along(merged$starts)) {
    results_bed[[length(results_bed) + 1]] <- data.frame(
      chrom  = chrom,
      start  = merged$starts[i] - 1L,
      end    = merged$ends[i],
      name   = sym,
      score  = 0L,
      strand = strand,
      stringsAsFactors = FALSE
    )
  }
  
  results_size[[sym]] <- data.frame(
    gene    = sym,
    entrez  = as.character(entrez_ids[[sym]]),
    chrom   = chrom,
    strand  = strand,
    cds_bp  = total_bp,
    cds_kb  = round(total_bp / 1000, 1),
    n_exons = length(merged$starts),
    stringsAsFactors = FALSE
  )
}

message("Loop complete: ", length(results_bed), " BED entries, ",
        length(results_size), " genes processed")

## ─── Write outputs ────────────────────────────────────────────────────────────

if (length(results_bed) == 0)
  stop("results_bed is empty — see messages above")
if (length(results_size) == 0)
  stop("results_size is empty — see messages above")

bed_df  <- do.call(rbind, results_bed)
size_df <- do.call(rbind, results_size)
size_df <- size_df[order(size_df$gene), ]

write.table(bed_df,
            "cds_coordinates_Kmt2d_Pkd1.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.csv(size_df, "cds_sizes_Kmt2d_Pkd1.csv", row.names = FALSE)

## ─── Console summary ──────────────────────────────────────────────────────────

message("\n", strrep("-", 55))
message(sprintf("%-12s  %8s  %7s  %s", "Gene", "CDS (bp)", "CDS (kb)", "Exons"))
message(strrep("-", 55))
for (i in seq_len(nrow(size_df))) {
  message(sprintf("%-12s  %8s  %7.1f  %d",
                  size_df$gene[i],
                  formatC(size_df$cds_bp[i], format = "d", big.mark = ","),
                  size_df$cds_kb[i],
                  size_df$n_exons[i]))
}
message(strrep("-", 55))
message(sprintf("TOTAL  %s bp  (%.1f kb)  across %d genes",
                formatC(sum(size_df$cds_bp), format = "d", big.mark = ","),
                sum(size_df$cds_bp) / 1000,
                nrow(size_df)))
message(strrep("-", 55))
message("\nOutput files written:")
message("  cds_sizes_Kmt2d_Pkd1.csv")
message("  cds_coordinates_Kmt2d_Pkd1.bed")

