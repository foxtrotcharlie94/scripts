## ============================================================
##  Lentivirus Vector × Mouse Transcriptome Overlap Pipeline
## ============================================================

suppressPackageStartupMessages({
  library(Biostrings)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(GenomicRanges)
  library(IRanges)
  library(GenomicFeatures)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(ensembldb)
  library(EnsDb.Mmusculus.v79)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(patchwork)
  library(scales)
})

## ─── USER SETTINGS ────────────────────────────────────────────────────────────

VECTOR_FASTA  <- "my_vector.fasta"   # path to your vector

OUT_DIR <- "overlap_results"

INCLUDE_MRNA       <- TRUE
INCLUDE_LNCRNA     <- TRUE
INCLUDE_PREMRNA    <- TRUE

WORD_SIZE     <- 11
EVALUE        <- 1e-2
PERC_IDENTITY <- 90
MIN_ALIGN_LEN <- 20
NUM_THREADS   <- 4

# BLAST+ binary directory
BLAST_BIN <- "C:/Program Files/NCBI/blast-2.17.0+/bin"

## ─── SET BLAST PATH ───────────────────────────────────────────────────────────

Sys.setenv(PATH = paste(BLAST_BIN, Sys.getenv("PATH"), sep = ";"))
MAKEBLASTDB <- file.path(BLAST_BIN, "makeblastdb")
BLASTN      <- file.path(BLAST_BIN, "blastn")

## ─── 1. READ VECTOR ───────────────────────────────────────────────────────────

message("\n[1/6] Reading vector sequence ...")

if (!file.exists(VECTOR_FASTA))
  stop("Vector FASTA not found: ", VECTOR_FASTA)

vector_seqs <- readDNAStringSet(VECTOR_FASTA)
message("  Loaded ", length(vector_seqs), " sequence(s), total ",
        sum(width(vector_seqs)), " bp")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

## ─── 2. BUILD TRANSCRIPTOME BLAST DATABASE ────────────────────────────────────

message("\n[2/6] Building transcriptome database ...")

txdb         <- TxDb.Mmusculus.UCSC.mm10.knownGene
edb          <- EnsDb.Mmusculus.v79
mouse_genome <- BSgenome.Mmusculus.UCSC.mm10
db_path      <- file.path(OUT_DIR, "mm10_transcriptome_blastdb")
txome_fasta  <- file.path(OUT_DIR, "mm10_transcriptome.fa")
db_flag      <- paste0(db_path, ".nhr")

if (!file.exists(db_flag)) {
  
  all_seqs <- DNAStringSet()
  
  if (INCLUDE_MRNA) {
    message("  Extracting spliced mRNA sequences ...")
    mrna_seqs <- extractTranscriptSeqs(mouse_genome, exonsBy(txdb, by = "tx"))
    names(mrna_seqs) <- paste0("mRNA|tx", names(mrna_seqs))
    all_seqs <- c(all_seqs, mrna_seqs)
    message("    ", length(mrna_seqs), " mRNA transcript sequences")
  }
  
  if (INCLUDE_LNCRNA) {
    message("  Extracting lncRNA sequences ...")
    lnc_tx <- transcripts(
      edb,
      filter = TxBiotypeFilter(c("lincRNA", "lncRNA",
                                 "antisense", "processed_transcript",
                                 "retained_intron", "sense_intronic",
                                 "sense_overlapping")),
      columns = c("tx_id", "tx_biotype", "gene_name")
    )
    seqlevelsStyle(lnc_tx) <- "UCSC"
    lnc_tx <- lnc_tx[seqnames(lnc_tx) %in% seqnames(mouse_genome)]
    
    if (length(lnc_tx) > 0) {
      lnc_exons <- exonsBy(
        edb,
        by = "tx",
        filter = TxBiotypeFilter(c("lincRNA", "lncRNA",
                                   "antisense", "processed_transcript",
                                   "retained_intron", "sense_intronic",
                                   "sense_overlapping"))
      )
      seqlevelsStyle(lnc_exons) <- "UCSC"
      lnc_exons <- lnc_exons[
        sapply(lnc_exons, function(g)
          all(as.character(seqnames(g)) %in% seqnames(mouse_genome)))
      ]
      lnc_seqs        <- extractTranscriptSeqs(mouse_genome, lnc_exons)
      names(lnc_seqs) <- paste0("lncRNA|tx", names(lnc_seqs))
      all_seqs <- c(all_seqs, lnc_seqs)
      message("    ", length(lnc_seqs), " lncRNA transcript sequences")
    } else {
      message("    No lncRNA transcripts found")
    }
  }
  
  if (INCLUDE_PREMRNA) {
    message("  Extracting pre-mRNA gene body sequences ...")
    gene_ranges <- genes(txdb)
    gene_ranges <- gene_ranges[seqnames(gene_ranges) %in% seqnames(mouse_genome)]
    premrna_seqs        <- getSeq(mouse_genome, gene_ranges)
    names(premrna_seqs) <- paste0("pre-mRNA|gene", names(gene_ranges))
    all_seqs <- c(all_seqs, premrna_seqs)
    message("    ", length(premrna_seqs), " gene body sequences")
  }
  
  message("  Removing duplicate sequences ...")
  all_seqs <- unique(all_seqs)
  message("  Total unique sequences for BLAST DB: ", length(all_seqs))
  
  writeXStringSet(all_seqs, filepath = txome_fasta)
  
  message("  Running makeblastdb ...")
  system2(MAKEBLASTDB,
          args = c("-in",    txome_fasta,
                   "-dbtype", "nucl",
                   "-out",   db_path,
                   "-title", "mm10_transcriptome"),
          stdout = TRUE, stderr = TRUE)
  
} else {
  message("  Found existing transcriptome database at ", db_path)
}

## ─── 3. RUN BLASTN ────────────────────────────────────────────────────────────

message("\n[3/6] Running BLASTn against transcriptome ...")

query_fasta <- file.path(OUT_DIR, "query_vector.fa")
blast_out   <- file.path(OUT_DIR, "blast_results.tsv")

writeXStringSet(vector_seqs, filepath = query_fasta)

blast_fmt <- paste(
  "6 qseqid sseqid pident length mismatch gapopen",
  "qstart qend sstart send evalue bitscore sstrand"
)

system2(BLASTN,
        args = c("-query",         query_fasta,
                 "-db",            db_path,
                 "-out",           blast_out,
                 "-outfmt",        shQuote(blast_fmt),
                 "-word_size",     WORD_SIZE,
                 "-evalue",        EVALUE,
                 "-perc_identity", PERC_IDENTITY,
                 "-num_threads",   NUM_THREADS,
                 "-dust",          "yes",
                 "-soft_masking",  "true"),
        stdout = TRUE, stderr = TRUE)

## ─── 4. PARSE & FILTER ────────────────────────────────────────────────────────

message("\n[4/6] Parsing BLAST results ...")

col_names <- c("query_id", "subject_tx", "pct_identity", "align_len",
               "mismatches", "gap_opens", "q_start", "q_end",
               "s_start", "s_end", "evalue", "bitscore", "strand")

hits_raw <- tryCatch(
  read.table(blast_out, sep = "\t", header = FALSE,
             col.names = col_names, stringsAsFactors = FALSE),
  error = function(e) NULL
)

if (is.null(hits_raw) || nrow(hits_raw) == 0) {
  message("\n  No hits — your vector has no detectable homology to the mouse transcriptome.")
  quit(save = "no")
}

hits <- hits_raw %>%
  dplyr::filter(align_len >= MIN_ALIGN_LEN) %>%
  mutate(
    strand    = ifelse(strand == "plus", "+", "-"),
    biotype   = sub("\\|.*", "", subject_tx),
    tx_id     = sub(".*\\|tx", "", subject_tx),
    query_len = width(vector_seqs)[match(query_id, names(vector_seqs))]
  ) %>%
  arrange(evalue, desc(bitscore))

message("  ", nrow(hits), " hits pass filters  (>=", MIN_ALIGN_LEN,
        " bp | >=", PERC_IDENTITY, "% id | E <= ", EVALUE, ")")

## ─── 5. ANNOTATE HITS ─────────────────────────────────────────────────────────

message("\n[5/6] Annotating hits ...")

entrez_to_sym <- function(ids) {
  tryCatch(
    mapIds(org.Mm.eg.db, keys = as.character(ids),
           column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"),
    error = function(e) setNames(rep(NA_character_, length(ids)), ids)
  )
}

hits$gene_symbol <- NA_character_

mrna_idx <- which(hits$biotype %in% c("mRNA", "pre-mRNA"))
if (length(mrna_idx) > 0) {
  tx_tbl <- select(txdb,
                   keys    = unique(hits$tx_id[mrna_idx]),
                   columns = c("TXID", "GENEID"),
                   keytype = "TXID") %>% distinct()
  entrez_map <- entrez_to_sym(tx_tbl$GENEID)
  sym_map    <- setNames(entrez_map[as.character(tx_tbl$GENEID)], tx_tbl$TXID)
  hits$gene_symbol[mrna_idx] <- sym_map[hits$tx_id[mrna_idx]]
}

lnc_idx <- which(hits$biotype == "lncRNA")
if (length(lnc_idx) > 0) {
  lnc_info <- tryCatch(
    transcripts(edb,
                filter  = TxIdFilter(unique(hits$tx_id[lnc_idx])),
                columns = c("tx_id", "gene_name")) %>%
      as.data.frame() %>% select(tx_id, gene_name) %>% distinct(),
    error = function(e) data.frame(tx_id = character(), gene_name = character())
  )
  lnc_sym_map <- setNames(lnc_info$gene_name, lnc_info$tx_id)
  hits$gene_symbol[lnc_idx] <- lnc_sym_map[hits$tx_id[lnc_idx]]
}

hits$tx_region <- NA_character_
hits$tx_region[hits$biotype == "pre-mRNA"] <- "intron/exon (pre-mRNA)"
hits$tx_region[hits$biotype == "lncRNA"]   <- "lncRNA exon"

mrna_only <- hits %>% dplyr::filter(biotype == "mRNA")
if (nrow(mrna_only) > 0) {
  cds_ranges <- cdsBy(txdb, by = "tx")
  utr5       <- fiveUTRsByTranscript(txdb)
  utr3       <- threeUTRsByTranscript(txdb)
  
  classify_tx_pos <- function(tx_id, s_start, s_end) {
    tid <- as.character(tx_id)
    pos <- IRanges(start = s_start, end = s_end)
    
    check_overlap <- function(feature_list, tid, pos) {
      if (!tid %in% names(feature_list)) return(FALSE)
      feat_ir <- IRanges(start = start(feature_list[[tid]]),
                         end   = end(feature_list[[tid]]))
      any(overlapsAny(pos, feat_ir))
    }
    
    in_cds  <- check_overlap(cds_ranges, tid, pos)
    in_5utr <- check_overlap(utr5,       tid, pos)
    in_3utr <- check_overlap(utr3,       tid, pos)
    
    dplyr::case_when(in_cds  ~ "CDS",
                     in_5utr ~ "5' UTR",
                     in_3utr ~ "3' UTR",
                     TRUE    ~ "mRNA exon (unclassified)")
  }
  
  hits$tx_region[hits$biotype == "mRNA"] <- mapply(
    classify_tx_pos,
    tx_id   = hits$tx_id[hits$biotype == "mRNA"],
    s_start = hits$s_start[hits$biotype == "mRNA"],
    s_end   = hits$s_end[hits$biotype == "mRNA"]
  )
}

## ─── 6. OUTPUTS ───────────────────────────────────────────────────────────────

message("\n[6/6] Writing outputs ...")

# Flatten any list columns to character
hits <- hits %>%
  mutate(across(where(is.list), ~ sapply(., function(x)
    paste(unlist(x), collapse = ";"))))

write.csv(hits,
          file.path(OUT_DIR, "vector_transcriptome_overlaps_full.csv"),
          row.names = FALSE)

write.csv(hits,
          file.path(OUT_DIR, "vector_transcriptome_overlaps_full.csv"),
          row.names = FALSE)

summary_tbl <- hits %>%
  group_by(query_id) %>%
  summarise(
    query_length_bp   = first(query_len),
    n_hits            = n(),
    total_bp_hit      = sum(align_len),
    pct_query_hit     = round(100 * sum(align_len) / first(query_len), 1),
    mean_pct_identity = round(mean(pct_identity), 1),
    best_evalue       = min(evalue),
    n_CDS             = sum(tx_region == "CDS",     na.rm = TRUE),
    n_UTR             = sum(grepl("UTR", tx_region), na.rm = TRUE),
    n_lncRNA          = sum(biotype == "lncRNA"),
    n_premrna         = sum(biotype == "pre-mRNA"),
    genes_hit         = paste(sort(unique(na.omit(gene_symbol))), collapse = ";"),
    .groups = "drop"
  )

write.csv(summary_tbl,
          file.path(OUT_DIR, "vector_transcriptome_overlaps_summary.csv"),
          row.names = FALSE)

region_colors <- c(
  "CDS"                      = "#e63946",
  "5' UTR"                   = "#f4a261",
  "3' UTR"                   = "#e9c46a",
  "mRNA exon (unclassified)" = "#a8dadc",
  "lncRNA exon"              = "#2a9d8f",
  "intron/exon (pre-mRNA)"   = "#adb5bd"
)

p_cov <- ggplot(hits,
                aes(xmin = q_start, xmax = q_end,
                    ymin = 0, ymax = 1, fill = tx_region)) +
  geom_rect(alpha = 0.75, colour = NA) +
  facet_wrap(~ query_id, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = region_colors, na.value = "#cccccc") +
  labs(title = "Transcriptome hits mapped onto vector",
       x = "Position in vector (bp)", fill = "Transcript region", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(), legend.position = "top")

p_bio <- hits %>%
  count(biotype) %>%
  ggplot(aes(x = reorder(biotype, n), y = n, fill = biotype)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = n), hjust = -0.2, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  scale_fill_manual(values = c(mRNA = "#457b9d", lncRNA = "#2a9d8f",
                               "pre-mRNA" = "#adb5bd")) +
  coord_flip() +
  labs(title = "Hits by RNA biotype", x = NULL, y = "Count") +
  theme_minimal(base_size = 11)

p_reg <- hits %>%
  dplyr::filter(!is.na(tx_region)) %>%
  count(tx_region) %>%
  ggplot(aes(x = reorder(tx_region, n), y = n, fill = tx_region)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = n), hjust = -0.2, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  scale_fill_manual(values = region_colors) +
  coord_flip() +
  labs(title = "Hits by transcript region", x = NULL, y = "Count") +
  theme_minimal(base_size = 11)

p_id <- ggplot(hits, aes(x = pct_identity, fill = biotype)) +
  geom_histogram(binwidth = 1, colour = "white", linewidth = 0.2, position = "stack") +
  scale_fill_manual(values = c(mRNA = "#457b9d", lncRNA = "#2a9d8f",
                               "pre-mRNA" = "#adb5bd")) +
  labs(title = "% identity distribution",
       x = "% identity", y = "Count", fill = "Biotype") +
  theme_minimal(base_size = 11)

combined <- p_cov /
  (p_bio | p_reg | p_id) +
  plot_layout(heights = c(1.6, 1)) +
  plot_annotation(
    title    = "Lentivirus Vector - Mouse Transcriptome (mm10) Overlap",
    subtitle = sprintf("%d hits  |  E <= %.0e  |  >=%d%% identity  |  >=%d bp",
                       nrow(hits), EVALUE, PERC_IDENTITY, MIN_ALIGN_LEN),
    theme    = theme(plot.title    = element_text(face = "bold", size = 14),
                     plot.subtitle = element_text(size = 10, colour = "grey40"))
  )

ggsave(file.path(OUT_DIR, "overlap_summary_plot.pdf"),
       combined, width = 14, height = 10)

## ─── CONSOLE SUMMARY ──────────────────────────────────────────────────────────

message("\n", strrep("-", 65))
message("  PIPELINE COMPLETE")
message(strrep("-", 65))
message("\nPer-element summary:")
print(as.data.frame(summary_tbl), row.names = FALSE)
message("\nHits by transcript region:")
print(table(hits$tx_region))
message("\nHits by RNA biotype:")
print(table(hits$biotype))
all_genes <- sort(unique(na.omit(hits$gene_symbol)))
if (length(all_genes))
  message("\nGenes with transcriptome overlap (n=", length(all_genes), "):\n  ",
          paste(all_genes, collapse = ", "))
message(strrep("-", 65))
message("Output files in: ", normalizePath(OUT_DIR))
message("  vector_transcriptome_overlaps_full.csv")
message("  vector_transcriptome_overlaps_summary.csv")
message("  overlap_summary_plot.pdf")
message(strrep("-", 65))

make_mask <- function(seq_name) {
  seq_len  <- length(vector_seqs[[seq_name]])
  seq_hits <- dplyr::filter(hits, query_id == seq_name)
  mask     <- logical(seq_len)
  if (nrow(seq_hits) == 0) return(mask)
  for (i in seq_len(nrow(seq_hits)))
    mask[seq_hits$q_start[i]:seq_hits$q_end[i]] <- TRUE
  mask
}

# Deleted
deleted_seqs <- DNAStringSet(lapply(names(vector_seqs), function(seq_name) {
  vector_seqs[[seq_name]][!make_mask(seq_name)]
}))
names(deleted_seqs) <- paste0(names(vector_seqs), "_overlaps_deleted")
writeXStringSet(deleted_seqs, file.path(OUT_DIR, "vector_overlaps_deleted.fasta"))

# N-masked
masked_seqs <- DNAStringSet(lapply(names(vector_seqs), function(seq_name) {
  seq  <- vector_seqs[[seq_name]]
  mask <- make_mask(seq_name)
  if (!any(mask)) return(seq)
  seq_chars       <- strsplit(as.character(seq), "")[[1]]
  seq_chars[mask] <- "N"
  DNAString(paste(seq_chars, collapse = ""))
}))
names(masked_seqs) <- paste0(names(vector_seqs), "_overlaps_Nmasked")
writeXStringSet(masked_seqs, file.path(OUT_DIR, "vector_overlaps_Nmasked.fasta"))

