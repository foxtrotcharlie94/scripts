# ============================================================
# MOUSE LIVER MUTAGENESIS ASSAY — CANDIDATE REGION FINDER
# Genome: mm10 | Tissue: liver (adult, 2 months)
# ------------------------------------------------------------
# 0. PACKAGES
# ------------------------------------------------------------
if (!require("BiocManager")) install.packages("BiocManager")

BiocManager::install(c(
  "GenomicRanges",
  "rtracklayer",
  "AnnotationHub",
  "BSgenome.Mmusculus.UCSC.mm10",
  "Biostrings",
  "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "biomaRt"
), update = FALSE)

install.packages(c("tidyverse", "httr", "jsonlite"))

library(GenomicRanges)
library(rtracklayer)
library(AnnotationHub)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(biomaRt)
library(tidyverse)
library(httr)
library(jsonlite)

# ------------------------------------------------------------
# 0b. DECOMPRESS ALL INPUT FILES
# ------------------------------------------------------------
if (!require("R.utils")) install.packages("R.utils")
library(R.utils)

gz_files <- c(
  "mm10-blacklist.v2.bed.gz",
  "mm10-cCREs.bed.gz",
  "rmsk.txt.gz",
  "genomicSuperDups.txt.gz",
  "mm10.k100.umap.bedgraph.gz",
  "liver_mm10_H3K27ac_peaks.bed.gz",
  "liver_mm10_H3K4me1_peaks.bed.gz",
  "liver_mm10_H3K4me3_peaks.bed.gz",
  "liver_mm10_CTCF_peaks.bed.gz",
  "liver_mm10_H3K27me3_peaks.bed.gz",
  "liver_mm10_H3K36me3_peaks.bed.gz",
  "liver_mm10_H3K9me3_peaks.bed.gz"
)

for (gz in gz_files) {
  if (file.exists(gz)) {
    cat("Decompressing", gz, "...\n")
    R.utils::gunzip(gz, remove = TRUE, overwrite = TRUE)
  }
}


genome  <- BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm10")
txdb    <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Autosomes only
autosomes   <- paste0("chr", 1:19)
chrom_sizes <- seqlengths(genome)[autosomes]

whole_genome <- GRanges(
  seqnames = names(chrom_sizes),
  ranges   = IRanges(start = 1, end = chrom_sizes)
)

# ============================================================
# 1. GENERIC EXCLUSION TRACKS
# ============================================================

# ------------------------------------------------------------
# 1a. ENCODE Blacklist v2
# ------------------------------------------------------------
blacklist <- import("mm10-blacklist.v2.bed")
blacklist <- keepSeqlevels(blacklist, autosomes, pruning.mode = "coarse")

# ------------------------------------------------------------
# 1b. Gene bodies + 10kb flanking buffer
# ------------------------------------------------------------
genes <- suppressMessages(genes(txdb))
seqlevelsStyle(genes) <- "UCSC"
genes <- keepSeqlevels(genes, autosomes, pruning.mode = "coarse")

genes_buffered <- IRanges::reduce(GenomicRanges::trim(GenomicRanges::resize(genes, width = GenomicRanges::width(genes) + 20000, fix = "center")))

transcripts_all <- transcripts(txdb)
seqlevelsStyle(transcripts_all) <- "UCSC"
transcripts_all <- keepSeqlevels(transcripts_all, autosomes, pruning.mode = "coarse")

transcripts_buffered <- IRanges::reduce(GenomicRanges::trim(GenomicRanges::resize(transcripts_all, width = GenomicRanges::width(transcripts_all) + 20000, fix = "center")))

# ------------------------------------------------------------
# 1c. ENCODE cCREs
# ------------------------------------------------------------
ccres_raw <- read_tsv("mm10-cCREs.bed", col_names = FALSE, col_select = 1:3,
                      show_col_types = FALSE)
ccres <- GRanges(seqnames = ccres_raw$X1,
                 ranges   = IRanges(start = ccres_raw$X2 + 1, end = ccres_raw$X3))
ccres <- keepSeqlevels(ccres, autosomes, pruning.mode = "coarse")

ccres_buffered <- GenomicRanges::resize(ccres, width = GenomicRanges::width(ccres) + 4000, fix = "center") |> GenomicRanges::trim()

# ------------------------------------------------------------
# 1d. RepeatMasker (raw UCSC table format)
# ------------------------------------------------------------
rmsk_raw <- read_tsv(
  "rmsk.txt",
  col_names = c("bin", "swScore", "milliDiv", "milliDel", "milliIns",
                "genoName", "genoStart", "genoEnd", "genoLeft",
                "strand", "repName", "repClass", "repFamily",
                "repStart", "repEnd", "repLeft", "id"),
  comment   = "#",
  show_col_types = FALSE
)

repeats <- GRanges(
  seqnames = rmsk_raw$genoName,
  ranges   = IRanges(start = rmsk_raw$genoStart + 1,  # 0-based to 1-based
                     end   = rmsk_raw$genoEnd)
)
repeats <- keepSeqlevels(repeats, autosomes, pruning.mode = "coarse")
repeats <- IRanges::reduce(repeats)

# ------------------------------------------------------------
# 1e. Segmental Duplications (raw UCSC table format)
# ------------------------------------------------------------
segdups_raw <- read_tsv(
  "genomicSuperDups.txt",
  col_names = FALSE,
  show_col_types = FALSE
)

segdups <- GRanges(
  seqnames = segdups_raw$X2,
  ranges   = IRanges(start = segdups_raw$X3 + 1,  # 0-based to 1-based
                     end   = segdups_raw$X4)
)
segdups <- keepSeqlevels(segdups, autosomes, pruning.mode = "coarse")
segdups <- IRanges::reduce(segdups)

# ------------------------------------------------------------
# 1f. Imprinted regions (generated from Geneimprint mouse database)
# ------------------------------------------------------------
imprinted_raw <- read_tsv("mm10_imprinted_regions.bed",
                          col_names = c("chrom", "start", "end"),
                          show_col_types = FALSE)
imprinted_raw$chrom <- ifelse(grepl("^chr", imprinted_raw$chrom),
                              imprinted_raw$chrom,
                              paste0("chr", imprinted_raw$chrom))
imprinted_buffered <- GRanges(
  seqnames = imprinted_raw$chrom,
  ranges   = IRanges(start = imprinted_raw$start, end = imprinted_raw$end)
)
imprinted_buffered <- imprinted_buffered[as.character(seqnames(imprinted_buffered)) %in% autosomes]

# ============================================================
# 2. LIVER-SPECIFIC EXCLUSION TRACKS
# ============================================================

read_bed3 <- function(path) {
  raw <- read_tsv(path, col_names = FALSE, col_select = 1:3,
                  show_col_types = FALSE, comment = "#")
  GRanges(seqnames = raw$X1,
          ranges   = IRanges(start = raw$X2 + 1, end = raw$X3))
}

liver_h3k27ac  <- read_bed3("liver_mm10_H3K27ac_peaks.bed")
liver_h3k4me1  <- read_bed3("liver_mm10_H3K4me1_peaks.bed")
liver_h3k4me3  <- read_bed3("liver_mm10_H3K4me3_peaks.bed")
liver_ctcf     <- read_bed3("liver_mm10_CTCF_peaks.bed")
liver_h3k27me3 <- read_bed3("liver_mm10_H3K27me3_peaks.bed")
liver_h3k36me3 <- read_bed3("liver_mm10_H3K36me3_peaks.bed")
liver_h3k9me3  <- read_bed3("liver_mm10_H3K9me3_peaks.bed")

# Keep only autosomes
for (obj_name in c("liver_h3k27ac", "liver_h3k4me1", "liver_h3k4me3",
                   "liver_ctcf", "liver_h3k27me3", "liver_h3k36me3",
                   "liver_h3k9me3")) {
  obj <- get(obj_name)
  obj <- keepSeqlevels(obj, autosomes, pruning.mode = "coarse")
  assign(obj_name, obj)
}

# Add buffers around active marks
liver_h3k27ac_buffered  <- GenomicRanges::resize(liver_h3k27ac,  width = GenomicRanges::width(liver_h3k27ac)  + 2000, fix = "center") |> GenomicRanges::trim()
liver_h3k4me1_buffered  <- GenomicRanges::resize(liver_h3k4me1,  width = GenomicRanges::width(liver_h3k4me1)  + 2000, fix = "center") |> GenomicRanges::trim()
liver_h3k4me3_buffered  <- GenomicRanges::resize(liver_h3k4me3,  width = GenomicRanges::width(liver_h3k4me3)  + 2000, fix = "center") |> GenomicRanges::trim()
liver_ctcf_buffered     <- GenomicRanges::resize(liver_ctcf,     width = GenomicRanges::width(liver_ctcf)     + 1000, fix = "center") |> GenomicRanges::trim()
liver_h3k27me3_buffered <- GenomicRanges::resize(liver_h3k27me3, width = GenomicRanges::width(liver_h3k27me3) + 1000, fix = "center") |> GenomicRanges::trim()
liver_h3k36me3_buffered <- GenomicRanges::resize(liver_h3k36me3, width = GenomicRanges::width(liver_h3k36me3) + 1000, fix = "center") |> GenomicRanges::trim()
liver_h3k9me3_buffered  <- GenomicRanges::resize(liver_h3k9me3,  width = GenomicRanges::width(liver_h3k9me3)  + 1000, fix = "center") |> GenomicRanges::trim()

# ============================================================
# 3. BUILD FULL EXCLUSION MASK
# ============================================================
exclusion_combined <- unlist(GRangesList(
  blacklist,
  genes_buffered,
  transcripts_buffered,
  ccres_buffered,
  repeats,
  segdups,
  imprinted_buffered,
  liver_h3k27ac_buffered,
  liver_h3k4me1_buffered,
  liver_h3k4me3_buffered,
  liver_ctcf_buffered,
  liver_h3k27me3_buffered,
  liver_h3k36me3_buffered,
  liver_h3k9me3_buffered
))

# Intersect with whole_genome to hard-clip all intervals to chromosome boundaries
exclusion_combined <- exclusion_combined[as.character(seqnames(exclusion_combined)) %in% autosomes]
seqlevels(exclusion_combined, pruning.mode = "coarse") <- autosomes
exclusion_mask <- IRanges::reduce(GenomicRanges::intersect(exclusion_combined, whole_genome))

# Save mask
rtracklayer::export(exclusion_mask, "exclusion_mask_mm10_liver.bed", format = "BED")

genome_size <- sum(as.numeric(chrom_sizes))
cat("Exclusion mask covers",
    round(sum(as.numeric(GenomicRanges::width(exclusion_mask))) / genome_size * 100, 1),
    "% of autosomes\n")

# ============================================================
# 4. FIND CANDIDATE WINDOWS
# ============================================================
# Strategy: tile genome into 200kb windows, keep those with
# >= 60kb of total unmasked sequence (need not be contiguous).
# For each passing window, usable sequence = setdiff with mask.

TILE_SIZE      <- 200000L
MIN_USABLE     <- 60000L

tiles <- tileGenome(chrom_sizes, tilewidth = TILE_SIZE, cut.last.tile.in.chrom = TRUE)
seqlevels(tiles) <- autosomes

# For each tile, calculate total unmasked bases
usable <- GenomicRanges::setdiff(tiles, exclusion_mask)
seqlevels(usable)  <- autosomes
seqlengths(usable) <- chrom_sizes

# Sum usable bases per tile
usable_per_tile <- sapply(seq_along(tiles), function(i) {
  hits <- subsetByOverlaps(usable, tiles[i])
  sum(GenomicRanges::width(hits))
})

tiles$usable_bp <- usable_per_tile
candidates <- tiles[tiles$usable_bp >= MIN_USABLE]

# Set seqlengths so getSeq() works correctly
seqlengths(candidates) <- chrom_sizes[as.character(seqlevels(candidates))]

cat("Candidate 200kb windows with >= 60kb usable sequence:", length(candidates), "\n")

# ============================================================
# 5. SCORE CANDIDATES
# ============================================================

# ------------------------------------------------------------
# 5a & 5b. GC content and CpG O/E — computed over unmasked
#           sequence only within each candidate window
# ------------------------------------------------------------
cpg_oe_fn <- function(seq) {
  n   <- length(seq)
  cpg <- dinucleotideFrequency(seq)["CG"]
  c_f <- letterFrequency(seq, "C")
  g_f <- letterFrequency(seq, "G")
  if (c_f == 0 | g_f == 0) return(NA)
  (cpg * n) / (c_f * g_f)
}

cat("Computing GC and CpG over unmasked sequence per window...\n")
gc_vals  <- numeric(length(candidates))
cpg_vals <- numeric(length(candidates))

for (i in seq_along(candidates)) {
  # Get unmasked sub-intervals within this tile
  unmasked_i <- subsetByOverlaps(usable, candidates[i])
  if (length(unmasked_i) == 0) { gc_vals[i] <- NA; cpg_vals[i] <- NA; next }
  seqs_i <- getSeq(genome, unmasked_i)
  # Concatenate into one sequence for accurate dinucleotide counting
  combined <- do.call(xscat, as.list(seqs_i))
  gc_vals[i]  <- letterFrequency(combined, letters = "GC", as.prob = TRUE)
  cpg_vals[i] <- cpg_oe_fn(combined)
}

candidates$gc_content <- gc_vals
candidates$cpg_oe     <- cpg_vals

# ------------------------------------------------------------
# 5c. Mappability (100-mer, Umap bedgraph)
# ------------------------------------------------------------
mappability_gr <- import(
  "mm10.k100.umap.bedgraph",
  format = "bedGraph"
)
# Keep only autosomes
mappability_gr <- keepSeqlevels(mappability_gr, autosomes, pruning.mode = "coarse")

cat("Scoring mappability (vectorized)...\n")

# Find all overlaps between mappability intervals and candidate tiles directly
map_hits   <- findOverlaps(mappability_gr, candidates)
map_idx    <- queryHits(map_hits)
cand_idx   <- subjectHits(map_hits)

# Compute actual overlap width (mappability interval clipped to tile boundary)
map_starts <- pmax(start(mappability_gr)[map_idx], start(candidates)[cand_idx])
map_ends   <- pmin(end(mappability_gr)[map_idx],   end(candidates)[cand_idx])
overlap_bp <- pmax(0L, map_ends - map_starts + 1L)

# Weighted mean mappability per candidate tile
map_df <- data.frame(
  cand_idx   = cand_idx,
  score      = mappability_gr$score[map_idx],
  overlap_bp = overlap_bp
)

map_by_tile <- map_df |>
  group_by(cand_idx) |>
  summarise(mappability = weighted.mean(score, overlap_bp), .groups = "drop")

# Assign back — tiles with no mappability data get 0
candidates$mappability <- 0
candidates$mappability[map_by_tile$cand_idx] <- map_by_tile$mappability

# ============================================================
# 6. FILTER AND RANK
# ============================================================
scored <- as.data.frame(candidates) |>
  as_tibble() |>
  filter(
    gc_content  >= 0.40,
    gc_content  <= 0.60,
    mappability >= 0.90,
    !is.na(cpg_oe),
    cpg_oe      <= 0.30
  ) |>
  mutate(
    # Lower score = better candidate
    # Favor: high mappability, low CpG, GC close to 50%
    score = rank(-mappability) +
      rank(cpg_oe) +
      rank(abs(gc_content - 0.50))
  ) |>
  arrange(score) |>
  mutate(rank = row_number())

cat("Candidate windows passing all filters:", nrow(scored), "\n")

# ============================================================
# 7. SELECT FINAL CANDIDATES — top 2 per chromosome
# ============================================================
top_per_chrom <- scored |>
  group_by(seqnames) |>
  slice_min(score, n = 2) |>
  ungroup() |>
  arrange(rank)

cat("\nTop candidates per chromosome:\n")
print(
  top_per_chrom |>
    dplyr::select(seqnames, start, end, width, gc_content,
                  mappability, cpg_oe, rank),
  n = 40
)

# ============================================================
# 8. EXPORT
# ============================================================
write_tsv(scored, "candidates_ranked_full.tsv")

top_gr <- makeGRangesFromDataFrame(top_per_chrom, keep.extra.columns = TRUE)
rtracklayer::export(top_gr, "top_candidates_liver.bed", format = "BED")

# Export the actual usable fragment coordinates within the top candidate tiles.
# These are the sequences to target with capture probes or amplicons.
top_usable <- subsetByOverlaps(usable, top_gr)
# Annotate each fragment with its parent tile rank
top_usable$tile_rank <- NA_integer_
for (i in seq_along(top_gr)) {
  hits <- which(countOverlaps(top_usable, top_gr[i]) > 0)
  top_usable$tile_rank[hits] <- top_per_chrom$rank[i]
}
top_usable <- top_usable[order(top_usable$tile_rank)]
rtracklayer::export(top_usable, "top_candidates_usable_fragments.bed", format = "BED")

cat("\nDone. Output files written:\n")
cat("  exclusion_mask_mm10_liver.bed          — full exclusion mask\n")
cat("  candidates_ranked_full.tsv             — all passing windows, ranked\n")
cat("  top_candidates_liver.bed               — top 2 per chromosome tile boundaries\n")
cat("  top_candidates_usable_fragments.bed    — actual fragments to capture within top tiles\n")

# ============================================================
# Fragment selector for mutagenesis assay
# ------------------------------------------------------------
# Requires these objects already in memory from the main pipeline:
#   candidates     — GRanges of 1969 scored 200kb windows
#   usable         — GRanges of all unmasked genome fragments
#   chrom_sizes    — named integer vector of chromosome lengths
#   autosomes      — character vector of autosome names
#
# Writes to working/outputs/:
#   assay_fragments_225bp.bed       — final selected 225bp tiles
#   assay_fragments_225bp_full.tsv  — same with parent-tile metadata
# ============================================================

candidates <- GRanges(
  seqnames = scored$seqnames,
  ranges   = IRanges(start = scored$start, end = scored$end),
  usable_bp   = scored$usable_bp,
  gc_content  = scored$gc_content,
  cpg_oe      = scored$cpg_oe,
  mappability = scored$mappability,
  score       = scored$score,
  rank        = scored$rank
)
# ============================================================
# Fragment selector for mutagenesis assay
# ------------------------------------------------------------
# Requires these objects already in memory from the main pipeline:
#   candidates     — GRanges of 1969 scored 200kb windows
#   usable         — GRanges of all unmasked genome fragments
#   chrom_sizes    — named integer vector of chromosome lengths
#   autosomes      — character vector of autosome names
#
# Writes to working/outputs/:
#   assay_fragments_225bp.bed       — final selected 225bp tiles
#   assay_fragments_225bp_full.tsv  — same with parent-tile metadata
# ============================================================

# ============================================================
# Fragment selector for mutagenesis assay
# ------------------------------------------------------------
# Requires these objects already in memory from the main pipeline:
#   candidates     — GRanges of 1969 scored 200kb windows
#   usable         — GRanges of all unmasked genome fragments
#   chrom_sizes    — named integer vector of chromosome lengths
#   autosomes      — character vector of autosome names
#
# Writes to working/outputs/:
#   assay_fragments_225bp.bed       — final selected 225bp tiles
#   assay_fragments_225bp_full.tsv  — same with parent-tile metadata
# ============================================================

candidates <- GRanges(
  seqnames = scored$seqnames,
  ranges   = IRanges(start = scored$start, end = scored$end),
  usable_bp   = scored$usable_bp,
  gc_content  = scored$gc_content,
  cpg_oe      = scored$cpg_oe,
  mappability = scored$mappability,
  score       = scored$score,
  rank        = scored$rank
)

library(GenomicRanges)
library(dplyr)
library(readr)

TILE_SIZE    <- 225L          # target fragment size (bp)
MIN_FRAG     <- 200L          # minimum acceptable fragment (bp)
TARGET_TOTAL <- 60000L        # target total bp
OUT_DIR      <- "working/outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Attach parent-tile rank/score to every usable fragment ────────────────
cat("Overlapping usable fragments with candidate tiles...\n")

# Find which candidate tile each usable fragment falls in
hits        <- findOverlaps(usable, candidates, type = "within")
frag_idx    <- queryHits(hits)
tile_idx    <- subjectHits(hits)

# Build a data frame of usable fragments annotated with parent-tile metadata
frag_df <- data.frame(
  chr        = as.character(seqnames(usable))[frag_idx],
  frag_start = start(usable)[frag_idx],
  frag_end   = end(usable)[frag_idx],
  tile_rank  = candidates$rank[tile_idx],
  tile_score = candidates$score[tile_idx],
  gc         = candidates$gc_content[tile_idx],
  cpg_oe     = candidates$cpg_oe[tile_idx],
  mappability= candidates$mappability[tile_idx],
  tile_start = start(candidates)[tile_idx],
  tile_end   = end(candidates)[tile_idx]
)

# Each usable fragment should belong to exactly one tile, but guard against
# fragments spanning tile boundaries by keeping only the highest-ranked hit
frag_df <- frag_df |>
  group_by(chr, frag_start, frag_end) |>
  slice_min(tile_rank, n = 1, with_ties = FALSE) |>
  ungroup()

cat("  Total usable fragments across all candidates:", nrow(frag_df), "\n")

# ── 2. Chop each usable fragment into 225bp tiles ────────────────────────────
cat("Chopping into", TILE_SIZE, "bp tiles...\n")

chop_fragment <- function(chr, frag_start, frag_end,
                          tile_rank, tile_score, gc, cpg_oe, mappability,
                          tile_start, tile_end) {
  frag_len <- frag_end - frag_start + 1L
  if (frag_len < MIN_FRAG) return(NULL)
  
  n_tiles  <- frag_len %/% TILE_SIZE
  remainder <- frag_len %% TILE_SIZE
  
  # Build tile boundaries
  starts <- frag_start + (seq_len(n_tiles) - 1L) * TILE_SIZE
  ends   <- starts + TILE_SIZE - 1L
  
  # Keep remainder as a final tile only if it meets minimum length
  if (remainder >= MIN_FRAG) {
    last_start <- frag_start + n_tiles * TILE_SIZE
    starts <- c(starts, last_start)
    ends   <- c(ends, frag_end)
  } else if (n_tiles == 0L) {
    return(NULL)   # fragment too short even for one tile
  }
  
  data.frame(
    chr         = chr,
    start       = starts,
    end         = ends,
    width       = ends - starts + 1L,
    tile_rank   = tile_rank,
    tile_score  = tile_score,
    gc          = gc,
    cpg_oe      = cpg_oe,
    mappability = mappability,
    tile_start  = tile_start,
    tile_end    = tile_end
  )
}

mini_tiles <- mapply(
  chop_fragment,
  frag_df$chr, frag_df$frag_start, frag_df$frag_end,
  frag_df$tile_rank, frag_df$tile_score,
  frag_df$gc, frag_df$cpg_oe, frag_df$mappability,
  frag_df$tile_start, frag_df$tile_end,
  SIMPLIFY = FALSE
)
mini_tiles <- do.call(rbind, Filter(Negate(is.null), mini_tiles))

cat("  Total 225bp candidate tiles:", nrow(mini_tiles), "\n")
cat("  Total bp available:", sum(mini_tiles$width), "\n")
cat("  Chromosomes represented:", length(unique(mini_tiles$chr)), "\n")

# ── 3. Select ~60kb spread evenly across chromosomes ─────────────────────────
cat("\nSelecting ~", TARGET_TOTAL, "bp distributed across chromosomes...\n")

chrs_present <- sort(unique(mini_tiles$chr))
n_chrs       <- length(chrs_present)
target_per_chr <- ceiling(TARGET_TOTAL / n_chrs)

cat("  Chromosomes:", n_chrs, "| Target per chromosome:", target_per_chr, "bp\n")

selected <- lapply(chrs_present, function(chr) {
  chr_tiles <- mini_tiles |>
    filter(.data$chr == .env$chr) |>
    arrange(tile_rank, tile_score)      # best-ranked parent tiles first
  
  # Greedy selection: accumulate until target_per_chr reached
  cumulative <- cumsum(chr_tiles$width)
  keep       <- which(cumulative <= target_per_chr)
  
  # If no single tile fits (shouldn't happen), take at least one
  if (length(keep) == 0L) keep <- 1L
  
  chr_tiles[keep, ]
})

selected_df <- do.call(rbind, selected)

# Summary
total_bp  <- sum(selected_df$width)
n_frags   <- nrow(selected_df)
n_chrs_out <- length(unique(selected_df$chr))

cat("\n── Selection summary ────────────────────────────────────\n")
cat("  Fragments selected :", n_frags, "\n")
cat("  Total bp selected  :", total_bp, "\n")
cat("  Chromosomes covered:", n_chrs_out, "\n")
cat("  Mean fragment size :", round(mean(selected_df$width), 1), "bp\n")

cat("\nPer-chromosome breakdown:\n")
selected_df |>
  group_by(chr) |>
  summarise(n_frags = n(), total_bp = sum(width), .groups = "drop") |>
  arrange(chr) |>
  as.data.frame() |>
  print()

# ── 4. Export ─────────────────────────────────────────────────────────────────
cat("\nWriting outputs...\n")

# BED file (0-based)
bed_out <- data.frame(
  chrom      = selected_df$chr,
  chromStart = selected_df$start - 1L,   # convert to 0-based
  chromEnd   = selected_df$end,
  name       = paste0("frag_rank", selected_df$tile_rank,
                      "_", selected_df$chr,
                      ":", selected_df$start - 1L,
                      "-", selected_df$end),
  score      = round(selected_df$tile_score),
  strand     = "."
)
write_tsv(bed_out, file.path(OUT_DIR, "assay_fragments_225bp.bed"),
          col_names = FALSE)

# Full TSV with all metadata
full_out <- selected_df |>
  mutate(
    parent_tile = paste0(chr, ":", tile_start - 1L, "-", tile_end),
    fragment    = paste0(chr, ":", start - 1L, "-", end)
  ) |>
  select(fragment, chr, start, end, width,
         parent_tile, tile_rank, tile_score,
         gc, cpg_oe, mappability)
write_tsv(full_out, file.path(OUT_DIR, "assay_fragments_225bp_full.tsv"))

cat("  Wrote: working/outputs/assay_fragments_225bp.bed\n")
cat("  Wrote: working/outputs/assay_fragments_225bp_full.tsv\n")
cat("\nDone.\n")

