################################################################################
## Generate _cross_comparison for every pseudobulk_DE_GSEA cell-type folder that
## HAS the GSEA inputs but is MISSING the cross-comparison output.
##
## Robust to the analysis script overwriting common globals: all driver state
## uses dotted names (.cc_*) that DE_GSEA_analysis.R does not use, and the loop
## iterates over a precomputed index so a clobbered list can't break it.
## Safe to re-run: folders that already have _cross_comparison are skipped.
##
## Run:  source("C:/Users/fc809/Downloads/run_missing_cross_comparison.R")
################################################################################

.cc_base   <- "C:/Users/fc809/Downloads/pseudobulk_DE_GSEA"
.cc_script <- "C:/Users/fc809/Downloads/DE_GSEA_analysis.R"
.cc_comps  <- c("HB_PosVsNeg", "LB_PosVsNeg", "LentiNeg_HBvsLB", "LentiPos_HBvsLB")

.cc_dirs <- list.dirs(.cc_base, recursive = FALSE)
.cc_todo <- .cc_dirs[ vapply(.cc_dirs, function(x)
  all(dir.exists(file.path(x, .cc_comps))) && !dir.exists(file.path(x, "_cross_comparison")),
  logical(1)) ]

if (length(.cc_todo) == 0) {
  message("Nothing to do — every folder with inputs already has _cross_comparison.")
} else {
  message("Will generate cross-comparison for ", length(.cc_todo), " folder(s):\n  ",
          paste(basename(.cc_todo), collapse = "\n  "))
  .cc_names <- basename(.cc_todo)          # saved up front; not re-derived after source()
  for (.cc_i in seq_along(.cc_todo)) {
    .cc_nm  <- .cc_names[.cc_i]
    pop_dir <- .cc_todo[.cc_i]             # picked up by DE_GSEA_analysis.R (guarded)
    message("\n================ ", .cc_nm, " ================")
    .cc_ok <- tryCatch({ source(.cc_script); TRUE },
                       error = function(e) { message("  FAILED: ", conditionMessage(e)); FALSE })
    message(if (isTRUE(.cc_ok)) paste0("  done -> ", .cc_nm, "/_cross_comparison")
            else                "  skipped (see error above)")
  }
  message("\nFinished.")
}
