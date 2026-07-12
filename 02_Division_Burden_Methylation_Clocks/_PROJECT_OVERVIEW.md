# Project overview — Division-Burden / Methylation-Clock Validation (GSE44117)

**Biological question:** does forced proliferation (via serial transplantation
of just 10 HSCs) produce DNA-methylation changes resembling natural
chronological aging, as claimed by Figure 2E of a 2024 Stem Cell Reports paper?
This project independently reanalyzes a public RRBS dataset (GSE44117,
Beerman et al. 2013) to test that claim.

**Design:** mouse HSC RRBS methylation data across young/old baseline and
young/old 10-cell-transplant-reconstitution groups. Four published mouse
epigenetic clocks are applied (Blood/Petkovich, WLMT/Meer, YOMT, Liver), and
the key test is whether a CpG's natural-aging methylation delta (old minus
young baseline) correlates with its forced-proliferation delta (reconstituted
minus young baseline).

**Pipeline order:**
1. `analyze_gse44117.py` — coarse global-methylation proxy (superseded).
2. `apply_clocks.py` / `apply_clocks2.py` — compute real clock scores/deltas.
3. `cpg_delta_correlation.py` — the core hypothesis test (Blood + WLMT).
4. `cpg_delta_overlay.py` → `overlay3.py` → `overlay4.py` — presentation
   iterations adding the YOMT clock and refining figure style.
5. `plot_clock_scores.py` — group-level bar charts of clock scores.

**Note:** a final R script (`run_analysis.R` in the connected
`clock_correlation_analysis/` folder, with its own `REPORT.md`) reproduces
this whole Python pipeline end-to-end; that R script is out of scope for this
Python-only portfolio.
