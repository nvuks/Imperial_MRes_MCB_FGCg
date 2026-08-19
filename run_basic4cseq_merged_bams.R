#!/usr/bin/env Rscript
#
# run_basic4cseq_merged.R
#
# Runs Basic4Cseq independently for each barcode/integration site produced
# by split_bam_by_barcode.py, using each barcode's known whitelist
# coordinate as the 4C viewpoint. Unlike FourCSeq, Basic4Cseq is natively
# single-viewpoint-per-run, which matches this experiment's actual design
# (9 distinct integration sites, no real replicate structure) -- see chat
# history for the full reasoning behind switching from FourCSeq.
#
# IMPORTANT: this script expects the virtual restriction-fragment library
# to already exist -- build it ONCE with build_fragment_library.R before
# running this script (and re-run that script if you ever change read
# length, genome build, or restriction enzymes):
#   Rscript build_fragment_library.R <read_length> <library_csv_path>
#
# Background / design notes:
#   - firstCutter  = "CATG"  (NlaIII, primary cutter -- the site
#     trimmed-and-restored by trim_forward.sh, so R1 reads start exactly here)
#   - secondCutter = "GATC"  (DpnII, secondary cutter)
#   - Each barcode's viewpoint is supplied directly from the whitelist
#     coordinate (parsed from the manifest.tsv position label), as a small
#     interval centered on that coordinate -- there is no PCR primer to
#     align, since the integration site is already known exactly
#   - The literal string "NA" as a barcode/viewpoint label is renamed to
#     "none_unmapped" defensively (R treats the bare string "NA" specially
#     in many contexts) and that barcode is skipped entirely, since it has
#     no defined viewpoint coordinate to analyze against
#
# Usage:
#   Rscript run_basic4cseq.R <manifest_tsv> <bam_dir> <library_csv_path> <read_length> <output_dir>
#
# Arguments:
#   manifest_tsv      manifest.tsv produced by split_bam_by_barcode.py
#                      (columns: position, sanitized_filename)
#   bam_dir           directory containing the per-barcode BAMs
#                      (split_bam_by_barcode.py's output_dir)
#   library_csv_path  path to the virtual fragment library produced by
#                      build_fragment_library.R
#   read_length       the (post-trimming) read length used in your
#                      sequencing run -- MUST MATCH the read length used
#                      when building the library
#   output_dir        directory to write per-barcode results and plots into
#
# Requires: Basic4Cseq, GenomicAlignments
# (Does NOT require BSgenome.Mmusculus.UCSC.mm10 -- that is only needed by
# build_fragment_library.R, since this script reuses an already-built
# library file rather than cutting the genome itself.)

suppressMessages({
  library(Basic4Cseq)
  library(GenomicAlignments)
})

# --- Command-line arguments --------------------------------------------------
# PATCHED: the script originally hardcoded these five paths (leftover from
# ad hoc interactive use) and ignored any arguments passed to it, even though
# run_basic4cseq_pbs.sh already called it with 5 positional arguments as if
# it did. Replaced with real argv parsing so it can actually be driven by a
# PBS wrapper / the combined pipeline script.
#
# Usage:
#   Rscript run_basic4cseq_merged_bams.R <manifest_tsv> <bam_dir> <library_csv_path> <read_length> <output_dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(
    "Usage: Rscript run_basic4cseq_merged_bams.R <manifest_tsv> <bam_dir> <library_csv_path> <read_length> <output_dir>\n",
    "Got ", length(args), " argument(s): ", paste(args, collapse = " ")
  )
}

manifest_path <- path.expand(args[1])
bam_dir       <- path.expand(args[2])
library_path  <- path.expand(args[3])
read_length   <- as.integer(args[4])
output_dir    <- args[5]

# Width of the viewpoint interval to build around each barcode's single
# integration-site coordinate. We use a fixed window (see chat history)
# since we don't have an actual viewpoint FRAGMENT boundary, only a point
# coordinate; this is only used to define the viewpoint region for near-cis
# plotting/exclusion, not for fragment assignment itself (which uses the
# real restriction map already baked into the fragment library).
VIEWPOINT_HALF_WIDTH <- 500

# --- Sanity checks on inputs -------------------------------------------------

if (!file.exists(manifest_path)) {
  stop("Manifest file not found: ", manifest_path)
}
if (!dir.exists(bam_dir)) {
  stop("BAM directory not found: ", bam_dir)
}
if (!file.exists(library_path)) {
  stop(
    "Fragment library not found: ", library_path, "\n",
    "Build it first with: Rscript build_fragment_library.R <read_length> ", library_path
  )
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load manifest -----------------------------------------------------------

manifest <- read.delim(manifest_path, header = TRUE, stringsAsFactors = FALSE, na.strings = character(0))
if (!all(c("position", "sanitized_filename") %in% colnames(manifest))) {
  stop(
    "Manifest must have columns 'position' and 'sanitized_filename', got: ",
    paste(colnames(manifest), collapse = ", ")
  )
}
cat(sprintf("Loaded manifest with %d barcode/viewpoint entries\n", nrow(manifest)))

# --- Defend against the literal string "NA" ---------------------------------
na_label_mask <- manifest$position == "NA" & !is.na(manifest$position)
if (any(na_label_mask)) {
  cat(
    "NOTE: renaming barcode/viewpoint label literal \"NA\" to ",
    "\"none_unmapped\" to avoid R's NA-sentinel handling issues.\n",
    sep = ""
  )
  manifest$position[na_label_mask] <- "none_unmapped"
}

# --- Parse viewpoint coordinates from position labels ------------------------
#
# assign_barcodes.py's whitelist format produces labels like:
#   "chr18:77596512(+)"   -> chrom=chr18, coord=77596512
#   "none_unmapped"        -> no known coordinate, skip this barcode entirely

parse_position_label <- function(label) {
  match_result <- regmatches(label, regexec("^(chr[A-Za-z0-9]+):([0-9]+)\\(([+-])\\)$", label))
  if (length(match_result[[1]]) == 4) {
    list(chr = match_result[[1]][2], pos = as.integer(match_result[[1]][3]))
  } else {
    list(chr = NA_character_, pos = NA_integer_)
  }
}

# --- Fragment library (already built by build_fragment_library.R) ----------

cat("Using pre-built virtual fragment library:", library_path, "\n")

# --- Per-barcode loop ---------------------------------------------------------

results_summary <- data.frame(
  viewpoint = character(0),
  bamFile = character(0),
  status = character(0),
  nReads = integer(0),
  nNearCisFragments = integer(0),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(manifest))) {
  viewpoint_label <- manifest$position[i]
  bam_filename    <- manifest$sanitized_filename[i]
  bam_path        <- file.path(bam_dir, bam_filename)

  cat("\n==============================\n")
  cat("Processing barcode/viewpoint:", viewpoint_label, "\n")
  cat("==============================\n")

  if (viewpoint_label == "none_unmapped") {
    cat("Skipping: no known integration site coordinate for this barcode.\n")
    results_summary <- rbind(results_summary, data.frame(
      viewpoint = viewpoint_label, bamFile = bam_filename,
      status = "skipped_no_coordinate", nReads = NA_integer_, nNearCisFragments = NA_integer_
    ))
    next
  }

  if (!file.exists(bam_path)) {
    warning("BAM file not found, skipping: ", bam_path)
    results_summary <- rbind(results_summary, data.frame(
      viewpoint = viewpoint_label, bamFile = bam_filename,
      status = "skipped_missing_bam", nReads = NA_integer_, nNearCisFragments = NA_integer_
    ))
    next
  }

  parsed <- parse_position_label(viewpoint_label)
  if (is.na(parsed$chr)) {
    warning("Could not parse viewpoint coordinate from label, skipping: ", viewpoint_label)
    results_summary <- rbind(results_summary, data.frame(
      viewpoint = viewpoint_label, bamFile = bam_filename,
      status = "skipped_unparseable_label", nReads = NA_integer_, nNearCisFragments = NA_integer_
    ))
    next
  }

  vp_chr <- parsed$chr
  vp_pos <- parsed$pos
  vp_interval <- c(vp_pos - VIEWPOINT_HALF_WIDTH, vp_pos + VIEWPOINT_HALF_WIDTH)

  # --- Read alignments for this barcode ---
  cat("Reading alignments from:", bam_path, "\n")
  reads <- tryCatch(
    readGAlignments(bam_path),
    error = function(e) {
      warning("Failed to read BAM for ", viewpoint_label, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(reads)) {
    results_summary <- rbind(results_summary, data.frame(
      viewpoint = viewpoint_label, bamFile = bam_filename,
      status = "failed_reading_bam", nReads = NA_integer_, nNearCisFragments = NA_integer_
    ))
    next
  }
  n_reads <- length(reads)
  cat("  Read", n_reads, "alignments.\n")

  # --- Build the Data4Cseq object for this viewpoint ---
  #
  # NOTE on pointsOfInterest: every confirmed example we could verify in
  # Basic4Cseq's documentation passes a POPULATED points-of-interest
  # data.frame (loaded from a BED file via readPointsOfInterestFile). We
  # could not find a confirmed example of omitting it or passing an empty
  # one, so this is the one part of this script that is NOT verified
  # against a real documented example -- it is our best inference from the
  # documented column format (chr, start, end, name, colour). If this
  # specific call fails on your system, the likely fix is to supply a
  # real one-row points-of-interest data.frame marking the viewpoint
  # itself (chr=vp_chr, start=vp_interval[1], end=vp_interval[2],
  # name="VP", colour="black"), matching the documented vignette example.
  empty_poi <- data.frame(
    chr = character(0), start = integer(0), end = integer(0),
    name = character(0), colour = character(0)
  )
  vp_data <- tryCatch(
    Data4Cseq(
      viewpointChromosome = vp_chr,
      viewpointInterval = vp_interval,
      readLength = read_length,
      pointsOfInterest = empty_poi,
      rawReads = reads
    ),
    error = function(e) {
      cat("  Empty pointsOfInterest data.frame was rejected (", conditionMessage(e),
          "), retrying with the viewpoint itself as the single point of interest ...\n", sep = "")
      fallback_poi <- data.frame(
        chr = vp_chr, start = vp_interval[1], end = vp_interval[2],
        name = "VP", colour = "black", stringsAsFactors = FALSE
      )
      Data4Cseq(
        viewpointChromosome = vp_chr,
        viewpointInterval = vp_interval,
        readLength = read_length,
        pointsOfInterest = fallback_poi,
        rawReads = reads
      )
    }
  )

  # --- Map reads to the (shared) virtual fragment library ---
  rawFragments(vp_data) <- readsToFragments(vp_data, library_path)

  # --- Choose near-cis fragments around the viewpoint ----------------------
  #
  # This step is REQUIRED before normalizeFragmentData -- it selects the
  # subset of fragments in the near-cis region around the viewpoint and
  # removes fragments immediately adjacent to it (self-ligation bias).
  # Without this step, normalizeFragmentData receives no valid input and
  # returns empty results. Default region is 1Mb upstream and downstream
  # of the viewpoint position.
  near_cis_start <- max(1, vp_pos - 1e6)
  near_cis_end   <- vp_pos + 1e6
  cat("  Choosing near-cis fragments in",
      vp_chr, ":", near_cis_start, "-", near_cis_end, "\n")
  nearCisFragments(vp_data) <- chooseNearCisFragments(
    vp_data,
    regionCoordinates = c(near_cis_start, near_cis_end)
  )
  n_chosen <- nrow(nearCisFragments(vp_data))
  cat("  Near-cis fragments chosen:", n_chosen, "\n")

  if (n_chosen == 0) {
    warning("No near-cis fragments found for ", viewpoint_label,
            " -- check that vp_chr/vp_interval overlap the fragment library, ",
            "and that the library chromosome names match the BAM.")
    results_summary <- rbind(results_summary, data.frame(
      viewpoint = viewpoint_label, bamFile = bam_filename,
      status = "failed_no_near_cis_fragments", nReads = n_reads,
      nNearCisFragments = 0L
    ))
    next
  }

  # --- Normalize near-cis fragment data ---
  nearCisFragments(vp_data) <- normalizeFragmentData(vp_data)
  n_near_cis <- nrow(nearCisFragments(vp_data))
  cat("  Near-cis fragments after normalization:", n_near_cis,
      "(started with", n_chosen, "chosen fragments)\n")

  # --- Export results ---
  barcode_out_dir <- file.path(output_dir, viewpoint_label)
  dir.create(barcode_out_dir, showWarnings = FALSE, recursive = TRUE)

  wig_path <- file.path(barcode_out_dir, paste0(viewpoint_label, ".wig"))
  cat("  Writing wiggle track to:", wig_path, "\n")
  printWigFile(vp_data, wigFileName = wig_path)

  fragment_table_path <- file.path(barcode_out_dir, paste0(viewpoint_label, "_near_cis_fragments.tsv"))
  write.table(
    nearCisFragments(vp_data),
    file = fragment_table_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  cat("  Wrote near-cis fragment table to:", fragment_table_path, "\n")

  coverage_plot_path <- file.path(barcode_out_dir, paste0(viewpoint_label, "_coverage.pdf"))
  tryCatch({
   visualizeViewpoint(vp_data, plotFileName = coverage_plot_path)
    cat("  Wrote coverage plot to:", coverage_plot_path, "\n")
  }, error = function(e) {
    warning("Coverage plot failed for ", viewpoint_label, ": ", conditionMessage(e))
  })

  results_summary <- rbind(results_summary, data.frame(
    viewpoint = viewpoint_label, bamFile = bam_filename,
    status = "ok", nReads = n_reads, nNearCisFragments = n_near_cis
  ))
}

# --- Write overall summary ---------------------------------------------------

summary_path <- file.path(output_dir, "run_summary.tsv")
write.table(results_summary, file = summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== Done ===\n")
cat("Summary written to:", summary_path, "\n")
print(results_summary)
