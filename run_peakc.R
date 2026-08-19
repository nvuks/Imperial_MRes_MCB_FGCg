#!/usr/bin/env Rscript
#
# run_peakc.R
#
# Runs peakC on the wiggle tracks produced by run_basic4cseq_merged_bams.R,
# one viewpoint at a time (single-replicate analysis).
#
# For each viewpoint in the manifest, this script:
#
#   1. Reads the Basic4Cseq wig for that viewpoint
#      (<basic4cseq_dir>/<viewpoint>/<viewpoint>.wig), which contains
#      variableStep blocks for ALL chromosomes.
#   2. Filters it down to only the block for the viewpoint's OWN chromosome,
#      writing <viewpoint>_<chr>.wig alongside the original. peakC's
#      readqWig() expects a single-chromosome track -- feeding it a
#      multi-chromosome wig would mix positions from different chromosomes
#      into one coordinate axis.
#   3. Runs readqWig() + single.analysis() on the filtered wig.
#   4. Writes plot_C() output to a PDF, plus the called peaks as a TSV.
#
# Equivalent, for one viewpoint, to:
#
#   library(peakC)
#   data <- readqWig("<...>/chr4:101390053(+)_chr4.wig",
#                    vp.pos = 101390053, window = 700e3)
#   res  <- single.analysis(data$data, vp.pos = 101390053, qWd = 2.5)
#   plot_C(res)
#
# Usage:
#   Rscript run_peakc.R <manifest_tsv> <basic4cseq_dir> <output_dir> [window] [qWd]
#
# Arguments:
#   manifest_tsv     manifest.tsv from split_bam_by_barcode.py (columns:
#                    position, sanitized_filename). Only the `position`
#                    column (the viewpoint label) is used.
#   basic4cseq_dir   the output_dir given to run_basic4cseq_merged_bams.R,
#                    i.e. the directory containing one subdirectory per
#                    viewpoint, each holding <viewpoint>.wig
#   output_dir       directory to write peakC plots, peak tables and the
#                    run summary into
#   window           optional; half-width (bp) around the viewpoint passed to
#                    readqWig(). Default 700e3, as in the worked example.
#   qWd              optional; quantile width parameter passed to
#                    single.analysis(). Default 2.5, as in the worked example.
#
# Notes:
#   - Viewpoints whose label can't be parsed as chr:pos(strand) (e.g. the
#     "none_unmapped" bucket) are skipped, as in run_basic4cseq_merged_bams.R.
#   - A viewpoint is also skipped, with a warning rather than an error, if its
#     wig is missing (Basic4Cseq may have skipped it) or if the wig contains no
#     data on the viewpoint's own chromosome. One bad viewpoint should not
#     abort the whole run.
#   - Single-replicate analysis only (single.analysis). If replicates are added
#     later, peakC's combined.analysis() would be the entry point instead.

suppressPackageStartupMessages(library(peakC))

# --- Command-line arguments --------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || length(args) > 5) {
  stop(
    "Usage: Rscript run_peakc.R <manifest_tsv> <basic4cseq_dir> <output_dir> [window] [qWd]\n",
    "Got ", length(args), " argument(s): ", paste(args, collapse = " ")
  )
}

manifest_path  <- path.expand(args[1])
basic4cseq_dir <- path.expand(args[2])
output_dir     <- args[3]
window         <- if (length(args) >= 4) as.numeric(args[4]) else 700e3
qWd            <- if (length(args) >= 5) as.numeric(args[5]) else 2.5

if (is.na(window) || window <= 0) stop("window must be a positive number, got: ", args[4])
if (is.na(qWd)    || qWd    <= 0) stop("qWd must be a positive number, got: ",    args[5])

if (!file.exists(manifest_path))  stop("Manifest not found: ", manifest_path)
if (!dir.exists(basic4cseq_dir))  stop("Basic4Cseq output directory not found: ", basic4cseq_dir)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Manifest:        ", manifest_path,  "\n")
cat("Basic4Cseq dir:  ", basic4cseq_dir, "\n")
cat("Output dir:      ", output_dir,     "\n")
cat("window:          ", format(window, scientific = FALSE), "\n")
cat("qWd:             ", qWd, "\n")

# --- Helpers ------------------------------------------------------------------

# Same label convention as run_basic4cseq_merged_bams.R: "chr4:101390053(+)"
parse_position_label <- function(label) {
  match_result <- regmatches(label, regexec("^(chr[A-Za-z0-9]+):([0-9]+)\\(([+-])\\)$", label))
  if (length(match_result[[1]]) == 4) {
    list(chr = match_result[[1]][2], pos = as.integer(match_result[[1]][3]))
  } else {
    NULL
  }
}

# Basic4Cseq's printWigFile() emits standard wiggle: a "variableStep chrom=<chr>"
# declaration line, then position/value rows, repeating per chromosome. To get a
# single-chromosome track we keep the declaration line for the target chromosome
# plus its data rows, and drop everything belonging to other chromosomes. Any
# leading "track" line is preserved so the output stays a valid wig.
filter_wig_to_chromosome <- function(wig_in, wig_out, target_chr) {
  lines <- readLines(wig_in, warn = FALSE)

  keep          <- logical(length(lines))
  current_chr   <- NA_character_
  n_data_kept   <- 0L

  for (i in seq_along(lines)) {
    line <- lines[i]

    if (grepl("^track", line)) {
      keep[i] <- TRUE
      next
    }

    # A step declaration switches which chromosome subsequent rows belong to.
    if (grepl("^(variableStep|fixedStep)", line)) {
      chrom_match <- regmatches(line, regexec("chrom=([^[:space:]]+)", line))
      current_chr <- if (length(chrom_match[[1]]) == 2) chrom_match[[1]][2] else NA_character_
      keep[i] <- !is.na(current_chr) && current_chr == target_chr
      next
    }

    if (!nzchar(trimws(line))) next  # drop blank lines

    # Data row: keep only if we're inside the target chromosome's block.
    if (!is.na(current_chr) && current_chr == target_chr) {
      keep[i]     <- TRUE
      n_data_kept <- n_data_kept + 1L
    }
  }

  writeLines(lines[keep], wig_out)
  n_data_kept
}

# --- Read manifest ------------------------------------------------------------

manifest <- read.delim(manifest_path, header = TRUE, stringsAsFactors = FALSE)
if (!"position" %in% colnames(manifest)) {
  stop("Manifest is missing the required 'position' column: ", manifest_path)
}
cat("Viewpoints in manifest:", nrow(manifest), "\n")

# --- Per-viewpoint loop -------------------------------------------------------

results_summary <- data.frame(
  viewpoint    = character(0),
  chromosome   = character(0),
  vp_pos       = integer(0),
  status       = character(0),
  nWigDataRows = integer(0),
  nPeaks       = integer(0),
  stringsAsFactors = FALSE
)

add_result <- function(df, viewpoint, chr, pos, status, n_wig = NA_integer_, n_peaks = NA_integer_) {
  rbind(df, data.frame(
    viewpoint = viewpoint, chromosome = chr, vp_pos = pos, status = status,
    nWigDataRows = n_wig, nPeaks = n_peaks, stringsAsFactors = FALSE
  ))
}

for (i in seq_len(nrow(manifest))) {
  viewpoint_label <- manifest$position[i]

  cat("\n==============================\n")
  cat("peakC for viewpoint:", viewpoint_label, "\n")
  cat("==============================\n")

  parsed <- parse_position_label(viewpoint_label)
  if (is.null(parsed)) {
    # Covers "none_unmapped" and any other non-coordinate label.
    cat("Skipping: label is not a parseable chr:pos(strand) coordinate.\n")
    results_summary <- add_result(results_summary, viewpoint_label,
                                  NA_character_, NA_integer_, "skipped_no_coordinate")
    next
  }

  vp_chr <- parsed$chr
  vp_pos <- parsed$pos
  cat("  Chromosome:", vp_chr, " Position:", vp_pos, "\n")

  wig_path <- file.path(basic4cseq_dir, viewpoint_label,
                        paste0(viewpoint_label, ".wig"))

  if (!file.exists(wig_path)) {
    warning("Wig not found, skipping: ", wig_path)
    results_summary <- add_result(results_summary, viewpoint_label, vp_chr, vp_pos,
                                  "skipped_wig_not_found")
    next
  }

  # --- 1. Filter the wig to the viewpoint's own chromosome --------------------
  filtered_wig_path <- file.path(basic4cseq_dir, viewpoint_label,
                                 paste0(viewpoint_label, "_", vp_chr, ".wig"))

  n_data_kept <- tryCatch(
    filter_wig_to_chromosome(wig_path, filtered_wig_path, vp_chr),
    error = function(e) {
      warning("Failed to filter wig for ", viewpoint_label, ": ", conditionMessage(e))
      NA_integer_
    }
  )

  if (is.na(n_data_kept)) {
    results_summary <- add_result(results_summary, viewpoint_label, vp_chr, vp_pos,
                                  "error_filtering_wig")
    next
  }

  cat("  Filtered wig ->", filtered_wig_path, "(", n_data_kept, "data rows on", vp_chr, ")\n")

  if (n_data_kept == 0L) {
    warning("No wig data on the viewpoint's own chromosome (", vp_chr, ") for ",
            viewpoint_label, " -- skipping. This is worth investigating: it means ",
            "no cis signal at all was recorded for this viewpoint.")
    results_summary <- add_result(results_summary, viewpoint_label, vp_chr, vp_pos,
                                  "skipped_no_cis_signal", n_wig = 0L)
    next
  }

  # --- 2. Run peakC ------------------------------------------------------------
  peakc_result <- tryCatch({
    data_vp <- readqWig(filtered_wig_path, vp.pos = vp_pos, window = window)
    res_vp  <- single.analysis(data_vp$data, vp.pos = vp_pos, qWd = qWd)
    res_vp
  }, error = function(e) {
    warning("peakC failed for ", viewpoint_label, ": ", conditionMessage(e))
    NULL
  })

  if (is.null(peakc_result)) {
    results_summary <- add_result(results_summary, viewpoint_label, vp_chr, vp_pos,
                                  "error_peakc_failed", n_wig = n_data_kept)
    next
  }

  # --- 3. Plot to PDF ----------------------------------------------------------
  # Filenames use the sanitized label: ':', '(' , ')' and '+' are awkward in
  # shell/filesystem contexts, even though Basic4Cseq itself uses them as-is
  # for its directory names.
  safe_label <- gsub("[^A-Za-z0-9._-]", "_", viewpoint_label)
  plot_path  <- file.path(output_dir, paste0(safe_label, "_peakC.pdf"))

  tryCatch({
    pdf(plot_path, width = 10, height = 6)
    on.exit(if (!is.null(dev.list())) dev.off(), add = TRUE)
    plot_C(peakc_result)
    title(main = paste0("peakC: ", viewpoint_label,
                        "  (window=", format(window, scientific = FALSE),
                        ", qWd=", qWd, ")"))
    dev.off()
    cat("  Wrote peakC plot to:", plot_path, "\n")
  }, error = function(e) {
    if (!is.null(dev.list())) dev.off()
    warning("Failed to write peakC plot for ", viewpoint_label, ": ", conditionMessage(e))
  })

  # --- 4. Write the called peaks -----------------------------------------------
  # single.analysis() returns a list; $peak holds the significant positions.
  n_peaks <- NA_integer_
  tryCatch({
    peaks <- peakc_result$peak
    if (!is.null(peaks) && length(peaks) > 0) {
      peak_df <- data.frame(
        viewpoint  = viewpoint_label,
        chromosome = vp_chr,
        vp_pos     = vp_pos,
        peak_pos   = as.integer(peaks),
        stringsAsFactors = FALSE
      )
      peak_df$distance_from_vp <- peak_df$peak_pos - vp_pos
      peaks_path <- file.path(output_dir, paste0(safe_label, "_peaks.tsv"))
      write.table(peak_df, peaks_path, sep = "\t", quote = FALSE, row.names = FALSE)
      n_peaks <- nrow(peak_df)
      cat("  Wrote", n_peaks, "peak(s) to:", peaks_path, "\n")
    } else {
      n_peaks <- 0L
      cat("  No significant peaks called for this viewpoint.\n")
    }
  }, error = function(e) {
    warning("Failed to extract peaks for ", viewpoint_label, ": ", conditionMessage(e))
  })

  results_summary <- add_result(results_summary, viewpoint_label, vp_chr, vp_pos,
                                "ok", n_wig = n_data_kept, n_peaks = n_peaks)
}

# --- Summary ------------------------------------------------------------------

summary_path <- file.path(output_dir, "peakc_run_summary.tsv")
write.table(results_summary, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n==============================\n")
cat("peakC run complete.\n")
cat("Summary written to:", summary_path, "\n")
print(results_summary)

n_ok <- sum(results_summary$status == "ok")
cat("\nViewpoints analysed successfully:", n_ok, "/", nrow(results_summary), "\n")
if (n_ok == 0 && nrow(results_summary) > 0) {
  warning("peakC produced no successful analyses for any viewpoint.")
}
