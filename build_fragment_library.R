#!/usr/bin/env Rscript
#
# build_fragment_library.R
#
# Builds the Basic4Cseq virtual restriction-fragment library for mm10 with
# the NlaIII (primary) / DpnII (secondary) enzyme pair, for a given read
# length. This only needs to be run ONCE per (genome, enzyme pair, read
# length) combination -- the resulting file is reused by run_basic4cseq.R
# for every barcode/viewpoint in the experiment.
#
# Re-run this script (producing a new output file, or overwriting the old
# one deliberately) if you ever change read length, genome build, or
# restriction enzymes.
#
# Usage:
#   Rscript build_fragment_library.R <read_length> <output_csv_path>
#
# Requires: Basic4Cseq, BSgenome.Mmusculus.UCSC.mm10
# (BSgenome.Mmusculus.UCSC.mm10 specifically -- UCSC-style "chr1", "chr2"...
# contig naming, matching your BAM's contig naming. See run_basic4cseq.R
# header comments for why this matters.)

suppressMessages({
  library(Basic4Cseq)
  library(BSgenome.Mmusculus.UCSC.mm10)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript build_fragment_library.R <read_length> <output_csv_path>")
}

read_length <- as.numeric(args[1])
output_path <- args[2]

FIRST_CUTTER  <- "CATG"  # NlaIII, primary cutter
SECOND_CUTTER <- "GATC"  # DpnII, secondary cutter

if (is.na(read_length) || read_length <= 0) {
  stop("read_length must be a positive number, got: ", args[1])
}

if (file.exists(output_path)) {
  stop(
    "Output file already exists: ", output_path, "\n",
    "Refusing to silently overwrite an existing fragment library -- ",
    "delete it first or choose a different output path if you intend ",
    "to rebuild it (e.g. after changing read length, genome, or enzymes)."
  )
}

dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)

cat(
  "Building virtual fragment library:\n",
  "  Genome:        mm10 (BSgenome.Mmusculus.UCSC.mm10)\n",
  "  First cutter:  ", FIRST_CUTTER, " (NlaIII)\n",
  "  Second cutter: ", SECOND_CUTTER, " (DpnII)\n",
  "  Read length:   ", read_length, "\n",
  "  Output:        ", output_path, "\n",
  "This can take a while and use significant memory for a full genome ...\n",
  sep = ""
)

createVirtualFragmentLibrary(
  chosenGenome = Mmusculus,
  firstCutter = FIRST_CUTTER,
  secondCutter = SECOND_CUTTER,
  readLength = read_length,
  libraryName = output_path
)

cat("Done. Fragment library written to:", output_path, "\n")
