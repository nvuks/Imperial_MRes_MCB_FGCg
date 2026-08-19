#!/bin/bash
#PBS -N build_fragment_library
#PBS -l select=1:ncpus=1:mem=96gb
#PBS -l walltime=04:00:00
#PBS -j oe
#PBS -o build_fragment_library.log

#
# build_fragment_library.pbs.sh
#
# PBS submission script: builds the Basic4Cseq virtual restriction-fragment
# library for mm10 (NlaIII / DpnII) at a given read length. This is a
# ONE-TIME step per (genome, enzyme pair, read length) combination --
# rerun only if any of those change.
#
# EDIT THE VARIABLES BELOW (or override at submission time), then submit:
#   qsub -v READ_LENGTH=50,LIBRARY_PATH="$HOME/analysis/TRIP-4C_pilot/fragment_library/mm10_NlaIII_DpnII_50bp.csv" build_fragment_library.pbs.sh
#
# Requires build_fragment_library.R to be in the same directory as this
# script, or on PATH.

set -euo pipefail

# --- User-configurable variables -------------------------------------------
# These can be overridden at submission time with:
#   qsub -v READ_LENGTH=...,LIBRARY_PATH=... build_fragment_library.pbs.sh
# The values below are only used as fallbacks if a variable isn't passed in.

READ_LENGTH="${READ_LENGTH:=50}"
LIBRARY_PATH="${LIBRARY_PATH:=$HOME/analysis/TRIP-4C_pilot/fragment_library/mm10_NlaIII_DpnII_${READ_LENGTH}bp.csv}"


# --- Environment setup ------------------------------------------------------

eval "$(~/anaconda3/bin/conda shell.bash hook)"
conda activate chicago
cd "$PBS_O_WORKDIR"

echo "Job started: $(date)"
echo "Running on host: $(hostname)"
echo "R version: $(Rscript --version 2>&1)"
echo "Read length: $READ_LENGTH"
echo "Library output path: $LIBRARY_PATH"

# --- Locate build_fragment_library.R ----------------------------------------

#SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#BUILD_SCRIPT="${SCRIPT_DIR}/build_fragment_library.R"
#if command -v build_fragment_library.R >/dev/null 2>&1; then
#        BUILD_SCRIPT="build_fragment_library.R"
#else
#        echo "Error: could not find build_fragment_library.R" >&2
#        exit 1
#fi

# --- Run -------------------------------------------------------------------

echo "Building virtual fragment library ..."
echo "Current working directory: $(pwd)"
Rscript build_fragment_library.R "$READ_LENGTH" "$LIBRARY_PATH"

echo "Done: $(date)"
echo "  Library written to: $LIBRARY_PATH"
