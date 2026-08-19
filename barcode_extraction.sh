#!/usr/bin/env bash
#
# run_barcode_pipeline.sh
#
# Full pipeline: extract 16bp barcode after a constant adapter from a
# gzipped FASTQ file, then assign each read's barcode to its nearest
# whitelist barcode (within 1 mismatch by default), attaching integration
# site coordinates.
#
# Usage:
#   ./run_barcode_pipeline.sh <reads.fastq.gz> <whitelist.tsv> [output_prefix] [max_mismatches]
#
# Arguments:
#   reads.fastq.gz   - gzipped FASTQ file of reads
#   whitelist.tsv    - tab-separated: label  barcode  start  end
#   output_prefix    - optional. Default: "barcode_pipeline_output"
#   max_mismatches   - optional. Default: 1
#
# Outputs (using default prefix):
#   barcode_pipeline_output.barcodes_with_ids.tsv   ReadID + extracted barcode (all reads w/ adapter match)
#   barcode_pipeline_output.reads_with_sites.tsv    ReadID + barcode + assigned site (final result)
#
# Requires: zcat, awk, python3 (with assign_barcodes.py in the same directory
# as this script, or on PATH)

set -euo pipefail

ADAPTER="GGCCGGCCACAACTCGAG"
BARCODE_LENGTH=16

# --- Argument parsing -------------------------------------------------------

if [[ $# -lt 2 || $# -gt 4 ]]; then
    echo "Usage: $0 <reads.fastq.gz> <whitelist.tsv> [output_prefix] [max_mismatches]" >&2
    exit 1
fi

FASTQ_GZ="$1"
WHITELIST="$2"
OUTPUT_PREFIX="${3:-barcode_pipeline_output}"
MAX_MISMATCHES="${4:-1}"

if [[ ! -f "$FASTQ_GZ" ]]; then
    echo "Error: FASTQ file not found: $FASTQ_GZ" >&2
    exit 1
fi

if [[ ! -f "$WHITELIST" ]]; then
    echo "Error: whitelist file not found: $WHITELIST" >&2
    exit 1
fi

# Locate assign_barcodes.py: prefer the same directory as this script,
# fall back to PATH.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ASSIGN_SCRIPT="${SCRIPT_DIR}/assign_barcodes.py"
if [[ ! -f "$ASSIGN_SCRIPT" ]]; then
    if command -v assign_barcodes.py >/dev/null 2>&1; then
        ASSIGN_SCRIPT="assign_barcodes.py"
    else
        echo "Error: could not find assign_barcodes.py (looked in ${SCRIPT_DIR} and PATH)" >&2
        exit 1
    fi
fi

BARCODES_FILE="${OUTPUT_PREFIX}.barcodes_with_ids.tsv"
FINAL_FILE="${OUTPUT_PREFIX}.reads_with_sites.tsv"

echo "[1/2] Extracting ${BARCODE_LENGTH}bp barcodes after adapter from ${FASTQ_GZ} ..." >&2

zcat "$FASTQ_GZ" | awk -v adapter="$ADAPTER" -v bclen="$BARCODE_LENGTH" '
NR % 4 == 1 { id = $0 }
NR % 4 == 2 {
    seq = $0
    pos = index(seq, adapter)
    if (pos > 0) {
        start = pos + length(adapter)
        frag = substr(seq, start, bclen)
        if (length(frag) == bclen) {
            print id "\t" frag
        }
    }
}
' > "$BARCODES_FILE"

N_EXTRACTED=$(wc -l < "$BARCODES_FILE")
echo "      Extracted barcodes for ${N_EXTRACTED} reads -> ${BARCODES_FILE}" >&2

echo "[2/2] Assigning barcodes to whitelist (max ${MAX_MISMATCHES} mismatch(es)) ..." >&2

python3 "$ASSIGN_SCRIPT" --max-mismatches "$MAX_MISMATCHES" "$WHITELIST" "$BARCODES_FILE" > "$FINAL_FILE"

echo "Done." >&2
echo "  Per-read barcodes:        $BARCODES_FILE" >&2
echo "  Final assigned output:    $FINAL_FILE" >&2
