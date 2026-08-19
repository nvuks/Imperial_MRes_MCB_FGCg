#!/bin/bash
#PBS -N trip4c_full_pipeline
#PBS -l select=1:ncpus=8:mem=96gb
#PBS -l walltime=24:00:00
#PBS -j oe
#
# run_full_pipeline_pbs.sh
#
# Single-job PBS submission script running the full TRIP-4C SOX2 pilot
# pipeline end to end, from raw paired-end FASTQ to per-viewpoint
# Basic4Cseq results:
#
#   1. Barcode extraction + assignment (from R2)      -> barcode_extraction.sh (+ assign_barcodes.py)
#   2. Read trimming (longest non-insert CATG/GATC)   -> trim_longest_reference_catg_gatc_stretch.py
#   3. Alignment (bowtie2, single-end, R1 and R2)      -> align_pbs.sh logic (inlined below)
#   4. Split aligned BAMs by barcode/viewpoint         -> split_bam_by_barcode.py (once per read)
#   5. Merge R1+R2 per-barcode BAMs                    -> merge_paired_bams.py
#   6. Build virtual restriction-fragment library      -> build_fragment_library.R (skipped if already built)
#   7. Run Basic4Cseq once per barcode/viewpoint, on the merged BAMs -> run_basic4cseq_merged_bams.R
#
# Diagnostic-only steps from the README (unaligned-vs-aligned read
# comparison, integration-site density peak checks, exploratory barcode
# clustering) are intentionally NOT included here -- see the README,
# section 1.6, for how to run them separately if needed.
#
# NOTE on run_basic4cseq_merged_bams.R and run_basic4cseq.R: as supplied,
# both hardcoded their five input paths internally and ignored any
# command-line arguments -- even though run_basic4cseq_pbs.sh already
# called run_basic4cseq.R with 5 positional arguments as if it read them.
# Both scripts have been patched (commandArgs() parsing added) so they can
# actually be driven by this pipeline / by run_basic4cseq_pbs.sh. Use the
# patched copies alongside this script, not the originals as uploaded.
#
# REQUIRED FILES (same directory as this script, or set SCRIPT_DIR below):
#   assign_barcodes.py
#   whitelist.tsv
#   trim_longest_reference_catg_gatc_stretch.py
#   split_bam_by_barcode.py
#   merge_paired_bams.py
#   build_fragment_library.R
#   run_basic4cseq_merged_bams.R   (patched -- see note above)
#
# Submit with (edit paths as needed, or override at submission time):
#   qsub -v R1_FASTQ_GZ="/path/4C_SOX2_R1_S1_L001_R1_001.fastq.gz",\
#R2_FASTQ_GZ="/path/4C_SOX2_R1_S1_L001_R2_001.fastq.gz",\
#WHITELIST="/path/whitelist.tsv",\
#BOWTIE2_INDEX="/rds/general/project/lms-spivakov-raw/live/bowtie_index/mm10",\
#OUTPUT_DIR="/path/to/results" \
#run_full_pipeline_pbs.sh

eval "$(~/anaconda3/bin/conda shell.bash hook)"
conda activate chicago
module load SAMtools/1.16.1-GCC-11.3.0
module load Bowtie2/2.4.5-GCC-11.3.0

python3 -m pip install pysam
python3 -m pip install edlib

set -euo pipefail

# --- User-configurable variables --------------------------------------------
# Override at submission time with -v; values below are fallbacks only.

SCRIPT_DIR="${SCRIPT_DIR:=$PBS_O_WORKDIR}"

R1_FASTQ_GZ="${R1_FASTQ_GZ:?ERROR: R1_FASTQ_GZ must be set via -v}"
R2_FASTQ_GZ="${R2_FASTQ_GZ:?ERROR: R2_FASTQ_GZ must be set via -v}"
WHITELIST="${WHITELIST:=$SCRIPT_DIR/whitelist.tsv}"
BOWTIE2_INDEX="${BOWTIE2_INDEX:=/rds/general/project/lms-spivakov-raw/live/bowtie_index/mm10}"

THREADS="${THREADS:=8}"
MAX_BARCODE_MISMATCHES="${MAX_BARCODE_MISMATCHES:=1}"
READ_LENGTH="${READ_LENGTH:=50}"

# peakC parameters (Step 8). Defaults match the worked example.
PEAKC_WINDOW="${PEAKC_WINDOW:=700e3}"   # half-width around the viewpoint, bp
PEAKC_QWD="${PEAKC_QWD:=2.5}"           # quantile width for single.analysis()

OUTPUT_DIR="${OUTPUT_DIR:=$PBS_O_WORKDIR/trip4c_pipeline_output_${PBS_JOBID}}"
FRAGMENT_LIBRARY_DIR="${FRAGMENT_LIBRARY_DIR:=/rds/general/project/lms-spivakov-analysis/live/TRIP-4C_pilot/fragment_library}"
FRAG_LIB="${FRAG_LIB:=${FRAGMENT_LIBRARY_DIR}/mm10_NlaIII_DpnII_${READ_LENGTH}bp.csv}"

mkdir -p "$OUTPUT_DIR"

# --- Environment setup -------------------------------------------------------

cd "$PBS_O_WORKDIR" 2>/dev/null || true

echo "Job started: $(date)"
echo "Running on host: $(hostname)"
echo "Script dir:      $SCRIPT_DIR"
echo "R1 FASTQ:         $R1_FASTQ_GZ"
echo "R2 FASTQ:         $R2_FASTQ_GZ"
echo "Whitelist:        $WHITELIST"
echo "Bowtie2 index:    $BOWTIE2_INDEX"
echo "Output dir:       $OUTPUT_DIR"

# --- Sanity checks ------------------------------------------------------------

for f in "$R1_FASTQ_GZ" "$R2_FASTQ_GZ" "$WHITELIST"; do
    if [[ ! -f "$f" ]]; then
        echo "Error: required input not found: $f" >&2
        exit 1
    fi
done

for f in assign_barcodes.py trim_longest_reference_catg_gatc_stretch.py \
         split_bam_by_barcode.py merge_paired_bams.py build_fragment_library.R \
         run_basic4cseq_merged_bams.R; do
    if [[ ! -f "${SCRIPT_DIR}/${f}" ]]; then
        echo "Error: required script not found: ${SCRIPT_DIR}/${f}" >&2
        exit 1
    fi
done

if ! python3 -c "import edlib" 2>/dev/null; then
    echo "Warning: edlib not found -- trimming will fall back to slower pure-Python matching." >&2
    echo "Install with: pip install edlib --break-system-packages" >&2
fi

if ! python3 -c "import pysam" 2>/dev/null; then
    echo "Error: pysam is required by split_bam_by_barcode.py but is not installed." >&2
    exit 1
fi

# ==============================================================================
# Step 1 -- Barcode extraction + assignment (from R2)
# ==============================================================================

echo ""
echo "=== Step 1/8: barcode extraction + assignment ==="

ADAPTER="GGCCGGCCACAACTCGAG"
BARCODE_LENGTH=16
BARCODE_PREFIX="${OUTPUT_DIR}/barcode_pipeline_output"
BARCODES_FILE="${BARCODE_PREFIX}.barcodes_with_ids.tsv"
READS_WITH_SITES="${BARCODE_PREFIX}.reads_with_sites.tsv"

echo "[1/2] Extracting ${BARCODE_LENGTH}bp barcodes after adapter from ${R2_FASTQ_GZ} ..."

zcat "$R2_FASTQ_GZ" | awk -v adapter="$ADAPTER" -v bclen="$BARCODE_LENGTH" '
NR % 4 == 1 { id = $0 }
NR % 4 == 2 {
    seq = $0
    pos = index(seq, adapter)
    if (pos > 0) {
        start = pos + length(adapter)
        frag = substr(seq, start, bclen)
        if (length(frag) == bclen) print id "\t" frag
    }
}' > "$BARCODES_FILE"

N_EXTRACTED=$(wc -l < "$BARCODES_FILE")
echo "      Extracted barcodes for ${N_EXTRACTED} reads -> ${BARCODES_FILE}"

echo "[2/2] Assigning barcodes to whitelist (max ${MAX_BARCODE_MISMATCHES} mismatch(es)) ..."

python3 "${SCRIPT_DIR}/assign_barcodes.py" \
    --max-mismatches "$MAX_BARCODE_MISMATCHES" \
    "$WHITELIST" "$BARCODES_FILE" > "$READS_WITH_SITES"

echo "  Barcode assignment written to: $READS_WITH_SITES"

# ==============================================================================
# Step 2 -- Trim reads to the longest non-insert restriction fragment
# ==============================================================================

echo ""
echo "=== Step 2/8: read trimming ==="

trim_reads () {
    local input_fastq_gz="$1"
    local anchor="$2"
    local base output_fastq_gz
    base="$(basename "$input_fastq_gz")"
    base="${base%.fastq.gz}"
    base="${base%.fq.gz}"
    output_fastq_gz="${OUTPUT_DIR}/${base}.longest_reference_catg_gatc_trimmed.fastq.gz"

    echo "Trimming $input_fastq_gz (anchor=$anchor) -> $output_fastq_gz"
    # shellcheck disable=SC2086
    python3 "${SCRIPT_DIR}/trim_longest_reference_catg_gatc_stretch.py" \
        --anchor "$anchor" \
        "$input_fastq_gz" "$output_fastq_gz"
    echo "$output_fastq_gz"
}

# R1: NlaIII/CATG-anchored
R1_TRIMMED=$(trim_reads "$R1_FASTQ_GZ" "CATG" | tail -1)

# R2: DpnII/GATC-anchored
R2_TRIMMED=$(trim_reads "$R2_FASTQ_GZ" "GATC" | tail -1)

echo "  R1 trimmed: $R1_TRIMMED"
echo "  R2 trimmed: $R2_TRIMMED"

# ==============================================================================
# Step 3 -- Alignment (bowtie2, single-end, R1 and R2 independently)
# ==============================================================================

echo ""
echo "=== Step 3/8: alignment ==="

#eval "$(~/anaconda3/bin/conda shell.bash hook)"
#conda activate chicago
#module load SAMtools/1.16.1-GCC-11.3.0
#module load Bowtie2/2.4.5-GCC-11.3.0

echo "bowtie2 version: $(bowtie2 --version | head -1)"
echo "samtools version: $(samtools --version | head -1)"

if [[ ! -f "${BOWTIE2_INDEX}.1.bt2" && ! -f "${BOWTIE2_INDEX}.1.bt2l" ]]; then
    echo "Error: bowtie2 index not found at prefix: $BOWTIE2_INDEX" >&2
    exit 1
fi

align_reads () {
    # Aligns a single-end FASTQ, splits into aligned-only / unaligned-only
    # BAMs. Equivalent to align_pbs.sh, inlined here for a single-job run.
    local input_fastq_gz="$1"
    local label="$2"
    local output_bam="${OUTPUT_DIR}/${label}_aligned.sorted.bam"
    local unaligned_bam="${OUTPUT_DIR}/${label}_unaligned.bam"
    local align_log="${OUTPUT_DIR}/${label}_align.log"
    local intermediate_bam="${OUTPUT_DIR}/${label}.all.sorted.bam"

    echo "Aligning $input_fastq_gz -> $output_bam"

    bowtie2 \
        --very-sensitive \
        -x "$BOWTIE2_INDEX" \
        --threads "$THREADS" \
        --reorder \
        -U "$input_fastq_gz" \
        2> "$align_log" \
      | samtools sort -@ "$THREADS" -o "$intermediate_bam" -

    echo "bowtie2 finished for $label. Summary:"
    cat "$align_log"

    if [[ ! -s "$intermediate_bam" ]]; then
        echo "Error: bowtie2 produced an empty or missing BAM: $intermediate_bam" >&2
        exit 1
    fi

    samtools index "$intermediate_bam"

    samtools view -@ "$THREADS" -b -F 4 -o "$output_bam" "$intermediate_bam"
    samtools index "$output_bam"

    samtools view -@ "$THREADS" -b -f 4 -o "$unaligned_bam" "$intermediate_bam"
    # Deliberately not indexed -- every read is unmapped, a coordinate
    # index provides no benefit (see align_pbs.sh header for details).

    local n_aligned n_unaligned
    n_aligned=$(samtools view -c "$output_bam")
    n_unaligned=$(samtools view -c "$unaligned_bam")
    echo "  Aligned reads:   $n_aligned"
    echo "  Unaligned reads: $n_unaligned"

    rm -f "$intermediate_bam" "${intermediate_bam}.bai"

    echo "$output_bam"
}

R1_ALIGNED_BAM=$(align_reads "$R1_TRIMMED" "4C_SOX2_R1" | tail -1)
R2_ALIGNED_BAM=$(align_reads "$R2_TRIMMED" "4C_SOX2_R2" | tail -1)

echo "  R1 aligned BAM: $R1_ALIGNED_BAM"
echo "  R2 aligned BAM: $R2_ALIGNED_BAM"

# ==============================================================================
# Step 4 -- Split aligned BAMs by barcode/viewpoint
# ==============================================================================

eval "$(~/anaconda3/bin/conda shell.bash hook)"
conda activate chicago

echo ""
echo "=== Step 4/8: split BAMs by barcode ==="

SPLIT_R1_DIR="${OUTPUT_DIR}/split_bams_R1"
SPLIT_R2_DIR="${OUTPUT_DIR}/split_bams_R2"

python3 "${SCRIPT_DIR}/split_bam_by_barcode.py" \
    "$R1_ALIGNED_BAM" "$READS_WITH_SITES" "$SPLIT_R1_DIR"

python3 "${SCRIPT_DIR}/split_bam_by_barcode.py" \
    "$R2_ALIGNED_BAM" "$READS_WITH_SITES" "$SPLIT_R2_DIR"

echo "  R1 per-barcode BAMs: $SPLIT_R1_DIR"
echo "  R2 per-barcode BAMs: $SPLIT_R2_DIR"

# ==============================================================================
# Step 5 -- Merge R1+R2 per-barcode BAMs
# ==============================================================================
#
# For each barcode/viewpoint, per read pair: keep whichever of R1/R2 has the
# longer mapped fragment if both mapped, otherwise keep whichever side
# mapped. The R1 manifest is used as the canonical list of viewpoints (same
# convention as run_basic4cseq_merged_bams.R's original header comment).

echo ""
echo "=== Step 5/8: merge R1+R2 per-barcode BAMs ==="

MERGED_DIR="${OUTPUT_DIR}/merged_bams"

python3 "${SCRIPT_DIR}/merge_paired_bams.py" \
    "$SPLIT_R1_DIR" "$SPLIT_R2_DIR" "$MERGED_DIR" \
    --manifest "${SPLIT_R1_DIR}/manifest.tsv"

echo "  Merged per-barcode BAMs: $MERGED_DIR"
echo "  Merge summary:           ${MERGED_DIR}/merge_summary.tsv"

# ==============================================================================
# Step 6 -- Build virtual restriction-fragment library (skip if already built)
# ==============================================================================

echo ""
echo "=== Step 6/8: fragment library ==="

#eval "$(~/anaconda3/bin/conda shell.bash hook)"
#conda activate chicago

#if [[ -f "$FRAG_LIB" ]]; then
   echo "Fragment library already exists, skipping build: $FRAG_LIB"
#else
#    mkdir -p "$FRAGMENT_LIBRARY_DIR"
#    echo "Building fragment library (read length ${READ_LENGTH}bp) -> $FRAG_LIB"
#    Rscript "${SCRIPT_DIR}/build_fragment_library.R" "$READ_LENGTH" "$FRAG_LIB"
#fi

# ==============================================================================
# Step 7 -- Run Basic4Cseq once per barcode/viewpoint, on the merged BAMs
# ==============================================================================

conda deactivate
conda activate chicago

echo ""
echo "=== Step 7/8: Basic4Cseq (merged R1+R2) ==="

RESULTS_DIR="${OUTPUT_DIR}/basic4cseq_results_merged"

Rscript "${SCRIPT_DIR}/run_basic4cseq_merged_bams.R" \
    "${SPLIT_R1_DIR}/manifest.tsv" "$MERGED_DIR" "$FRAG_LIB" "$READ_LENGTH" "$RESULTS_DIR"


# ==============================================================================
# Step 8 -- peakC on the Basic4Cseq wiggle tracks (single-replicate analysis)
# ==============================================================================
#
# For each viewpoint, run_peakc.R filters that viewpoint's Basic4Cseq wig down
# to its own chromosome (peakC's readqWig() expects a single-chromosome track;
# Basic4Cseq writes all chromosomes into one wig), then runs
# readqWig() + single.analysis() and writes a plot_C() PDF plus a peaks TSV.
 
echo ""
echo "=== Step 8/8: peakC ==="
 
PEAKC_DIR="${OUTPUT_DIR}/peakc_results"
 
Rscript "${SCRIPT_DIR}/run_peakc.R" \
    "${SPLIT_R1_DIR}/manifest.tsv" "$RESULTS_DIR" "$PEAKC_DIR" \
    "$PEAKC_WINDOW" "$PEAKC_QWD"

echo ""
echo "Job finished: $(date)"
echo "  Barcode assignment: $READS_WITH_SITES"
echo "  R1 aligned BAM:      $R1_ALIGNED_BAM"
echo "  R2 aligned BAM:      $R2_ALIGNED_BAM"
echo "  Merged per-barcode BAMs: $MERGED_DIR"
echo "  Basic4Cseq results:      $RESULTS_DIR"
echo "  peakC results:           $PEAKC_DIR"
echo "    (plots: <viewpoint>_peakC.pdf, peaks: <viewpoint>_peaks.tsv,"
echo "     summary: ${PEAKC_DIR}/peakc_run_summary.tsv)"
