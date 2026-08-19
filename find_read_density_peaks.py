#!/usr/bin/env python3
"""
find_read_density_peaks.py

For each per-barcode BAM (from split_bam_by_barcode.py), progressively
finds the region with the highest read density at decreasing window sizes:
    1. Which chromosome has the most reads?
    2. Within that chromosome, which 10 Mb window has the most reads?
    3. Within that 10 Mb window, which 1 Mb window has the most reads?
    4. Within that 1 Mb window, which 100 kb window has the most reads?

This helps diagnose whether per-barcode BAMs have reads concentrated near
the declared integration site (expected) or somewhere unexpected (suggests
a barcode allocation or alignment problem).

Usage:
    python3 find_read_density_peaks.py <bam_dir> [--manifest manifest.tsv]
    python3 find_read_density_peaks.py <bam_dir> --manifest manifest.tsv --viewpoints whitelist.tsv

Arguments:
    bam_dir      directory containing the per-barcode BAMs
    --manifest   optional: manifest.tsv from split_bam_by_barcode.py
                 (columns: position, sanitized_filename). If given, only
                 the BAMs listed in the manifest are processed, and the
                 output labels each BAM with its position/barcode label.
    --viewpoints optional: whitelist TSV (columns: position, barcode)
                 e.g. from assign_barcodes.py. If given, the declared
                 integration site is shown alongside the peak so you can
                 see at a glance whether the read peak matches expectation.
"""

import sys
import os
import argparse
from collections import defaultdict
import pysam


WINDOW_SIZES = [10_000_000, 1_000_000, 100_000]
WINDOW_SIZE_LABELS = ["10 Mb", "1 Mb", "100 kb"]


def count_reads_by_chrom(bam):
    """Return a dict of {chrom: read_count} from the BAM index stats."""
    counts = defaultdict(int)
    for stat in bam.get_index_statistics():
        counts[stat.contig] = stat.mapped
    return counts


def count_reads_in_windows(bam, chrom, chrom_length, window_size):
    """
    Return (best_start, best_end, best_count) for the non-overlapping
    window of `window_size` on `chrom` with the most aligned reads.
    Uses pysam.fetch() to count reads starting within each window.
    """
    best_start = 0
    best_count = 0

    start = 0
    while start < chrom_length:
        end = min(start + window_size, chrom_length)
        count = bam.count(chrom, start, end)
        if count > best_count:
            best_count = count
            best_start = start
        start += window_size

    return best_start, min(best_start + window_size, chrom_length), best_count


def process_bam(bam_path, label, declared_position=None):
    """
    Run the progressive window narrowing for one BAM file.
    Returns a dict with the results, prints a formatted summary.
    """
    print(f"\n{'='*60}")
    print(f"BAM:      {os.path.basename(bam_path)}")
    print(f"Label:    {label}")
    if declared_position:
        print(f"Declared: {declared_position}")
    print(f"{'='*60}")

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        print(f"  ERROR opening BAM: {e}")
        return None

    # Check the BAM is indexed
    try:
        bam.check_index()
    except Exception:
        print("  WARNING: BAM is not indexed -- attempting to index now ...")
        bam.close()
        try:
            pysam.index(bam_path)
            bam = pysam.AlignmentFile(bam_path, "rb")
        except Exception as e:
            print(f"  ERROR: could not index BAM: {e}")
            return None

    total_reads = bam.mapped
    if total_reads == 0:
        print("  No mapped reads in this BAM.")
        bam.close()
        return {"label": label, "total_reads": 0}

    print(f"  Total mapped reads: {total_reads:,}")

    # Build chrom -> length map from BAM header
    chrom_lengths = {sq["SN"]: sq["LN"] for sq in bam.header.to_dict()["SQ"]}

    # Step 1: best chromosome
    chrom_counts = count_reads_by_chrom(bam)
    best_chrom = max(chrom_counts, key=lambda c: chrom_counts[c])
    best_chrom_count = chrom_counts[best_chrom]
    best_chrom_pct = 100 * best_chrom_count / total_reads
    print(f"\n  Step 1 — best chromosome:")
    print(f"    {best_chrom}: {best_chrom_count:,} reads ({best_chrom_pct:.1f}%)")

    # Print top 5 chromosomes for context
    top_chroms = sorted(chrom_counts.items(), key=lambda x: -x[1])[:5]
    if len(top_chroms) > 1:
        print("    (top 5 chromosomes: " +
              ", ".join(f"{c}={n:,}" for c, n in top_chroms) + ")")

    chrom_len = chrom_lengths.get(best_chrom, 0)
    if chrom_len == 0:
        print(f"  WARNING: chromosome {best_chrom} not found in BAM header, cannot proceed.")
        bam.close()
        return None

    # Steps 2-4: progressive window narrowing within the best region
    current_region_start = 0
    current_region_end = chrom_len
    best_windows = []

    for window_size, window_label in zip(WINDOW_SIZES, WINDOW_SIZE_LABELS):
        # Only search within the window found at the previous level
        # Build windows within the current region
        best_start_in_region = current_region_start
        best_count_in_region = 0

        start = current_region_start
        while start < current_region_end:
            end = min(start + window_size, current_region_end)
            count = bam.count(best_chrom, start, end)
            if count > best_count_in_region:
                best_count_in_region = count
                best_start_in_region = start
            start += window_size

        best_end_in_region = min(best_start_in_region + window_size, current_region_end)
        pct = 100 * best_count_in_region / total_reads
        step_num = WINDOW_SIZES.index(window_size) + 2

        print(f"\n  Step {step_num} — best {window_label} window:")
        print(f"    {best_chrom}:{best_start_in_region:,}-{best_end_in_region:,}: "
              f"{best_count_in_region:,} reads ({pct:.1f}% of total)")

        best_windows.append((best_chrom, best_start_in_region, best_end_in_region,
                              best_count_in_region, window_label))

        # Narrow search region to this window for the next level
        current_region_start = best_start_in_region
        current_region_end = best_end_in_region

    bam.close()

    # Summary line: compare the 100kb peak to the declared viewpoint
    final_chrom, final_start, final_end, final_count, _ = best_windows[-1]
    if declared_position:
        print(f"\n  Declared integration site: {declared_position}")
        # Try to parse declared position for distance check
        import re
        m = re.match(r"^(chr[A-Za-z0-9]+):([0-9]+)\([+-]\)$", declared_position)
        if m:
            decl_chr = m.group(1)
            decl_pos = int(m.group(2))
            if decl_chr != final_chrom:
                print(f"  *** MISMATCH: reads peak on {final_chrom}, "
                      f"declared site is on {decl_chr} ***")
            elif decl_pos < final_start or decl_pos > final_end:
                dist = min(abs(decl_pos - final_start), abs(decl_pos - final_end))
                print(f"  *** MISMATCH: declared site at {decl_pos:,} is "
                      f"{dist:,} bp outside the 100kb peak window ***")
            else:
                print(f"  OK: declared site at {decl_pos:,} falls within "
                      f"the 100kb peak window.")

    return {
        "label": label,
        "total_reads": total_reads,
        "best_chrom": best_chrom,
        "best_chrom_pct": best_chrom_pct,
        "best_100kb": f"{final_chrom}:{final_start}-{final_end}",
        "best_100kb_count": final_count,
    }


def load_manifest(manifest_path):
    """Return list of (position_label, sanitized_filename) tuples."""
    entries = []
    with open(manifest_path) as f:
        header = f.readline()  # skip header
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) >= 2:
                entries.append((fields[0], fields[1]))
    return entries


def load_viewpoints(viewpoints_path):
    """Return dict: position_label -> declared_position (same thing, just
    confirms the label is in the whitelist). Also useful if the whitelist
    format has extra info you want to cross-reference."""
    vp = {}
    with open(viewpoints_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) >= 2:
                position, barcode = fields[0], fields[1]
                vp[position] = position  # label maps to itself for display
    return vp


def main():
    parser = argparse.ArgumentParser(
        description="Progressively find the highest-read-density region in each per-barcode BAM."
    )
    parser.add_argument("bam_dir", help="Directory containing per-barcode BAM files")
    parser.add_argument(
        "--manifest", type=str, default=None,
        help="manifest.tsv from split_bam_by_barcode.py (position, sanitized_filename)"
    )
    parser.add_argument(
        "--viewpoints", type=str, default=None,
        help="Whitelist TSV (position, barcode) -- used to show declared integration site "
             "alongside the peak for each BAM"
    )
    args = parser.parse_args()

    if not os.path.isdir(args.bam_dir):
        sys.exit(f"Error: BAM directory not found: {args.bam_dir}")

    viewpoints = {}
    if args.viewpoints:
        viewpoints = load_viewpoints(args.viewpoints)
        print(f"Loaded {len(viewpoints)} viewpoints from {args.viewpoints}")

    if args.manifest:
        entries = load_manifest(args.manifest)
        print(f"Processing {len(entries)} BAMs from manifest: {args.manifest}")
        bam_entries = [
            (os.path.join(args.bam_dir, fname), label)
            for label, fname in entries
        ]
    else:
        # Fall back to all *.bam files in the directory
        bam_files = sorted(f for f in os.listdir(args.bam_dir) if f.endswith(".bam"))
        bam_entries = [
            (os.path.join(args.bam_dir, f), f.replace(".bam", ""))
            for f in bam_files
        ]
        print(f"Processing {len(bam_entries)} BAM files found in {args.bam_dir}")

    all_results = []
    for bam_path, label in bam_entries:
        # When using a manifest, the label IS the position/barcode label
        # (e.g. "chr18:77596512(+)"), so use it directly as declared_position
        # unless a separate --viewpoints file overrides it.
        declared = viewpoints.get(label, label)
        result = process_bam(bam_path, label, declared_position=declared)
        if result:
            all_results.append(result)

    # Print a final summary table
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"{'Label':<30} {'TotalReads':>12} {'BestChrom%':>12} {'Best100kb':>35} {'Count':>8}")
    print("-" * 105)
    for r in all_results:
        if r["total_reads"] == 0:
            print(f"{r['label']:<30} {'0':>12} {'-':>12} {'-':>35} {'-':>8}")
        else:
            print(
                f"{r['label']:<30} "
                f"{r['total_reads']:>12,} "
                f"{r['best_chrom_pct']:>11.1f}% "
                f"{r['best_100kb']:>35} "
                f"{r['best_100kb_count']:>8,}"
            )


if __name__ == "__main__":
    main()
