#!/usr/bin/env python3
"""
merge_paired_bams.py

Merges per-barcode R1 and R2 BAM directories produced by two separate runs
of split_bam_by_barcode.py (one on the aligned R1 FASTQ, one on the aligned
R2 FASTQ). For each viewpoint/barcode:

  - Reads present in BOTH R1 and R2 BAMs (matched on base ReadID, i.e.
    query_name before the first space): keep whichever has the longer
    query_length. One BAM record per pair in the output.
  - Reads present in only ONE of R1 or R2: keep that read as-is.
  - Reads present in neither: not possible by construction.

Output: one merged BAM per viewpoint (sorted + indexed), plus a summary
TSV reporting per-viewpoint read counts broken down by origin:
    both_r1_won, both_r2_won, both_tied_r1, r1_only, r2_only

This breakdown tells you how much is gained by running the analysis in a
paired-end setting (r1_only + r2_only are reads that would have been missed
if you'd only processed one side).

The unassigned.bam files written by split_bam_by_barcode.py (reads whose
barcode could not be assigned to a viewpoint) are paired the same way and
written to <output_dir>/unassigned.bam, so the paired-end merge is applied
to unassigned reads too. This is on by default; pass --no-include-unassigned
to skip it. Note the manifest does NOT list unassigned, so without this the
manifest-driven path would silently ignore those files. The unassigned row
appears in merge_summary.tsv but is kept out of the "aggregate across all
viewpoints" totals, since it is not a real viewpoint.

Usage:
    python3 merge_paired_bams.py <r1_bam_dir> <r2_bam_dir> <output_dir> [--manifest manifest.tsv]

Arguments:
    r1_bam_dir             Directory of per-barcode R1 BAMs (from split_bam_by_barcode.py)
    r2_bam_dir             Directory of per-barcode R2 BAMs (from split_bam_by_barcode.py)
    output_dir             Directory to write merged per-barcode BAMs and summary into
    --manifest             manifest.tsv from split_bam_by_barcode.py (preferred -- ensures
                           both dirs are processed consistently). If omitted, all *.bam
                           files found in r1_bam_dir are used as the reference list.
    --no-include-unassigned  Do not pair/merge the unassigned.bam files (default: include).

Notes:
    - Matching is done on the BASE ReadID (query_name before the first space),
      since R1 and R2 BAM records for the same pair share this base ID but
      may differ in the trailing mate-pair suffix.
    - Both BAMs for a viewpoint are read entirely into memory (one dict per
      side). For very large per-barcode BAMs this may use significant memory;
      if that becomes an issue, sort both BAMs by query name and stream them
      instead.
    - Requires pysam.
"""

import sys
import os
import argparse
from collections import defaultdict
import pysam


def base_read_id(query_name):
    """Strip mate-pair suffix (everything after first space) from a BAM
    query_name, giving the base ID shared between R1 and R2."""
    return query_name.split()[0] if query_name else query_name


def load_bam_by_readid(bam_path):
    """Read all mapped reads from a BAM into a dict: base_read_id -> read.
    Only mapped reads are included (unmapped reads have no useful alignment
    to contribute to the merged output)."""
    reads = {}
    if not os.path.exists(bam_path):
        return reads
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            rid = base_read_id(read.query_name)
            reads[rid] = read
    return reads


def load_manifest(manifest_path):
    """Return list of (position_label, sanitized_filename) tuples."""
    entries = []
    with open(manifest_path) as f:
        f.readline()  # skip header
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) >= 2:
                entries.append((fields[0], fields[1]))
    return entries


def merge_viewpoint(r1_bam_path, r2_bam_path, out_bam_path):
    """
    Merge R1 and R2 BAMs for one viewpoint. Returns a dict of counts:
        both_r1_won, both_r2_won, both_tied_r1, r1_only, r2_only, total
    """
    r1_reads = load_bam_by_readid(r1_bam_path)
    r2_reads = load_bam_by_readid(r2_bam_path)

    counts = defaultdict(int)

    # We need the BAM header from whichever file exists, to write output.
    # Prefer R1; fall back to R2.
    header_source = r1_bam_path if os.path.exists(r1_bam_path) else r2_bam_path
    if not os.path.exists(header_source):
        return counts

    with pysam.AlignmentFile(header_source, "rb") as src:
        header = src.header.to_dict()

    # Build the merged set of reads to write
    output_reads = []
    all_ids = set(r1_reads) | set(r2_reads)

    for rid in all_ids:
        r1 = r1_reads.get(rid)
        r2 = r2_reads.get(rid)

        if r1 is not None and r2 is not None:
            # Both sides aligned -- keep the longer fragment by query_length.
            # Ties go to R2.
            l1 = r1.query_length or 0
            l2 = r2.query_length or 0
            if l2 >= l1:
                output_reads.append(r2)
                counts["both_r2_won" if l2 > l1 else "both_tied_r2"] += 1
            else:
                output_reads.append(r1)
                counts["both_r1_won"] += 1
        elif r1 is not None:
            output_reads.append(r1)
            counts["r1_only"] += 1
        else:
            output_reads.append(r2)
            counts["r2_only"] += 1

    counts["total"] = len(output_reads)

    # Write the merged BAM
    os.makedirs(os.path.dirname(os.path.abspath(out_bam_path)), exist_ok=True)
    tmp_path = out_bam_path + ".tmp.bam"
    with pysam.AlignmentFile(tmp_path, "wb", header=header) as out:
        for read in output_reads:
            out.write(read)

    # Sort and index
    pysam.sort("-o", out_bam_path, tmp_path)
    os.remove(tmp_path)
    pysam.index(out_bam_path)

    return counts


def main():
    parser = argparse.ArgumentParser(
        description="Merge per-barcode R1 and R2 BAMs, keeping longer fragment per pair."
    )
    parser.add_argument("r1_bam_dir", help="Directory of per-barcode R1 BAMs")
    parser.add_argument("r2_bam_dir", help="Directory of per-barcode R2 BAMs")
    parser.add_argument("output_dir", help="Directory for merged BAMs and summary")
    parser.add_argument(
        "--manifest", type=str, default=None,
        help="manifest.tsv from split_bam_by_barcode.py (columns: position, "
             "sanitized_filename). If omitted, all *.bam files in r1_bam_dir "
             "are used."
    )
    parser.add_argument(
        "--include-unassigned", dest="include_unassigned", action="store_true",
        help="Also pair/merge unassigned.bam from each dir (default: on)."
    )
    parser.add_argument(
        "--no-include-unassigned", dest="include_unassigned", action="store_false",
        help="Do not pair/merge unassigned.bam."
    )
    parser.set_defaults(include_unassigned=True)
    args = parser.parse_args()

    for d in [args.r1_bam_dir, args.r2_bam_dir]:
        if not os.path.isdir(d):
            sys.exit(f"Error: directory not found: {d}")
    os.makedirs(args.output_dir, exist_ok=True)

    # Build list of (label, filename) to process
    if args.manifest:
        entries = load_manifest(args.manifest)
        print(f"Processing {len(entries)} viewpoint(s) from manifest: {args.manifest}")
    else:
        bam_files = sorted(f for f in os.listdir(args.r1_bam_dir)
                           if f.endswith(".bam") and not f.endswith(".bai"))
        entries = [(f.replace(".bam", ""), f) for f in bam_files]
        print(f"Processing {len(entries)} BAM file(s) found in {args.r1_bam_dir}")

    # Ensure the unassigned reads get paired too. The manifest never lists
    # "unassigned", so on the manifest-driven path these files would otherwise
    # be ignored. On the glob path unassigned.bam may already be in `entries`,
    # so guard against adding it twice.
    UNASSIGNED_LABEL = "unassigned"
    UNASSIGNED_FILENAME = "unassigned.bam"
    if args.include_unassigned:
        already_listed = any(fn == UNASSIGNED_FILENAME for _, fn in entries)
        if not already_listed:
            in_r1 = os.path.exists(os.path.join(args.r1_bam_dir, UNASSIGNED_FILENAME))
            in_r2 = os.path.exists(os.path.join(args.r2_bam_dir, UNASSIGNED_FILENAME))
            if in_r1 or in_r2:
                entries.append((UNASSIGNED_LABEL, UNASSIGNED_FILENAME))
                print(f"Including unassigned reads: {UNASSIGNED_FILENAME}")
            else:
                print(f"No {UNASSIGNED_FILENAME} found in either dir; "
                      f"nothing to merge for unassigned reads")

    summary_rows = []

    for label, filename in entries:
        r1_path = os.path.join(args.r1_bam_dir, filename)
        r2_path = os.path.join(args.r2_bam_dir, filename)
        out_path = os.path.join(args.output_dir, filename)

        r1_exists = os.path.exists(r1_path)
        r2_exists = os.path.exists(r2_path)

        if not r1_exists and not r2_exists:
            print(f"  SKIP {label}: BAM missing from both R1 and R2 dirs")
            continue

        if not r1_exists:
            print(f"  WARNING {label}: R1 BAM missing, using R2 only")
        if not r2_exists:
            print(f"  WARNING {label}: R2 BAM missing, using R1 only")

        print(f"  Merging: {label} ...", end=" ", flush=True)
        counts = merge_viewpoint(r1_path, r2_path, out_path)

        total = counts["total"]
        both = counts["both_r1_won"] + counts["both_r2_won"] + counts["both_tied_r2"]
        r1_only = counts["r1_only"]
        r2_only = counts["r2_only"]
        paired_pct = 100 * both / total if total else 0
        r1_only_pct = 100 * r1_only / total if total else 0
        r2_only_pct = 100 * r2_only / total if total else 0

        print(
            f"total={total:,}  "
            f"both={both:,} ({paired_pct:.0f}%)  "
            f"r1_only={r1_only:,} ({r1_only_pct:.0f}%)  "
            f"r2_only={r2_only:,} ({r2_only_pct:.0f}%)"
        )

        summary_rows.append({
            "viewpoint": label,
            "total": total,
            "both_r1_won": counts["both_r1_won"],
            "both_r2_won": counts["both_r2_won"],
            "both_tied_r2": counts["both_tied_r2"],
            "both_total": both,
            "r1_only": r1_only,
            "r2_only": r2_only,
            "paired_pct": f"{paired_pct:.1f}",
            "r1_only_pct": f"{r1_only_pct:.1f}",
            "r2_only_pct": f"{r2_only_pct:.1f}",
        })

    # Write summary TSV
    summary_path = os.path.join(args.output_dir, "merge_summary.tsv")
    cols = ["viewpoint", "total", "both_total", "both_r1_won", "both_r2_won",
            "both_tied_r2", "r1_only", "r2_only",
            "paired_pct", "r1_only_pct", "r2_only_pct"]
    with open(summary_path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for row in summary_rows:
            f.write("\t".join(str(row[c]) for c in cols) + "\n")

    print(f"\nSummary written to: {summary_path}")
    print(f"Merged BAMs written to: {args.output_dir}")

    # Print a quick aggregate to stdout. The unassigned row is reported
    # separately so it does not distort the per-viewpoint recovery stats.
    viewpoint_rows = [r for r in summary_rows if r["viewpoint"] != UNASSIGNED_LABEL]
    if viewpoint_rows:
        total_all = sum(r["total"] for r in viewpoint_rows)
        both_all = sum(r["both_total"] for r in viewpoint_rows)
        r1_only_all = sum(r["r1_only"] for r in viewpoint_rows)
        r2_only_all = sum(r["r2_only"] for r in viewpoint_rows)
        print(f"\nAggregate across all viewpoints:")
        print(f"  Total reads:   {total_all:,}")
        if total_all:
            print(f"  Both sides:    {both_all:,} ({100*both_all/total_all:.1f}%)")
            print(f"  R1 only:       {r1_only_all:,} ({100*r1_only_all/total_all:.1f}%)")
            print(f"  R2 only:       {r2_only_all:,} ({100*r2_only_all/total_all:.1f}%)")
            print(f"  Gained from R1 (missed with R2 alone): "
                  f"{r1_only_all:,} ({100*r1_only_all/total_all:.1f}%)")
            print(f"  Gained from R2 (missed with R1 alone): "
                  f"{r2_only_all:,} ({100*r2_only_all/total_all:.1f}%)")

    unassigned_rows = [r for r in summary_rows if r["viewpoint"] == UNASSIGNED_LABEL]
    if unassigned_rows:
        ur = unassigned_rows[0]
        u_total = ur["total"]
        print(f"\nUnassigned reads (not included in the viewpoint aggregate above):")
        print(f"  Total reads:   {u_total:,}")
        if u_total:
            print(f"  Both sides:    {ur['both_total']:,} "
                  f"({100*ur['both_total']/u_total:.1f}%)")
            print(f"  R1 only:       {ur['r1_only']:,} "
                  f"({100*ur['r1_only']/u_total:.1f}%)")
            print(f"  R2 only:       {ur['r2_only']:,} "
                  f"({100*ur['r2_only']/u_total:.1f}%)")


if __name__ == "__main__":
    main()
