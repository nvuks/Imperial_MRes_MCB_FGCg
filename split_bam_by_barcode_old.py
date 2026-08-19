#!/usr/bin/env python3
"""
Split an aligned R1 BAM into one BAM per whitelist barcode / integration
site, using the R2-derived barcode assignment TSV to look up, for each
ReadID, which barcode (viewpoint) it belongs to.

This produces one BAM per viewpoint, ready to feed into a multi-sample 4C
analysis package (e.g. FourCSeq), where each viewpoint's known integration
site acts as a synthetic 4C "anchor".

Inputs
------
R1 BAM (aligned, e.g. from bowtie2 + samtools sort):
    Standard BAM, indexed or not. Only ReadID, chrom, pos, strand are used.

Barcode assignment TSV (output of assign_barcodes.py):
    ReadID  whitelist_barcode  position  orientation  match_type
e.g.
    read1  TGTATCGGGACAGGGA  NA                Forward  Exact
    read2  TTGTTGCGTTCGATCA  chr18:77596512(+)  Forward  Mismatch
Reads with position == "unassigned" are skipped (not written to any
per-barcode BAM).

Output
------
One BAM per distinct, assigned `position` value, named:
    <output_dir>/<sanitized_position>.bam
e.g. chr18:77596512(+) becomes chr18_77596512_plus.bam (see
sanitize_label() for exact rules). A manifest TSV mapping sanitized
filenames back to their original position labels is also written, since
filenames can't contain all the original characters safely.

A read whose ReadID appears in the BAM but not in the assignment TSV (or
vice versa) is simply not written to any output BAM; counts of such
mismatches are reported at the end for sanity-checking.

Usage:
    python3 split_bam_by_barcode.py aligned.sorted.bam barcode_assignments.tsv output_dir/
"""

import sys
import os
import argparse
import re
from collections import defaultdict

import pysam


def sanitize_label(label):
    """Turn a position label like 'chr18:77596512(+)' into a safe filename
    component, e.g. 'chr18_77596512_plus'. Keeps things human-readable
    rather than hashing, since there are only a handful of barcodes."""
    label = label.replace(":", "_").replace("(", "_").replace(")", "")
    label = label.replace("+", "plus").replace("-", "minus")
    label = re.sub(r"[^A-Za-z0-9_]", "_", label)
    label = re.sub(r"_+", "_", label).strip("_")
    return label


def load_assignments(path):
    """Return dict: read_id -> position label (only for assigned reads)."""
    assignments = {}
    n_total = 0
    n_unassigned = 0
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) < 3:
                sys.exit(f"Assignment line does not have enough fields: {line!r}")
            read_id, _discarded_part_of_read_id, _barcode, position = fields[0], fields[1], fields[2], fields[3]
            #sys.stderr.write(f"{read_id}, {_discarded_part_of_read_id}, {_barcode}, {position}\n")
            n_total += 1
            if position == "unassigned":
                n_unassigned += 1
                continue
            assignments[read_id] = position
    sys.stderr.write(
        f"Loaded {n_total} read assignments ({n_unassigned} unassigned, "
        f"{len(assignments)} assigned to a barcode/viewpoint)\n"
    )
    return assignments


def main():
    parser = argparse.ArgumentParser(
        description="Split an aligned R1 BAM into one BAM per whitelist barcode/viewpoint."
    )
    parser.add_argument("input_bam", help="Aligned R1 BAM file")
    parser.add_argument("assignment_tsv", help="ReadID-to-barcode assignment TSV (output of assign_barcodes.py)")
    parser.add_argument("output_dir", help="Directory to write per-barcode BAM files into")
    args = parser.parse_args()

    if not os.path.isfile(args.input_bam):
        sys.exit(f"Error: input BAM not found: {args.input_bam}")
    if not os.path.isfile(args.assignment_tsv):
        sys.exit(f"Error: assignment TSV not found: {args.assignment_tsv}")

    os.makedirs(args.output_dir, exist_ok=True)

    assignments = load_assignments(args.assignment_tsv)

    # Determine the set of distinct position labels up front, so we can
    # open one output BAM writer per label.
    distinct_positions = sorted(set(assignments.values()))
    sanitized_map = {}
    used_sanitized = {}
    for position in distinct_positions:
        sanitized = sanitize_label(position)
        if sanitized in used_sanitized and used_sanitized[sanitized] != position:
            sys.exit(
                f"Error: position labels '{position}' and "
                f"'{used_sanitized[sanitized]}' both sanitize to the same "
                f"filename '{sanitized}'. Please rename one in the whitelist."
            )
        used_sanitized[sanitized] = position
        sanitized_map[position] = sanitized

    manifest_path = os.path.join(args.output_dir, "manifest.tsv")
    with open(manifest_path, "w") as manifest:
        manifest.write("position\tsanitized_filename\n")
        for position, sanitized in sanitized_map.items():
            manifest.write(f"{position}\t{sanitized}.bam\n")

    in_bam = pysam.AlignmentFile(args.input_bam, "rb")

    writers = {}
    for position, sanitized in sanitized_map.items():
        out_path = os.path.join(args.output_dir, f"{sanitized}.bam")
        writers[position] = pysam.AlignmentFile(out_path, "wb", template=in_bam)

    n_reads_total = 0
    n_written = defaultdict(int)
    n_no_assignment = 0  # aligned read whose ID isn't in the assignment TSV at all

    for read in in_bam.fetch(until_eof=True):
        n_reads_total += 1
        read_id = read.query_name
        #sys.stderr.write(f"{read_id}\n")
        position = assignments.get("@" + read_id)
        if position is None:
            n_no_assignment += 1
            continue
        writers[position].write(read)
        n_written[position] += 1

    in_bam.close()
    for writer in writers.values():
        writer.close()

    # Sort + index each output BAM (FourCSeq/Basic4Cseq and most downstream
    # tools expect coordinate-sorted, indexed BAMs).
    for position, sanitized in sanitized_map.items():
        out_path = os.path.join(args.output_dir, f"{sanitized}.bam")
        sorted_path = os.path.join(args.output_dir, f"{sanitized}.sorted.bam")
        pysam.sort("-o", sorted_path, out_path)
        os.replace(sorted_path, out_path)
        pysam.index(out_path)

    sys.stderr.write(f"\nTotal aligned reads in input BAM: {n_reads_total}\n")
    sys.stderr.write(f"Reads with no entry in assignment TSV (skipped): {n_no_assignment}\n")
    sys.stderr.write("Reads written per barcode/viewpoint:\n")
    for position in distinct_positions:
        sys.stderr.write(f"  {position}: {n_written[position]}\n")
    sys.stderr.write(f"\nManifest written to: {manifest_path}\n")
    sys.stderr.write(f"Per-barcode BAMs (sorted + indexed) written to: {args.output_dir}\n")


if __name__ == "__main__":
    main()
