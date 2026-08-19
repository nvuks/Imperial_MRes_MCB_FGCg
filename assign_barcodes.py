#!/usr/bin/env python3
"""
Assign each read's extracted barcode to its nearest whitelist barcode
(within a configurable number of mismatches), and attach that whitelist
entry's integration-site coordinates.

Whitelist input format (tab-separated, no header):
    label   barcode   start   end
e.g.
    NA                       TGTATCGGGACAGGGA   1716448   1956464
    chr18:77596512(+)        TTGTTGCGTTCGATCA   1428500   2958692

Reads input format (tab-separated, no header):
    ReadID   barcode
(i.e. the output of the earlier extraction pipeline)

Output (tab-separated):
    ReadID   whitelist_barcode   position   orientation   match_type
"whitelist_barcode" is the matched whitelist entry's own (forward-orientation)
barcode sequence -- not the read's raw extracted barcode. "orientation" is
"Forward" if the read's barcode matched the whitelist barcode as-is, or "RC"
if it matched the whitelist barcode's reverse complement. "match_type" is
"Exact" or "Mismatch" depending on whether the read's barcode matched
perfectly or required allowed mismatch(es) to assign. A read whose barcode
does not match any whitelist barcode within the allowed mismatches (in
either orientation), or matches more than one whitelist barcode ambiguously,
is reported with whitelist_barcode unchanged (the read's own barcode) and
position/orientation/match_type = "unassigned".

Usage:
    python3 assign_barcodes.py whitelist.tsv barcodes_with_ids.tsv > reads_with_sites.tsv
    python3 assign_barcodes.py --max-mismatches 1 whitelist.tsv barcodes_with_ids.tsv > out.tsv
"""

import sys
import argparse
from itertools import combinations, product


def hamming_distance_n_variants(barcode, n, alphabet="ACGT"):
    """Yield all sequences within Hamming distance exactly `n` of `barcode`."""
    length = len(barcode)
    for positions in combinations(range(length), n):
        original_bases = [barcode[p] for p in positions]
        choices_per_position = [
            [b for b in alphabet if b != original_bases[idx]]
            for idx in range(n)
        ]
        for replacement_bases in product(*choices_per_position):
            variant = list(barcode)
            for pos, new_base in zip(positions, replacement_bases):
                variant[pos] = new_base
            yield "".join(variant)


def variants_by_distance_up_to_n(barcode, n, alphabet="ACGT"):
    """Yield (distance, variant) pairs for distance 0..n of `barcode`,
    where distance 0 is the barcode itself."""
    yield 0, barcode
    for d in range(1, n + 1):
        for variant in hamming_distance_n_variants(barcode, d, alphabet):
            yield d, variant


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(complement.get(base, "N") for base in reversed(seq))


def load_whitelist(path):
    """Return list of (position, barcode) tuples."""
    entries = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                sys.exit(f"Whitelist line does not have 2 tab-separated fields: {line!r}")
            position, barcode = fields
            entries.append((position, barcode.strip()))
    return entries


def build_lookup(whitelist_entries, max_mismatches):
    """
    Build a dict mapping every barcode variant (within max_mismatches of a
    whitelist barcode, in EITHER forward or reverse-complement orientation)
    to a (position, whitelist_barcode, orientation, match_type) tuple, where:
      - orientation is "Forward" or "RC", depending on which orientation of
        the whitelist barcode the variant was generated from.
      - match_type is "Exact" (distance 0) or "Mismatch" (distance 1..n).

    If a variant is reachable from more than one whitelist barcode -- whether
    via a different barcode's forward sequence or its reverse complement --
    it is mapped to None so it is treated as unassigned rather than silently
    picking one.
    """
    lookup = {}

    def register(variant, position, barcode, orientation, match_type, source_barcode):
        existing = lookup.get(variant)
        if variant in lookup and existing is not None and existing[1] != source_barcode:
            lookup[variant] = None  # ambiguous: collides with a different whitelist barcode
        elif variant not in lookup:
            lookup[variant] = (position, barcode, orientation, match_type)

    for entry in whitelist_entries:
        position, barcode = entry
        rc_barcode = reverse_complement(barcode)

        for distance, variant in variants_by_distance_up_to_n(barcode, max_mismatches):
            match_type = "Exact" if distance == 0 else "Mismatch"
            register(variant, position, barcode, "Forward", match_type, barcode)

        for distance, variant in variants_by_distance_up_to_n(rc_barcode, max_mismatches):
            match_type = "Exact" if distance == 0 else "Mismatch"
            register(variant, position, barcode, "RC", match_type, barcode)

    return lookup


def main():
    parser = argparse.ArgumentParser(
        description="Assign reads to nearest whitelist barcode within N mismatches."
    )
    parser.add_argument("whitelist_file", help="Tab-separated: position, barcode")
    parser.add_argument("reads_file", help="Tab-separated: ReadID, barcode")
    parser.add_argument(
        "--max-mismatches", type=int, default=1,
        help="Maximum Hamming distance allowed for assignment (default: 1)"
    )
    args = parser.parse_args()

    whitelist_entries = load_whitelist(args.whitelist_file)

    # Sanity check: warn if whitelist barcodes vary in length (Hamming distance
    # comparisons assume equal length)
    lengths = {len(b) for _, b in whitelist_entries}
    if len(lengths) > 1:
        sys.stderr.write(f"Warning: whitelist barcodes have differing lengths: {lengths}\n")

    lookup = build_lookup(whitelist_entries, args.max_mismatches)

    n_total = 0
    n_assigned = 0
    n_unassigned = 0
    stats = {
        ("Forward", "Exact"): 0,
        ("Forward", "Mismatch"): 0,
        ("RC", "Exact"): 0,
        ("RC", "Mismatch"): 0,
    }

    with open(args.reads_file) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            read_id, barcode = line.split("\t", 1)
            n_total += 1

            entry = lookup.get(barcode)
            if entry is not None:
                position, wl_barcode, orientation, match_type = entry
                print(f"{read_id}\t{wl_barcode}\t{position}\t{orientation}\t{match_type}")
                n_assigned += 1
                stats[(orientation, match_type)] += 1
            else:
                print(f"{read_id}\t{barcode}\tunassigned\tunassigned\tunassigned")
                n_unassigned += 1

    sys.stderr.write(
        f"Total reads: {n_total}, assigned: {n_assigned}, unassigned: {n_unassigned}\n"
    )
    sys.stderr.write(
        f"  Forward, exact match:    {stats[('Forward', 'Exact')]}\n"
        f"  Forward, with mismatch:  {stats[('Forward', 'Mismatch')]}\n"
        f"  RC, exact match:         {stats[('RC', 'Exact')]}\n"
        f"  RC, with mismatch:       {stats[('RC', 'Mismatch')]}\n"
    )


if __name__ == "__main__":
    main()
