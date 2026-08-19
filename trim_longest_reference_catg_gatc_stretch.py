#!/usr/bin/env python3
"""
trim_longest_reference_catg_gatc_stretch_v6.py

Replaces v5's insert / barcode-splice / dimer / chimera-trim logic with a
single rule, per discussion:

    Any undigested-chromatin read (pure insert, insert+barcode+insert,
    insert+insert dimer, or insert/genomic chimera) is of no interest to
    4C -- there's no need to classify *which* of those it is, or to
    rescue a genomic remainder from it. All that matters is: does this
    cut-site-flanked fragment contain enough insert-derived sequence to
    call it undigested chromatin? If yes, drop the whole fragment. If
    no, keep it whole, untouched.

WHAT THIS REPLACES
-------------------
- matches_insert_edlib (whole-fragment / prefix-vs-anchor-tail check)
- matches_barcode_spanning (motif location + splice + re-test)
- the 4-orientation dimer reference set + --check-dimers
- find_chimera_trim + --trim-chimeras (prefix/suffix window rescue)
- --max-insert-mismatches, --motif-max-mismatch, --r2-barcode-skip

...with one function, insert_coverage_fraction(), that greedily tiles a
fragment with local insert-matching windows (both orientations, no need
to handle the barcode's N-run specially -- it just becomes an uncovered
gap between two independently-matching windows, and a dimer junction is
just two independently-matching windows from two copies) and returns
what fraction of the fragment those windows cover. One threshold
(--insert-coverage-threshold) decides keep-whole vs discard-whole.

Adapter trim, poly-G trim, cut-site segmentation, and min-fragment-length
are unchanged from v4/v5.
"""

import sys
import gzip
import argparse

try:
    import edlib
    EDLIB_AVAILABLE = True
except ImportError:
    EDLIB_AVAILABLE = False

CUTTERS = ["CATG", "GATC"]
CUTTER_LEN = 4

DEFAULT_ADAPTER = "AGATCGGAAGAGC"  # standard Illumina universal adapter (Read 1)
DEFAULT_POLYG_MIN_RUN = 8

INSERT_SEQUENCE = (
    "ttaaccctagaaagatagtctgcgtaaaattgacgcatgcattcttgaaatattgctctctctttctaaatag"
    "cgcgaatccgtcgctgtgcatttaggacatctcagtcgccgcttggagctcccgtgaggcgtgcttgtcaatg"
    "cggtaagtgtcactgattttgaactataacgaccgcgtgagtcaaaatgacgcatgattatcttttacgtgact"
    "tttaagatttaactcatacgataattatattgttatttcatgttctacttacgtgataacttattatatatatat"
    "tttcttgttatagatatcaactagaatgctagcatgggcccatctcgaggatccaccggtctagaaAGCTTGAAC"
    "AATGGGAGCGGTGCAGAACAATGGCCGTACTCCTGAACAATGGGGAGAGTATAGAACAATGGGCTCGAAATCTA"
    "GTGGGTTTGGGCGTTTACTATGGGAGGTCTATATAAAGCAGAGCTCGTTTAGTGAACCGTCAGATCCATCGTCGC"
    "ATCCAAGAGCcatgGGTAAGCCTATCCCTAACCCTCTCCTCGGTCTCGATTCTACGtaagaattcgcggccgca"
    "tacgatttaggtgacactgcaGGATCANNNNNNNNNNNNNNNNCTCGAGTTGTGGCCGGCCCTTGTGACTGgga"
    "aaaccctggcgtaaataaaatacgaaatgactagttaaaagttttgttactttatagaagaaattttgagttttt"
    "gtttttttttaataaataaataaacataaataaattgtttgttgaatttattattagtatgtaagtgtaaatata"
    "ataaaacttaatatctattcaaattaataaataaacctcgatatacagaccgataaaacacatgcgtcaattttac"
    "gcatgattatctttaacgtacgtcacaatatgattatctttctagggttaa"
).upper()


def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def build_insert_references(insert):
    """
    Split the insert on its N-run (the barcode placeholder) so we never
    hand edlib a run of Ns to align against -- the two flanking pieces
    are matched independently, which is also exactly what lets a
    barcode-spanning read get covered by two separate windows below.
    """
    n_start = insert.find("N")
    if n_start == -1:
        anchors = [insert]
    else:
        n_end = n_start
        while n_end < len(insert) and insert[n_end] == "N":
            n_end += 1
        anchors = [a for a in [insert[:n_start], insert[n_end:]] if a]
    anchors_rc = [reverse_complement(a) for a in anchors]
    return anchors + anchors_rc


# --------------------------------------------------------------------
# adapter / poly-G cleaning (unchanged from v5)
# --------------------------------------------------------------------

def find_adapter_trim_pos(seq, adapter, max_edits_full=2, min_overlap=6,
                           max_mismatch_rate=0.10):
    if not EDLIB_AVAILABLE:
        return None
    s = seq.upper()
    candidates = []

    if len(s) >= len(adapter):
        r = edlib.align(adapter, s, mode="HW", task="locations", k=max_edits_full)
        if r["editDistance"] != -1 and r["locations"]:
            candidates.append(min(loc[0] for loc in r["locations"]))

    max_ov = min(len(adapter), len(s))
    for ov in range(max_ov, min_overlap - 1, -1):
        read_tail = s[-ov:]
        adapter_prefix = adapter[:ov]
        max_mm = int(ov * max_mismatch_rate)
        mismatches = sum(a != b for a, b in zip(read_tail, adapter_prefix))
        if mismatches <= max_mm:
            candidates.append(len(s) - ov)
            break

    return min(candidates) if candidates else None


def find_polyg_trim_pos(seq, min_run):
    s = seq.upper()
    run_start = None
    run_len = 0
    for i, c in enumerate(s):
        if c == "G":
            if run_len == 0:
                run_start = i
            run_len += 1
            if run_len >= min_run:
                return run_start
        else:
            run_len = 0
    return None


def clean_read(seq, qual, adapter=None, polyg_min_run=None):
    candidates = []
    if adapter:
        p = find_adapter_trim_pos(seq, adapter)
        if p is not None:
            candidates.append((p, "adapter"))
    if polyg_min_run:
        p = find_polyg_trim_pos(seq, polyg_min_run)
        if p is not None:
            candidates.append((p, "polyg"))

    if not candidates:
        return seq, qual, None

    candidates.sort(key=lambda x: (x[0], x[1] != "adapter"))
    cut, reason = candidates[0]
    return seq[:cut], qual[:cut], reason


# --------------------------------------------------------------------
# NEW: single insert-coverage-fraction rule
# --------------------------------------------------------------------

def build_insert_kmer_set(refs, k=15):
    kmers = set()
    for ref in refs:
        for i in range(len(ref) - k + 1):
            kmers.add(ref[i:i + k])
    return kmers


def insert_coverage_fraction(fragment, kmer_set, k=15, gap_fill=10):
    """
    v7: O(n) k-mer-seeded replacement for v6's O(n^2 * edlib) window-
    tiling scan. Slide an exact k-mer window (step 1) across the
    fragment and mark positions covered by any k-mer that also occurs
    in the insert (either orientation, precomputed once at startup).
    Overlapping k-mers (step 1) already tolerate isolated single-base
    errors -- a mismatch only kills the ~k kmers touching it, kmers
    further away still hit exactly. A short post-pass fills small
    uncovered gaps flanked by covered bases on both sides (an isolated
    mismatch/short indel inside an otherwise-matching block), so we
    don't undercount those the way pure exact-kmer coverage would.
    Zero edlib calls.
    """
    n = len(fragment)
    if n < k:
        return 0.0
    s = fragment.upper()
    covered = bytearray(n)
    for i in range(n - k + 1):
        if s[i:i + k] in kmer_set:
            for j in range(i, i + k):
                covered[j] = 1
    i = 0
    while i < n:
        if covered[i] == 0:
            j = i
            while j < n and covered[j] == 0:
                j += 1
            if i > 0 and j < n and (j - i) <= gap_fill:
                for x in range(i, j):
                    covered[x] = 1
            i = j
        else:
            i += 1
    return sum(covered) / n


def find_cut_positions(seq, anchor):

    first_anchor = seq.find(anchor)
    if first_anchor == -1:
        return None

    cuts = [(first_anchor, anchor)]
    pos = first_anchor + CUTTER_LEN

    while pos < len(seq):
        next_pos = len(seq)
        next_site = None
        for cutter in CUTTERS:
            idx = seq.find(cutter, pos)
            if idx != -1 and idx < next_pos:
                next_pos = idx
                next_site = cutter
        if next_site is None:
            break
        cuts.append((next_pos, next_site))
        pos = next_pos + CUTTER_LEN

    return cuts


def best_kept_stretch(seq, kmer_set, anchor, coverage_threshold,
                       k=15, gap_fill=10, min_fragment_length=0):
    """
    Segment `seq` on cut sites, drop any stretch whose insert-coverage
    fraction clears `coverage_threshold` (whole stretch discarded, no
    partial rescue), drop any surviving stretch shorter than
    `min_fragment_length`, and return the longest stretch left -- same
    "longest surviving stretch per read" policy as v4/v5, just with a
    single discard rule instead of several.
    """
    cuts = find_cut_positions(seq, anchor)
    if cuts is None:
        return None, f"no_{anchor.lower()}"

    best_start, best_end = None, None
    best_length = -1
    n_insert = 0
    n_too_short = 0

    for i, (cut_start, _site) in enumerate(cuts):
        if i + 1 < len(cuts):
            stretch_end = cuts[i + 1][0] + CUTTER_LEN
        else:
            stretch_end = len(seq)

        stretch = seq[cut_start:stretch_end]

        frac = insert_coverage_fraction(stretch, kmer_set, k=k, gap_fill=gap_fill)

        if frac >= coverage_threshold:
            n_insert += 1
            continue

        stretch_length = stretch_end - cut_start
        if stretch_length < min_fragment_length:
            n_too_short += 1
            continue

        if stretch_length > best_length:
            best_length = stretch_length
            best_start, best_end = cut_start, stretch_end

    if best_start is None:
        if n_insert > 0 and n_too_short > 0:
            reason = "insert_or_too_short"
        elif n_insert > 0:
            reason = "all_insert"
        elif n_too_short > 0:
            reason = "all_too_short"
        else:
            reason = f"no_{anchor.lower()}"
        return None, reason

    return (best_start, best_end), None


def open_maybe_gzip_read(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def open_maybe_gzip_write(path):
    if path == "-":
        return sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "w")


def read_fastq_records(f):
    while True:
        header = f.readline().rstrip("\n")
        if not header:
            break
        seq = f.readline().rstrip("\n")
        plus = f.readline()
        qual = f.readline().rstrip("\n")
        if not plus or not qual:
            break
        yield header, seq, qual


def main():
    parser = argparse.ArgumentParser(
        description="Trim reads to the longest cut-site-flanked stretch that "
                     "is NOT mostly insert-derived (undigested chromatin). "
                     "Stretches above --insert-coverage-threshold are dropped "
                     "whole; survivors are kept whole, untouched."
    )
    parser.add_argument("input_fastq")
    parser.add_argument("output_fastq")
    parser.add_argument("--anchor", choices=["CATG", "GATC"], default="CATG")
    parser.add_argument(
        "--insert-coverage-threshold", type=float, default=0.2,
        help="Fraction (0-1) of a cut-site-flanked stretch's length that "
             "must be covered by exact-kmer insert matches before the whole "
             "stretch is discarded as undigested chromatin. Default: 0.2, "
             "tuned for R1/CATG on real aligned-vs-unaligned samples: no "
             "real aligned R1 read exceeded 0.195 coverage by chance "
             "(n=3000), so 0.2 gives ~0% false-positive discards of real "
             "reads while still catching ~3.6% of the unaligned pool. "
             "USE --insert-coverage-threshold 0.3 FOR R2/GATC: R2's "
             "barcode-spanning chromatin population is much larger and "
             "separates cleanly at 0.3 (0% FP, ~45% of the unaligned pool "
             "caught, vs no real aligned R2 read exceeding 0.283)."
    )
    parser.add_argument(
        "--kmer-size", type=int, default=15,
        help="Length of the exact k-mer seed used to build insert coverage "
             "(default: 15bp). Overlapping k-mers (step 1) already tolerate "
             "isolated single-base errors -- only kmers touching a mismatch "
             "fail, kmers further away still hit exactly."
    )
    parser.add_argument(
        "--gap-fill", type=int, default=10,
        help="Fill uncovered gaps up to this length (bp) when flanked by "
             "covered bases on both sides -- treats a short run of missed "
             "kmers around an isolated mismatch/indel as still covered, "
             "instead of undercounting it (default: 10bp)."
    )
    parser.add_argument(
        "--adapter-seq", default=DEFAULT_ADAPTER,
        help=f"Adapter sequence to scan for and trim (default: standard "
             f"Illumina universal adapter, {DEFAULT_ADAPTER})."
    )
    parser.add_argument(
        "--no-adapter-trim", action="store_true",
        help="Disable adapter read-through trimming."
    )
    parser.add_argument(
        "--polyg-min-run", type=int, default=DEFAULT_POLYG_MIN_RUN,
        help=f"Minimum length of a G-run to trim as a dark-cycle artefact "
             f"(default: {DEFAULT_POLYG_MIN_RUN}bp)."
    )
    parser.add_argument(
        "--no-polyg-trim", action="store_true",
        help="Disable poly-G tail trimming."
    )
    parser.add_argument(
        "--min-fragment-length", type=int, default=50,
        help="Minimum length (bp, including flanking cut site(s)) for a "
             "surviving stretch to be kept. Default: 50bp. Set to 0 to "
             "disable."
    )
    args = parser.parse_args()

    if not EDLIB_AVAILABLE:
        sys.stderr.write(
            "WARNING: edlib not found -- adapter trimming will be SKIPPED "
            "(insert-coverage discard no longer needs edlib at all in v7). "
            "Install with: pip install edlib\n"
        )

    adapter = None if args.no_adapter_trim else args.adapter_seq.upper()
    polyg_min_run = None if args.no_polyg_trim else args.polyg_min_run

    refs = build_insert_references(INSERT_SEQUENCE)
    kmer_set = build_insert_kmer_set(refs, k=args.kmer_size)

    mode = "R1/forward" if args.anchor == "CATG" else "R2/reverse"
    sys.stderr.write(
        f"Anchor site: {args.anchor}  ({mode} mode)\n"
        f"Insert-coverage discard threshold: {args.insert_coverage_threshold:.0%} "
        f"(kmer size {args.kmer_size}bp, gap-fill {args.gap_fill}bp, "
        f"{len(kmer_set):,} reference kmers)\n"
        f"Adapter trim: {adapter if adapter else 'off'}\n"
        f"Poly-G trim: {'min run ' + str(polyg_min_run) + 'bp' if polyg_min_run else 'off'}\n"
        f"Min fragment length: "
        f"{str(args.min_fragment_length) + 'bp' if args.min_fragment_length else 'off'}\n"
    )

    n_total = 0
    n_kept = 0
    n_no_anchor = 0
    n_all_insert = 0
    n_all_too_short = 0
    n_insert_or_too_short = 0
    n_adapter_trimmed = 0
    n_polyg_trimmed = 0
    bases_removed_adapter = 0
    bases_removed_polyg = 0

    infile = open_maybe_gzip_read(args.input_fastq)
    outfile = open_maybe_gzip_write(args.output_fastq)
    try:
        for header, seq, qual in read_fastq_records(infile):
            n_total += 1

            clean_seq, clean_qual, reason = clean_read(
                seq, qual, adapter=adapter, polyg_min_run=polyg_min_run
            )
            if reason == "adapter":
                n_adapter_trimmed += 1
                bases_removed_adapter += len(seq) - len(clean_seq)
            elif reason == "polyg":
                n_polyg_trimmed += 1
                bases_removed_polyg += len(seq) - len(clean_seq)

            result, cut_reason = best_kept_stretch(
                clean_seq, kmer_set, args.anchor, args.insert_coverage_threshold,
                k=args.kmer_size, gap_fill=args.gap_fill,
                min_fragment_length=args.min_fragment_length,
            )

            if result is None:
                if cut_reason == "all_insert":
                    n_all_insert += 1
                elif cut_reason == "all_too_short":
                    n_all_too_short += 1
                elif cut_reason == "insert_or_too_short":
                    n_insert_or_too_short += 1
                else:
                    n_no_anchor += 1
                continue

            start, end = result
            outfile.write(f"{header}\n{clean_seq[start:end]}\n+\n{clean_qual[start:end]}\n")
            n_kept += 1

    finally:
        if infile is not sys.stdin:
            infile.close()
        if outfile is not sys.stdout:
            outfile.close()

    if n_total:
        a = args.anchor
        pad = " " * (14 - len(a))
        sys.stderr.write(
            f"\nTotal reads:                              {n_total:,}\n"
            f"Kept (longest non-insert stretch):        {n_kept:,} ({100*n_kept/n_total:.1f}%)\n"
            f"Discarded (no {a}):{pad}{n_no_anchor:,} ({100*n_no_anchor/n_total:.1f}%)\n"
            f"Discarded (all stretches = insert):       {n_all_insert:,} ({100*n_all_insert/n_total:.1f}%)\n"
            f"Discarded (all stretches < {args.min_fragment_length}bp):        {n_all_too_short:,} ({100*n_all_too_short/n_total:.1f}%)\n"
            f"Discarded (insert or too-short mix):      {n_insert_or_too_short:,} ({100*n_insert_or_too_short/n_total:.1f}%)\n"
            f"\n"
            f"Adapter-trimmed reads:                    {n_adapter_trimmed:,} "
            f"({100*n_adapter_trimmed/n_total:.1f}%), "
            f"{bases_removed_adapter:,} bases removed\n"
            f"Poly-G-trimmed reads:                     {n_polyg_trimmed:,} "
            f"({100*n_polyg_trimmed/n_total:.1f}%), "
            f"{bases_removed_polyg:,} bases removed\n"
        )
    else:
        sys.stderr.write("Total reads: 0\n")


if __name__ == "__main__":
    main()
