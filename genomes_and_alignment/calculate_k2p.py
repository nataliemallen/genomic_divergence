#!/usr/bin/env python3

# calculate genome-wide pairwise K2P nucleotide divergence from minimap2

# input:  Directory of sorted PAF files (*.srt.paf)
# output: Tab-delimited results file with one row per pair
#
# to use:
#   python3 calc_k2p.py /path/to/paf_dir results_k2p.tsv

import argparse
import glob
import math
import os
import sys

DEFAULT_PAF_DIR = "/scratch/gautschi/allen715/divergence/batch2_alignments"
DEFAULT_OUT_FILE = "results_python_k2p_batch2alignments.tsv"
DEFAULT_MIN_MAPQ = 10
DEFAULT_MIN_ALEN = 500


TRANSITIONS = {"AG", "GA", "CT", "TC"}

def parse_cs(cs, counts):
    """Parse a short-form cs string, updating counts in place.

    counts is a dict with integer keys 'matches', 'ts', 'tv'.
    Mirrors the awk while-loop character by character:
      :N   -> N exact matches           (add N to matches)
      *XY  -> single-base substitution  (classify ts / tv), 3 chars
      +... -> insertion in query        (skip the acgt run)
      -... -> deletion in query         (skip the acgt run)
      other-> unknown operator          (advance one char)
    """
    n = len(cs)
    j = 0
    while j < n:
        op = cs[j]

        if op == ":":
            j += 1
            num = ""
            while j < n and cs[j].isdigit():
                num += cs[j]
                j += 1
            if num:
                counts["matches"] += int(num)

        elif op == "*":
            ref_base = cs[j + 1:j + 2].upper()
            alt_base = cs[j + 2:j + 3].upper()
            pair_bases = ref_base + alt_base
            if pair_bases in TRANSITIONS:
                counts["ts"] += 1
            else:
                counts["tv"] += 1
            j += 3

        elif op == "+" or op == "-":
            j += 1
            while j < n and cs[j] in "acgt":
                j += 1

        else:
            j += 1


def process_paf(paf_path, min_mapq, min_alen):

    with open(paf_path) as fh:
        for line in fh:
            if "tp:A:P" not in line:            
                continue
            fields = line.split()
            if len(fields) < 12:               
                continue
            try:
                mapq = int(fields[11])          
                alen = int(fields[10])          
            except ValueError:
                continue
            if mapq < min_mapq:                 
                continue
            if alen < min_alen:                
                continue

            for i in range(11, len(fields)):
                if fields[i].startswith("cs:Z:"):
                    parse_cs(fields[i][5:], counts)

    matches, ts, tv = counts["matches"], counts["ts"], counts["tv"]
    aligned = matches + ts + tv

    if aligned == 0:
        return (0, 0, 0, None, None, None, "NO_ALIGNMENTS_PASSED_FILTERS")

    P = ts / aligned
    Q = tv / aligned
    arg1 = 1 - 2 * P - Q
    arg2 = 1 - 2 * Q

    if arg1 <= 0 or arg2 <= 0:
        return (aligned, ts, tv, P, Q, None, "K2P_UNDEFINED_SATURATION")

    k2p = -0.5 * math.log(arg1) - 0.25 * math.log(arg2)
    return (aligned, ts, tv, P, Q, k2p, "OK")


def format_row(pair, aligned, ts, tv, P, Q, k2p, status):
    if status == "NO_ALIGNMENTS_PASSED_FILTERS":
        return f"{pair}\t0\t0\t0\tNA\tNA\tNA\t{status}"
    if status == "K2P_UNDEFINED_SATURATION":
        return f"{pair}\t{aligned}\t{ts}\t{tv}\t{P:.6f}\t{Q:.6f}\tNA\t{status}"
    return (f"{pair}\t{aligned}\t{ts}\t{tv}\t"
            f"{P:.6f}\t{Q:.6f}\t{k2p:.6f}\t{status}")


def main():
    ap = argparse.ArgumentParser(
        description="Genome-wide pairwise K2P divergence from cs-tagged PAF files.")
    ap.add_argument("paf_dir", nargs="?", default=DEFAULT_PAF_DIR,
                    help="Directory containing *.srt.paf files")
    ap.add_argument("output_file", nargs="?", default=DEFAULT_OUT_FILE,
                    help="Output TSV path")
    ap.add_argument("--min-mapq", type=int, default=DEFAULT_MIN_MAPQ,
                    help=f"Minimum mapping quality (default {DEFAULT_MIN_MAPQ})")
    ap.add_argument("--min-alen", type=int, default=DEFAULT_MIN_ALEN,
                    help=f"Minimum alignment block length (default {DEFAULT_MIN_ALEN})")
    args = ap.parse_args()

    if not os.path.isdir(args.paf_dir):
        sys.exit(f"Error: directory not found: {args.paf_dir}")

    paf_files = sorted(glob.glob(os.path.join(args.paf_dir, "*.srt.paf")))
    if not paf_files:
        sys.exit(f"Error: no *.srt.paf files found in {args.paf_dir}")

    print(f"Found {len(paf_files)} PAF files in {args.paf_dir}")
    print(f"Filters: primary alignments | MAPQ >= {args.min_mapq} | "
          f"alignment length >= {args.min_alen} bp")
    print()

    n_ok = 0
    with open(args.output_file, "w") as out:
        out.write("pair\taligned_bases\ttransitions\ttransversions\t"
                  "P\tQ\tK2P\tstatus\n")

        for paf in paf_files:
            pair = os.path.basename(paf)[:-len(".srt.paf")]
            print(f"Processing {pair} ... ", end="")

            aligned, ts, tv, P, Q, k2p, status = process_paf(
                paf, args.min_mapq, args.min_alen)
            out.write(format_row(pair, aligned, ts, tv, P, Q, k2p, status) + "\n")

            if status == "OK":
                n_ok += 1
            k2p_str = f"{k2p:.6f}" if k2p is not None else "NA"
            print(f"K2P = {k2p_str}  [{status}]")

    total = len(paf_files)
    failed = total - n_ok
    print()
    print(f"Done. Results written to: {args.output_file}")
    print(f"  Pairs processed : {total}")
    print(f"  Successful (OK) : {n_ok}")
    print(f"  Failed/flagged  : {failed}")


if __name__ == "__main__":
    main()
