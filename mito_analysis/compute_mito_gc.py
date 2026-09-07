#!/usr/bin/env python3
# compute mitochondrial GC content

import argparse, csv, glob, os, sys

def read_seq(path):
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()

def gc_percent(seq):
    a = seq.count("A"); c = seq.count("C"); g = seq.count("G"); t = seq.count("T")
    n = a + c + g + t
    return (100.0 * (g + c) / n) if n else None, n

def species_to_accession(divergence_csv):
    m = {}
    with open(divergence_csv) as f:
        for r in csv.DictReader(f):
            for sp_col, acc_col in (("species1", "genome1"), ("species2", "genome2")):
                sp = r.get(sp_col, "").strip()
                acc = r.get(acc_col, "").strip()
                if sp and acc:
                    m.setdefault(sp, acc)
    return m

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta-dir", required=True,
                    help="Directory of <Species_name>.fasta mitogenomes")
    ap.add_argument("--divergence", required=True,
                    help="Complete_divergence_with_mito_*.csv")
    ap.add_argument("--out", default="mito_gc.csv")
    args = ap.parse_args()

    sp2acc = species_to_accession(args.divergence)
    files = sorted(glob.glob(os.path.join(args.fasta_dir, "*.fasta")))
    if not files:
        sys.exit(f"No .fasta files found")

    rows, no_acc = [], []
    for fp in files:
        species = os.path.splitext(os.path.basename(fp))[0].replace("_", " ")
        seq = read_seq(fp)
        gc, length = gc_percent(seq)
        if gc is None:
            continue
        acc = sp2acc.get(species)
        if not acc:
            no_acc.append(species)
            continue
        rows.append({"accession": acc, "species": species,
                     "mito_gc": round(gc, 4), "mito_length": length})

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["accession", "species", "mito_gc", "mito_length"])
        w.writeheader()
        w.writerows(rows)

    gcs = [r["mito_gc"] for r in rows]
    print(f"Wrote {len(rows)} species -> {args.out}")
    if gcs:
        gcs_sorted = sorted(gcs)
        print(f"mito GC%: min={gcs_sorted[0]:.1f}  median={gcs_sorted[len(gcs)//2]:.1f}  max={gcs_sorted[-1]:.1f}")
    if no_acc:
        print(f"NOTE: {len(no_acc)} FASTA species had no nuclear accession match)

if __name__ == "__main__":
    main()
