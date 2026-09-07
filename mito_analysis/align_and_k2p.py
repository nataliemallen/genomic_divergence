#!/usr/bin/env python3
# align 13 mitochondrial genes and compute k2p 
import argparse, csv, math, os, subprocess, sys, tempfile

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


CANON = ["ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3",
         "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB"]
GENE_SYN = {
    "COI": "COX1", "CO1": "COX1", "COXI": "COX1", "MTCO1": "COX1",
    "COII": "COX2", "CO2": "COX2", "COXII": "COX2", "MTCO2": "COX2",
    "COIII": "COX3", "CO3": "COX3", "COXIII": "COX3", "MTCO3": "COX3",
    "COB": "CYTB", "CYB": "CYTB", "MTCYB": "CYTB",
    "ATPASE6": "ATP6", "ATP6": "ATP6", "ATPASE8": "ATP8", "ATP8": "ATP8",
    "NAD1": "ND1", "NAD2": "ND2", "NAD3": "ND3", "NAD4": "ND4",
    "NAD4L": "ND4L", "NAD5": "ND5", "NAD6": "ND6",
    "NADH1": "ND1", "NADH2": "ND2", "NADH3": "ND3", "NADH4": "ND4",
    "NADH4L": "ND4L", "NADH5": "ND5", "NADH6": "ND6",
}

def norm_gene(raw):
    if not raw:
        return None
    g = "".join(ch for ch in raw.upper() if ch.isalnum())
    g = g.replace("MT", "", 1) if g.startswith("MT") and len(g) > 2 else g
    if g in CANON:
        return g
    if g in GENE_SYN:
        return GENE_SYN[g]
    # last resort: keyword-match a protein description
    r = raw.upper()
    if "CYTOCHROME B" in r or r.strip() in ("COB", "CYTB"):
        return "CYTB"
    for n, name in ((1, "COX1"), (2, "COX2"), (3, "COX3")):
        if f"OXIDASE SUBUNIT {n}" in r or f"OXIDASE SUBUNIT I{'I'*(n-1)}" in r:
            return name
    if "ATP SYNTHASE" in r and "8" in r:
        return "ATP8"
    if "ATP SYNTHASE" in r and "6" in r:
        return "ATP6"
    if "NADH" in r:
        for tag, name in (("4L", "ND4L"), ("1", "ND1"), ("2", "ND2"), ("3", "ND3"),
                          ("4", "ND4"), ("5", "ND5"), ("6", "ND6")):
            if f"SUBUNIT {tag}" in r:
                return name
    return None

def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()

def parse_cds_fasta(path):
    genes, cur_gene, cur_seq = {}, None, []

    def flush():
        if cur_gene and cur_seq:
            s = "".join(cur_seq).upper()
            if cur_gene not in genes or len(s) > len(genes[cur_gene]):
                genes[cur_gene] = s

    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                flush(); cur_seq = []
                gene = None
                for tag in ("gene", "protein"):
                    key = f"[{tag}="
                    if key in line:
                        val = line.split(key, 1)[1].split("]", 1)[0]
                        gene = norm_gene(val)
                        if gene:
                            break
                cur_gene = gene
            else:
                cur_seq.append(line.strip())
        flush()
    return genes

def revcomp(s):
    return s.translate(COMP)[::-1]

def rotate_to_anchor(ref, qry, k=21, tries=8):
    doubled_variants = {"fwd": qry + qry, "rev": (rc := revcomp(qry)) + rc}
    for start in range(0, min(len(ref) - k, tries * k), k):
        anchor = ref[start:start + k]
        if "N" in anchor:
            continue
        for strand, dbl in doubled_variants.items():
            pos = dbl.find(anchor)
            if 0 <= pos < len(qry):
                base = qry if strand == "fwd" else revcomp(qry)
                origin = (pos - start) % len(base)
                rot = base[origin:] + base[:origin]
                return rot, origin != 0, strand == "rev", True
    return qry, False, False, False

def mafft_align(seq1, seq2, name1, name2, threads, mafft_cmd, tmpdir):
    os.makedirs(tmpdir, exist_ok=True)
    fd, inp = tempfile.mkstemp(suffix=".fa", dir=tmpdir)
    with os.fdopen(fd, "w") as tf:
        tf.write(f">{name1}\n{seq1}\n>{name2}\n{seq2}\n")
    cmd = mafft_cmd + ["--auto", "--thread", str(threads), "--quiet", inp]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            err = (proc.stderr or "").strip().replace("\n", " | ")
            raise RuntimeError(err[:500] or f"exit {proc.returncode}, no stderr")
        out = proc.stdout
    finally:
        os.unlink(inp)
    seqs, cur = {}, None
    for line in out.splitlines():
        if line.startswith(">"):
            cur = line[1:].strip(); seqs[cur] = []
        elif cur is not None:
            seqs[cur].append(line.strip())
    a = "".join(seqs[name1]).upper(); b = "".join(seqs[name2]).upper()
    return a, b

def count_ts_tv(a, b):
    ts = tv = valid = 0
    purines = set("AG")
    for x, y in zip(a, b):
        if x in "ACGT" and y in "ACGT":
            valid += 1
            if x != y:
                ts += 1 if (x in purines) == (y in purines) else 0
                tv += 0 if (x in purines) == (y in purines) else 1
    return ts, tv, valid

def k2p_from_counts(ts, tv, valid):
    if valid == 0:
        return None, None
    p, q = ts / valid, tv / valid
    try:
        d = -0.5 * math.log(1 - 2 * p - q) - 0.25 * math.log(1 - 2 * q)
        if d < 0 or math.isinf(d) or math.isnan(d):
            d = None
    except ValueError:
        d = None  # 1-2p-q <= 0 : saturation, undefined
    return d, (ts + tv) / valid

def whole_k2p(seq1, seq2, k, threads, mafft_cmd, tmpdir):
    rot, rotated, rc, found = rotate_to_anchor(seq1, seq2, k=k)
    a, b = mafft_align(seq1, rot, "s1", "s2", threads, mafft_cmd, tmpdir)
    ts, tv, valid = count_ts_tv(a, b)
    d, p = k2p_from_counts(ts, tv, valid)
    return d, valid, ts, tv, p, int(rotated), int(rc), found

def coding_k2p(genes1, genes2, threads, mafft_cmd, tmpdir):
    shared = [g for g in CANON if g in genes1 and g in genes2]
    ts = tv = valid = 0
    used = []
    for g in shared:
        try:
            a, b = mafft_align(genes1[g], genes2[g], "a", "b", threads, mafft_cmd, tmpdir)
        except Exception:
            continue
        t, v, val = count_ts_tv(a, b)
        ts += t; tv += v; valid += val; used.append(g)
    d, p = k2p_from_counts(ts, tv, valid)
    return d, valid, ts, tv, p, len(used), used

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coverage", default="results/pair_coverage.csv")
    ap.add_argument("--fasta", default="results/fasta")
    ap.add_argument("--cds", default="results/cds")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--k", type=int, default=21, help="anchor k-mer for whole-mito rotation")
    ap.add_argument("--min-genes", type=int, default=8,
                    help="min shared PCGs required to trust a coding K2P (of 13)")
    ap.add_argument("--min-coding-sites", type=int, default=3000,
                    help="min aligned coding sites required to trust a coding K2P")
    ap.add_argument("--high-div-note", type=float, default=0.30,
                    help="coding K2P above this is NOTED (not dropped) as high/near-saturation")
    ap.add_argument("--mafft-cmd", default="mafft",
                    help='MAFFT command; quote if it has args (Singularity wrapper).')
    ap.add_argument("--tmpdir", default=None)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    mafft_cmd = args.mafft_cmd.split()
    tmpdir = args.tmpdir or os.path.join(args.outdir, "_aln_tmp")

    if args.selftest:
        os.makedirs(args.outdir, exist_ok=True)
        try:
            a, b = mafft_align("ACGTACGTACGTACGT", "ACGTACGTACGTACGT",
                               "x", "y", args.threads, mafft_cmd, tmpdir)
            print(f"MAFFT selftest OK — aligned length {len(a)}.")
        except Exception as e:
            print(f"MAFFT selftest FAILED:\n  {e}")
        return

    def fpath(base, sp):
        return os.path.join(base, sp.replace(" ", "_").replace("/", "_") + ".fasta")

    rows = list(csv.DictReader(open(args.coverage)))
    todo = [r for r in rows if r.get("both_available") == "1"]

    compare_cols = ["pair_id", "taxonomic_group", "species1", "species2",
                    "k2p_whole", "sites_whole", "k2p_coding", "sites_coding",
                    "n_genes", "ts_coding", "tv_coding", "p_dist_coding",
                    "rotated", "revcomp", "coding_available", "flag", "keep", "note"]
    down_cols = ["pair_id", "taxonomic_group", "species1", "species2", "k2p",
                 "aligned_sites", "ts", "tv", "p_dist", "n_genes", "note"]
    whole_cols = ["pair_id", "taxonomic_group", "species1", "species2", "k2p",
                  "aligned_sites", "ts", "tv", "p_dist", "rotated", "revcomp", "note"]

    cmp_f = open(os.path.join(args.outdir, "mito_k2p_compare.csv"), "w", newline="")
    dn_f = open(os.path.join(args.outdir, "mito_k2p_coding.csv"), "w", newline="")
    wh_f = open(os.path.join(args.outdir, "mito_k2p_whole.csv"), "w", newline="")
    cmp_w = csv.DictWriter(cmp_f, fieldnames=compare_cols); cmp_w.writeheader()
    dn_w = csv.DictWriter(dn_f, fieldnames=down_cols); dn_w.writeheader()
    wh_w = csv.DictWriter(wh_f, fieldnames=whole_cols); wh_w.writeheader()

    stats = {"whole": [], "coding": []}  # (group, k2p) for kept/defined values
    n_lost_no_cds = 0

    for i, r in enumerate(todo, 1):
        s1, s2, grp, pid = r["species1"], r["species2"], r["taxonomic_group"], r["pair_id"]
        sys.stderr.write(f"[{i}/{len(todo)}] {s1} vs {s2}\n")

        ### whole mitogenome
        kw = sw = tsw = tvw = pw = None; rotated = rc = 0; wnote = ""
        try:
            seq1, seq2 = read_fasta(fpath(args.fasta, s1)), read_fasta(fpath(args.fasta, s2))
            kw, sw, tsw, tvw, pw, rotated, rc, found = whole_k2p(
                seq1, seq2, args.k, args.threads, mafft_cmd, tmpdir)
            if not found:
                wnote = "no_anchor"
            if kw is None:
                wnote = (wnote + ";" if wnote else "") + "whole_saturated"
        except Exception as e:
            wnote = f"whole_error:{str(e)[:100]}"
        wh_w.writerow(dict(pair_id=pid, taxonomic_group=grp, species1=s1, species2=s2,
                           k2p=f"{kw:.6f}" if kw is not None else "", aligned_sites=sw or "",
                           ts=tsw or "", tv=tvw or "", p_dist=f"{pw:.6f}" if pw is not None else "",
                           rotated=rotated, revcomp=rc, note=wnote))
        if kw is not None:
            stats["whole"].append((grp, kw))

        ### coding (13 PCGs)
        coding_avail = r.get("both_cds_available", "0") == "1"
        kc = sc = tsc = tvc = pc = None; n_genes = 0; cnote = ""
        if coding_avail:
            try:
                g1 = parse_cds_fasta(fpath(args.cds, s1))
                g2 = parse_cds_fasta(fpath(args.cds, s2))
                kc, sc, tsc, tvc, pc, n_genes, _ = coding_k2p(
                    g1, g2, args.threads, mafft_cmd, tmpdir)
            except Exception as e:
                cnote = f"coding_error:{str(e)[:100]}"
        else:
            n_lost_no_cds += 1

        ### quality flag 
        flag = ""
        if not coding_avail:
            flag = "no_cds"
        elif cnote.startswith("coding_error"):
            flag = "coding_error"
        elif kc is None:
            flag = "saturated"  # K2P undefined -> unreliable
        elif n_genes < args.min_genes or (sc or 0) < args.min_coding_sites:
            flag = "low_coverage"  # too little sequence to trust
        keep = 1 if flag == "" else 0
        if keep and kc is not None and kc > args.high_div_note:
            cnote = (cnote + ";" if cnote else "") + f"high_div({kc:.2f})"  # kept, just noted

        cmp_w.writerow(dict(
            pair_id=pid, taxonomic_group=grp, species1=s1, species2=s2,
            k2p_whole=f"{kw:.6f}" if kw is not None else "",
            sites_whole=sw or "",
            k2p_coding=f"{kc:.6f}" if kc is not None else "",
            sites_coding=sc or "", n_genes=n_genes,
            ts_coding=tsc or "", tv_coding=tvc or "",
            p_dist_coding=f"{pc:.6f}" if pc is not None else "",
            rotated=rotated, revcomp=rc, coding_available=int(coding_avail),
            flag=flag, keep=keep, note=(wnote + (";" if wnote and cnote else "") + cnote)))

        if keep:
            dn_w.writerow(dict(pair_id=pid, taxonomic_group=grp, species1=s1, species2=s2,
                               k2p=f"{kc:.6f}", aligned_sites=sc, ts=tsc, tv=tvc,
                               p_dist=f"{pc:.6f}", n_genes=n_genes, note=cnote))
            stats["coding"].append((grp, kc))

    for fh in (cmp_f, dn_f, wh_f):
        fh.close()

    write_summary(args.outdir, stats, len(todo), n_lost_no_cds, args)
    print(f"Done. Downstream file -> {os.path.join(args.outdir, 'mito_k2p_coding.csv')}")

def write_summary(outdir, stats, n_pairs, n_lost_no_cds, args):
    import statistics as st
    from collections import defaultdict

    def block(name, data):
        lines = [f"[{name}]  n = {len(data)}"]
        if not data:
            return lines + ["  (none)"]
        vals = [v for _, v in data]
        lines.append(f"  overall: min {min(vals):.4f}  median {st.median(vals):.4f}  "
                     f"p95 {sorted(vals)[int(0.95*len(vals))-1]:.4f}  max {max(vals):.4f}")
        by = defaultdict(list)
        for g, v in data:
            by[g].append(v)
        for g in sorted(by):
            gv = by[g]
            lines.append(f"    {g:9s} n={len(gv):4d}  median {st.median(gv):.4f}  max {max(gv):.4f}")
        return lines

    wd = {(g, round(v, 9)) for g, v in stats["whole"]}  
    compare_path = os.path.join(outdir, "mito_k2p_compare.csv")
    xs, ys = [], []
    ratios = []
    with open(compare_path) as f:
        for r in csv.DictReader(f):
            try:
                w = float(r["k2p_whole"]); c = float(r["k2p_coding"])
            except (ValueError, KeyError):
                continue
            if r["keep"] == "1":
                xs.append(w); ys.append(c);
                if w > 0:
                    ratios.append(c / w)

    def pearson(x, y):
        if len(x) < 3:
            return float("nan")
        mx, my = sum(x)/len(x), sum(y)/len(y)
        num = sum((a-mx)*(b-my) for a, b in zip(x, y))
        dx = math.sqrt(sum((a-mx)**2 for a in x)); dy = math.sqrt(sum((b-my)**2 for b in y))
        return num/(dx*dy) if dx and dy else float("nan")

    out = []
    path = os.path.join(outdir, "mito_k2p_summary.txt")
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n")
    print("\n".join(out))

if __name__ == "__main__":
    main()
