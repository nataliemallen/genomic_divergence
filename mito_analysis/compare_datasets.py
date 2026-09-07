#!/usr/bin/env python3
# used to compare nuclear and mitochondrial yardstick datasets

import argparse, csv
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CLASSES = ["Aves", "Mammalia", "Reptilia", "Amphibia"]
C_NUC, C_MITO, C_BOTH = "#2C7FB8", "#D95F0E", "#6A51A3"
CLASS_COLORS = {"Aves": "#4C9BD4", "Mammalia": "#C0392B",
                "Reptilia": "#27AE60", "Amphibia": "#8E7CC3"}

def fnum(x):
    try: return float(x)
    except: return None

def in_nuclear(r):
    q, ref, k, t = (fnum(r["query_alignment_percent"]), fnum(r["ref_alignment_percent"]),
                    fnum(r["k2p"]), fnum(r["timetree_div"]))
    return (None not in (q, ref, k, t) and q >= 80 and ref >= 80 and k > 0 and t > 0
            and r["taxonomic_group"] in CLASSES)

def in_mito(r):
    q, ref, k, t = (fnum(r["mito_query_align_pct"]), fnum(r["mito_ref_align_pct"]),
                    fnum(r["mito_k2p"]), fnum(r["timetree_div"]))
    return (r["mito_k2p_flag"] != "artifact" and None not in (q, ref, k, t)
            and q >= 80 and ref >= 80 and k > 0 and t > 0 and r["taxonomic_group"] in CLASSES)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--divergence", default="Complete_divergence_with_mito_05-14-26.csv")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.divergence)))
    counts = defaultdict(lambda: dict(nuc_only=0, both=0, mito_only=0))
    comp = {"nuclear": defaultdict(int), "mito": defaultdict(int)}
    tdiv = {"nuclear": [], "mito": [], "both": []}

    with open(f"{args.outdir}/dataset_membership.csv", "w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["pair_id", "taxonomic_group", "species1", "species2", "timetree_div",
                    "in_nuclear_yardstick", "in_mito_yardstick", "membership"])
        for r in rows:
            n, m = in_nuclear(r), in_mito(r)
            if not (n or m):
                continue
            g = r["taxonomic_group"]; t = fnum(r["timetree_div"])
            mem = "both" if (n and m) else ("nuclear_only" if n else "mito_only")
            counts[g]["both" if (n and m) else ("nuc_only" if n else "mito_only")] += 1
            if n: comp["nuclear"][g] += 1; tdiv["nuclear"].append(t)
            if m: comp["mito"][g] += 1;    tdiv["mito"].append(t)
            if n and m: tdiv["both"].append(t)
            w.writerow([r["ID"], g, r["species1"], r["species2"], r["timetree_div"],
                        int(n), int(m), mem])

    tot = dict(nuc_only=sum(c["nuc_only"] for c in counts.values()),
               both=sum(c["both"] for c in counts.values()),
               mito_only=sum(c["mito_only"] for c in counts.values()))

    ### figure
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle("Nuclear vs mitochondrial yardstick datasets: pair overlap and composition",
                 fontsize=14, fontweight="bold")

    a = ax[0, 0]
    y = np.arange(len(CLASSES))[::-1]
    no = [counts[g]["nuc_only"] for g in CLASSES]
    bo = [counts[g]["both"] for g in CLASSES]
    mo = [counts[g]["mito_only"] for g in CLASSES]
    a.barh(y, no, color=C_NUC, label="nuclear only")
    a.barh(y, bo, left=no, color=C_BOTH, label="both")
    a.barh(y, mo, left=np.array(no) + np.array(bo), color=C_MITO, label="mito only")
    for i, g in enumerate(CLASSES):
        yy = y[i]; row = counts[g]
        segs = [(row["nuc_only"], 0), (row["both"], row["nuc_only"]),
                (row["mito_only"], row["nuc_only"] + row["both"])]
        for val, left in segs:
            if val > 0:
                a.text(left + val / 2, yy, str(val), va="center", ha="center",
                       fontsize=8, color="white", fontweight="bold")
    a.set_yticks(y); a.set_yticklabels(CLASSES)
    a.set_xlabel("number of species pairs")
    a.set_title("A. Overlap by taxonomic class", loc="left", fontweight="bold")
    a.legend(fontsize=8, loc="lower right", frameon=False)

    b = ax[0, 1]
    for i, ds in enumerate(["nuclear", "mito"]):
        bottom = 0
        for g in CLASSES:
            v = comp[ds][g]
            b.bar(i, v, bottom=bottom, color=CLASS_COLORS[g], width=0.6,
                  label=g if i == 0 else None)
            if v > 0:
                b.text(i, bottom + v / 2, f"{g[:4]}\n{v}", ha="center", va="center",
                       fontsize=7.5, color="white", fontweight="bold")
            bottom += v
    b.set_xticks([0, 1]); b.set_xticklabels(["nuclear\nyardstick", "mito\nyardstick"])
    b.set_ylabel("number of species pairs")
    b.set_title("B. Composition differs (nuclear = bird-heavy, mito = mammal-heavy)",
                loc="left", fontweight="bold", fontsize=9.5)

    c = ax[1, 0]
    bins = np.logspace(np.log10(0.05), np.log10(60), 26)
    c.hist(tdiv["nuclear"], bins=bins, alpha=0.55, color=C_NUC, label=f"nuclear (n={len(tdiv['nuclear'])})")
    c.hist(tdiv["mito"], bins=bins, alpha=0.55, color=C_MITO, label=f"mito (n={len(tdiv['mito'])})")
    c.set_xscale("log")
    c.set_xlabel("divergence time (My, log scale)")
    c.set_ylabel("number of pairs")
    c.set_title("C. Both datasets span the same divergence-time range", loc="left", fontweight="bold", fontsize=9.5)
    c.legend(fontsize=8, frameon=False)
    for ds, col in [("nuclear", C_NUC), ("mito", C_MITO)]:
        c.axvline(np.median(tdiv[ds]), color=col, ls="--", lw=1.2)

    d = ax[1, 1]
    d.barh([0], [tot["nuc_only"]], color=C_NUC)
    d.barh([0], [tot["both"]], left=tot["nuc_only"], color=C_BOTH)
    d.barh([0], [tot["mito_only"]], left=tot["nuc_only"] + tot["both"], color=C_MITO)
    jac = tot["both"] / (tot["nuc_only"] + tot["both"] + tot["mito_only"])
    n_nuc = tot["nuc_only"] + tot["both"]; n_mito = tot["mito_only"] + tot["both"]
    d.text(tot["nuc_only"] / 2, 0, f"nuclear only\n{tot['nuc_only']}", ha="center", va="center",
           color="white", fontweight="bold", fontsize=9)
    d.text(tot["nuc_only"] + tot["both"] / 2, 0, f"both\n{tot['both']}", ha="center", va="center",
           color="white", fontweight="bold", fontsize=9)
    d.text(tot["nuc_only"] + tot["both"] + tot["mito_only"] / 2, 0, f"mito only\n{tot['mito_only']}",
           ha="center", va="center", color="white", fontweight="bold", fontsize=9)
    d.set_yticks([])
    d.set_xlabel("number of species pairs")
    d.set_title("D. Overall overlap", loc="left", fontweight="bold")
    d.text(0.5, -0.9,
           f"nuclear yardstick: {n_nuc} pairs   |   mito yardstick: {n_mito} pairs\n"
           f"shared: {tot['both']}   →   Jaccard overlap = {jac:.2f}   "
           f"({100*tot['both']/n_nuc:.0f}% of nuclear, {100*tot['both']/n_mito:.0f}% of mito)",
           transform=d.transAxes, ha="center", va="center", fontsize=9)
    d.set_ylim(-1.2, 0.6)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out_png = f"{args.outdir}/yardstick_dataset_comparison.png"
    plt.savefig(out_png, dpi=150, bbox_inches="tight")

if __name__ == "__main__":
    main()
