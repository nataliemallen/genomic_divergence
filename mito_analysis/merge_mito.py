#!/usr/bin/env python3
# add mito k2p to divergence table
import argparse, csv

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--divergence", required=True,
                    help="nuclear divergence table (must have an ID / pair_id column)")
    ap.add_argument("--compare", default="results/mito_k2p_compare.csv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--id-col", default=None,
                    help="pair-id column in the divergence table (auto: ID or pair_id)")
    args = ap.parse_args()

    # index the mito results by pair_id
    mito = {}
    with open(args.compare) as f:
        for r in csv.DictReader(f):
            mito[str(r["pair_id"])] = r

    div = list(csv.DictReader(open(args.divergence)))
    cols = list(div[0].keys())
    id_col = args.id_col or ("ID" if "ID" in cols else ("pair_id" if "pair_id" in cols else None))
    if id_col is None:
        raise SystemExit("Could not find an ID/pair_id column in the divergence table; "
                         "pass --id-col.")

    add = ["mito_k2p", "mito_k2p_whole", "mito_k2p_flag",
           "mito_query_align_pct", "mito_ref_align_pct"]
    out_cols = cols + [c for c in add if c not in cols]

    n_kept = n_flagged = n_nomito = 0
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_cols); w.writeheader()
        for row in div:
            m = mito.get(str(row[id_col]))
            if m is None:
                row.update(mito_k2p="", mito_k2p_whole="", mito_k2p_flag="",
                           mito_query_align_pct="", mito_ref_align_pct="")
                n_nomito += 1
            elif m["keep"] == "1":
                row.update(mito_k2p=m["k2p_coding"], mito_k2p_whole=m["k2p_whole"],
                           mito_k2p_flag="", mito_query_align_pct="100",
                           mito_ref_align_pct="100")
                n_kept += 1
            else:
                row.update(mito_k2p="", mito_k2p_whole=m["k2p_whole"],
                           mito_k2p_flag=m["flag"] or "dropped",
                           mito_query_align_pct="", mito_ref_align_pct="")
                n_flagged += 1
            w.writerow(row)

    print(f"Wrote {args.out}\n  kept (coding): {n_kept}\n  dropped by QC: {n_flagged}"
          f"\n  no mito data:  {n_nomito}")

if __name__ == "__main__":
    main()
