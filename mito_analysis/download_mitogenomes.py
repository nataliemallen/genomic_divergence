#!/usr/bin/env python3
# used to download mitochondrial genomes from NCBI

import argparse, csv, json, os, sys, time, urllib.parse, urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MIN_LEN, MAX_LEN = 14000, 20000  # plausible vertebrate mitogenome size window

def _params(extra):
    p = {"email": os.environ.get("NCBI_EMAIL", ""),
         "tool": "genomic_divergence_mito"}
    key = os.environ.get("NCBI_API_KEY")
    if key:
        p["api_key"] = key
    p.update(extra)
    return p

def _get(endpoint, extra, retries=4):
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(_params(extra))
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=40) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            wait = 2 ** attempt
            sys.stderr.write(f"  [retry {attempt+1}] {endpoint}: {e} (sleep {wait}s)\n")
            time.sleep(wait)
    return None

def resolve_taxid(species):
    txt = _get("esearch.fcgi", {"db": "taxonomy",
                                "term": f'"{species}"[Scientific Name]', "retmode": "json"})
    if not txt:
        return None
    try:
        ids = json.loads(txt)["esearchresult"]["idlist"]
        return ids[0] if ids else None
    except Exception:
        return None

def search_mito(species, taxid):
    org = f"txid{taxid}[Organism:exp]" if taxid else f'"{species}"[Organism]'
    term = (f'{org} AND mitochondrion[filter] AND '
            f'("complete genome"[Title] OR "complete sequence"[Title]) AND '
            f'{MIN_LEN}:{MAX_LEN}[SLEN]')
    hits = []
    for refseq_only in (True, False):
        t = term + (' AND refseq[filter]' if refseq_only else '')
        txt = _get("esearch.fcgi", {"db": "nuccore", "term": t, "retmode": "json",
                                    "retmax": "50", "sort": "slen"})
        if txt:
            try:
                for i in json.loads(txt)["esearchresult"]["idlist"]:
                    if i not in hits:
                        hits.append(i)
            except Exception:
                pass
        if refseq_only and hits:
            return hits, "RefSeq"
    return hits, ("GenBank" if hits else None)

def summarize(uids):
    if not uids:
        return []
    txt = _get("esummary.fcgi", {"db": "nuccore", "id": ",".join(uids), "retmode": "json"})
    out = []
    if not txt:
        return out
    try:
        res = json.loads(txt)["result"]
        for uid in res.get("uids", []):
            d = res[uid]
            acc = d.get("accessionversion", d.get("caption", ""))
            length = int(d.get("slen", 0) or 0)
            out.append((uid, acc, length, acc.startswith("NC_")))
    except Exception:
        pass
    return out

def fetch_fasta(uid):
    return _get("efetch.fcgi", {"db": "nuccore", "id": uid,
                                "rettype": "fasta", "retmode": "text"})

def fetch_cds(uid):
    return _get("efetch.fcgi", {"db": "nuccore", "id": uid,
                                "rettype": "fasta_cds_na", "retmode": "text"})

def count_cds(text):
    return text.count(">") if (text and text.lstrip().startswith(">")) else 0

def load_done(path):
    done = {}
    if os.path.exists(path):
        with open(path) as f:
            for row in csv.DictReader(f):
                done[row["species"]] = row
    return done

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="inputs/species_list.csv")
    ap.add_argument("--pairs", default="inputs/pairs.csv")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--sleep", type=float, default=None,
                    help="Seconds between requests. Default 0.11 with API key, 0.34 without.")
    args = ap.parse_args()

    if not os.environ.get("NCBI_EMAIL"):
        sys.exit("ERROR: set NCBI_EMAIL environment variable")
    sleep = args.sleep if args.sleep is not None else (0.11 if os.environ.get("NCBI_API_KEY") else 0.34)

    fasta_dir = os.path.join(args.outdir, "fasta")
    cds_dir = os.path.join(args.outdir, "cds")
    os.makedirs(fasta_dir, exist_ok=True)
    os.makedirs(cds_dir, exist_ok=True)
    acc_path = os.path.join(args.outdir, "mito_accessions.csv")

    species = list(csv.DictReader(open(args.species)))
    records = dict(load_done(acc_path))
    fieldnames = ["species", "taxonomic_group", "taxid", "accession", "length",
                  "source", "status", "n_cds", "cds_status"]

    for n, row in enumerate(species, 1):
        sp = row["species"].strip()
        grp = row["taxonomic_group"]
        prev = records.get(sp)
        if prev and prev.get("status") == "ok":
            if prev.get("cds_status"):
                continue  # fully done (genome + CDS info)
            # backfill CDS for species downloaded before annotation support existed
            acc = prev.get("accession")
            if acc:
                sys.stderr.write(f"[{n}/{len(species)}] {sp} (CDS backfill) ...\n")
                cds = fetch_cds(acc); time.sleep(sleep)
                n_cds = count_cds(cds)
                if n_cds > 0:
                    safe = sp.replace(" ", "_").replace("/", "_")
                    with open(os.path.join(cds_dir, f"{safe}.fasta"), "w") as f:
                        f.write(cds)
                    prev.update(n_cds=n_cds, cds_status="ok")
                else:
                    prev.update(n_cds=0, cds_status="no_annotation")
                records[sp] = prev
                if n % 25 == 0:
                    _write_acc(acc_path, fieldnames, records, species)
            continue
        sys.stderr.write(f"[{n}/{len(species)}] {sp} ...\n")
        rec = {"species": sp, "taxonomic_group": grp, "taxid": "",
               "accession": "", "length": "", "source": "", "status": "not_found",
               "n_cds": "", "cds_status": ""}
        taxid = resolve_taxid(sp); time.sleep(sleep)
        rec["taxid"] = taxid or ""
        uids, _ = search_mito(sp, taxid); time.sleep(sleep)
        summ = [s for s in summarize(uids) if MIN_LEN <= s[2] <= MAX_LEN]; time.sleep(sleep)
        summ.sort(key=lambda s: (s[3], s[2]), reverse=True)  # RefSeq first, then longest
        if summ:
            uid, acc, length, is_ref = summ[0]
            fa = fetch_fasta(uid); time.sleep(sleep)
            if fa and fa.startswith(">"):
                safe = sp.replace(" ", "_").replace("/", "_")
                with open(os.path.join(fasta_dir, f"{safe}.fasta"), "w") as f:
                    f.write(fa)
                rec.update(accession=acc, length=length,
                           source=("RefSeq" if is_ref else "GenBank"), status="ok")
                # also fetch the annotated CDS (coding-only, literature-comparable)
                cds = fetch_cds(uid); time.sleep(sleep)
                n_cds = count_cds(cds)
                if n_cds > 0:
                    with open(os.path.join(cds_dir, f"{safe}.fasta"), "w") as f:
                        f.write(cds)
                    rec.update(n_cds=n_cds, cds_status="ok")
                else:
                    rec.update(n_cds=0, cds_status="no_annotation")
        records[sp] = rec
        if n % 25 == 0:
            _write_acc(acc_path, fieldnames, records, species)

    _write_acc(acc_path, fieldnames, records, species)
    write_coverage(args.pairs, records, args.outdir)
    print(f"\nDone. Accessions -> {acc_path}")

def _write_acc(path, fieldnames, records, species):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames); w.writeheader()
        for row in species:
            sp = row["species"].strip()
            if sp in records:
                w.writerow({k: records[sp].get(k, "") for k in fieldnames})

def write_coverage(pairs_path, records, outdir):
    from collections import Counter
    pairs = list(csv.DictReader(open(pairs_path)))
    cov_path = os.path.join(outdir, "pair_coverage.csv")
    tot, both, both_cds = Counter(), Counter(), Counter()

    def has_cds(sp):
        return str(records.get(sp, {}).get("cds_status", "")) == "ok"

    with open(cov_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pair_id", "taxonomic_group", "species1", "species2",
                    "sp1_available", "sp2_available", "both_available",
                    "sp1_cds", "sp2_cds", "both_cds_available", "missing"])
        for p in pairs:
            s1, s2, grp = p["species1"].strip(), p["species2"].strip(), p["taxonomic_group"]
            a1 = records.get(s1, {}).get("status") == "ok"
            a2 = records.get(s2, {}).get("status") == "ok"
            c1, c2 = has_cds(s1), has_cds(s2)
            miss = ";".join([s for s, a in ((s1, a1), (s2, a2)) if not a])
            w.writerow([p["pair_id"], grp, s1, s2, int(a1), int(a2), int(a1 and a2),
                        int(c1), int(c2), int(c1 and c2), miss])
            tot[grp] += 1
            if a1 and a2:
                both[grp] += 1
            if c1 and c2:
                both_cds[grp] += 1
    with open(os.path.join(outdir, "coverage_summary.txt"), "w") as f:
        f.write("Per-pair mitogenome coverage\n")
        f.write("=" * 60 + "\n")
        f.write(f"{'GROUP':12} {'whole-mito (both)':>20}   {'coding/CDS (both)':>20}\n")
        for g in sorted(tot):
            f.write(f"{g:12} {both[g]:>8}/{tot[g]:<5} ({100*both[g]/tot[g]:4.1f}%)   "
                    f"{both_cds[g]:>8}/{tot[g]:<5} ({100*both_cds[g]/tot[g]:4.1f}%)\n")
        tb, tc, tt = sum(both.values()), sum(both_cds.values()), sum(tot.values())
        f.write(f"{'TOTAL':12} {tb:>8}/{tt:<5} ({100*tb/tt:4.1f}%)   "
                f"{tc:>8}/{tt:<5} ({100*tc/tt:4.1f}%)\n")
        f.write(f"\nPairs lost by requiring annotation (both_available minus both_cds): "
                f"{tb - tc}\n")
    print(open(os.path.join(outdir, "coverage_summary.txt")).read())

if __name__ == "__main__":
    main()
