#!/usr/bin/env python3
"""Compare WBPS genome V4 amplicons against the PR2 amplicons.

    bin/compare_wbps_pr2.py <pr2 pcr_results.json> <wbps outdir> <outdir>

pcr_results.json is pr2_insilico_pcr.py's cache. "Nearest" for a variant
not in PR2 is the PR2 sequence sharing the most 12-mers, and
kmer_containment is the fraction of the variant's 12-mers it shares: a
rough proxy for identity, not an alignment.
"""
import glob, json, os, sys
from collections import Counter, defaultdict

PR2_JSON, WBPS, OUT = sys.argv[1:4]
K = 12
W1, W2 = 132, 128


def kmers(s):
    return {s[i:i + K] for i in range(len(s) - K + 1)}


def ill(s):
    return s if len(s) <= W1 + W2 else (s[:W1], s[-W2:])


res = json.load(open(PR2_JSON))
pr2 = [(r[0] + [""] * 9, r[2]) for r in res if r[1] == "amplicon"]
by_ins = defaultdict(list)
by_ill = defaultdict(list)
for t, s in pr2:
    by_ins[s].append(t)
    by_ill[ill(s)].append(t)
pr2_species = {t[8].lower() for t, _ in pr2}
pr2_all_species = {(r[0] + [""] * 9)[8].lower() for r in res}
pr2_genera = {(r[0] + [""] * 9)[7].lower() for r in res}
# k-mer index over unique PR2 inserts for nearest-neighbour lookups
uniq = list(by_ins)
index = defaultdict(list)
for j, s in enumerate(uniq):
    for k in kmers(s):
        index[k].append(j)


def nearest(q):
    qk = kmers(q)
    hits = Counter()
    for k in qk:
        for j in index.get(k, ()):
            hits[j] += 1
    if not hits:
        return None, 0.0
    j, n = hits.most_common(1)[0]
    # containment of the query's k-mers: a rough identity proxy
    return j, n / max(1, len(qk))


rows_out = []
genomes = sorted(glob.glob(os.path.join(WBPS, "*.done")))
failed = sorted(glob.glob(os.path.join(WBPS, "*.failed")))
summary = []
for d in genomes:
    stem = os.path.basename(d)[:-5]
    sp_dir = stem.split(".")[0]                     # e.g. haemonchus_contortus
    genus = sp_dir.split("_")[0]
    species = sp_dir.capitalize()
    bases, n, secs = open(d).read().split()
    tsv = os.path.join(WBPS, stem + ".tsv")
    inserts = []
    if os.path.exists(tsv):
        for line in open(tsv).read().splitlines()[1:]:
            inserts.append(line.split("\t")[-1])
    variants = Counter(inserts)
    rec = {"genome": stem, "species": species, "bases_mb": int(bases) / 1e6,
           "copies": len(inserts), "variants": len(variants),
           "species_in_pr2_any": sp_dir in pr2_all_species,
           "species_amplicon_in_pr2": sp_dir in pr2_species,
           "genus_in_pr2": genus in pr2_genera, "calls": []}
    for ins, cnt in variants.most_common():
        exact = by_ins.get(ins, [])
        if exact:
            labs = Counter(f"{t[4]}|{t[7]}|{t[8]}" for t in exact)
            j, idn = None, 1.0
            nearest_lab = labs.most_common(1)[0][0]
        else:
            j, idn = nearest(ins)
            nearest_lab = (Counter(f"{t[4]}|{t[7]}|{t[8]}" for t in by_ins[uniq[j]])
                           .most_common(1)[0][0] if j is not None else "")
        same_ill = by_ill.get(ill(ins), [])
        rec["calls"].append({
            "count": cnt, "len": len(ins), "exact_in_pr2": bool(exact),
            "exact_labels": sorted({f"{t[4]}|{t[7]}|{t[8]}" for t in exact})[:5],
            "nearest": nearest_lab, "kmer_containment": round(idn, 3),
            "illumina_window_shared_with": sorted({t[8] for t in same_ill})[:6],
            "insert": ins})
    summary.append(rec)

os.makedirs(OUT, exist_ok=True)
json.dump(summary, open(os.path.join(OUT, "wbps_vs_pr2.json"), "w"), indent=1)

# Human-readable digest
L = []
P = L.append
n_g = len(summary)
with_amp = [r for r in summary if r["copies"]]
P(f"Genomes scanned: {n_g}  (failed downloads: {len(failed)})")
P(f"Genomes yielding >=1 V4 amplicon: {len(with_amp)} ({100*len(with_amp)/max(1,n_g):.0f}%)")
sp_all = {r["species"] for r in summary}
sp_amp = {r["species"] for r in with_amp}
P(f"Species: {len(sp_all)}; with an amplicon from at least one assembly: {len(sp_amp)}")
new = {r["species"] for r in with_amp if not r["species_amplicon_in_pr2"]}
P(f"Species with a genome amplicon but NO amplicon-spanning PR2 entry: {len(new)}")
P("")
# Does the genome's main variant agree with PR2 for species PR2 does have?
agree = disagree = 0
dis = []
for r in with_amp:
    if not r["species_amplicon_in_pr2"]:
        continue
    top = r["calls"][0]
    labs = top["exact_labels"] or [top["nearest"]]
    sp_hit = any(l.split("|")[-1].lower() == r["species"].lower() for l in labs)
    if sp_hit:
        agree += 1
    else:
        disagree += 1
        dis.append((r["species"], top["exact_in_pr2"], top["kmer_containment"], labs[:3]))
P(f"Species in both: main genome variant matches PR2's own label for that species: {agree}; does not: {disagree}")
for s, ex, idn, labs in dis[:25]:
    P(f"  {s:36s} {'exact' if ex else f'~{idn:.2f}':>6} -> {'; '.join(labs)}")
P("")
# Possible contaminants: genome variant whose nearest PR2 hit is a different phylum-level group
P("Possible contaminant hits (genome is a helminth, variant nearest a non-helminth group):")
helm = {"Nematoda", "Platyhelminthes", "Acanthocephala"}
for r in with_amp:
    for c in r["calls"]:
        grp = (c["exact_labels"][0] if c["exact_labels"] else c["nearest"]).split("|")[0]
        if grp and grp not in helm:
            P(f"  {r['genome']:55s} x{c['count']:<3} {c['len']}bp -> {c['exact_labels'][:1] or c['nearest']} ({c['kmer_containment']})")
P("")
P("Gap taxa, genome result:")
for g in ["teladorsagia", "cooperia", "oesophagostomum", "trichostrongylus", "haemonchus",
          "necator", "ancylostoma", "ascaris", "trichuris", "strongyloides", "toxocara",
          "schistosoma", "fasciola", "taenia", "echinococcus", "brugia", "onchocerca"]:
    for r in summary:
        if r["species"].lower().startswith(g):
            top = r["calls"][0] if r["calls"] else None
            desc = ("no amplicon" if not top else
                    f"{r['copies']} copies, {r['variants']} variants, {top['len']}bp, "
                    + ("exact PR2: " + ";".join(x.split('|')[-1] for x in top["exact_labels"][:3])
                       if top["exact_in_pr2"] else f"novel, nearest {top['nearest'].split('|')[-1]} ({top['kmer_containment']})"))
            P(f"  {r['genome']:55s} {desc}")
P("")
nocopy = sorted(r["genome"] for r in summary if not r["copies"])
P(f"Genomes with no amplicon ({len(nocopy)}): " + ", ".join(nocopy))
text = "\n".join(L)
open(os.path.join(OUT, "wbps_vs_pr2.txt"), "w").write(text + "\n")
print(text)
