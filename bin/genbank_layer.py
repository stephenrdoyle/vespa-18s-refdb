#!/usr/bin/env python3
"""Layer: GenBank 18S records PR2 does not already hold, trimmed to V4.

    bin/genbank_layer.py work/pr2_layer.tsv.gz work/genbank_layer.tsv.gz \
        work/genbank_review.tsv <genbank_fetch.py fasta.gz> ...

PR2 is curated from GenBank, so most GenBank 18S is already in it; this
layer is what is not. Each record is
  - dropped if its GenBank accession is already a PR2 entry;
  - put through the in-silico PCR (and its reverse complement, since
    GenBank does not guarantee orientation); dropped unless it spans V4;
  - dropped if its insert is more than 1% N, or its title says
    UNVERIFIED;
  - given PR2 ranks from its NCBI taxid (vespa.PR2Index.ranks_for);
  - sent to review instead of the layer if its insert's nearest PR2
    sequence is in another PR2 class - the same contamination check the
    genome layer uses.

NCBI names become PR2-style labels: "Genus species" -> Genus_species,
anything less ("Genus sp. isolate 7") -> Genus_sp., a placeholder that
names the genus and nothing below.

TRUST. This layer is the least trusted: PR2 is GenBank curated, so a raw
GenBank label must not overrule PR2's. Pass it to collapse.py first.
"""

from __future__ import annotations

import gzip
import json
import os
import re
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vespa import RANKS, PR2Index, pcr  # noqa: E402

SOURCE = "GenBank"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
_RC = str.maketrans("ACGTRYSWKMBDHVNacgtn", "TGCAYRSWMKVHDBNtgcan")


def read_fetch(path: str):
    with gzip.open(path, "rt") as fh:
        head = None
        for line in fh:
            if line.startswith(">"):
                head = line[1:].rstrip("\n")
            elif head:
                acc, tx, title = (head.split(" ", 2) + ["", ""])[:3]
                yield acc, tx.replace("taxid=", ""), title, line.strip()


def _pcr_both(rec):
    acc, tx, title, seq = rec
    status, ins, fc, rc = pcr(seq)
    if status != "amplicon":
        s2, i2, f2, r2 = pcr(seq.translate(_RC)[::-1])
        if s2 == "amplicon":
            status, ins, fc, rc = s2, i2, f2, r2
    return acc, tx, title, status, ins, fc, rc


def lineages(taxids, cache_path: str) -> dict:
    cache = json.load(open(cache_path)) if os.path.exists(cache_path) else {}
    # Entries cached before ranks were kept are fetched again.
    todo = sorted({t for t in taxids if t and "rank" not in cache.get(t, {})})
    for i in range(0, len(todo), 200):
        url = EUTILS + "efetch.fcgi?db=taxonomy&id=" + ",".join(todo[i:i + 200])
        for attempt in range(5):
            try:
                root = ET.fromstring(urllib.request.urlopen(url, timeout=120).read())
                break
            except Exception:
                time.sleep(5 * (attempt + 1))
        else:
            raise RuntimeError("NCBI taxonomy fetch failed five times")
        for t in root.findall("Taxon"):
            tid = t.findtext("TaxId")
            lin = t.findall("LineageEx/Taxon")
            cache[tid] = {"name": t.findtext("ScientificName"), "rank": t.findtext("Rank"),
                          "lineage": [x.findtext("ScientificName") for x in lin],
                          "ranks": [x.findtext("Rank") for x in lin]}
            # A merged taxid answers under its new id; keep the asked one too.
            for old in t.findall("AkaTaxIds/TaxId"):
                cache[old.text] = cache[tid]
        time.sleep(0.4)
        print(f"  taxonomy {min(i + 200, len(todo)):,}/{len(todo):,}", file=sys.stderr)
    json.dump(cache, open(cache_path, "w"))
    return cache


def pr2_label(t: dict):
    """NCBI taxon -> (genus, PR2-style species label).

    The genus is NCBI's genus, from the ranked lineage - never simply the
    first word of the name: "Nematoda sp." and "Dicrocoeliidae sp." have
    no genus, and taking the first word made Nematoda a genus. A taxon
    with no genus gets an empty genus and a placeholder species named
    after its lowest named rank. "Genus species" -> Genus_species;
    anything less specific -> Genus_sp., a placeholder.
    """
    names, ranks = t["lineage"] + [t["name"]], t["ranks"] + [t["rank"]]
    genus = next((n for n, r in zip(names, ranks) if r == "genus"), "")
    if not genus:
        return "", t["name"].split()[0] + "_sp."
    w = t["name"].split()
    if (t["rank"] in ("species", "subspecies", "varietas", "forma") and len(w) >= 2
            and w[0] == genus and re.fullmatch(r"[a-z][a-z-]+", w[1])
            and w[1] not in ("sp", "cf", "aff")):
        return genus, f"{genus}_{w[1]}"
    return genus, f"{genus}_sp."


def main() -> int:
    pr2_path, out_path, review_path, *fetched = sys.argv[1:]
    workdir = os.path.dirname(out_path) or "."
    pr2 = PR2Index(pr2_path)

    stats = Counter()
    records = []
    for path in fetched:
        for rec in read_fetch(path):
            stats["fetched"] += 1
            acc = rec[0].split(".")[0]
            if acc in pr2.accessions:
                stats["already_in_PR2"] += 1
                continue
            if "UNVERIFIED" in rec[2]:
                stats["unverified"] += 1
                continue
            records.append(rec)
    # The in-silico PCR is the slow part and depends only on the record,
    # so it is cached by accession across runs.
    pcr_cache_path = os.path.join(workdir, "genbank_pcr_cache.json.gz")
    pcr_cache = (json.load(gzip.open(pcr_cache_path, "rt"))
                 if os.path.exists(pcr_cache_path) else {})
    todo = [r for r in records if r[0] not in pcr_cache]
    if todo:
        with Pool(os.cpu_count()) as pool:
            for acc, tx, title, status, ins, fc, rc in pool.map(_pcr_both, todo, chunksize=200):
                pcr_cache[acc] = [status, ins, fc, rc]
        json.dump(pcr_cache, gzip.open(pcr_cache_path, "wt"))
    amplified = [(r[0], r[1], r[2], *pcr_cache[r[0]]) for r in records]

    keep = []
    for acc, tx, title, status, ins, fc, rc in amplified:
        if status != "amplicon":
            stats["no_amplicon"] += 1
            continue
        if ins.count("N") > 0.01 * len(ins):
            stats["too_many_N"] += 1
            continue
        keep.append((acc, tx, title, ins, fc, rc))
    tax = lineages([k[1] for k in keep], os.path.join(workdir, "ncbi_taxids.json"))

    rows, review, how = [], [], Counter()
    for acc, tx, title, ins, fc, rc in keep:
        t = tax.get(tx)
        if not t:
            stats["no_taxonomy"] += 1
            continue
        genus, species = pr2_label(t)
        cls, share = pr2.nearest_class(ins)
        ranks, h = pr2.ranks_for(genus, species, t["lineage"], cls)
        how[h.split()[0]] += 1
        if ranks[4] and cls and cls != ranks[4]:
            stats["review_other_class"] += 1
            review.append([acc, t["name"], ranks[4], cls, f"{share:.2f}", title[:120], ins])
            continue
        rows.append([SOURCE, acc, "18S_rRNA", "nucleus", *ranks, "amplicon",
                     str(fc), str(rc), ins])
        stats["entered"] += 1

    with gzip.open(out_path, "wt") as fh:
        fh.write("\t".join(["source", "accession", "gene", "organelle", *RANKS,
                            "status", "f_edits", "r_edits", "insert"]) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    with open(review_path, "w") as fh:
        fh.write("accession\tncbi_name\texpected_class\tnearest_pr2_class\t"
                 "kmer_share\ttitle\tinsert\n")
        for r in review:
            fh.write("\t".join(r) + "\n")

    for k in ("fetched", "already_in_PR2", "unverified", "no_amplicon", "too_many_N",
              "no_taxonomy", "review_other_class", "entered"):
        print(f"  {k:20s} {stats[k]:>8,}")
    print("  ranks from: " + ", ".join(f"{k} {v:,}" for k, v in how.items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
