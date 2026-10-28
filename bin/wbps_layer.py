#!/usr/bin/env python3
"""Layer 2: V4 from WormBase ParaSite genomes, filtered and given PR2 ranks.

    bin/wbps_layer.py results/wbps19_v4_amplicons.tsv.gz work/pr2_layer.tsv.gz \
        work/wbps_layer.tsv.gz work/wbps_review.tsv

Input is wbps_insilico_pcr.py's hits: every V4 copy found in every
assembly. Three rules (docs/DESIGN.md, genome layer) decide what enters:

1. SUPPORT. A genome's rDNA array holds many copies; a variant seen once
   is as likely an assembly error as a real variant, and old drafts are
   full of them. Per species, pooled across its assemblies, a variant
   enters if it is seen in at least two copies. A species with no such
   variant contributes its most common one, marked single_copy.
2. PHYLUM. Assemblies carry other organisms' rDNA - host, vector, gut
   contents, lab contaminants. A variant enters only if its nearest PR2
   sequence is in the same PR2 class (Nematoda, Platyhelminthes, ...) as
   the species itself. The rest go to the review file, not the layer.
3. LABEL. Entering variants carry the genome's species name. collapse.py
   lists this layer after PR2, so where they disagree the genome wins.

Ranks come from PR2 where PR2 knows the genus. Otherwise the NCBI
lineage is walked upward to the first name PR2 has at any rank; PR2's
ranks are used down to there and the rest filled with PR2-style
placeholders. NCBI lineages are cached in work/ncbi_lineages.json.
"""

from __future__ import annotations

import gzip
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vespa import RANKS, PR2Index, informative, open_text  # noqa: E402

SOURCE = "WBPS19"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def ncbi_lineage(name: str, cache: dict) -> list:
    """Scientific names from root to the species, or [] if NCBI has none."""
    if name in cache:
        return cache[name]
    q = urllib.parse.quote(name.replace("_", " "))
    lineage = []
    try:
        with urllib.request.urlopen(f"{EUTILS}esearch.fcgi?db=taxonomy&term={q}", timeout=30) as r:
            ids = re.findall(r"<Id>(\d+)</Id>", r.read().decode())
        time.sleep(0.4)
        if ids:
            with urllib.request.urlopen(f"{EUTILS}efetch.fcgi?db=taxonomy&id={ids[0]}", timeout=30) as r:
                t = ET.fromstring(r.read()).find("Taxon")
            lineage = [x.text for x in t.findall("LineageEx/Taxon/ScientificName")]
            lineage.append(t.findtext("ScientificName"))
            time.sleep(0.4)
    except Exception as e:  # leave uncached so a rerun tries again
        print(f"  NCBI lookup failed for {name}: {e}", file=sys.stderr)
        return []
    cache[name] = lineage
    return lineage


def main() -> int:
    hits_path, pr2_path, out_path, review_path = sys.argv[1:5]
    workdir = os.path.dirname(out_path) or "."
    cache_path = os.path.join(workdir, "ncbi_lineages.json")
    cache = json.load(open(cache_path)) if os.path.exists(cache_path) else {}

    pr2 = PR2Index(pr2_path)
    nearest_class = pr2.nearest_class

    def ranks_for(species: str) -> tuple:
        genus = species.split("_")[0].capitalize()
        sp = species[0].upper() + species[1:]
        # Strain-coded names (Auanema_ju1783) are unknown to NCBI by that
        # name; the genus is not.
        return pr2.ranks_for(genus, sp, lambda: ncbi_lineage(sp, cache)[:-1]
                             or ncbi_lineage(genus, cache))

    # --- genome hits, pooled per species
    copies = defaultdict(Counter)
    where = {}
    with open_text(hits_path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(head)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            species = f[col["genome"]].split(".")[0]
            ins = f[col["insert"]]
            copies[species][ins] += 1
            where.setdefault((species, ins), f"{f[col['genome']]}:{f[col['contig']]}:"
                                            f"{f[col['start']]}-{f[col['end']]}")

    out_rows, review = [], []
    stats = Counter()
    for species in sorted(copies):
        ranks, how = ranks_for(species)
        stats[how.split()[0]] += 1
        variants = copies[species]
        supported = [v for v, n in variants.items() if n >= 2]
        single = not supported
        if single:
            supported = [variants.most_common(1)[0][0]]
        stats["variants_seen"] += len(variants)
        for v in variants:
            if v not in supported:
                stats["dropped_unsupported"] += 1
                continue
            cls, score = nearest_class(v)
            row = [SOURCE, where[(species, v)], "18S_rRNA", "nucleus", *ranks,
                   "amplicon", "", "", v]
            if ranks[4] and informative(ranks[4]) and cls and cls != ranks[4]:
                stats["review_other_class"] += 1
                review.append([species, where[(species, v)], str(variants[v]),
                               ranks[4], cls, f"{score:.2f}", how, v])
                continue
            stats["entered"] += 1
            stats["single_copy"] += single
            out_rows.append(row)
    json.dump(cache, open(cache_path, "w"), indent=1)

    with gzip.open(out_path, "wt") as fh:
        fh.write("\t".join(["source", "accession", "gene", "organelle", *RANKS,
                            "status", "f_edits", "r_edits", "insert"]) + "\n")
        for r in out_rows:
            fh.write("\t".join(r) + "\n")
    with open(review_path, "w") as fh:
        fh.write("species\tcopy\tcopies\texpected_class\tnearest_pr2_class\t"
                 "kmer_share\tranks_from\tinsert\n")
        for r in review:
            fh.write("\t".join(r) + "\n")

    print(f"species with V4: {len(copies)}  (ranks from PR2 genus "
          f"{stats['PR2']}, NCBI {stats['NCBI']}, unplaced {stats['unplaced']})")
    print(f"variants seen {stats['variants_seen']:,}; dropped as single-copy "
          f"{stats['dropped_unsupported']:,}; sent to review as another class "
          f"{stats['review_other_class']}; entered {stats['entered']:,} "
          f"({stats['single_copy']} from single-copy species)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
