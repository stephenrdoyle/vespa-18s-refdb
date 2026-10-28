#!/usr/bin/env python3
"""Layer: V4 from a set of genomes, under the genome-layer rules.

    bin/pararef_layer.py <hits dir> <taxa.json> work/pr2_layer.tsv.gz \
        <layer.tsv.gz> <review.tsv> [source]

Built for ParaRef (source "ParaRef", the default); also builds the layer
for NCBI genome sets from ncbi_genomes_insilico_pcr.py, whose taxa.json
has the same shape - e.g. source "NCBI_genomes".

An assembly listed under both its GenBank (GCA_) and RefSeq (GCF_)
accession is the same genome twice; only the GCA_ copy is counted, or a
single-copy variant would pass the two-copy rule.

<pararef_out> is pararef_insilico_pcr.py's per-genome hits. Each genome
is tied to its taxon through ParaRef's own map files (assembly -> NCBI
taxid), resolved to names and lineages once and kept in
pararef_taxa.json. Species are pooled across assemblies, and the rules
are the genome layer's (vespa.genome_layer_rows): two copies, own class.

ParaRef's genomes are decontaminated, and decontamination can mask rDNA:
18S is conserved enough across kingdoms to look like contamination. So
an assembly here may hold less V4 than the same assembly in WBPS. That
is measured, not fixed, here: see docs/DESIGN.md.
"""

from __future__ import annotations

import glob
import json
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vespa import (GENOME_REVIEW_COLUMNS, LAYER_COLUMNS, PR2Index,  # noqa: E402
                   genome_layer_rows, write_table)

SOURCE = "ParaRef"


def assembly_of(stem: str, known: set) -> str | None:
    """GCA_900618025.1_B_timori_..._genomic -> GCA_900618025.1."""
    m = re.match(r"(GC[AF]_\d+\.\d+)", stem)
    if m and m.group(1) in known:
        return m.group(1)
    return next((k for k in known if stem.startswith(k)), None)


def main() -> int:
    hits_dir, taxa_path, pr2_path, out_path, review_path = sys.argv[1:6]
    source = sys.argv[6] if len(sys.argv) > 6 else SOURCE
    taxa = json.load(open(taxa_path))
    asm, info = taxa["assemblies"], taxa["taxa"]
    pr2 = PR2Index(pr2_path)

    copies = defaultdict(Counter)
    where, lineage_of = {}, {}
    unmatched, untaxed = [], []
    tsvs = sorted(glob.glob(os.path.join(hits_dir, "*.tsv")))
    stems = {os.path.basename(t)[:-4] for t in tsvs}
    twins = 0
    for tsv in tsvs:
        stem = os.path.basename(tsv)[:-4]
        if stem.startswith("GCF_") and "GCA_" + stem[4:] in stems:
            twins += 1
            continue
        if stem in ("assemblies",):
            continue
        a = assembly_of(stem, set(asm))
        if a is None:
            unmatched.append(stem)
            continue
        t = info.get(asm[a])
        if not t:
            untaxed.append(a)
            continue
        w = t["name"].split()
        species = (f"{w[0]}_{w[1]}" if len(w) > 1 and re.fullmatch(r"[a-z][a-z-]+", w[1])
                   and w[1] not in ("sp", "cf", "aff") else f"{w[0]}_sp.")
        lineage_of[species] = t["lineage"]
        with open(tsv) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                ins = f[-1]
                copies[species][ins] += 1
                where.setdefault((species, ins), f"{a}:{f[1]}:{f[3]}-{f[4]}")

    def ranks_for(species: str, seq_class: str = ""):
        return pr2.ranks_for(species.split("_")[0], species, lineage_of[species], seq_class)

    rows, review, stats = genome_layer_rows(copies, where, ranks_for, pr2, source)
    write_table(out_path, LAYER_COLUMNS, rows)
    write_table(review_path, GENOME_REVIEW_COLUMNS, review)

    print(f"genomes with hits files: {len(tsvs)} ({twins} RefSeq twins of GenBank assemblies skipped); "
          f"not matched to a map: {len(unmatched)}; no taxonomy: {len(untaxed)}")
    if unmatched:
        print("  unmatched: " + ", ".join(unmatched[:10]))
    print(f"species with V4: {len(copies)}  (ranks from PR2 genus {stats['PR2']}, "
          f"NCBI {stats['NCBI']}, unplaced {stats['unplaced']})")
    print(f"variants seen {stats['variants_seen']:,}; dropped as single-copy "
          f"{stats['dropped_unsupported']:,}; sent to review as another class "
          f"{stats['review_other_class']}; entered {stats['entered']:,} "
          f"({stats['single_copy']} from single-copy species)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
