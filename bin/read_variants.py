#!/usr/bin/env python3
"""Turn long-read V4 scans into genome-layer input.

    bin/read_variants.py curation/longread_runs.tsv <scans dir> <out dir>

curation/longread_runs.tsv lists, per run: run accession, NCBI organism
name, and the scan's file stem in <scans dir> (from reads_insilico_pcr.py).
For each run, a variant is kept if it is at least MIN_SHARE of the run's
V4 amplicons and seen in at least MIN_READS reads: HiFi errors recur at
homopolymers, so a handful of reads is not evidence, but a real minor
rDNA copy holds a steady share. Kept variants are written as genome hit
files (one row per supporting read), named for the run, with a taxa.json
in pararef_layer.py's format - so the same genome rules (two copies, the
nearest-PR2-class check that sends contaminants to review) build the
layer:

    bin/pararef_layer.py <out dir> <out dir>/taxa.json work/pr2_layer.tsv.gz \
        work/reads_layer.tsv.gz work/reads_review.tsv ENA_reads
"""

from __future__ import annotations

import collections
import csv
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

MIN_SHARE, MIN_READS = 0.02, 5
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
COLS = ["genome", "contig", "strand", "start", "end", "f_edits", "r_edits", "insert_len", "insert"]


def lineage(name):
    q = urllib.parse.quote(name)
    ids = re.findall(r"<Id>(\d+)</Id>", urllib.request.urlopen(
        f"{EUTILS}esearch.fcgi?db=taxonomy&term={q}", timeout=60).read().decode())
    time.sleep(0.4)
    t = ET.fromstring(urllib.request.urlopen(
        f"{EUTILS}efetch.fcgi?db=taxonomy&id={ids[0]}", timeout=60).read()).find("Taxon")
    time.sleep(0.4)
    return ids[0], t.findtext("ScientificName"), [x.text for x in t.findall("LineageEx/Taxon/ScientificName")]


def main() -> int:
    runs_path, scans, out = sys.argv[1:4]
    os.makedirs(out, exist_ok=True)
    taxa = {"assemblies": {}, "taxa": {}}
    for r in csv.DictReader((l for l in open(runs_path) if not l.startswith("#")), delimiter="\t"):
        rows = list(csv.DictReader(open(os.path.join(scans, r["stem"] + ".tsv")), delimiter="\t"))
        n = collections.Counter(x["insert"] for x in rows)
        total = sum(n.values())
        keep = {s for s, k in n.items() if k >= MIN_READS and k >= MIN_SHARE * total}
        with open(os.path.join(out, r["run"] + ".tsv"), "w") as fh:
            fh.write("\t".join(COLS) + "\n")
            for x in rows:
                if x["insert"] in keep:
                    fh.write("\t".join([r["run"], x["read"], x["strand"], "0", "0", x["f_edits"],
                                        x["r_edits"], x["insert_len"], x["insert"]]) + "\n")
        open(os.path.join(out, r["run"] + ".done"), "w").write(f"{total}\n")
        taxid, name, lin = lineage(r["organism"])
        taxa["assemblies"][r["run"]] = taxid
        taxa["taxa"][taxid] = {"name": name, "lineage": lin}
        print(f"  {r['run']} {name}: {total} amplicons, {len(keep)} variants kept "
              f"({sum(n[s] for s in keep)} reads)", file=sys.stderr)
    json.dump(taxa, open(os.path.join(out, "taxa.json"), "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
