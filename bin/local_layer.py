#!/usr/bin/env python3
"""Layer: our own V4 sequences, from local/sequences.tsv. The most trusted.

    bin/local_layer.py local/sequences.tsv work/pr2_layer.tsv.gz \
        work/local_layer.tsv.gz work/local_review.tsv

local/sequences.tsv is curated by hand and committed: one row per
sequence, with
  id         a stable name for it
  ncbi_name  the taxon, as NCBI writes it ("Teladorsagia circumcincta")
  sequence   either a read or region that spans both primers - it is
             trimmed to the insert here - or the insert itself, already
             trimmed (as the genome scans write it)
  evidence   where it came from, specifically enough to check again
  added_by, date
  f_edits, r_edits  optional: for a bare insert, the primer edits the scan
             that found it recorded (0 and 0 = both primers exact)

Ranks come from the NCBI lineage of ncbi_name, mapped onto PR2's as for
every other layer. A sequence whose nearest PR2 insert is in another PR2
class goes to review rather than the layer: being ours does not make a
contaminated read right.

Pass this layer to collapse.py last.
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vespa import LAYER_COLUMNS, PR2Index, pcr_full  # noqa: E402

SOURCE = "local"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def lineage(name: str, cache: dict) -> list:
    if name not in cache:
        q = urllib.parse.quote(name)
        ids = re.findall(r"<Id>(\d+)</Id>", urllib.request.urlopen(
            f"{EUTILS}esearch.fcgi?db=taxonomy&term={q}", timeout=30).read().decode())
        time.sleep(0.4)
        names = []
        if ids:
            t = ET.fromstring(urllib.request.urlopen(
                f"{EUTILS}efetch.fcgi?db=taxonomy&id={ids[0]}", timeout=30).read()).find("Taxon")
            names = [x.text for x in t.findall("LineageEx/Taxon/ScientificName")]
            time.sleep(0.4)
        cache[name] = names
    return cache[name]


def main() -> int:
    table, pr2_path, out_path, review_path = sys.argv[1:5]
    cache_path = os.path.join(os.path.dirname(out_path) or ".", "ncbi_lineages.json")
    cache = json.load(open(cache_path)) if os.path.exists(cache_path) else {}
    pr2 = PR2Index(pr2_path)

    rows, review = [], []
    with open(table) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            seq = r["sequence"].strip().upper()
            status, ins, fc, rc, fsite, rsite = pcr_full(seq)
            if status != "amplicon":
                # Already the insert, as genome scans write it: the primer
                # sites are not known, but the scan's edit counts may be
                # given (f_edits, r_edits) - both 0 is an exact fit.
                ins, fsite, rsite = seq, "", ""
                fc, rc = (r.get("f_edits") or "").strip(), (r.get("r_edits") or "").strip()
            w = r["ncbi_name"].split()
            genus = w[0]
            species = f"{w[0]}_{w[1]}" if len(w) > 1 and w[1] not in ("sp.", "sp") else f"{genus}_sp."
            cls, share = pr2.nearest_class(ins)
            ranks, how = pr2.ranks_for(genus, species, lineage(r["ncbi_name"], cache), cls)
            if ranks[4] and cls and cls != ranks[4]:
                review.append([r["id"], r["ncbi_name"], ranks[4], cls, f"{share:.2f}"])
                continue
            rows.append([SOURCE, r["id"], "18S_rRNA", "nucleus", *ranks, "amplicon",
                         str(fc), str(rc), fsite or "", rsite or "", ins])
            print(f"  {r['id']}: {r['ncbi_name']} -> {';'.join(ranks[4:])} "
                  f"({how}; {len(ins)}bp{', trimmed from a longer sequence' if status == 'amplicon' else ''})")
    json.dump(cache, open(cache_path, "w"), indent=1)

    with gzip.open(out_path, "wt") as fh:
        fh.write("\t".join(LAYER_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    with open(review_path, "w") as fh:
        fh.write("id\tncbi_name\texpected_class\tnearest_pr2_class\tkmer_share\n")
        for r in review:
            fh.write("\t".join(r) + "\n")
    print(f"local layer: {len(rows)} entered, {len(review)} to review")
    return 0


if __name__ == "__main__":
    sys.exit(main())
