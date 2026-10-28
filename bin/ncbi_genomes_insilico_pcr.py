#!/usr/bin/env python3
"""Stream every NCBI genome assembly for a taxon through the V4 scan.

    bin/ncbi_genomes_insilico_pcr.py <taxon> <outdir> [workers]

Lists the taxon's assemblies with the NCBI Datasets API, streams each
*_genomic.fna.gz from the NCBI FTP site through wbps_insilico_pcr.py's
scan (nothing is stored but the hits), and writes:
  <outdir>/<accession>.tsv / .done   per-assembly hits, as for WBPS
  <outdir>/taxa.json                 assembly -> taxid, taxid -> name and
                                     lineage: the format pararef_layer.py
                                     reads, so the same layer code builds it
  <outdir>/assemblies.tsv            what was listed: level, size, N50, date

Raw NCBI assemblies, not decontaminated ones: ParaRef showed that
decontamination masks rDNA.
"""

from __future__ import annotations

import gzip
import io
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wbps_insilico_pcr import amplicons_in, write_hits  # noqa: E402

DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{}/dataset_report?page_size=1000"
FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def assemblies(taxon: str) -> list:
    out, token = [], ""
    while True:
        url = DATASETS.format(urllib.parse.quote(taxon)) + (f"&page_token={token}" if token else "")
        d = json.load(urllib.request.urlopen(url, timeout=120))
        for r in d.get("reports", []):
            a, st = r["assembly_info"], r.get("assembly_stats", {})
            out.append({"accession": r["accession"],
                        "name": a.get("assembly_name", ""),
                        "organism": r["organism"]["organism_name"],
                        "taxid": str(r["organism"]["tax_id"]),
                        "level": a.get("assembly_level", ""),
                        "length": st.get("total_sequence_length", ""),
                        "contig_n50": st.get("contig_n50", ""),
                        "released": a.get("release_date", "")})
        token = d.get("next_page_token")
        if not token:
            return out


def ftp_url(acc: str, name: str) -> str:
    stem = f"{acc}_{re.sub(r'[^A-Za-z0-9._]', '_', name)}"
    pre, num = acc[:3], acc[4:13]
    return f"{FTP}/{pre}/{num[0:3]}/{num[3:6]}/{num[6:9]}/{stem}/{stem}_genomic.fna.gz"


def one(acc: str, url: str, outdir: str):
    t0 = time.time()
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=300) as resp:
                text = io.TextIOWrapper(gzip.GzipFile(fileobj=resp), "ascii", errors="replace")
                bases, rows = amplicons_in(text)
            write_hits(outdir, acc, rows, bases, time.time() - t0)
            return acc, bases, len(rows)
        except Exception as e:
            err = e
            time.sleep(10 * (attempt + 1))
    with open(os.path.join(outdir, acc + ".failed"), "w") as fh:
        fh.write(f"{url}\n{err!r}\n")
    return acc, 0, -1


def lineages(taxids) -> dict:
    out = {}
    ids = sorted(set(taxids))
    for i in range(0, len(ids), 200):
        root = ET.fromstring(urllib.request.urlopen(
            EUTILS + "efetch.fcgi?db=taxonomy&id=" + ",".join(ids[i:i + 200]), timeout=120).read())
        for t in root.findall("Taxon"):
            rec = {"name": t.findtext("ScientificName"),
                   "lineage": [x.findtext("ScientificName") for x in t.findall("LineageEx/Taxon")]}
            out[t.findtext("TaxId")] = rec
            for old in t.findall("AkaTaxIds/TaxId"):
                out[old.text] = rec
        time.sleep(0.4)
    return out


def main() -> int:
    taxon, outdir = sys.argv[1], sys.argv[2]
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    os.makedirs(outdir, exist_ok=True)
    asms = assemblies(taxon)
    with open(os.path.join(outdir, "assemblies.tsv"), "w") as fh:
        cols = ["accession", "name", "organism", "taxid", "level", "length", "contig_n50", "released"]
        fh.write("\t".join(cols) + "\n")
        for a in asms:
            fh.write("\t".join(str(a[c]) for c in cols) + "\n")
    json.dump({"assemblies": {a["accession"]: a["taxid"] for a in asms},
               "taxa": lineages(a["taxid"] for a in asms)},
              open(os.path.join(outdir, "taxa.json"), "w"), indent=1)
    todo = [a for a in asms if not os.path.exists(os.path.join(outdir, a["accession"] + ".done"))]
    print(f"{taxon}: {len(asms)} assemblies, {len(todo)} to scan", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(workers) as pool:
        futs = [pool.submit(one, a["accession"], ftp_url(a["accession"], a["name"]), outdir)
                for a in todo]
        for k, f in enumerate(as_completed(futs), 1):
            acc, b, n = f.result()
            print(f"[{k}/{len(todo)} {time.time() - t0:5.0f}s] {acc}  {b / 1e6:6.1f} Mb  "
                  + ("FAILED" if n < 0 else f"{n} amplicons"), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
