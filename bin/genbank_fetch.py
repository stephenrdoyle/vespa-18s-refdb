#!/usr/bin/env python3
"""Fetch 18S records from GenBank for a set of taxa, as FASTA plus taxids.

    bin/genbank_fetch.py <outdir> <taxon> [<taxon> ...]

For each taxon, searches nuccore for 18S / small subunit records of
500-4000bp, excluding uncultured and environmental ones, and fetches them
in batches with their NCBI taxid. Writes <outdir>/<taxon>.fasta.gz with
headers "accession.version taxid=<id> <title>". Resumable: a taxon whose
file exists is skipped.

Uses NCBI E-utilities at their documented rate (3 requests a second
without an API key); set NCBI_API_KEY to go faster.
"""

from __future__ import annotations

import gzip
import os
import re
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
KEY = os.environ.get("NCBI_API_KEY", "")
PAUSE = 0.12 if KEY else 0.4
BATCH = 300
QUERY = ('{taxon}[Organism] AND (18S OR "small subunit ribosomal") '
         'AND 500:4000[SLEN] NOT uncultured NOT environmental')


def get(endpoint: str, **params) -> str:
    if KEY:
        params["api_key"] = KEY
    # POST, so a batch of a few hundred ids fits.
    url = EUTILS + endpoint
    data = urllib.parse.urlencode(params).encode()
    # NCBI turns away bursts (HTTP 429), worse when something else on the
    # same address is fetching too; a long, growing back-off rides it out
    # where five quick retries did not.
    for attempt in range(8):
        try:
            with urllib.request.urlopen(url, data=data, timeout=120) as r:
                body = r.read().decode()
            time.sleep(PAUSE)
            return body
        except Exception as e:
            err = e
            time.sleep(15 * (attempt + 1))
    raise RuntimeError(f"NCBI request failed eight times: {endpoint}: {err!r}")


def fetch_taxon(taxon: str, outdir: str) -> int:
    out = os.path.join(outdir, taxon.replace(" ", "_") + ".fasta.gz")
    if os.path.exists(out):
        return -1
    # Ids first, then batches by id. NCBI's stored-search route (WebEnv,
    # query_key) intermittently answers "Cannot retrieve history data" for
    # a whole taxon, and no amount of retrying fixes it; ids always work.
    body = get("esearch.fcgi", db="nuccore", term=QUERY.format(taxon=taxon),
               retmax=100000)
    count = int(re.search(r"<Count>(\d+)</Count>", body).group(1))
    ids = re.findall(r"<Id>(\d+)</Id>", body)
    if len(ids) != count:
        raise RuntimeError(f"{taxon}: esearch listed {len(ids)} of {count} ids")
    tmp = out + ".partial"
    n = 0
    with gzip.open(tmp, "wt") as fh:
        for start in range(0, count, BATCH):
            batch = ",".join(ids[start:start + BATCH])
            # docsum for taxids, fasta for sequence; the same ids.
            summ = get("esummary.fcgi", db="nuccore", id=batch)
            # Per DocSum: TaxId comes before AccessionVersion inside one,
            # so pairing them across the whole response would shift them.
            taxid = {}
            for doc in summ.split("<DocSum>")[1:]:
                a = re.search(r'Name="AccessionVersion" Type="String">([^<]+)<', doc)
                t = re.search(r'Name="TaxId" Type="Integer">(\d+)<', doc)
                if a and t:
                    taxid[a.group(1)] = t.group(1)
            fa = get("efetch.fcgi", db="nuccore", id=batch, rettype="fasta", retmode="text")
            for rec in fa.split(">")[1:]:
                head, _, seq = rec.partition("\n")
                acc = head.split()[0]
                title = head[len(acc):].strip()
                fh.write(f">{acc} taxid={taxid.get(acc, '')} {title}\n"
                         f"{seq.replace(chr(10), '')}\n")
                n += 1
            print(f"  {taxon}: {min(start + BATCH, count):,}/{count:,}",
                  file=sys.stderr, flush=True)
    os.replace(tmp, out)
    return n


def main() -> int:
    outdir, taxa = sys.argv[1], sys.argv[2:]
    os.makedirs(outdir, exist_ok=True)
    for t in taxa:
        n = fetch_taxon(t, outdir)
        print(f"{t}: " + ("already fetched" if n < 0 else f"{n:,} records"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
