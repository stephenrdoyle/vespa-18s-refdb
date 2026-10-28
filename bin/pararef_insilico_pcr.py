#!/usr/bin/env python3
"""Stream the ParaRef genome archive and find VESPA V4 amplicons in it.

    bin/pararef_insilico_pcr.py <outdir> [workers]

ParaRef (Niemann et al. 2025, Genome Biology; zenodo.org/records/13744644,
CC-BY 4.0) is 645 decontaminated parasite assemblies in one 21.8GB
tar.gz. It is read as a stream - nothing is written to disk but the
hits - and each genome is handed to a worker pool as it arrives. The
scan is wbps_insilico_pcr.py's, so the two sources are judged alike.

Output per genome is the same as wbps_insilico_pcr.py's (<stem>.tsv and
<stem>.done), stem being the assembly id. A tar stream cannot seek, so a
rerun downloads everything again; it does skip scanning genomes already
done.
"""

from __future__ import annotations

import gzip
import io
import os
import sys
import tarfile
import time
import urllib.request
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wbps_insilico_pcr import amplicons_in, write_hits  # noqa: E402

URL = ("https://zenodo.org/api/records/13744644/files/"
       "Decontaminated_Parasite_Genomes.tar.gz/content")


def stem_of(member: str) -> str:
    """The assembly id: FungiDB-61-Acandida2VRR_Genome.decontaminated.fna.gz
    -> FungiDB-61-Acandida2VRR, matching the names in ParaRef's map files."""
    base, changed = os.path.basename(member), True
    while changed:
        changed = False
        for suffix in (".gz", ".fna", ".fasta", ".fa", ".decontaminated", "_Genome"):
            if base.endswith(suffix):
                base, changed = base[: -len(suffix)], True
    return base


def one(stem: str, data: bytes, gz: bool, outdir: str):
    try:
        return _one(stem, data, gz, outdir)
    except Exception as e:  # one bad genome must not end a 22GB stream
        with open(os.path.join(outdir, stem + ".failed"), "w") as fh:
            fh.write(repr(e) + "\n")
        return stem, 0, -1


def _one(stem: str, data: bytes, gz: bool, outdir: str):
    t0 = time.time()
    raw = io.BytesIO(data)
    text = io.TextIOWrapper(gzip.GzipFile(fileobj=raw) if gz else raw,
                            "ascii", errors="replace")
    bases, rows = amplicons_in(text)
    write_hits(outdir, stem, rows, bases, time.time() - t0)
    return stem, bases, len(rows)


def main() -> int:
    outdir = sys.argv[1]
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()
    n = 0
    pending = set()
    # Each queued genome is held in memory compressed; bounding the queue
    # bounds that, and the download simply waits for a free worker.
    limit = workers + 2
    with urllib.request.urlopen(URL, timeout=300) as resp, \
         tarfile.open(fileobj=resp, mode="r|gz") as tar, \
         ProcessPoolExecutor(workers) as pool:
        for m in tar:
            if not m.isfile() or not any(x in m.name for x in (".fna", ".fa", ".fasta")):
                continue
            stem = stem_of(m.name)
            if os.path.exists(os.path.join(outdir, stem + ".done")):
                continue
            data = tar.extractfile(m).read()
            pending.add(pool.submit(one, stem, data, m.name.endswith(".gz"), outdir))
            while len(pending) >= limit:
                done, pending = wait(pending, return_when=FIRST_COMPLETED)
                for f in done:
                    s, b, k = f.result()
                    n += 1
                    print(f"[{n} {time.time() - t0:6.0f}s] {s}  {b / 1e6:8.1f} Mb  "
                          f"{k} amplicons", flush=True)
        for f in pending:
            s, b, k = f.result()
            n += 1
            print(f"[{n} {time.time() - t0:6.0f}s] {s}  {b / 1e6:8.1f} Mb  {k} amplicons",
                  flush=True)
    print(f"done: {n} genomes in {time.time() - t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
