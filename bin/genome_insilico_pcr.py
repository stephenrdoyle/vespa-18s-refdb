#!/usr/bin/env python3
"""Find VESPA V4 amplicons in one genome, from a URL or a local file.

    bin/genome_insilico_pcr.py <url or path> <stem> <outdir>

For assemblies outside WormBase ParaSite and ParaRef - our own, for a
start. Plain or gzipped FASTA; a URL is streamed, not stored. The scan
and the output (<stem>.tsv, <stem>.done) are wbps_insilico_pcr.py's, so
the hits feed the same layer-building as any other genome.
"""

from __future__ import annotations

import gzip
import io
import os
import sys
import time
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wbps_insilico_pcr import amplicons_in, write_hits  # noqa: E402


def main() -> int:
    src, stem, outdir = sys.argv[1:4]
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()
    raw = (urllib.request.urlopen(src, timeout=300)
           if src.startswith(("http://", "https://")) else open(src, "rb"))
    with raw:
        stream = gzip.GzipFile(fileobj=raw) if src.endswith(".gz") else raw
        bases, rows = amplicons_in(io.TextIOWrapper(stream, "ascii", errors="replace"))
    write_hits(outdir, stem, rows, bases, time.time() - t0)
    print(f"{stem}: {bases / 1e6:.1f} Mb, {len(rows)} amplicons, "
          f"{time.time() - t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
