#!/usr/bin/env python3
"""Find V4 amplicons in raw long reads, streamed from ENA.

    bin/reads_insilico_pcr.py <fastq.gz URL or path> <out.tsv> [max reads] [workers]

For species whose genome assembly lacks rDNA, or that have reads but no
assembly yet. rDNA is present in tens to hundreds of copies, so the
first few hundred thousand long reads (PacBio HiFi, nanopore) already
cross V4 hundreds of times; nothing but the hits is kept. Each read is
scanned as a contig, with the same scan as the genome layers
(wbps_insilico_pcr.amplicons_in: both strands, forward primer <= 1 edit
by seed, reverse <= 2).

Writes one row per amplicon: read, strand, primer edits, insert length,
insert. A variant seen in several independent reads is a real rDNA copy;
HiFi errors are random, so a single-read variant is not evidence.
"""

from __future__ import annotations

import gzip
import io
import os
import sys
import urllib.request
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wbps_insilico_pcr import amplicons_in  # noqa: E402

CHUNK = 2000


def open_source(src):
    """A local file, or a URL. EBI's HTTPS sometimes refuses (403) under
    load: retry, then fall back to its FTP server."""
    if not src.startswith(("ftp.", "http://", "https://", "ftp://")):
        return open(src, "rb")
    import time
    path = src.split("://", 1)[-1]
    for attempt, url in enumerate(["https://" + path, "https://" + path, "ftp://" + path,
                                   "ftp://" + path]):
        try:
            return urllib.request.urlopen(url, timeout=300)
        except Exception as e:  # noqa: BLE001
            err = e
            time.sleep(10 * (attempt + 1))
    raise err


def fastq(stream, limit):
    n = 0
    while n < limit:
        h = stream.readline()
        if not h:
            return
        s = stream.readline().strip()
        stream.readline()
        stream.readline()
        n += 1
        yield h[1:].split()[0], s


def scan(chunk):
    text = "".join(f">{rid}\n{seq}\n" for rid, seq in chunk)
    _, rows = amplicons_in(io.StringIO(text))
    return rows


def main() -> int:
    src, out = sys.argv[1:3]
    limit = int(sys.argv[3]) if len(sys.argv) > 3 else 300_000
    workers = int(sys.argv[4]) if len(sys.argv) > 4 else os.cpu_count()
    raw = open_source(src)
    stream = io.TextIOWrapper(gzip.GzipFile(fileobj=raw), "ascii", errors="replace")

    def chunks():
        buf = []
        for rec in fastq(stream, limit):
            buf.append(rec)
            if len(buf) == CHUNK:
                yield buf
                buf = []
        if buf:
            yield buf

    n_hits = 0
    with Pool(workers) as pool, open(out, "w") as fh:
        fh.write("read\tstrand\tf_edits\tr_edits\tinsert_len\tinsert\n")
        for rows in pool.imap(scan, chunks()):
            for ctg, strand, a, b, fc, rc, ins in rows:
                fh.write(f"{ctg}\t{strand}\t{fc}\t{rc}\t{len(ins)}\t{ins}\n")
                n_hits += 1
    raw.close()
    print(f"{src}: first {limit:,} reads, {n_hits} V4 amplicons", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
