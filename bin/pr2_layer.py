#!/usr/bin/env python3
"""Layer 1: PR2 trimmed to the VESPA V4 insert, one row per PR2 entry.

    bin/pr2_layer.py sources/pr2_version_5.1.1_SSU_taxo_long.fasta.gz \
        work/pr2_layer.tsv.gz [workers]

Every entry is kept in the output with its PCR status, so what was set
aside and why stays countable; only status "amplicon" carries an insert.
About ten minutes on 8 workers.
"""

from __future__ import annotations

import gzip
import os
import sys
from collections import Counter
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vespa import LAYER_COLUMNS, pcr_full, read_taxo_long  # noqa: E402


def _one(rec):
    acc, gene, organelle, ranks, seq = rec
    return (acc, gene, organelle, ranks, *pcr_full(seq))


def main() -> int:
    src, out = sys.argv[1], sys.argv[2]
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else os.cpu_count()
    version = os.path.basename(src).split("_SSU")[0].replace("pr2_version_", "PR2_")
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    counts: Counter = Counter()
    tmp = out + ".partial"
    with Pool(workers) as pool, gzip.open(tmp, "wt") as fh:
        fh.write("\t".join(LAYER_COLUMNS) + "\n")
        for acc, gene, org, ranks, status, ins, fc, rc, fsite, rsite in pool.imap(
                _one, read_taxo_long(src), chunksize=500):
            counts[status] += 1
            ranks = (list(ranks) + [""] * 9)[:9]
            fh.write("\t".join([version, acc, gene, org, *ranks, status,
                                "" if fc is None else str(fc),
                                "" if rc is None else str(rc),
                                fsite or "", rsite or "", ins or ""]) + "\n")
    os.replace(tmp, out)
    total = sum(counts.values())
    print(f"{version}: {total:,} entries")
    for k in ("amplicon", "F_only", "R_only", "neither"):
        print(f"  {k:9s} {counts[k]:>8,}  {100 * counts[k] / total:5.1f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
