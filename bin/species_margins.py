#!/usr/bin/env python3
"""How firmly does each species-level record name its species?

    bin/species_margins.py <build dir> [W1xW2 ...] [--workers N]

    bin/species_margins.py work/build 132x128 250x200

A resolution table says whether a record resolves to species, not by how
much. *Trichinella nelsoni* resolves through 2x301 windows, but by one
base from the record *T. britovi*, *murrelli* and *nativa* share: a read
from a lineage the release lacks, a base away, is named *T. nelsoni* just
as confidently (v0.4.0 feedback §5e). This writes, for every record
resolved to species, the edits to the nearest record of another species
- read whole, and through each Illumina geometry's two windows.

"Another species": any record whose species label differs, including a
higher-rank record that several species share. Candidates come from
minimap2 (whole insert, all against all; the 30 most similar records of
another species per record), then exact edit distances; read whole and
through the windows the nearest is taken among them, so a margin is an
upper bound if the true nearest is not among the 30 - with V4's shared
conserved stretches, it is in practice.

Writes <build dir>/species_margins.tsv: id, margin_whole, then
margin_<W1>x<W2> per geometry, nearest_whole (the id it was measured to).
bin/release.py adds the matching margin to each resolution table.
"""

from __future__ import annotations

import collections
import csv
import os
import subprocess
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vespa import edit_distance  # noqa: E402

TOP = 30
_S = _T = _G = None


def _init(seqs, tax, geoms):
    global _S, _T, _G
    _S, _T, _G = seqs, tax, geoms


def _one(job):
    rid, cands = job
    s = _S[rid]
    whole = sorted((edit_distance(s, _S[c]), c) for c in cands)
    out = [rid, str(whole[0][0]) if whole else ""]
    for w1, w2 in _G:
        if len(s) <= w1 + w2:            # the windows cover the whole insert
            out.append(out[1])
            continue
        best = None
        for _, c in whole:
            t = _S[c]
            d = edit_distance(s[:w1], t[:w1]) + edit_distance(s[-w2:], t[-w2:])
            best = d if best is None else min(best, d)
        out.append("" if best is None else str(best))
    out.append(whole[0][1] if whole else "")
    return out


def main() -> int:
    args = sys.argv[1:]
    workers = os.cpu_count()
    if "--workers" in args:
        k = args.index("--workers")
        workers = int(args[k + 1])
        del args[k:k + 2]
    build, geoms = args[0], [tuple(int(x) for x in g.split("x")) for g in args[1:]]
    with open(os.path.join(build, "taxonomy.tsv")) as fh:
        tax = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    seqs, cur = {}, None
    with open(os.path.join(build, "refdb.fasta")) as fh:
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].strip()
                seqs[cur] = ""
            else:
                seqs[cur] += line.strip()

    fa = os.path.join(build, "refdb.fasta")
    p = subprocess.Popen(["minimap2", "-t", str(workers), "-X", "-k", "13", "-w", "5",
                          "-N", "2000", "-p", "0.5", fa, fa],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    score = collections.defaultdict(dict)
    for line in p.stdout:
        f = line.split("\t", 11)
        q, t, m = f[0], f[5], int(f[9])
        # Only near neighbours matter: chain matches over half the longer
        # insert (as in neighbour_check.py), and at most 4 x TOP kept per
        # record, so dense clusters do not fill memory.
        if m < 0.5 * max(int(f[1]), int(f[6])):
            continue
        for a, b in ((q, t), (t, q)):
            if tax[a]["deepest_rank"] == "species" and tax[b]["species"] != tax[a]["species"]:
                d = score[a]
                if m > d.get(b, -1):
                    d[b] = m
                    if len(d) > 8 * TOP:
                        for c in sorted(d, key=d.get)[:4 * TOP]:
                            del d[c]
    if p.wait():
        sys.exit("minimap2 failed")
    jobs = []
    for rid, t in tax.items():
        if t["deepest_rank"] != "species":
            continue
        top = sorted(score[rid], key=lambda c: -score[rid][c])[:TOP]
        jobs.append((rid, top))
    with Pool(workers, initializer=_init, initargs=(seqs, tax, geoms)) as pool:
        rows = pool.map(_one, jobs, chunksize=200)

    path = os.path.join(build, "species_margins.tsv")
    with open(path, "w") as fh:
        fh.write("\t".join(["id", "margin_whole"]
                           + [f"margin_{w1}x{w2}" for w1, w2 in geoms] + ["nearest_whole"]) + "\n")
        for r in sorted(rows):
            fh.write("\t".join(r) + "\n")
    c = collections.Counter(r[1] for r in rows)
    thin = sum(v for k, v in c.items() if k and int(k) <= 2)
    print(f"{len(rows):,} species-level records; {thin:,} within 2 edits (read whole) "
          f"of another species; no candidate: {c['']:,}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
