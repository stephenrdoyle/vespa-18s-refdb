#!/usr/bin/env python3
"""Records whose near-identical neighbours disagree at a high rank.

    minimap2 -t 8 -X -k 21 -w 11 -N 500 -p 0.8 refdb.fasta refdb.fasta \
        | awk '{m=($2>$7?$2:$7); if ($10/m >= 0.85) print $1"\t"$6}' \
        | gzip > candidates.tsv.gz
    bin/neighbour_check.py RELEASE_DIR candidates.tsv.gz OUT.tsv [MIN_IDENTITY]

minimap2 only proposes pairs (long shared exact stretches, chain covering
>= 85% of the longer insert - generous for a 99% cut); every pair is then
measured exactly. A full base-level all-against-all of 97,708 inserts
runs for hours.

Two V4 inserts >= 99% identical (default) that sit in different domains,
supergroups, divisions or classes are almost always a mislabelled or
contaminated member on one side, not biology: the same evidence the
v0.4.0 feedback asked for (section 1a). Identity is global: 1 - exact
edit distance / the longer insert.

One row per pair, ranked by how high they disagree, with the members of
each side, and whether a side rests only on assembly contigs (WGS/TSA),
which are the likelier culprit.
"""

from __future__ import annotations

import csv
import collections
import gzip
import re
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blastocystis_subtypes import edit_distance  # noqa: E402

RANKS = ["domain", "supergroup", "division", "subdivision", "class"]
ASSEMBLY = re.compile(r"^[A-Z]{4}(?:[A-Z]{2})?\d{8,}")


def clean(v: str) -> str:
    return v.split(":")[0]


def main() -> int:
    rel, paf, out = sys.argv[1:4]
    min_id = float(sys.argv[4]) if len(sys.argv) > 4 else 0.99
    with open(f"{rel}/taxonomy.tsv") as fh:
        tax = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    members = collections.defaultdict(list)
    with open(f"{rel}/members.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["role"] == "kept":
                members[r["id"]].append(r)

    seqs, name = {}, None
    for line in open(f"{rel}/refdb.fasta"):
        if line.startswith(">"):
            name = line[1:].split()[0]
        else:
            seqs[name] = seqs.get(name, "") + line.strip()
    # Only pairs whose labels disagree at class or above are reported, so
    # only they are measured: the taxonomy test is cheap, the distance is
    # not, and candidates run to tens of millions.
    def disagreement(a, b):
        ta, tb = tax[a], tax[b]
        for level, rank in enumerate(RANKS):
            x, y = clean(ta[rank]), clean(tb[rank])
            if x and y and x != y:
                return level, rank
        return None

    seen, rows, n_cand = set(), [], 0
    with gzip.open(paf, "rt") as fh:
        for line in fh:
            q, t = line.rstrip("\n").split("\t")[:2]
            n_cand += 1
            if q == t:
                continue
            d = disagreement(q, t)
            key = (q, t) if q < t else (t, q)
            if d is None or key in seen:
                continue
            seen.add(key)
            a, b = seqs[q], seqs[t]
            longer = max(len(a), len(b))
            if abs(len(a) - len(b)) > (1 - min_id) * longer:
                continue
            edits = edit_distance(a, b)
            if 1 - edits / longer >= min_id:
                rows.append((d[0], d[1], *key, edits))

    def side(i):
        t = tax[i]
        ms = members[i]
        contig_only = all(ASSEMBLY.match(m["accession"]) for m in ms)
        label = ";".join(clean(t[r]) for r in RANKS + ["genus", "species"] if t[r])
        accs = ",".join(m["accession"] for m in ms[:4]) + (",..." if len(ms) > 4 else "")
        return [label, str(len(ms)), "yes" if contig_only else "no", accs]

    rows.sort()
    with open(out, "w") as fh:
        fh.write("disagree_at\tedits\tid_a\ttaxonomy_a\tmembers_a\tcontigs_only_a\taccessions_a"
                 "\tid_b\ttaxonomy_b\tmembers_b\tcontigs_only_b\taccessions_b\n")
        for level, rank, a, b, edits in rows:
            fh.write("\t".join([rank, str(edits), a, *side(a), b, *side(b)]) + "\n")
    print(f"{n_cand:,} candidate pairs, {len(seen):,} labelled differently at class or "
          f"above; {len(rows):,} of those >= {min_id:.0%} identical:",
          dict(collections.Counter(r[1] for r in rows)), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
