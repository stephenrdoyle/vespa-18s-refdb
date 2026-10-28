#!/usr/bin/env python3
"""Cut a release from a collapsed build (docs/DESIGN.md, release format).

    bin/release.py <build dir> <releases dir> <version> <layer.tsv.gz> ... \
        [--read-length 150] [--workers 8]

<build dir> is collapse.py's output; the layers are the ones it was built
from, in the same order, so the manifest can fingerprint them. Writes
releases/vespa18s_<version>/ and a .tar.gz of it next to it.

Adds to the collapsed files:
  resolution_nanopore.tsv      the rank each record reaches read whole
  resolution_illumina_RxR.tsv  the rank it reaches through the two read
                               windows: records whose windows are identical
                               are indistinguishable, and share the merge
                               of their taxonomies
  anchor_lookalikes.tsv        adapter-anchor near-copies, in pre-flight's
                               own columns, so the pipeline reads it as is
  manifest.json                versions, inputs, counts, checksums, commit
  summary.txt                  the same, for a person
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
from collections import Counter, defaultdict
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from collapse import merge  # noqa: E402
from seqmatch import ANCHORS, IUPAC, LOOKALIKE_MAX_EDITS, find_anchor_lookalikes  # noqa: E402
from vespa import FWD, PRIMER_EDITS, RANKS, REV, REV_RC, informative  # noqa: E402

COPIED = ["refdb.fasta", "taxonomy.tsv", "members.tsv", "conflicts.tsv"]
HELMINTHS = ("Nematoda", "Platyhelminthes", "Acanthocephala")


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_fasta(path: str) -> dict:
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip()
            else:
                seqs[name] = line.strip()
    return seqs


def deepest(ranks) -> tuple:
    i = max((k for k, r in enumerate(ranks) if informative(r)), default=-1)
    return (RANKS[i], ranks[i]) if i >= 0 else ("", "")


def concrete(primer: str) -> str:
    # One representative base per degenerate position. The look-alike
    # context only needs to sit where a read's primer sits.
    return "".join(IUPAC[c][0] for c in primer)


def lookalike_chunk(items):
    return find_anchor_lookalikes(items, ANCHORS, LOOKALIKE_MAX_EDITS)


def git_commit() -> str:
    try:
        rev = subprocess.run(["git", "-C", HERE, "rev-parse", "HEAD"],
                             capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", HERE, "status", "--porcelain", "--", "bin"],
                               capture_output=True, text=True).stdout.strip()
        return rev + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("build")
    ap.add_argument("releases")
    ap.add_argument("version")
    ap.add_argument("layers", nargs="+")
    ap.add_argument("--read-length", type=int, default=150)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args()

    name = f"vespa18s_{a.version}"
    out = os.path.join(a.releases, name)
    if os.path.exists(out):
        sys.exit(f"{out} exists: a release is never overwritten")
    os.makedirs(out)
    for f in COPIED:
        shutil.copy2(os.path.join(a.build, f), out)

    seqs = read_fasta(os.path.join(out, "refdb.fasta"))
    tax = {}
    with open(os.path.join(out, "taxonomy.tsv")) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            tax[f[0]] = [f[head.index(r)] for r in RANKS]

    # --- nanopore: the whole insert, so a record is exactly itself
    with open(os.path.join(out, "resolution_nanopore.tsv"), "w") as fh:
        fh.write("id\tdeepest_rank\tname\n")
        for rid in seqs:
            r, n = deepest(tax[rid])
            fh.write(f"{rid}\t{r}\t{n}\n")

    # --- Illumina: only the two read windows are seen
    w1, w2 = a.read_length - len(FWD), a.read_length - len(REV)
    groups = defaultdict(list)
    for rid, s in seqs.items():
        groups[s if len(s) <= w1 + w2 else (s[:w1], s[-w2:])].append(rid)
    ill_name = f"resolution_illumina_{a.read_length}x{a.read_length}.tsv"
    ill_rank = {}
    with open(os.path.join(out, ill_name), "w") as fh:
        fh.write("id\twindow_group\trecords_in_group\tdeepest_rank\tname\n")
        for key, ids in groups.items():
            gid = "W_" + hashlib.sha1(repr(key).encode()).hexdigest()[:12]
            merged = merge([{"ranks": tax[i]} for i in ids])
            r, n = deepest(merged)
            for rid in ids:
                ill_rank[rid] = r
                fh.write(f"{rid}\t{gid}\t{len(ids)}\t{r}\t{n}\n")

    # --- anchor look-alikes, scanned on primer + insert + primer as read
    f5, r3 = concrete(FWD), concrete(REV_RC)
    items = [(rid, f5 + s + r3) for rid, s in seqs.items()]
    chunks = [items[i:i + 250] for i in range(0, len(items), 250)]
    with Pool(a.workers) as pool:
        rows = [r for part in pool.map(lookalike_chunk, chunks) for r in part]
    rows.sort(key=lambda r: (r["seq_id"], r["anchor"], r["strand"], r["position"]))
    cols = ["seq_id", "anchor", "strand", "position", "edits", "context"]
    with open(os.path.join(out, "anchor_lookalikes.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # --- counts
    n_members = Counter()
    with open(os.path.join(out, "members.tsv")) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            n_members[(f[1], f[4])] += 1
    conflicts = Counter()
    with open(os.path.join(out, "conflicts.tsv")) as fh:
        fh.readline()
        for line in fh:
            conflicts[line.split("\t")[1]] += 1
    nano = Counter(deepest(tax[i])[0] for i in seqs)
    ill = Counter(ill_rank.values())
    helm = [i for i in seqs if tax[i][4] in HELMINTHS]
    helm_nano = Counter(deepest(tax[i])[0] for i in helm)
    helm_ill = Counter(ill_rank[i] for i in helm)
    look = Counter(r["edits"] for r in rows)

    files = sorted(os.listdir(out))
    manifest = {
        "name": name,
        "version": a.version,
        "built": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "code_commit": git_commit(),
        "marker": "18S rRNA V4",
        "primers": {"forward": {"name": "AMPSEQ_VESPA18S_29F", "sequence": FWD},
                    "reverse": {"name": "AMPSEQ_VESPA18S_21b8_I_R_Nmod", "sequence": REV},
                    "max_edits_in_silico": PRIMER_EDITS},
        "insert": "between the primers; primers not included",
        "illumina_read_length": a.read_length,
        "illumina_windows": [w1, w2],
        "anchors": ANCHORS,
        "lookalike_max_edits": LOOKALIKE_MAX_EDITS,
        "layers": [{"file": os.path.basename(p), "sha256": sha256(p),
                    "trust_order": k} for k, p in enumerate(a.layers)],
        "sources": sorted({s for s, _ in n_members}),
        "counts": {
            "records": len(seqs),
            "members": {f"{s} {role}": n for (s, role), n in sorted(n_members.items())},
            "conflicts": dict(conflicts),
            "deepest_rank_nanopore": dict(nano),
            "deepest_rank_illumina": dict(ill),
            "helminth_records": len(helm),
            "illumina_window_groups": len(groups),
            "anchor_lookalikes_by_edits": {str(k): v for k, v in sorted(look.items())},
        },
        "files": {f: sha256(os.path.join(out, f)) for f in files},
    }
    with open(os.path.join(out, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=1)

    def table(c: Counter, n: int) -> list:
        return [f"    {r or '(none)':12s} {c.get(r, 0):>8,}  {100 * c.get(r, 0) / max(1, n):5.1f}%"
                for r in ["species", "genus", "family", "order", "class",
                          "subdivision", "division", "supergroup", "domain", ""]]

    lines = [f"{name}  built {manifest['built']}  code {manifest['code_commit']}",
             f"primers {FWD} / {REV}",
             "layers: " + ", ".join(f"{l['file']}" for l in manifest["layers"]),
             f"records {len(seqs):,}; Illumina {a.read_length}x{a.read_length} window "
             f"groups {len(groups):,}",
             "members: " + ", ".join(f"{k} {v:,}" for k, v in manifest["counts"]["members"].items()),
             "conflicts: " + ", ".join(f"{k} {v}" for k, v in sorted(conflicts.items())),
             "", "Deepest rank, all records - nanopore:", *table(nano, len(seqs)),
             "  Illumina:", *table(ill, len(seqs)),
             "", f"Helminth records ({len(helm):,}) - nanopore:", *table(helm_nano, len(helm)),
             "  Illumina:", *table(helm_ill, len(helm)),
             "", "Anchor look-alikes by edits: "
             + ", ".join(f"{k}: {v:,}" for k, v in sorted(look.items())),
             "", "Files:", *[f"  {f}  {manifest['files'][f][:16]}" for f in files]]
    with open(os.path.join(out, "summary.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")

    tar = out + ".tar.gz"
    with tarfile.open(tar, "w:gz") as t:
        t.add(out, arcname=name)
    print("\n".join(lines))
    print(f"\n{tar}  {os.path.getsize(tar) / 1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
