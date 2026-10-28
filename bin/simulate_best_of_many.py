#!/usr/bin/env python3
"""Does a densely sampled species pull reads away from a sparse sister?

    bin/simulate_best_of_many.py <release dir> <pipeline bin dir> <outdir> [workers]

The "best of many" question: a species with hundreds of V4 variants
offers hundreds of chances to sit close to any read, so a read from a
sister species with one record may tie with, or lose to, one of them. The
nearest-hit-plus-LCA rule the pipeline will use is simulated here, read
by read, with the pipeline's own edit distance (ampseqlib.bounded_edit).

Test cases are genera where a dense species (>= DENSE species-level
records) has a sparse sister (<= SPARSE records) within NEAR edits of one
of its variants. Reads are simulated from the sparse species' records -
the ones at risk - and from a sample of the dense species' records.

  Illumina  the two 2x150 read windows only, 0.5% substitutions; a
            read pair is scored window by window against each record's
            windows; ties exact (tolerance 0)
  nanopore  the whole insert, 2% errors (0.8% substitution, 0.6%
            insertion, 0.6% deletion); tolerance 0 and 2

A read is scored against every record of the genus; its call is the rank
all best hits (within the tolerance) share. Each is run against the
release as it is and with each species' variants merged: a variant within
MERGE edits of a more common variant of the same species is folded into
it. Best case throughout - the record a read came from is in the
database.

Writes <outdir>/best_of_many.tsv (every case x platform x database) and
best_of_many.txt (the summary).
"""

from __future__ import annotations

import collections
import csv
import os
import random
import sys
import zlib
from concurrent.futures import ProcessPoolExecutor

DENSE, SPARSE, NEAR, MERGE = 20, 3, 8, 2
READS_SPARSE, DENSE_SAMPLE, READS_DENSE = 20, 10, 3  # per record
MAX_SPARSE_RECORDS, MAX_CASES, MIN_HELMINTH_CASES = 6, 12, 4
HELMINTHS = ("Nematoda", "Platyhelminthes", "Acanthocephala")
W1, W2 = 132, 128
ILLUMINA_SUB = 0.005
NANO = (0.008, 0.006, 0.006)  # substitution, insertion, deletion per base

_EDIT = None


def edit(a, b, k):
    return _EDIT(a, b, k)


def mutate_sub(s, rate, rng):
    return "".join(rng.choice([b for b in "ACGT" if b != c]) if rng.random() < rate else c
                   for c in s)


def mutate_nano(s, rng):
    sub, ins, dele = NANO
    out = []
    for c in s:
        r = rng.random()
        if r < dele:
            continue
        if r < dele + sub:
            out.append(rng.choice([b for b in "ACGT" if b != c]))
        else:
            out.append(c)
        if rng.random() < ins:
            out.append(rng.choice("ACGT"))
    return "".join(out)


def call(best_ids, tax):
    """The rank all best hits share: 'species', 'genus' or 'above'."""
    sp = {tax[i]["species"] for i in best_ids}
    if len(sp) == 1 and next(iter(sp)):
        return "species", next(iter(sp))
    ge = {tax[i]["genus"] for i in best_ids}
    if len(ge) == 1 and next(iter(ge)):
        return "genus", next(iter(ge))
    return "above", ""


def score_case(job):
    """One case x platform x database: simulate, assign, tally."""
    (case, platform, tol, dbname, pool, reads, tax, seqs, seed, edit_dir) = job
    global _EDIT
    if _EDIT is None:
        sys.path.insert(0, edit_dir)
        from ampseqlib import bounded_edit
        _EDIT = bounded_edit
    rng = random.Random(seed)
    tally = collections.Counter()
    for kind, src in reads:
        true_sp, true_ge = tax[src]["species"], tax[src]["genus"]
        s = seqs[src]
        if platform == "illumina":
            r1 = mutate_sub(s[:W1], ILLUMINA_SUB, rng)
            r2 = mutate_sub(s[-W2:], ILLUMINA_SUB, rng)
            budget = 20
            dist = {i: edit(r1, seqs[i][:W1], budget) + edit(r2, seqs[i][-W2:], budget)
                    for i in pool}
        else:
            r = mutate_nano(s, rng)
            budget = max(30, len(r) // 10)
            dist = {i: edit(r, seqs[i], budget) for i in pool}
        best = min(dist.values())
        rank, name = call([i for i, d in dist.items() if d <= best + tol], tax)
        if rank == "species":
            out = "correct species" if name == true_sp else "WRONG species"
        elif rank == "genus":
            out = "genus only" if name == true_ge else "WRONG genus"
        else:
            out = "above genus"
        tally[(kind, out)] += 1
    return case, platform, tol, dbname, tally


def merged_pool(pool, tax, seqs, members, edit_dir):
    """Fold each species' variants into a more common variant within MERGE."""
    global _EDIT
    if _EDIT is None:
        sys.path.insert(0, edit_dir)
        from ampseqlib import bounded_edit
        _EDIT = bounded_edit
    by_sp = collections.defaultdict(list)
    keep = []
    for i in pool:
        if tax[i]["deepest_rank"] == "species":
            by_sp[tax[i]["species"]].append(i)
        else:
            keep.append(i)
    for ids in by_sp.values():
        centroids = []
        for i in sorted(ids, key=lambda i: (-members[i], i)):
            if not any(edit(seqs[i], seqs[c], MERGE) <= MERGE for c in centroids):
                centroids.append(i)
        keep.extend(centroids)
    return keep


def main() -> int:
    rel, edit_dir, outdir = sys.argv[1:4]
    workers = int(sys.argv[4]) if len(sys.argv) > 4 else os.cpu_count()
    os.makedirs(outdir, exist_ok=True)
    sys.path.insert(0, edit_dir)
    from ampseqlib import bounded_edit
    global _EDIT
    _EDIT = bounded_edit

    with open(os.path.join(rel, "taxonomy.tsv")) as fh:
        tax = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    seqs, cur = {}, None
    with open(os.path.join(rel, "refdb.fasta")) as fh:
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].strip()
            else:
                seqs[cur] = line.strip()
    members = {i: int(t["members"]) for i, t in tax.items()}

    genus_ids = collections.defaultdict(list)
    for i, t in tax.items():
        if t["genus"]:
            genus_ids[t["genus"]].append(i)

    cases = []
    for g, ids in genus_ids.items():
        sp = collections.defaultdict(list)
        for i in ids:
            if tax[i]["deepest_rank"] == "species":
                sp[tax[i]["species"]].append(i)
        dense = [s for s, v in sp.items() if len(v) >= DENSE]
        sparse = [s for s, v in sp.items() if len(v) <= SPARSE]
        for d in dense:
            near = []
            for s in sparse:
                for i in sp[s]:
                    if any(edit(seqs[i], seqs[j], NEAR) <= NEAR for j in sp[d]):
                        near.append(i)
            if near:
                cases.append((g, d, near, sp[d], ids))
    # One case per genus (its densest species), the largest first, and
    # helminth genera always among them.
    best = {}
    for c in sorted(cases, key=lambda c: -len(c[3])):
        best.setdefault(c[0], c)
    ranked = list(best.values())
    helm = [c for c in ranked if tax[c[3][0]]["class"].split(":")[0] in HELMINTHS]
    chosen = helm[:MIN_HELMINTH_CASES]
    chosen += [c for c in ranked if c not in chosen][:MAX_CASES - len(chosen)]
    cases = sorted(chosen, key=lambda c: -len(c[3]))
    pick = random.Random(7)
    cases = [(g, d, sorted(pick.sample(n, min(MAX_SPARSE_RECORDS, len(n)))), dv, ids)
             for g, d, n, dv, ids in cases]
    print(f"{len(cases)} cases: " + ", ".join(f"{g}: {d} ({len(dv)}) vs {len(n)} sparse records"
                                             for g, d, n, dv, _ in cases), file=sys.stderr)

    jobs = []
    rng = random.Random(1)
    for g, d, near, dv, ids in cases:
        case = f"{g} ({d.split('_', 1)[1]} x{len(dv)})"
        reads = [("sparse", i) for i in near for _ in range(READS_SPARSE)]
        reads += [("dense", i) for i in rng.sample(dv, min(DENSE_SAMPLE, len(dv)))
                  for _ in range(READS_DENSE)]
        merged = merged_pool(ids, tax, seqs, members, edit_dir)
        for dbname, pool in (("as released", ids), ("variants merged", merged)):
            sub_tax = {i: tax[i] for i in set(pool) | {r for _, r in reads}}
            sub_seqs = {i: seqs[i] for i in sub_tax}
            for platform, tol in (("illumina", 0), ("nanopore", 0), ("nanopore", 2)):
                jobs.append((case, platform, tol, f"{dbname} ({len(pool)} records)",
                             pool, reads, sub_tax, sub_seqs,
                             zlib.crc32(f"{g}|{d}|{platform}".encode()),
                             edit_dir))

    results = []
    with ProcessPoolExecutor(workers) as ex:
        for case, platform, tol, dbname, tally in ex.map(score_case, jobs):
            results.append((case, platform, tol, dbname, tally))
            print(f"  done: {case} {platform} tol {tol} {dbname}", file=sys.stderr, flush=True)

    outcomes = ["correct species", "genus only", "WRONG species", "WRONG genus", "above genus"]
    with open(os.path.join(outdir, "best_of_many.tsv"), "w") as fh:
        fh.write("case\tplatform\ttolerance\tdatabase\treads_from\t" + "\t".join(outcomes) + "\n")
        for case, platform, tol, dbname, tally in results:
            for kind in ("sparse", "dense"):
                fh.write("\t".join([case, platform, str(tol), dbname, kind]
                                   + [str(tally[(kind, o)]) for o in outcomes]) + "\n")
    lines = []
    for case in dict.fromkeys(r[0] for r in results):
        lines.append(f"\n{case}")
        for case_, platform, tol, dbname, tally in results:
            if case_ != case:
                continue
            for kind in ("sparse", "dense"):
                n = sum(v for (k, _), v in tally.items() if k == kind)
                cells = "  ".join(f"{o.split()[0][:7]} {100 * tally[(kind, o)] / max(1, n):3.0f}%"
                                  for o in outcomes[:3])
                lines.append(f"  {platform:8s} tol {tol}  {dbname:30s} {kind:6s} n={n:4d}  {cells}")
    text = "\n".join(lines)
    open(os.path.join(outdir, "best_of_many.txt"), "w").write(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
