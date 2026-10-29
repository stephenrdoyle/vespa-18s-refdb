#!/usr/bin/env python3
"""What happens to reads from a species the database does not hold?

    bin/simulate_leave_one_out.py <release dir> <pipeline bin dir> <outdir> [workers] [--present]

With --present nothing is held out: the same reads, scored against the
release as it is, give the distances of correct calls - the other half
of any cut-off. Output files then end _present.

The best-of-many simulation assumed the read's own record is in the
database. Most sampled organisms are not: 8,200 named species have no
V4 record at all (docs/DESIGN.md §8). This removes each species in turn
and asks where its reads land - at the right genus (the best possible),
at the right family or above, or confidently somewhere wrong.

Held out: every species-level record of one species. Species whose V4
is also in a merged, higher-rank record are skipped (removing them would
need relabelling that record). Tested: helminths (Nematoda,
Platyhelminthes, Acanthocephala) and Apicomplexa, one record per
species, picked at random.

Reads, as in simulate_best_of_many.py, scored with the pipeline's edit
distance (ampseqlib.bounded_edit) against candidate records that
minimap2 proposes (whole insert, up to 500 hits) minus the held-out
species:
  nanopore exact  the whole insert, no errors; ties exact
  nanopore        2% errors; tolerance 2
  Illumina        the two 2x150 windows, 0.5% substitutions; ties exact
The call is the deepest rank all best hits share (the pipeline's LCA).

Outcomes, judged against the held-out species' own ranks:
  sister species   named to species - another species of its genus
  WRONG species    named to a species of another genus
  correct genus    stops at its genus (the best possible when the genus
                   has other species left)
  correct higher   stops at family or above, consistent with its ranks
  WRONG genus      named to another genus
  WRONG higher     stops above genus on a lineage it does not belong to
Each read also carries its distance to the best hit, so a cut-off on
identity can be judged.

Writes <outdir>/leave_one_out.tsv (one row per read) and
leave_one_out.txt (the summary).
"""

from __future__ import annotations

import collections
import csv
import os
import random
import subprocess
import sys
import tempfile
import zlib
from concurrent.futures import ProcessPoolExecutor

HELMINTHS = ("Nematoda", "Platyhelminthes", "Acanthocephala")
W1, W2 = 132, 128
ILLUMINA_SUB = 0.005
NANO = (0.008, 0.006, 0.006)
RANKS = ["domain", "supergroup", "division", "subdivision", "class",
         "order", "family", "genus", "species"]

_EDIT = None
_TAX = _SEQS = None


def group(t):
    c = t["class"].split(":")[0]
    if c in HELMINTHS:
        return c
    if t["subdivision"].split(":")[0] == "Apicomplexa":
        return "Apicomplexa"
    return None


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
        out.append(rng.choice([b for b in "ACGT" if b != c]) if r < dele + sub else c)
        if rng.random() < ins:
            out.append(rng.choice("ACGT"))
    return "".join(out)


def informative(name):
    return bool(name) and not name.endswith("_sp.") and "_X" not in name


def lca(ids):
    """(rank index, name) of the deepest informative rank all share."""
    deepest = (-1, "")
    for k, r in enumerate(RANKS):
        names = {_TAX[i][r].split(":")[0] for i in ids}
        if len(names) != 1:
            break
        n = next(iter(names))
        if not n:
            break
        if informative(n) or r != "species":
            deepest = (k, n)
    return deepest


def judge(truth, k, name):
    g, s = RANKS.index("genus"), RANKS.index("species")
    if k == s:
        if name == truth["species"].split(":")[0]:
            return "correct species"
        if truth["genus"].split(":")[0] == _TAX_GENUS_OF[name]:
            return "sister species"
        return "WRONG species"
    if k == g:
        return "correct genus" if name == truth["genus"].split(":")[0] else "WRONG genus"
    if k < 0:
        return "correct higher"
    return ("correct higher" if truth[RANKS[k]].split(":")[0] == name else "WRONG higher")


_TAX_GENUS_OF = {}


def init(edit_dir, tax, seqs):
    global _EDIT, _TAX, _SEQS, _TAX_GENUS_OF
    sys.path.insert(0, edit_dir)
    from ampseqlib import bounded_edit
    _EDIT, _TAX, _SEQS = bounded_edit, tax, seqs
    _TAX_GENUS_OF = {t["species"]: t["genus"].split(":")[0] for t in tax.values()}


def run(job):
    rid, sp, cands, seed = job
    rng = random.Random(seed)
    truth, s = _TAX[rid], _SEQS[rid]
    out = []
    plats = (("nanopore exact", 0, s), ("nanopore", 2, mutate_nano(s, rng)),
             ("illumina", 0, (mutate_sub(s[:W1], ILLUMINA_SUB, rng),
                              mutate_sub(s[-W2:], ILLUMINA_SUB, rng))))
    for plat, tol, read in plats:
        if plat == "illumina":
            r1, r2 = read
            dist = {i: _EDIT(r1, _SEQS[i][:W1], 40) + _EDIT(r2, _SEQS[i][-W2:], 40)
                    for i in cands}
            length = W1 + W2
        else:
            budget = max(60, len(read) // 5)
            dist = {i: _EDIT(read, _SEQS[i], budget) for i in cands}
            length = len(read)
        if not dist:
            out.append((rid, sp, group(truth), plat, "no candidate", "", "", ""))
            continue
        best = min(dist.values())
        hits = [i for i, d in dist.items() if d <= best + tol]
        k, name = lca(hits)
        out.append((rid, sp, group(truth), plat, judge(truth, k, name),
                    RANKS[k] if k >= 0 else "none", name, f"{best / length:.4f}"))
    return out


def main() -> int:
    args = [a for a in sys.argv[1:] if a != "--present"]
    present = "--present" in sys.argv
    rel, edit_dir, outdir = args[:3]
    workers = int(args[3]) if len(args) > 3 else os.cpu_count()
    tag = "_present" if present else ""
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(rel, "taxonomy.tsv")) as fh:
        tax = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    seqs, cur = {}, None
    with open(os.path.join(rel, "refdb.fasta")) as fh:
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].strip()
            else:
                seqs[cur] = seqs.get(cur, "") + line.strip()

    # Species whose sequence also sits in a merged record are skipped.
    in_merged = set()
    with open(os.path.join(rel, "members.tsv")) as fh:
        for m in csv.DictReader(fh, delimiter="\t"):
            if m["role"] == "kept" and tax[m["id"]]["deepest_rank"] != "species":
                in_merged.add(m["label"].split(";")[8].split(":")[0])
    by_sp = collections.defaultdict(list)
    for i, t in tax.items():
        if t["deepest_rank"] == "species" and group(t):
            by_sp[t["species"]].append(i)
    pick = random.Random(11)
    held = {sp: ids for sp, ids in by_sp.items() if sp not in in_merged}
    chosen = {sp: pick.choice(sorted(ids)) for sp, ids in held.items()}
    print(f"{len(chosen)} species held out one at a time "
          f"({len(by_sp) - len(held)} skipped: their V4 is also in a merged record)",
          file=sys.stderr)

    # Candidates: minimap2, whole insert against the whole release.
    with tempfile.TemporaryDirectory() as tmp:
        q = os.path.join(tmp, "q.fa")
        with open(q, "w") as fh:
            for sp, rid in chosen.items():
                fh.write(f">{rid}\n{seqs[rid]}\n")
        p = subprocess.run(["minimap2", "-t", str(workers), "-k", "13", "-w", "5",
                            "-N", "500", "-p", "0.3", os.path.join(rel, "refdb.fasta"), q],
                           capture_output=True, text=True, check=True)
    cands = collections.defaultdict(set)
    for line in p.stdout.splitlines():
        f = line.split("\t")
        cands[f[0]].add(f[5])
    jobs = []
    for sp, rid in chosen.items():
        own = set() if present else set(by_sp[sp])
        jobs.append((rid, sp, sorted(cands[rid] - own), zlib.crc32(sp.encode())))

    rows = []
    with ProcessPoolExecutor(workers, initializer=init, initargs=(edit_dir, tax, seqs)) as ex:
        for k, out in enumerate(ex.map(run, jobs, chunksize=20), 1):
            rows.extend(out)
            if k % 500 == 0:
                print(f"  {k}/{len(jobs)}", file=sys.stderr, flush=True)

    with open(os.path.join(outdir, f"leave_one_out{tag}.tsv"), "w") as fh:
        fh.write("record\tspecies\tgroup\tplatform\toutcome\tcall_rank\tcall\tbest_distance\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    outcomes = ["correct species", "sister species", "WRONG species", "correct genus", "correct higher",
                "WRONG genus", "WRONG higher"]
    genus_n = collections.Counter(t["genus"] for t in tax.values())
    genus_left = {sp: genus_n[tax[rid]["genus"]] > len(by_sp[sp]) for sp, rid in chosen.items()}
    lines = ["Leave-one-out: each species removed in turn; where do its reads land?",
             "(percent of species; 'genus left' = its genus still has other records)", ""]
    for plat in ("nanopore exact", "nanopore", "illumina"):
        lines.append(plat)
        for g in ("Nematoda", "Platyhelminthes", "Acanthocephala", "Apicomplexa"):
            for left in (True, False):
                sel = [r for r in rows if r[3] == plat and r[2] == g and genus_left[r[1]] == left]
                if not sel:
                    continue
                c = collections.Counter(r[4] for r in sel)
                cells = "  ".join(f"{o} {100 * c[o] / len(sel):4.1f}%" for o in outcomes)
                lines.append(f"  {g:16s} genus {'left ' if left else 'GONE '} n={len(sel):5d}  {cells}")
        lines.append("")
    text = "\n".join(lines)
    open(os.path.join(outdir, f"leave_one_out{tag}.txt"), "w").write(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
