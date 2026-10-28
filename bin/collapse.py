#!/usr/bin/env python3
"""Collapse layers into unique V4 inserts, merge their taxonomy, and flag
label conflicts.

    bin/collapse.py <outdir> <layer.tsv.gz> [<layer.tsv.gz> ...]

A layer is a TSV with at least: source, accession, the nine PR2 ranks,
status, insert. Rows whose status is not "amplicon" are ignored. Layers
are given in increasing order of trust - PR2 first, local sequences last -
and that order decides conflicts (below).

Writes refdb.fasta, taxonomy.tsv, members.tsv and conflicts.tsv to
<outdir>: the parts of a release that come from collapsing. Resolution
tables, the look-alike mask and the manifest are separate steps.

MERGING. Identical inserts become one record. At each rank the record
takes the single informative name its members agree on. Placeholders
(Metazoa_X, Genus_sp.) do not vote: an entry labelled only to class does
not stop its identical twin from naming the genus. Two different
informative names at a rank stop the record there - that rank and every
rank below it are left blank.

CONFLICTS. Members that disagree at class or above cannot all be right:
identical V4 across Nematoda and Platyhelminthes, or a nematode and
"bacteria", is a mislabel or a contaminant. Such a record is
  - resolved by trust: if the most trusted layer present agrees with
    itself, its label wins and the others are overruled;
  - otherwise resolved by majority: a label carried by at least two
    members and at least two thirds of them wins;
  - otherwise unresolved: kept, with every member merged, so the record
    carries only the ranks they all share (Nematoda versus
    Platyhelminthes becomes Metazoa; a plant versus a bacterium becomes
    nothing at all) and is listed for review.
Below class, a more trusted layer can also overrule - but narrowly. A
genome showing a sequence in *Haemonchus contortus* does not show it is
ONLY there: PR2's identical *H. placei* and *H. similis* are biology, and
the record should stop at genus. So species never overrule. At order,
family and genus, a less trusted label that conflicts with the trusted
one is overruled only when a single accession carries it AND that
accession either was never identified to species (PR2's lone
"*Haemonchus* sp." identical to the *Necator* genome) or is outvoted by
at least three members agreeing with the trusted label (one *Brugia*
genome-project contig among nine *C. elegans* records). Otherwise a
named species is believed, even alone: the *Teladorsagia* genome's V4 is
identical to PR2's single *Trichostrongylus colubriformis*, and
trichostrongylids really are one or two bases apart over V4 - so the
record stops at the rank they share rather than erasing a probably
genuine record.

Unresolved records are kept rather than dropped on purpose. A sequence
missing from the database does not stop reads arriving from it; they are
pushed onto the nearest remaining neighbour - a confident wrong call.
Kept at a high rank, they are assigned no deeper than the evidence goes.
Overruled members stay in members.tsv, marked, so nothing is silently
dropped.
"""

from __future__ import annotations

import hashlib
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vespa import CONFLICT_DEPTH, RANKS, informative, open_text  # noqa: E402


def record_id(insert: str) -> str:
    # From the sequence, so a record keeps its id across builds for as long
    # as its sequence exists.
    return "V4_" + hashlib.sha1(insert.encode()).hexdigest()[:12]


def high_label(ranks) -> tuple:
    """The informative names at class and above; placeholders blanked."""
    return tuple(r if informative(r) else "" for r in ranks[:CONFLICT_DEPTH + 1])


def compatible(a: tuple, b: tuple) -> bool:
    return all(x == y or not x or not y for x, y in zip(a, b))


def merge(members) -> list:
    out, stopped = [], False
    for i in range(len(RANKS)):
        if stopped:
            out.append("")
            continue
        names = {m["ranks"][i] for m in members}
        info = {n for n in names if informative(n)}
        if len(info) == 1:
            out.append(info.pop())
        elif len(info) > 1:
            out.append("")
            stopped = True
        else:
            out.append(names.pop() if len(names) == 1 else "")
    return out


def resolve(members, trust: dict):
    """(kept, overruled, decision) for one group of identical inserts."""
    labels = [high_label(m["ranks"]) for m in members]
    # Conflict = two labels that cannot both be true at class or above.
    distinct = list(dict.fromkeys(labels))
    conflict = any(not compatible(a, b)
                   for i, a in enumerate(distinct) for b in distinct[i + 1:])
    if not conflict:
        return members, [], "none"

    top = max(trust[m["source"]] for m in members)
    top_labels = {high_label(m["ranks"]) for m in members
                  if trust[m["source"]] == top}
    if top > 0 and all(compatible(a, b) for a in top_labels for b in top_labels):
        win = next(iter(top_labels))
        how = "trust"
    else:
        counts = Counter(labels)
        win, n = counts.most_common(1)[0]
        supporters = sum(c for lab, c in counts.items() if compatible(lab, win))
        if n >= 2 and supporters * 3 >= 2 * len(members):
            how = "majority"
        else:
            return members, [], "unresolved"
    kept = [m for m in members if compatible(high_label(m["ranks"]), win)]
    over = [m for m in members if not compatible(high_label(m["ranks"]), win)]
    return kept, over, how


TRUST_RANKS = range(RANKS.index("order"), RANKS.index("genus") + 1)


def trust_below_class(kept, trust: dict):
    """Overrule lone lower-trust labels conflicting at order..genus that are
    unnamed or outvoted three to one."""
    levels = {trust[m["source"]] for m in kept}
    if len(levels) < 2:
        return kept, []
    top = max(levels)
    ref = merge([m for m in kept if trust[m["source"]] == top])
    dissent = defaultdict(list)
    for m in kept:
        if trust[m["source"]] == top:
            continue
        for i in TRUST_RANKS:
            a, b = m["ranks"][i], ref[i]
            if informative(a) and informative(b) and a != b:
                dissent[(i, a)].append(m)
                break
    over = []
    for (i, _), ms in dissent.items():
        if len(ms) != 1:
            continue
        agree = sum(1 for m in kept if m["ranks"][i] == ref[i])
        if not informative(ms[0]["ranks"][-1]) or agree >= 3:
            over.append(ms[0])
    return [m for m in kept if m not in over], over


def main() -> int:
    outdir, layers = sys.argv[1], sys.argv[2:]
    os.makedirs(outdir, exist_ok=True)

    groups = defaultdict(list)
    trust, order = {}, []
    n_rows = Counter()
    for li, path in enumerate(layers):
        with open_text(path) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            col = {c: i for i, c in enumerate(head)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[col["status"]] != "amplicon" or not f[col["insert"]]:
                    continue
                src = f[col["source"]]
                if src not in trust:
                    # A layer earlier on the command line is trusted less.
                    trust[src] = li
                    order.append(src)
                n_rows[src] += 1
                groups[f[col["insert"]]].append({
                    "source": src, "accession": f[col["accession"]],
                    "ranks": [f[col[r]] for r in RANKS]})

    decisions = Counter()
    with open(os.path.join(outdir, "refdb.fasta"), "w") as fa, \
         open(os.path.join(outdir, "taxonomy.tsv"), "w") as tx, \
         open(os.path.join(outdir, "members.tsv"), "w") as mb, \
         open(os.path.join(outdir, "conflicts.tsv"), "w") as cf:
        tx.write("\t".join(["id", *RANKS, "deepest_rank", "members", "sources"]) + "\n")
        mb.write("id\tsource\taccession\tlabel\trole\n")
        cf.write("id\tdecision\tmembers\tlabels\tkept_label\n")
        for insert in sorted(groups, key=record_id):
            members = groups[insert]
            rid = record_id(insert)
            kept, over, how = resolve(members, trust)
            kept, over2 = trust_below_class(kept, trust)
            if over2:
                over += over2
                how = "trust" if how == "none" else how
            decisions[how] += 1
            for m in members:
                role = "kept" if m in kept else "overruled"
                mb.write(f"{rid}\t{m['source']}\t{m['accession']}\t"
                         f"{';'.join(m['ranks'])}\t{role}\n")
            if how != "none":
                labs = Counter(";".join(x for x in high_label(m["ranks"]))
                               for m in members)
                kept_lab = (";".join(high_label(kept[0]["ranks"])) if kept else "")
                cf.write(f"{rid}\t{how}\t{len(members)}\t"
                         + " | ".join(f"{k} x{v}" for k, v in labs.most_common())
                         + f"\t{kept_lab}\n")
            if not kept:
                continue
            ranks = merge(kept)
            deepest = max((i for i, r in enumerate(ranks) if informative(r)),
                          default=-1)
            srcs = sorted({m["source"] for m in kept}, key=order.index)
            tx.write("\t".join([rid, *ranks,
                                RANKS[deepest] if deepest >= 0 else "",
                                str(len(kept)), ",".join(srcs)]) + "\n")
            fa.write(f">{rid}\n{insert}\n")

    print(f"layers: " + ", ".join(f"{s} ({n_rows[s]:,} amplicons)" for s in order))
    print(f"unique inserts: {len(groups):,}")
    for k in ("none", "trust", "majority", "unresolved"):
        print(f"  conflict {k:9s} {decisions[k]:>7,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
