#!/usr/bin/env python3
"""Assign Blastocystis V4 inserts to subtypes, by the curated references.

    bin/blastocystis_subtypes.py curation/blastocystis_subtypes.tsv \
        work/genbank/Blastocystis.fasta.gz work/pr2_layer.tsv.gz work/blastocystis \
        curation/blastocystis_v4_subtypes.tsv

Subtypes (ST1, ST2 ...) are not species and PR2 and GenBank barely label
them, so they are set here from curation/blastocystis_subtypes.tsv: the
published full-length reference sequences for each accepted subtype.

1. The references are fetched from GenBank and trimmed to the V4 insert.
   They are also written as OUTDIR/references_genbank.fasta.gz, in
   genbank_fetch.py's format, so bin/genbank_layer.py can make them a
   layer: some never reach a title search, and every subtype should have
   at least its own reference as a record.
2. Every Blastocystis insert we hold - the references, GenBank records
   that span both primers, and PR2's - is compared with every reference and
   given its nearest subtype, the distance to it (d1), and the distance
   to the nearest *other* subtype (d2): exact global edit distances.
3. Where a GenBank title names a subtype, it is compared with the
   nearest one.

Writes OUTDIR/references_v4.fasta, OUTDIR/assignments.tsv (one row per
distinct insert, by record id) and prints a summary. The `assigned`
column is the subtype given (see MAX_D1, MIN_MARGIN); blank means
Blastocystis sp. The assigned rows alone go to the curation table (the
last argument), which bin/collapse.py applies: there they can be read and
corrected by hand before a release.
"""

from __future__ import annotations

import collections
import csv
import gzip
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collapse import record_id  # noqa: E402
from vespa import pcr_full  # noqa: E402

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
# A subtype is given when the nearest reference subtype is within MAX_D1
# of the insert's length and every other subtype is at least MIN_MARGIN
# edits further. V4 varies a lot inside some subtypes (ST3 titles run to
# 7.2% from the ST3 references) but subtypes stay well apart (the same
# ST3 inserts are >= 11% from any other), so a wide distance with a firm
# margin. Set from the v0.6.0 data, where it contradicts no GenBank title
# in the current numbering.
MAX_D1, MIN_MARGIN = 0.08, 10
TITLE_ST = re.compile(r"\b(?:subtype|ST)[\s:_-]*(\d{1,2})\b", re.I)


def fasta(text):
    name, buf = None, []
    for line in text:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield name, "".join(buf).upper()
            name, buf = line[1:], []
        elif line:
            buf.append(line)
    if name:
        yield name, "".join(buf).upper()


def eutils(endpoint: str, **params) -> str:
    data = urllib.parse.urlencode(params).encode()
    for attempt in range(6):
        try:
            return urllib.request.urlopen(EUTILS + endpoint, data, timeout=120).read().decode()
        except Exception:
            time.sleep(5 * (attempt + 1))
    sys.exit(f"could not reach NCBI ({endpoint})")


def fetch(accs: list, path: str) -> None:
    """The references as FASTA, headers "accession taxid=N title"."""
    if os.path.exists(path):
        return
    text = eutils("efetch.fcgi", db="nuccore", id=",".join(accs), rettype="fasta",
                  retmode="text")
    summ = json.loads(eutils("esummary.fcgi", db="nuccore", id=",".join(accs),
                             retmode="json"))["result"]
    taxid = {summ[u]["accessionversion"]: summ[u]["taxid"] for u in summ["uids"]}
    with gzip.open(path, "wt") as fh:
        for name, seq in fasta(text.splitlines()):
            acc, title = name.split(" ", 1)
            fh.write(f">{acc} taxid={taxid[acc]} {title}\n{seq}\n")


def edit_distance(a: str, b: str) -> int:
    """Global (Levenshtein) distance, Myers' bit-parallel algorithm with
    Python integers as bit vectors: one pass over b."""
    m = len(a)
    if not m:
        return len(b)
    peq = collections.defaultdict(int)
    for i, c in enumerate(a):
        peq[c] |= 1 << i
    full, top = (1 << m) - 1, 1 << (m - 1)
    pv, mv, score = full, 0, m
    for c in b:
        eq = peq[c]
        xv = eq | mv
        xh = (((eq & pv) + pv) ^ pv) | eq
        ph = mv | (~(xh | pv) & full)
        mh = pv & xh
        if ph & top:
            score += 1
        elif mh & top:
            score -= 1
        ph = ((ph << 1) | 1) & full
        mh = (mh << 1) & full
        pv = mh | (~(xv | ph) & full)
        mv = ph & xv
    return score


def title_subtype(title: str) -> str:
    sts = {f"ST{int(m)}" for m in TITLE_ST.findall(title)}
    return sts.pop() if len(sts) == 1 else ("ambiguous" if sts else "")


def main() -> int:
    table, gb_path, pr2_path, out, curated = sys.argv[1:6]
    os.makedirs(out, exist_ok=True)
    refs = [r for r in csv.DictReader(open(table), delimiter="\t") if r["status"] == "accepted"]
    ref_path = os.path.join(out, "references_genbank.fasta.gz")
    fetch([r["accession"] for r in refs], ref_path)
    st_of = {r["accession"].split(".")[0]: r["subtype"] for r in refs}
    ref_v4 = {}
    for name, seq in fasta(gzip.open(ref_path, "rt")):
        acc = name.split()[0].split(".")[0]
        status, ins, *_ = pcr_full(seq)
        if status == "amplicon":
            ref_v4[f"{st_of[acc]}|{acc}"] = ins
        else:
            print(f"  reference {acc} ({st_of[acc]}): no V4 insert ({status})", file=sys.stderr)
    with open(os.path.join(out, "references_v4.fasta"), "w") as fh:
        for k, s in ref_v4.items():
            fh.write(f">{k}\n{s}\n")

    # Every insert we hold, with where it came from and any subtype title.
    seen = collections.defaultdict(list)     # insert -> [(source, accession, title_st)]
    for k, ins in ref_v4.items():
        st, acc = k.split("|")
        seen[ins].append(("reference", acc, st))
    for name, seq in fasta(gzip.open(gb_path, "rt")):
        status, ins, *_ = pcr_full(seq)
        if status == "amplicon":
            acc = name.split()[0]
            seen[ins].append(("GenBank", acc, title_subtype(name)))
    with gzip.open(pr2_path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["genus"].startswith("Blastocystis") and r["insert"]:
                seen[r["insert"]].append(("PR2", r["accession"], ""))
    ids = {ins: record_id(ins) for ins in seen}

    # Exact global edit distance to every reference.
    dist = {i: {ref: edit_distance(ins, r) for ref, r in ref_v4.items()}
            for ins, i in ids.items()}

    rows = []
    for ins, i in ids.items():
        per_st = {}
        for ref, e in dist[i].items():
            st = ref.split("|")[0]
            per_st[st] = min(e, per_st.get(st, 10**6))
        ranked = sorted(per_st.items(), key=lambda x: x[1])
        st1, d1 = ranked[0] if ranked else ("", "")
        st2, d2 = ranked[1] if len(ranked) > 1 else ("", "")
        titles = sorted({t for _, _, t in seen[ins] if t})
        assigned = st1 if d1 <= MAX_D1 * len(ins) and d2 - d1 >= MIN_MARGIN else ""
        rows.append([i, str(len(ins)), assigned, st1, str(d1), st2, str(d2), ",".join(titles),
                     str(len(seen[ins])), ",".join(sorted({s for s, _, _ in seen[ins]})),
                     ",".join(a for _, a, _ in seen[ins][:6]) + (",..." if len(seen[ins]) > 6 else ""),
                     ins])
    with open(os.path.join(out, "assignments.tsv"), "w") as fh:
        fh.write("id\tlength\tassigned\tnearest_st\td1\tnext_st\td2\ttitle_st\tmembers\tsources"
                 "\taccessions\tinsert\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    with open(curated, "w") as fh:
        fh.write("# Blastocystis subtype of each V4 record, written by\n"
                 "# bin/blastocystis_subtypes.py from curation/blastocystis_subtypes.tsv\n"
                 f"# (nearest reference subtype within {MAX_D1:.0%} of the insert, every\n"
                 f"# other subtype >= {MIN_MARGIN} edits further). bin/collapse.py names\n"
                 "# the Blastocystis members of these records Blastocystis_STn. Edit by\n"
                 "# hand with care: rerunning the script rewrites it.\n")
        fh.write("id\tsubtype\td1\tnext_st\td2\tlength\n")
        for r in sorted(rows):
            if r[2]:
                fh.write("\t".join([r[0], r[2], r[4], r[5], r[6], r[1]]) + "\n")

    agree = sum(1 for r in rows if r[7] and "," not in r[7] and r[7] == r[3])
    titled = sum(1 for r in rows if r[7] and "," not in r[7] and r[7] != "ambiguous")
    print(f"{len(ref_v4)} reference V4 inserts for "
          f"{len({k.split('|')[0] for k in ref_v4})} subtypes; {len(rows):,} distinct "
          f"Blastocystis inserts; {titled:,} with one subtype in a title, "
          f"{agree:,} of them nearest that subtype; {sum(1 for r in rows if r[2]):,} "
          f"given a subtype ({sum(int(r[8]) for r in rows if r[2]):,} members)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
