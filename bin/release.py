#!/usr/bin/env python3
"""Cut a release from a collapsed build (docs/DESIGN.md, release format).

    bin/release.py <build dir> <releases dir> <version> <layer.tsv.gz> ... \
        [--illumina 150] [--illumina 300:250x200] [--workers 8]

<build dir> is collapse.py's output; the layers are the ones it was built
from, in the same order, so the manifest can fingerprint them. Run
bin/species_margins.py on the build first (same geometries) to fill the
margin columns; without it they are blank. Writes
releases/vespa18s_<version>/ and a .tar.gz of it next to it.

Adds to the collapsed files:
  resolution_nanopore.tsv      the rank each record reaches read whole, and
                               for species-level records the margin: edits
                               to the nearest record of another species
                               (blank on a species row: none near)
  resolution_illumina_RxR.tsv  the rank it reaches through the two read
                               windows: records whose windows are identical
                               are indistinguishable, and share the merge
                               of their taxonomies; margin as above, over
                               the windows
  rank_names.tsv               what the positional taxonomy columns hold, by
                               clade (curation/rank_names.tsv)
  anchor_lookalikes.tsv        adapter-anchor near-copies, in pre-flight's
                               own columns, so the pipeline reads it as is
  manifest.json                versions, inputs, counts, checksums, commit
(primer_fit.tsv comes from collapse.py with the other collapsed files.)
  summary.txt                  the same, for a person
"""

from __future__ import annotations

import argparse
import csv
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
from seqmatch import ANCHORS, IUPAC, LOOKALIKE_MAX_EDITS, find_anchor_lookalikes  # noqa: E402
from vespa import FWD, PRIMER_EDITS, RANKS, REV, REV_RC, informative  # noqa: E402

COPIED = ["refdb.fasta", "taxonomy.tsv", "members.tsv", "conflicts.tsv", "primer_fit.tsv"]
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


def merge_records(records) -> list:
    """The common ancestor of records the read windows cannot tell apart.

    Strict: a rank is kept only if every record has exactly the same value
    there. Not collapse.merge, where placeholders abstain - that merges the
    labels of ONE sequence. These are DIFFERENT sequences, and a read
    matching their shared windows may come from any of them: if one names
    only its family (Eurotiomycetes_X) or was left blank by a conflict, the
    read cannot be assigned the genus another names. Merged loosely, 3,512
    v0.1.0 records were reported deeper through the windows than whole.
    """
    out, stopped = [], False
    for i in range(len(RANKS)):
        vals = {r[i] for r in records}
        if stopped or len(vals) != 1:
            stopped = True
            out.append("")
        else:
            out.append(vals.pop())
    return out


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
    ap.add_argument("--illumina", action="append", metavar="READLEN[:W1xW2]",
                    help="an Illumina read geometry to tabulate: read length, and "
                         "optionally the insert bases each read covers (default: "
                         "read length minus the primer). Repeatable; the first is "
                         "the summary's. Default 150")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args()

    name = f"vespa18s_{a.version}"
    out = os.path.join(a.releases, name)
    if os.path.exists(out):
        sys.exit(f"{out} exists: a release is never overwritten")
    os.makedirs(out)
    for f in COPIED:
        shutil.copy2(os.path.join(a.build, f), out)
    # What the positional taxonomy columns hold, by clade (curated).
    shutil.copy2(os.path.join(HERE, os.pardir, "curation", "rank_names.tsv"), out)
    # How firmly each species-level record names its species: from
    # bin/species_margins.py, run on the same build. Blank when absent.
    margins = {}
    mpath = os.path.join(a.build, "species_margins.tsv")
    if os.path.exists(mpath):
        with open(mpath) as fh:
            margins = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    seqs = read_fasta(os.path.join(out, "refdb.fasta"))
    tax = {}
    with open(os.path.join(out, "taxonomy.tsv")) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            tax[f[0]] = [f[head.index(r)] for r in RANKS]

    # --- nanopore: the whole insert, so a record is exactly itself
    with open(os.path.join(out, "resolution_nanopore.tsv"), "w") as fh:
        # margin: edits to the nearest record of another species, read
        # whole - how many bases a species call rests on. Species rows only.
        fh.write("id\tdeepest_rank\tname\tmargin\n")
        for rid in seqs:
            r, n = deepest(tax[rid])
            m = margins.get(rid, {}).get("margin_whole", "") if r == "species" else ""
            fh.write(f"{rid}\t{r}\t{n}\t{m}\n")

    # --- Illumina: only the two read windows are seen, one table per read
    # geometry. The first geometry is the summary's.
    order = {r: k for k, r in enumerate(RANKS)}
    geometries = []
    for spec in a.illumina or ["150"]:
        read_len, _, win = spec.partition(":")
        read_len = int(read_len)
        if win:
            w1, w2 = (int(x) for x in win.split("x"))
        else:
            w1, w2 = read_len - len(FWD), read_len - len(REV)
        groups = defaultdict(list)
        for rid, s in seqs.items():
            groups[s if len(s) <= w1 + w2 else (s[:w1], s[-w2:])].append(rid)
        table = f"resolution_illumina_{read_len}x{read_len}.tsv"
        rank_of = {}
        with open(os.path.join(out, table), "w") as fh:
            fh.write("id\twindow_group\trecords_in_group\tdeepest_rank\tname\tmargin\n")
            mcol = f"margin_{w1}x{w2}"
            for key, ids in groups.items():
                gid = "W_" + hashlib.sha1(repr(key).encode()).hexdigest()[:12]
                r, n = deepest(merge_records([tax[i] for i in ids]))
                for rid in ids:
                    rank_of[rid] = r
                    m = margins.get(rid, {}).get(mcol, "") if r == "species" else ""
                    fh.write(f"{rid}\t{gid}\t{len(ids)}\t{r}\t{n}\t{m}\n")
        # Through the windows a record can never be resolved deeper than
        # read whole; if it is, the merge above is wrong. Refuse the release.
        worse = [i for i in seqs
                 if order.get(rank_of[i], -1) > order.get(deepest(tax[i])[0], -1)]
        if worse:
            sys.exit(f"{len(worse)} records resolve deeper through the {table} "
                     f"windows than whole, e.g. {worse[:3]}: not releasing")
        geometries.append({"read_length": read_len, "windows": [w1, w2], "file": table,
                           "window_groups": len(groups), "rank_of": rank_of})
    primary = geometries[0]
    ill_rank, groups = primary["rank_of"], range(primary["window_groups"])
    a.read_length, (w1, w2) = primary["read_length"], primary["windows"]

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
    with open(os.path.join(out, "primer_fit.tsv")) as fh:
        fh.readline()
        fit = Counter(line.split("\t")[1] for line in fh)
        fit_helm = Counter()
    # By domain too: prokaryotes amplify only at the 2-edit limit, 3' end
    # mismatched, so their share of each fit class shows the release's
    # bacterial exposure (v0.4.0 feedback, §3).
    fit_dom = defaultdict(Counter)
    with open(os.path.join(out, "primer_fit.tsv")) as fh:
        fh.readline()
        for line in fh:
            f = line.split("\t")
            if tax[f[0]][4] in HELMINTHS:
                fit_helm[f[1]] += 1
            fit_dom[tax[f[0]][0].split(":")[0] or "(none)"][f[1]] += 1

    # Margins: how many bases species calls rest on, read whole and
    # through the summary's windows.
    def margin_bins(col, ids):
        c = Counter()
        for i in ids:
            v = margins.get(i, {}).get(col, "")
            c["none near" if v == "" else "1" if int(v) <= 1 else "2" if int(v) == 2
              else "3-5" if int(v) <= 5 else ">5"] += 1
        return c
    sp_whole = [i for i in seqs if deepest(tax[i])[0] == "species"]
    sp_ill = [i for i in seqs if ill_rank[i] == "species"]
    mb_whole = margin_bins("margin_whole", sp_whole)
    mb_ill = margin_bins(f"margin_{w1}x{w2}", sp_ill)

    # Clades V4 cannot split: by the class column (see rank_names.tsv for
    # what it holds), the share of records that stop above genus.
    by_clade = defaultdict(list)
    for i in seqs:
        by_clade[tax[i][4].split(":")[0] or "(none)"].append(i)
    clade_rows = []
    for c, ids in by_clade.items():
        if len(ids) < 300:
            continue
        above = lambda r: r not in ("species", "genus")  # noqa: E731
        clade_rows.append((c, len(ids),
                           sum(deepest(tax[i])[0] == "species" for i in ids) / len(ids),
                           sum(above(deepest(tax[i])[0]) for i in ids) / len(ids),
                           sum(above(ill_rank[i]) for i in ids) / len(ids)))
    clade_rows.sort(key=lambda r: -r[3])

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
        "illumina_geometries": [{k: g[k] for k in ("read_length", "windows", "file",
                                                    "window_groups")} for g in geometries],
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
            "primer_fit": dict(fit),
            "primer_fit_helminth": dict(fit_helm),
            "primer_fit_by_domain": {d: dict(c) for d, c in sorted(fit_dom.items())},
            "species_margins": {"read_whole": dict(mb_whole), f"illumina_{w1}x{w2}": dict(mb_ill)}
            if margins else None,
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
             "", "Primer fit (best-fitting member of each record):",
             *[f"    {k:12s} {fit.get(k, 0):>8,}  {100 * fit.get(k, 0) / len(seqs):5.1f}%   "
               f"helminth {fit_helm.get(k, 0):>6,}"
               for k in ("exact", "away_3prime", "near_3prime", "indel", "edits_only", "unknown")],
             "  by domain (records: exact / away_3prime / near_3prime / other):",
             *[f"    {d:12s} {sum(c.values()):>8,}: " + " / ".join(
                 f"{c.get(k, 0):,}" for k in ("exact", "away_3prime", "near_3prime"))
               + f" / {sum(v for k, v in c.items() if k not in ('exact', 'away_3prime', 'near_3prime')):,}"
               for d, c in sorted(fit_dom.items(), key=lambda x: -sum(x[1].values()))],
             "", "Species margins - edits from a species-level record to the nearest record of",
             "another species (species_margins.tsv; the 'margin' column of each resolution",
             "table). A species call resting on 1-2 edits would also be made for a lineage",
             "the release lacks, a base or two away. 'none near': no record of another",
             "species shares even half its insert - a large margin:",
             *[f"    {label:28s} " + "  ".join(f"{k} {c.get(k, 0):>6,}" for k in
                                               ("1", "2", "3-5", ">5", "none near"))
               for label, c in ((f"read whole ({len(sp_whole):,})", mb_whole),
                                (f"Illumina {w1}x{w2} ({len(sp_ill):,})", mb_ill))],
             "", "Clades V4 splits least (class column, >= 300 records; what the column",
             "holds by clade: rank_names.tsv). A call stopping in one of these means",
             "'genus unknowable from V4', not 'badly classified':",
             f"    {'clade':22s} {'records':>8s}  {'species':>7s}  {'stop above genus':>16s}  {'(Illumina)':>10s}",
             *[f"    {c:22s} {n:>8,}  {100 * sp:6.0f}%  {100 * ab:15.0f}%  {100 * abi:9.0f}%"
               for c, n, sp, ab, abi in clade_rows[:15]],
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
