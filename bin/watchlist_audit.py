#!/usr/bin/env python3
"""Where do the helminths we care about stand in a release?

    bin/watchlist_audit.py <release dir> [curation/helminth_watchlist.tsv] > audit.txt

Per species on the watch list: records naming it to species, records it is
a member of (a shared V4 sits in a higher-rank record), the deepest rank it
reaches read whole and through each Illumina table, and the largest margin
among its species-level records (edits to the nearest other species;
'far' = none near). MISSING: no V4-spanning sequence in any layer.
"""
from __future__ import annotations
import collections, csv, glob, os, sys

ORDER = ["species", "genus", "family", "order", "class", "subdivision", "division",
         "supergroup", "domain", ""]


def main() -> int:
    rel = sys.argv[1]
    wl = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), os.pardir, "curation", "helminth_watchlist.tsv")
    tax = {t["id"]: t for t in csv.DictReader(open(os.path.join(rel, "taxonomy.tsv")), delimiter="\t")}
    tables = {"whole": "resolution_nanopore.tsv"}
    for f in sorted(glob.glob(os.path.join(rel, "resolution_illumina_*.tsv")), reverse=True):
        tables[os.path.basename(f)[len("resolution_illumina_"):-4]] = os.path.basename(f)
    res = {k: {r["id"]: r for r in csv.DictReader(open(os.path.join(rel, f)), delimiter="\t")}
           for k, f in tables.items()}
    members = collections.defaultdict(set)
    for m in csv.DictReader(open(os.path.join(rel, "members.tsv")), delimiter="\t"):
        if m["role"] == "kept":
            members[m["label"].split(";")[8].split(":")[0]].add(m["id"])
    for t in tax.values():
        if t["species"]:
            members[t["species"]]
    group = None
    n = collections.Counter()
    for r in csv.DictReader((l for l in open(wl) if not l.startswith("#")), delimiter="\t"):
        if r["group"] != group:
            group = r["group"]
            print(f"\n## {group}")
        sp = r["species"]
        own = [i for i, t in tax.items() if t["species"] == sp and t["deepest_rank"] == "species"]
        ids = members.get(sp, set()) | set(own)
        if not ids:
            g = sum(1 for t in tax.values() if t["genus"] == sp.split("_")[0])
            print(f"  {sp:32s} MISSING (genus records: {g})")
            n["missing"] += 1
            continue
        cells = []
        for k, tb in res.items():
            rank = ORDER[min(ORDER.index(tb[i]["deepest_rank"]) for i in ids if i in tb)]
            ms = [tb[i].get("margin", "") for i in own if tb.get(i, {}).get("deepest_rank") == "species"]
            m = max((int(x) for x in ms if x), default=None) if any(ms) else ("far" if ms else "-")
            cells.append(f"{k} {rank:7s} m={m}")
        top = ORDER[min(ORDER.index(res['whole'][i]['deepest_rank']) for i in ids if i in res['whole'])]
        n["species" if top == "species" else "shared"] += 1
        print(f"  {sp:32s} own {len(own):3d} in {len(ids):3d}  " + "  ".join(cells))
    print(f"\n{sum(n.values())} species: {n['species']} reach species read whole, "
          f"{n['shared']} only in shared records, {n['missing']} missing", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
