#!/usr/bin/env python3
"""Collapse layers into unique V4 inserts, merge their taxonomy, and flag
label conflicts.

    bin/collapse.py <outdir> <layer.tsv.gz> [<layer.tsv.gz> ...]

A layer is a TSV with at least: source, accession, the nine PR2 ranks,
status, insert. Rows whose status is not "amplicon" are ignored. Layers
are given in increasing order of trust - PR2 first, local sequences last -
and that order decides conflicts (below).

Writes refdb.fasta, taxonomy.tsv, members.tsv, conflicts.tsv and
primer_fit.tsv (see best_fit) to
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
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vespa import (CONFLICT_DEPTH, FWD, RANKS, REV, THREE_PRIME_WINDOW,  # noqa: E402
                   informative, normalise_ranks, open_text, primer_fit)


MIN_INSERT = 100

# A bacterial or archaeal 16S amplifies with these primers to about 370bp;
# a eukaryote's 18S to about 550. A prokaryote label on an insert this long
# is eukaryotic 18S inside a bacterial genome assembly - host plant, insect
# or human DNA in the sequenced culture.
EUKARYOTE_LENGTH = 450

# Organelle suffixes PR2 puts on labels. Not "nucl": a nucleomorph's 18S
# is an 18S.
ORGANELLES = ("mito", "plas", "apic", "chro")


# GenBank WGS and TSA accessions: four (or six) letters, then digits.
# Contigs from genome or transcriptome assemblies, which carry whatever was
# in the sample, so a directly sequenced, named isolate outranks them.
ASSEMBLY_ACCESSION = re.compile(r"^[A-Z]{4}(?:[A-Z]{2})?\d{8,}")


def load_overrides(path: str):
    """curation/overrides.tsv: judgements no rule can make, with evidence.

    action "drop": first column is a source accession (as members.tsv
    writes it, stable across builds where record ids follow the sequence);
    that member is dropped. action "relabel": the same, but the member is
    kept with the label in the fourth column (the nine ranks, ";"-joined,
    PR2 style) - for a sequence that is right but filed under the wrong
    name, such as a genus homonym. action "host_species": first column is
    a species label that names the host, not the organism; it becomes
    Genus_sp. Returns ({accession: reason}, {species label},
    {accession: (ranks, reason)}).
    """
    drop, host, relabel = {}, set(), {}
    if path and os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if f[0] == "accession" or len(f) < 3:
                    continue
                if f[1] == "drop":
                    drop[f[0]] = f[2]
                elif f[1] == "host_species":
                    host.add(f[0])
                elif f[1] == "relabel":
                    ranks = f[3].split(";") if len(f) > 3 else []
                    if len(ranks) != len(RANKS):
                        sys.exit(f"{path}: relabel of {f[0]} needs {len(RANKS)} ranks")
                    relabel[f[0]] = (ranks, f[2])
    return drop, host, relabel


def load_subtypes(path: str) -> dict:
    """curation/blastocystis_v4_subtypes.tsv, from
    bin/blastocystis_subtypes.py: {record id: subtype}."""
    out = {}
    if path and os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if line.startswith("#") or f[0] == "id" or len(f) < 2:
                    continue
                out[f[0]] = f[1]
    return out


def apply_subtype(members, subtype: str) -> int:
    """Blastocystis subtypes are not species, and PR2 and GenBank name few
    of them: "Blastocystis sp." and the host-based names (B. hominis,
    B. ratti) say nothing a subtype does not. Every Blastocystis member of
    a record with a curated subtype takes Blastocystis_STn as its species;
    members of other genera keep their labels, so the conflict rules still
    judge a record whose members disagree. Returns how many changed."""
    n = 0
    for m in members:
        r = m["ranks"]
        if r[7].split(":")[0].split("-")[0] == "Blastocystis":
            m["ranks"] = r[:8] + [f"Blastocystis_{subtype}"]
            n += 1
    return n


def screen(members, insert: str, overrides: dict):
    """(kept, dropped) before conflicts are judged: members whose label
    cannot belong to this sequence, each dropped with its reason."""
    kept, dropped = [], []
    nuclear = any(":" not in m["ranks"][0] for m in members)
    for m in members:
        dom = m["ranks"][0]
        if m["accession"] in overrides:
            dropped.append((m, "override: " + overrides[m["accession"]]))
        elif dom in ("Bacteria", "Archaea") and len(insert) >= EUKARYOTE_LENGTH:
            dropped.append((m, f"prokaryote label on a {len(insert)}bp insert: "
                               "eukaryotic 18S in a bacterial assembly"))
        elif dom.split(":")[-1] in ORGANELLES and nuclear:
            dropped.append((m, "organelle label beside nuclear ones on an 18S"))
        else:
            kept.append(m)
    return kept, dropped


def fix_host_species(members, genus_classes: dict, host_named: set = frozenset()) -> int:
    """A species name that is the host's, not the organism's (a gregarine
    labelled with its weevil, Leidyana / Trachyphloeus): the species
    becomes Genus_sp. Caught two ways: the name's genus is used in another
    class here, or the species is listed as host-named in the overrides
    (checked against NCBI taxonomy - a parasite lineage is full of species
    PR2 files under a clade's genus, which a name test cannot tell from a
    host). Returns how many were fixed."""
    n = 0
    for m in members:
        r = m["ranks"]
        genus, species, cls = r[7], r[8], r[4]
        if not informative(genus) or not informative(species):
            continue
        word = species.split(":")[0].split("_")[0]
        # The genus without PR2's organelle suffix or homonym tag
        # (Nemastoma_arachnid, Serpula-fungus): the bare name species use.
        base = genus.split(":")[0].split("-")[0].split("_")[0]
        if word == base:
            continue
        other = genus_classes.get(word, set())
        foreign = bool(other) and cls not in other or species in host_named
        if foreign:
            m["ranks"] = r[:8] + [genus.split(":")[0] + "_sp." + (":" + genus.split(":")[1] if ":" in genus else "")]
            n += 1
    return n


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

    # A directly sequenced isolate outranks an assembly contig: if the
    # isolates agree among themselves, conflicting contigs lose.
    direct = [m for m in members if not ASSEMBLY_ACCESSION.match(m["accession"])]
    direct_labels = {high_label(m["ranks"]) for m in direct}
    if direct and len(direct) < len(members) and \
            all(compatible(a, b) for a in direct_labels for b in direct_labels):
        win = next(iter(direct_labels))
        over = [m for m in members if not compatible(high_label(m["ranks"]), win)]
        if over and all(ASSEMBLY_ACCESSION.match(m["accession"]) for m in over):
            return [m for m in members if m not in over], over, "isolate"

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


# Sources whose members are an organism's own sequenced rDNA (genome
# assemblies, raw reads, our own sequences), not a deposited label.
GENOME_SOURCES = {"WBPS19", "ParaRef", "NCBI_genomes", "ENA_reads", "ToL_prerelease", "local"}


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
        # A genome's or a read set's own rDNA is the organism, not a label
        # someone attached: never outvoted here (it still faces the class
        # check and screen()). Found when cyathostomin genera sharing one V4
        # had Cyathostomum's and Coronocyclus's reads overruled.
        if ms[0]["source"] in GENOME_SOURCES:
            continue
        agree = sum(1 for m in kept if m["ranks"][i] == ref[i])
        if not informative(ms[0]["ranks"][-1]) or agree >= 3:
            over.append(ms[0])
    return [m for m in kept if m not in over], over


FIT_COLUMNS = ["id", "fit", "f_edits", "r_edits", "f_mismatch_from_3prime",
               "r_mismatch_from_3prime", "f_site", "r_site", "best_member"]


def best_fit(kept) -> list:
    """The primer fit of a record: its best-fitting member's.

    A record is one insert; its members may carry different primer sites
    (the primers sit outside the insert). The best is the one PCR is most
    likely to amplify: known sites over unknown, fewest edits, and the
    farthest any mismatch sits from a primer's 3' end. Classes:
      exact        both sites match their primer
      away_3prime  mismatches, none in the last THREE_PRIME_WINDOW bases
      near_3prime  a mismatch within the last THREE_PRIME_WINDOW bases,
                   where it most impairs extension
      indel        a site is not the primer's length
      edits_only   edit counts known, sites not (genome scans kept counts)
      unknown      nothing known (a local sequence given as a bare insert)
    """
    def score(m):
        fe, re_, fs, rs = m["fit"]
        fp, rp = primer_fit(fs, FWD), primer_fit(rs, REV)
        if fe == "" or re_ == "":
            return (3, 99, 0), "unknown", m, fp, rp
        edits = int(fe) + int(re_)
        if edits == 0:
            return (0, 0, 0), "exact", m, fp or [], rp or []
        if fp is None or rp is None:
            return (2, edits, 0), ("edits_only" if not fs or not rs else "indel"), m, fp, rp
        near = min(fp + rp)
        cls = "near_3prime" if near <= THREE_PRIME_WINDOW else "away_3prime"
        return (1 if cls == "away_3prime" else 1.5, edits, -near), cls, m, fp, rp
    _, cls, m, fp, rp = min((score(m) for m in kept), key=lambda t: t[0])
    fe, re_, fs, rs = m["fit"]
    pos = lambda p: "" if p is None else ",".join(map(str, sorted(p)))  # noqa: E731
    return [cls, fe, re_, pos(fp), pos(rp), fs, rs, f"{m['source']}:{m['accession']}"]


def main() -> int:
    args = sys.argv[1:]
    overrides_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  os.pardir, "curation", "overrides.tsv")
    if "--overrides" in args:
        k = args.index("--overrides")
        overrides_path = args[k + 1]
        del args[k:k + 2]
    subtypes_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 os.pardir, "curation", "blastocystis_v4_subtypes.tsv")
    if "--subtypes" in args:
        k = args.index("--subtypes")
        subtypes_path = args[k + 1]
        del args[k:k + 2]
    outdir, layers = args[0], args[1:]
    os.makedirs(outdir, exist_ok=True)
    overrides, host_named, relabels = load_overrides(overrides_path)
    relabelled = 0

    groups = defaultdict(list)
    trust, order = {}, []
    n_rows = Counter()
    n_short = Counter()
    for li, path in enumerate(layers):
        with open_text(path) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            col = {c: i for i, c in enumerate(head)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[col["status"]] != "amplicon":
                    continue
                if len(f[col["insert"]]) < MIN_INSERT:
                    # Primer matches that touch or overlap "amplify" nothing;
                    # the genome scans already apply the same floor.
                    n_short[f[col["source"]]] += 1
                    continue
                src = f[col["source"]]
                if src not in trust:
                    # A layer earlier on the command line is trusted less.
                    trust[src] = li
                    order.append(src)
                n_rows[src] += 1
                label = [f[col[r]] for r in RANKS]
                get = lambda c: f[col[c]] if c in col else ""  # noqa: E731
                groups[f[col["insert"]]].append({
                    "source": src, "accession": f[col["accession"]],
                    "ranks": normalise_ranks(label), "label": label,
                    "fit": (get("f_edits"), get("r_edits"), get("f_site"), get("r_site"))})
                if f[col["accession"]] in relabels:
                    ranks, why = relabels[f[col["accession"]]]
                    m = groups[f[col["insert"]]][-1]
                    m["ranks"], m["note"] = normalise_ranks(ranks), "relabelled: " + why
                    relabelled += 1

    # Which classes each genus name is used in, across every source - so a
    # species name borrowed from another class (a host) can be recognised.
    genus_classes = defaultdict(set)
    for members in groups.values():
        for m in members:
            g = m["ranks"][7].split(":")[0].split("-")[0].split("_")[0]
            if informative(m["ranks"][7]):
                genus_classes[g].add(m["ranks"][4])
    host_fixed = sum(fix_host_species(ms, genus_classes, host_named) for ms in groups.values())
    subtypes = load_subtypes(subtypes_path)
    subtyped = sum(apply_subtype(ms, subtypes[record_id(ins)])
                   for ins, ms in groups.items() if record_id(ins) in subtypes)
    # "B. hominis" is a host-based name, not an identification: with no
    # subtype to say more, it says only Blastocystis.
    for ins, ms in groups.items():
        if record_id(ins) not in subtypes:
            for m in ms:
                if m["ranks"][8] == "Blastocystis_hominis":
                    m["ranks"] = m["ranks"][:8] + ["Blastocystis_sp."]

    decisions = Counter()
    fits = Counter()
    drops = Counter()
    excluded = 0
    with open(os.path.join(outdir, "refdb.fasta"), "w") as fa, \
         open(os.path.join(outdir, "taxonomy.tsv"), "w") as tx, \
         open(os.path.join(outdir, "members.tsv"), "w") as mb, \
         open(os.path.join(outdir, "conflicts.tsv"), "w") as cf, \
         open(os.path.join(outdir, "primer_fit.tsv"), "w") as pf:
        pf.write("\t".join(FIT_COLUMNS) + "\n")
        tx.write("\t".join(["id", *RANKS, "deepest_rank", "members", "sources"]) + "\n")
        mb.write("id\tsource\taccession\tlabel\trole\tnote\n")
        cf.write("id\tdecision\tmembers\tlabels\tkept_label\n")
        for insert in sorted(groups, key=record_id):
            members = groups[insert]
            rid = record_id(insert)
            screened, dropped = screen(members, insert, overrides)
            for m, why in dropped:
                drops[why.split(":")[0].split(" on a ")[0]] += 1
                mb.write(f"{rid}\t{m['source']}\t{m['accession']}\t"
                         f"{';'.join(m['label'])}\tdropped\t{why}\n")
            if not screened:
                # Nothing left to say what this sequence is: it is not a
                # record. (Human 18S from a bacterial assembly, for one.)
                excluded += 1
                continue
            members = screened
            kept, over, how = resolve(members, trust)
            kept, over2 = trust_below_class(kept, trust)
            if over2:
                over += over2
                how = "trust" if how == "none" else how
            decisions[how] += 1
            for m in members:
                role = "kept" if m in kept else "overruled"
                mb.write(f"{rid}\t{m['source']}\t{m['accession']}\t"
                         f"{';'.join(m['label'])}\t{role}\t{m.get('note', '')}\n")
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
            fit = best_fit(kept)
            fits[fit[0]] += 1
            pf.write("\t".join([rid, *fit]) + "\n")

    print(f"layers: " + ", ".join(f"{s} ({n_rows[s]:,} amplicons)" for s in order))
    if n_short:
        print("  set aside as shorter than {}bp: ".format(MIN_INSERT)
              + ", ".join(f"{s} {n}" for s, n in n_short.items()))
    print(f"unique inserts: {len(groups):,}; excluded (no member left): {excluded}")
    print("members dropped: " + (", ".join(f"{k} {v}" for k, v in drops.items()) or "none")
          + f"; host-named species fixed: {host_fixed}"
          + f"; Blastocystis members given a subtype: {subtyped}"
          + f"; members relabelled: {relabelled}")
    for k in ("none", "isolate", "trust", "majority", "unresolved"):
        print(f"  conflict {k:9s} {decisions[k]:>7,}")
    print("primer fit: " + ", ".join(f"{k} {fits[k]:,}" for k in
          ("exact", "away_3prime", "near_3prime", "indel", "edits_only", "unknown")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
