"""What every build step shares: the primers, the in-silico PCR, and the
PR2 rank conventions."""

from __future__ import annotations

import gzip
import re

from seqmatch import IUPAC, find_anchor_cost

FWD = "AGCAGCCGCGGTAATTCC"                # AMPSEQ_VESPA18S_29F
REV = "TCCGTCAATTYCTTNAASTTTC"            # AMPSEQ_VESPA18S_21b8_I_R_Nmod
_COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
REV_RC = REV.translate(_COMP)[::-1]
PRIMER_EDITS = 2

RANKS = ["domain", "supergroup", "division", "subdivision", "class",
         "order", "family", "genus", "species"]
# Identical V4 across different PR2 classes (Nematoda, Platyhelminthes,
# Craniata, ...) is not biology; it is a mislabel or a contaminant.
CONFLICT_DEPTH = RANKS.index("class")


def _regex(p: str):
    return re.compile("".join(c if c in "ACGT" else f"[{IUPAC[c]}]" for c in p))


_RF, _RR = _regex(FWD), _regex(REV_RC)


def pcr(seq: str):
    """In-silico PCR on one sense-strand sequence.

    Returns (status, insert, f_edits, r_edits); status is "amplicon",
    "F_only", "R_only" or "neither". Exact match first, then up to
    PRIMER_EDITS edits.
    """
    seq = seq.upper()
    m = _RF.search(seq)
    if m:
        fe, fc = m.end(), 0
    else:
        h = find_anchor_cost(seq, FWD, PRIMER_EDITS)
        fe, fc = (h[1], h[2]) if h else (None, None)
    if fe is None:
        r = _RR.search(seq) or find_anchor_cost(seq, REV_RC, PRIMER_EDITS)
        return ("R_only" if r else "neither"), None, None, None
    rest = seq[fe:]
    m = _RR.search(rest)
    if m:
        rs, rc = m.start(), 0
    else:
        h = find_anchor_cost(rest, REV_RC, PRIMER_EDITS)
        rs, rc = (h[0], h[2]) if h else (None, None)
    if rs is None:
        return "F_only", None, fc, None
    return "amplicon", rest[:rs], fc, rc


def informative(name: str) -> bool:
    """A PR2 placeholder (Metazoa_X, Enoplea_XX, Genus_sp.) names nothing."""
    # Organelle entries carry a suffix (Diatomeae_X:plas); judge the name.
    base = name.split(":")[0] if name else ""
    return bool(base) and not re.search(r"(_X+(_|$)|_sp\.?$)", base)


def read_taxo_long(path: str):
    """PR2 taxo_long FASTA: yields (accession, gene, organelle, ranks, seq).

    Header: accession|gene|organelle|strain|domain|...|species
    """
    head, buf = None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if head is not None:
                    yield (*head, "".join(buf))
                f = line[1:].rstrip("\n").split("|")
                head, buf = (f[0], f[1], f[2], f[4:13]), []
            else:
                buf.append(line.strip())
    if head is not None:
        yield (*head, "".join(buf))


def open_text(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode[0])


class PR2Index:
    """What a non-PR2 layer needs from the PR2 layer: PR2's ranks for a
    genus or for any higher name, and the PR2 class of the nearest PR2
    insert. Built once from work/pr2_layer.tsv.gz.

    Mapping a taxon onto PR2's nine ranks: if PR2 knows the genus, its
    ranks are PR2's for that genus. Otherwise the NCBI lineage is walked
    upward to the first name PR2 has at some rank down to family; PR2's
    ranks are used to there and the rest filled with PR2-style
    placeholders (Rhabditida_X, Rhabditida_XX). Both layers map this way,
    so a species gets the same ranks whichever layer brings it.
    """

    K, STRIDE = 16, 4

    def __init__(self, pr2_layer: str):
        from collections import Counter, defaultdict
        self.by_genus = defaultdict(Counter)
        self.by_name = {}
        self.ins_class = {}
        self.accessions = set()
        with open_text(pr2_layer) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            col = {c: i for i, c in enumerate(head)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                ranks = tuple(f[col[r]] for r in RANKS)
                # PR2 accessions are AB353770.1.1740_U: GenBank id first.
                self.accessions.add(f[col["accession"]].split(".")[0])
                if informative(ranks[7]):
                    self.by_genus[ranks[7]][ranks[:8]] += 1
                for i, n in enumerate(ranks[:8]):
                    if informative(n) and n not in self.by_name:
                        self.by_name[n] = (i, ranks[:i + 1])
                if f[col["status"]] == "amplicon":
                    self.ins_class.setdefault(f[col["insert"]], ranks[4])
        self._uniq = None

    def ranks_for(self, genus: str, species: str, lineage) -> tuple:
        """(nine ranks, how). lineage: NCBI names root-first, or a callable
        returning them."""
        names = (lineage() if callable(lineage) else lineage) or []
        if self.by_genus.get(genus):
            # Genus names are not unique across kingdoms: Setaria is a
            # nematode and a grass, Labrys a nematode and a bacterium.
            # PR2's ranks for the name count only if their higher names
            # agree with the NCBI lineage; without a lineage, trust PR2.
            # The most specific of them must agree - Metazoa is shared by
            # every animal, so it cannot tell a nematode from a vertebrate.
            seen = set(names)
            for cand, _ in self.by_genus[genus].most_common():
                higher = [n.split(":")[0] for n in cand[1:5] if informative(n)]
                if not seen or not higher or higher[-1] in seen:
                    return (*cand, species), "PR2 genus"
        for name in reversed(names):
            if name in self.by_name and self.by_name[name][0] <= 6:
                depth, prefix = self.by_name[name]
                fill = [f"{prefix[-1]}_{'X' * k}" for k in range(1, 7 - depth)]
                return (*prefix, *fill, genus, species), f"NCBI via {name}"
        return ("",) * 7 + (genus, species), "unplaced"

    def nearest_class(self, q: str):
        """(PR2 class, share) of the PR2 insert sharing most sparse k-mers."""
        from collections import Counter, defaultdict
        if q in self.ins_class:
            return self.ins_class[q], 1.0
        if self._uniq is None:
            self._uniq = list(self.ins_class)
            self._index = defaultdict(list)
            for j, s in enumerate(self._uniq):
                for i in range(0, len(s) - self.K + 1, self.STRIDE):
                    self._index[s[i:i + self.K]].append(j)
        hits = Counter()
        for i in range(len(q) - self.K + 1):
            for j in self._index.get(q[i:i + self.K], ()):
                hits[j] += 1
        if not hits:
            return "", 0.0
        j, n = hits.most_common(1)[0]
        s = self._uniq[j]
        return self.ins_class[s], n / max(1, (len(s) - self.K) // self.STRIDE + 1)


LAYER_COLUMNS = ["source", "accession", "gene", "organelle", *RANKS,
                 "status", "f_edits", "r_edits", "insert"]
GENOME_REVIEW_COLUMNS = ["species", "copy", "copies", "expected_class",
                         "nearest_pr2_class", "kmer_share", "ranks_from", "insert"]


def genome_layer_rows(copies: dict, where: dict, ranks_for, pr2: "PR2Index",
                      source: str):
    """The genome-layer rules (docs/DESIGN.md), for any set of assemblies.

    copies:    species -> Counter(insert -> copies), pooled across assemblies
    where:     (species, insert) -> one copy's location, as its accession
    ranks_for: species -> (nine ranks, how they were found)

    SUPPORT: a variant enters if seen in at least two copies; a species
    with no such variant contributes its most common one. PHYLUM: it
    enters only if its nearest PR2 insert is in the species' own PR2
    class; otherwise it goes to review as another organism's rDNA.
    Returns (rows, review, stats).
    """
    from collections import Counter
    rows, review, stats = [], [], Counter()
    for species in sorted(copies):
        ranks, how = ranks_for(species)
        stats[how.split()[0]] += 1
        variants = copies[species]
        supported = [v for v, n in variants.items() if n >= 2]
        single = not supported
        if single:
            supported = [variants.most_common(1)[0][0]]
        stats["variants_seen"] += len(variants)
        for v in variants:
            if v not in supported:
                stats["dropped_unsupported"] += 1
                continue
            cls, score = pr2.nearest_class(v)
            if ranks[4] and informative(ranks[4]) and cls and cls != ranks[4]:
                stats["review_other_class"] += 1
                review.append([species, where[(species, v)], str(variants[v]),
                               ranks[4], cls, f"{score:.2f}", how, v])
                continue
            stats["entered"] += 1
            stats["single_copy"] += single
            rows.append([source, where[(species, v)], "18S_rRNA", "nucleus", *ranks,
                         "amplicon", "", "", v])
    return rows, review, stats


def write_table(path: str, columns, rows) -> None:
    with (gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")) as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
