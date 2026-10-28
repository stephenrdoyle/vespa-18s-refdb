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
