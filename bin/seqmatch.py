"""Sequence matching shared with nf-helminth-amplicon.

Copied verbatim from nf-helminth-amplicon bin/ampseqlib.py at commit
29362fa. The anchor look-alike scan in particular MUST stay identical to the
pipeline's: the pipeline reads the mask this repository writes, and a
mask built by a different matcher protects different sites. If either
copy changes, change both.
"""

from __future__ import annotations

# The inner Nextera anchors the pipeline's splitter and demultiplexer
# search for (assets/assay/ampseq_helminth.json in the pipeline).
ANCHORS = {
    "i5_inner_anchor": "TCGTCGGCAGCGTC",
    "i7_inner_anchor": "CCGAGCCCACGAGAC",
}


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


# IUPAC ambiguity, as the set of concrete bases each code stands for. Real
# primers are routinely degenerate, and a matcher that treats Y as a literal
# character fails to find the primer in every sequence that carries either
# base it was designed to cover.
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))


# --------------------------------------------------------------------------
# approximate matching
# --------------------------------------------------------------------------


def find_anchor_cost(
    text: str, pattern: str, max_edits: int
) -> tuple[int, int, int] | None:
    """As find_anchor, but also returns the edit cost.

    Primer trimming needs the cost: several primers compete for the same
    read start, and the one that fits best is the one that was used.
    """
    m = len(pattern)
    if not m or m > 64 or not text:
        return None

    # Bits per concrete base, expanding any degenerate code in the pattern,
    # then one mask per IUPAC code the text might use - so degeneracy on
    # either side matches, and an N in either matches anything.
    base: dict[str, int] = {"A": 0, "C": 0, "G": 0, "T": 0}
    for i, c in enumerate(pattern):
        for b in IUPAC.get(c, ""):
            base[b] |= 1 << i
    peq = {
        code: (base["A"] if "A" in bases else 0)
        | (base["C"] if "C" in bases else 0)
        | (base["G"] if "G" in bases else 0)
        | (base["T"] if "T" in bases else 0)
        for code, bases in IUPAC.items()
    }

    full = (1 << m) - 1
    top = 1 << (m - 1)
    pv, mv, score = full, 0, m
    best_end, best_cost = -1, max_edits + 1

    for j, c in enumerate(text, 1):
        eq = peq.get(c, 0)
        xv = eq | mv
        xh = ((((eq & pv) + pv) ^ pv) | eq) & full
        ph = (mv | ~(xh | pv)) & full
        mh = pv & xh
        if ph & top:
            score += 1
        elif mh & top:
            score -= 1
        ph = (ph << 1) & full
        mh = (mh << 1) & full
        pv = (mh | ~(xv | ph)) & full
        mv = ph & xv
        if score < best_cost:
            best_cost, best_end = score, j

    if best_end < 0 or best_cost > max_edits:
        return None
    return max(0, best_end - m), best_end, best_cost


def find_anchor_leftmost(
    text: str, pattern: str, max_edits: int
) -> tuple[int, int, int] | None:
    """The FIRST occurrence of `pattern` within tolerance, not the best.

    find_anchor_cost answers "where does this pattern fit best in this
    window", which is the right question when the window is known to
    hold one copy - a primer at a read start, an anchor at a read end.
    It is the wrong question for a caller that wants every copy, left to
    right. Such a caller searches the remainder after each hit, and if
    the remainder holds a one-edit copy followed by an exact one, the
    best hit is the exact one: the search jumps straight past the
    earlier copy and never comes back for it. In the splitter that loses
    exactly the anchor that marks a junction, and the read goes on to be
    demultiplexed on one molecule's i5 and another's i7.

    So: the first end position whose cost is within tolerance. Not that
    position alone, though - an exact copy ending at j also matches
    with k deletions ending at j - k, so the first position within
    tolerance is usually a worse-scoring early view of the same copy.
    The scan carries on for max_edits further columns and takes the
    lowest cost in that window (the earliest, on a tie), which is the
    end of the copy it found rather than a ragged prefix of it. Beyond
    that window any further hit is a different copy, and is left for
    the caller's next search.

    Returns (start, end, cost) like find_anchor_cost, with start taken
    the same way - end minus pattern length - so the two agree on any
    window holding a single copy.
    """
    m = len(pattern)
    if not m or m > 64 or not text:
        return None

    base: dict[str, int] = {"A": 0, "C": 0, "G": 0, "T": 0}
    for i, c in enumerate(pattern):
        for b in IUPAC.get(c, ""):
            base[b] |= 1 << i
    peq = {
        code: (base["A"] if "A" in bases else 0)
        | (base["C"] if "C" in bases else 0)
        | (base["G"] if "G" in bases else 0)
        | (base["T"] if "T" in bases else 0)
        for code, bases in IUPAC.items()
    }

    full = (1 << m) - 1
    top = 1 << (m - 1)
    pv, mv, score = full, 0, m
    best_end, best_cost, stop_at = -1, max_edits + 1, -1

    for j, c in enumerate(text, 1):
        eq = peq.get(c, 0)
        xv = eq | mv
        xh = ((((eq & pv) + pv) ^ pv) | eq) & full
        ph = (mv | ~(xh | pv)) & full
        mh = pv & xh
        if ph & top:
            score += 1
        elif mh & top:
            score -= 1
        ph = (ph << 1) & full
        mh = (mh << 1) & full
        pv = (mh | ~(xv | ph)) & full
        mv = ph & xv
        if score < best_cost:
            best_cost, best_end = score, j
            if stop_at < 0:
                stop_at = j + max_edits
        if 0 <= stop_at <= j:
            break

    if best_end < 0 or best_cost > max_edits:
        return None
    return max(0, best_end - m), best_end, best_cost





# --------------------------------------------------------------------------
# anchor look-alikes inside the panel's own amplicons
# --------------------------------------------------------------------------

# Sequence kept either side of a look-alike to recognise it again in a
# read. 24 + a 14-15bp anchor stays inside the 64bp the bit-parallel
# matcher handles, and is long enough that no adapter junction resembles
# it: a junction's flanks are index and Nextera overhang, not insert.
LOOKALIKE_FLANK = 24

# Found at one edit looser than the splitter and demultiplexer search at,
# so a look-alike that a sequencing error carries inside their tolerance
# is still in the mask. Masking a site that never triggers costs nothing:
# it is insert either way.
LOOKALIKE_MAX_EDITS = 3


def find_anchor_lookalikes(entries, anchors: dict, max_edits: int,
                           flank: int = LOOKALIKE_FLANK) -> list[dict]:
    """Every place an adapter anchor is nearly spelled out inside an amplicon.

    The splitter cuts a read wherever it finds an anchor, and the chimera
    check quarantines a read with an anchor in its interior - both on the
    assumption that an amplicon's insert contains no adapter. A panel
    entry whose insert happens to contain a near-copy of an anchor breaks
    that assumption for every one of its reads: each is cut in two, or set
    aside as a chimera, and the target is quietly under-counted. Found on
    a real panel, where a Trichuris incognita amplicon carried the i5
    anchor at two edits and lost most of its reads.

    entries: (seq_id, sequence) pairs; anchors: name -> sequence. Both
    strands are searched, since reads arrive in either orientation.
    """
    rows = []
    for seq_id, seq in entries:
        seq = (seq or "").upper()
        for name, anchor in anchors.items():
            if not anchor:
                continue
            for strand, s in (("+", seq), ("-", revcomp(seq))):
                off = 0
                while off < len(s):
                    h = find_anchor_leftmost(s[off:], anchor, max_edits)
                    if h is None:
                        break
                    st, en, cost = h[0] + off, h[1] + off, h[2]
                    rows.append({
                        "seq_id": seq_id, "anchor": name, "strand": strand,
                        "position": st, "edits": cost,
                        "context": s[max(0, st - flank):en + flank],
                    })
                    off = max(en, st + 1)
    return rows


