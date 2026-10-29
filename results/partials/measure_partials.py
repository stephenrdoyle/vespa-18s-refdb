"""Exploratory, v0.7.0 (docs/DESIGN.md §8): step 1 of partial references: how far does each near-complete partial
fall short of the V4 insert end it is missing?

For each species lost to the database whose best partial record has one
primer site and 350-699bp of insert after (F) or before (R) it, align the
fragment to full V4 records of the same genus (else family) in v0.7.0:
fragment fully aligned, reference free at the far end. missing = reference
bases beyond the fragment's end; overhang = fragment bases beyond the
reference's end (sequence continues past where the other primer should be).
"""
import collections, csv, glob, gzip, json, sys
sys.path.insert(0, "bin")
from vespa import _RF, _RR, FWD, REV_RC, PRIMER_EDITS, find_anchor_cost

# Run from the repository root: python3 results/partials/measure_partials.py results/partials
S = sys.argv[1]
rows = [l.rstrip("\n").split("\t") for l in open(f"{S}/partial_only_species.tsv")][1:]
rows = [(sp, g, side, int(n), int(k), acc) for sp, g, side, n, k, acc in rows]
want = {acc: (sp, g, side, n) for sp, g, side, n, k, acc in rows
        if side in ("F", "R") and 350 <= n < 700}

COMP = str.maketrans("ACGTRYKMBDHVN", "TGCAYRMKVHDBN")


def site(seq, rx, p):
    m = rx.search(seq)
    if m:
        return m.start(), m.end()
    h = find_anchor_cost(seq, p, PRIMER_EDITS)
    return (h[0], h[1]) if h else None


def fragment(seq):
    seq = seq.upper()
    for s in (seq, seq.translate(COMP)[::-1]):
        f, r = site(s, _RF, FWD), site(s, _RR, REV_RC)
        if f and not r:
            return "F", s[f[1]:]
        if r and not f:
            return "R", s[:r[0]]
    return None, None


def fasta(fh):
    h, b = None, []
    for l in fh:
        l = l.rstrip()
        if l.startswith(">"):
            if h:
                yield h, "".join(b)
            h, b = l[1:], []
        else:
            b.append(l)
    if h:
        yield h, "".join(b)


frag = {}
for h, s in fasta(gzip.open("sources/pr2_version_5.1.1_SSU_taxo_long.fasta.gz", "rt")):
    f = h.split("|")
    if f[0] in want:
        frag[f[0]] = (f[-2], f[-3], *fragment(s))      # genus, family
for p in glob.glob("work/genbank/*.fasta.gz"):
    for h, s in fasta(gzip.open(p, "rt")):
        acc = h.split()[0]
        if acc in want:
            t = h.split(None, 2)[2].split()
            frag[acc] = (t[0], "", *fragment(s))

R = "releases/vespa18s_v0.7.0/"
tax = {t["id"]: t for t in csv.DictReader(open(R + "taxonomy.tsv"), delimiter="\t")}
seqs, name = {}, None
for l in open(R + "refdb.fasta"):
    if l.startswith(">"):
        name = l[1:].strip()
    else:
        seqs[name] = seqs.get(name, "") + l.strip()
by_genus, by_family = collections.defaultdict(list), collections.defaultdict(list)
for i, t in tax.items():
    by_genus[t["genus"].split(":")[0].rstrip("0123456789")].append(i)
    by_family[t["family"].split(":")[0]].append(i)


def end_scores(p, t):
    """Edit distance of all of p against t[:j], for each j (Myers; both
    start at position 0, text end free). Returns list over j = 0..len(t)."""
    m = len(p)
    peq = collections.defaultdict(int)
    for i, c in enumerate(p):
        peq[c] |= 1 << i
    full, top = (1 << m) - 1, 1 << (m - 1)
    pv, mv, score, out = full, 0, m, [m]   # D[m][0] = m
    for c in t:
        eq = peq[c]
        xv = eq | mv
        xh = (((eq & pv) + pv) ^ pv) | eq
        ph = mv | (~(xh | pv) & full)
        mh = pv & xh
        if ph & top:
            score += 1
        elif mh & top:
            score -= 1
        ph = ((ph << 1) | 1) & full    # anchored start: D[0][j] = j
        mh = (mh << 1) & full
        pv = mh | (~(xv | ph) & full)
        mv = ph & xv
        out.append(score)
    return out


def measure(side, fr, ref):
    """(edits, missing, overhang) with the fragment anchored at the
    primer end it has and free at the other."""
    if side == "R":
        fr, ref = fr[::-1], ref[::-1]
    # Anchored start: score fragment prefixes against ref; we want the
    # whole fragment against a ref prefix, or a fragment prefix against
    # the whole ref (overhang). Use end_scores on (fragment vs ref) and
    # (ref vs fragment).
    a = end_scores(fr, ref)            # whole fragment, ref end free
    j = min(range(len(a)), key=lambda x: (a[x], -x))
    b = end_scores(ref, fr)            # whole ref, fragment end free
    k = min(range(len(b)), key=lambda x: (b[x], -x))
    if b[k] <= a[j] and k < len(fr):
        return b[k], 0, len(fr) - k
    return a[j], len(ref) - j, 0


out = []
for acc, (sp, g, side, n) in want.items():
    if acc not in frag or frag[acc][2] is None:
        continue
    genus, family, fside, fr = frag[acc]
    refs = by_genus.get(genus.rstrip("0123456789"), []) or by_family.get(family, [])
    refs = [i for i in refs if tax[i]["deepest_rank"] in ("species", "genus")][:40]
    if not refs:
        out.append([sp, g, acc, fside, str(len(fr)), "", "", "", "", "no full record in genus or family"])
        continue
    best = min((measure(fside, fr, seqs[i]) + (i,) for i in refs), key=lambda x: (x[0], x[1]))
    e, miss, over, ref = best
    out.append([sp, g, acc, fside, str(len(fr)), ref, tax[ref]["species"] or tax[ref]["genus"],
                str(e), str(miss), str(over)])

with open(f"{S}/near_complete_partials.tsv", "w") as fh:
    fh.write("species\tgroup\taccession\thas_primer\tfragment_len\tref\tref_taxon\tedits\tmissing_bp\toverhang_bp\n")
    for r in out:
        fh.write("\t".join(r) + "\n")
print(len(out), "measured")
