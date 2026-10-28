#!/usr/bin/env python3
"""What does PR2 look like through the VESPA primers?

In-silico PCR, insert lengths, resolution over the whole insert (nanopore)
and over the two 2x150 read windows (Illumina), and anchor look-alikes.
See docs/DESIGN.md for what it found.

    bin/pr2_insilico_pcr.py <pr2_..._dada2.fasta.gz> <outdir>

The in-silico PCR is cached as <outdir>/pcr_results.json (about 120MB,
not committed); a rerun reuses it and takes seconds rather than ten
minutes. SKIP_LOOKALIKES=1 skips the anchor scan.
"""
import gzip, json, os, re, sys
from collections import Counter, defaultdict
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seqmatch import ANCHORS, IUPAC, find_anchor_cost, find_anchor_lookalikes  # noqa

PR2 = sys.argv[1]
OUT = sys.argv[2]
FWD = "AGCAGCCGCGGTAATTCC"
REV = "TCCGTCAATTYCTTNAASTTTC"
COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
REV_RC = REV.translate(COMP)[::-1]
MAX_EDITS = 2
READ_LEN = 150
W1, W2 = READ_LEN - len(FWD), READ_LEN - len(REV)
RANKS = ["domain", "supergroup", "division", "subdivision", "class",
         "order", "family", "genus", "species"]


def to_regex(p):
    return re.compile("".join(c if c in "ACGT" else f"[{IUPAC[c]}]" for c in p))


RF, RR = to_regex(FWD), to_regex(REV_RC)


def informative(name):
    return bool(name) and not re.search(r"(_X+(_|$)|_sp\.?$|_X+_sp)", name)


def pcr(item):
    tax, seq = item
    seq = seq.upper()
    m = RF.search(seq)
    if m:
        fs, fe, fc = m.start(), m.end(), 0
    else:
        h = find_anchor_cost(seq, FWD, MAX_EDITS)
        fs, fe, fc = h if h else (None, None, None)
    if fe is None:
        # Is the reverse site there at all? (truncated at the 5' end)
        r = RR.search(seq) or find_anchor_cost(seq, REV_RC, MAX_EDITS)
        return (tax, "R_only" if r else "neither", None, None, None)
    rest = seq[fe:]
    m = RR.search(rest)
    if m:
        rs, re_, rc = m.start(), m.end(), 0
    else:
        h = find_anchor_cost(rest, REV_RC, MAX_EDITS)
        rs, re_, rc = h if h else (None, None, None)
    if rs is None:
        return (tax, "F_only", None, None, None)
    return (tax, "amplicon", rest[:rs], fc, rc)


def read_fasta(path):
    tax, buf = None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if tax is not None:
                    yield tax, "".join(buf)
                tax, buf = line[1:].rstrip(";").split(";"), []
            else:
                buf.append(line)
    if tax is not None:
        yield tax, "".join(buf)


def group_of(tax):
    t = tax + [""] * 9
    if ":plas" in t[0]: return "plastid"
    if ":mito" in t[0]: return "mitochondrion"
    if t[0] != "Eukaryota": return t[0]  # Bacteria / Archaea
    for name in ("Nematoda", "Platyhelminthes", "Acanthocephala"):
        if t[4] == name: return name
    if t[4] == "Craniata": return "Craniata (host)"
    if t[3] == "Metazoa": return "other Metazoa"
    if t[3] == "Fungi": return "Fungi"
    if t[4] == "Apicomplexa": return "Apicomplexa"
    if t[2] == "Metamonada": return "Metamonada"
    if t[7] == "Blastocystis": return "Blastocystis"
    if t[7] == "Entamoeba": return "Entamoeba"
    if t[3] == "Ciliophora": return "Ciliophora"
    return "other eukaryotes"


def resolved_depth(entries_by_key, key_of, amps):
    """Per amplicon entry: 9 if every entry sharing its key has the same
    species, 8 if the same genus, else 0. Names are compared directly, so
    a placeholder rank higher up (Chromadorea_X) does not stop the walk."""
    groups = defaultdict(list)
    for i, a in enumerate(amps):
        groups[key_of(a)].append(i)
    depth = [0] * len(amps)
    for idxs in groups.values():
        taxa = [amps[i]["tax"] + [""] * 9 for i in idxs]
        sp = {t[8] for t in taxa}; ge = {t[7] for t in taxa}
        d = 9 if len(sp) == 1 else (8 if len(ge) == 1 else 0)
        for i in idxs:
            depth[i] = d
    return depth, len(groups)


def lookalike_job(seqs):
    return find_anchor_lookalikes(seqs, ANCHORS, 3)


def pct(n, d):
    return f"{100 * n / d:5.1f}%" if d else "    -"


def main():
    cache = os.path.join(OUT, "pcr_results.json")
    if os.path.exists(cache):
        res = json.load(open(cache)); n_records = len(res)
    else:
        records = list(read_fasta(PR2))
        with Pool(os.cpu_count()) as pool:
            res = pool.map(pcr, records, chunksize=500)
        os.makedirs(OUT, exist_ok=True)
        json.dump(res, open(cache, "w"))
    n_records = len(res)

    status = Counter(r[1] for r in res)
    status_by_group = defaultdict(Counter)
    amps = []
    for tax, st, ins, fc, rc in res:
        g = group_of(tax)
        status_by_group[g][st] += 1
        if st == "amplicon":
            amps.append({"tax": tax, "group": g, "ins": ins,
                         "fc": fc, "rc": rc})

    # Resolution: nanopore = whole insert; Illumina = two read windows.
    nano_depth, n_nano = resolved_depth(None, lambda a: a["ins"], amps)
    def ill_key(a):
        s = a["ins"]
        return s if len(s) <= W1 + W2 else (s[:W1], s[-W2:])
    ill_depth, n_ill = resolved_depth(None, ill_key, amps)

    lines = []
    P = lines.append
    P(f"PR2 entries: {n_records:,}")
    P(f"primers: F {FWD}  R {REV}  (search: exact, then <= {MAX_EDITS} edits)")
    P("")
    P("In-silico PCR outcome:")
    for k in ("amplicon", "F_only", "R_only", "neither"):
        P(f"  {k:9s} {status[k]:>8,}  {pct(status[k], n_records)}")
    fuzzy = sum(1 for a in amps if a["fc"] or a["rc"])
    P(f"  amplicons needing >=1 primer edit: {fuzzy:,} ({pct(fuzzy, len(amps)).strip()})")
    P("")

    groups = sorted(status_by_group, key=lambda g: -sum(status_by_group[g].values()))
    P("By group: entries | amplified | spans-neither (partial seqs) | insert length p5/median/p95")
    by_group_amps = defaultdict(list)
    for i, a in enumerate(amps):
        by_group_amps[a["group"]].append(i)
    for g in groups:
        c = status_by_group[g]; tot = sum(c.values())
        L = sorted(len(amps[i]["ins"]) for i in by_group_amps[g])
        q = (lambda p: L[min(len(L) - 1, int(p * len(L)))]) if L else (lambda p: "-")
        P(f"  {g:18s} {tot:>7,} | {c['amplicon']:>6,} ({pct(c['amplicon'], tot).strip():>6}) "
          f"| {c['neither']:>6,} | {q(.05)}/{q(.5)}/{q(.95)}")
    P("")
    allL = sorted(len(a["ins"]) for a in amps)
    over = sum(1 for x in allL if x > W1 + W2)
    P(f"Insert length (all amplicons): min {allL[0]}, median {allL[len(allL)//2]}, max {allL[-1]}")
    P(f"Illumina windows: R1 {W1} + R2 {W2} = {W1 + W2} insert bases; "
      f"{over:,} of {len(allL):,} inserts ({pct(over, len(allL)).strip()}) are longer, i.e. leave an unread gap")
    P("")
    P(f"Unique keys: nanopore (whole insert) {n_nano:,}; Illumina (two windows) {n_ill:,}")
    P("")
    P("Resolution, per amplified entry with an informative species label:")
    P("  rank reached = deepest rank shared by every entry with an identical sequence (best case: exact sequence)")
    P(f"  {'group':18s} {'n':>6}  {'nanopore sp / gen / fam+':>26}   {'Illumina sp / gen / fam+':>26}")
    detail = {}
    for g in groups:
        idx = [i for i in by_group_amps[g] if informative((amps[i]["tax"] + [""] * 9)[8])]
        if not idx:
            continue
        def summ(depths):
            sp = sum(1 for i in idx if depths[i] >= 9)
            ge = sum(1 for i in idx if depths[i] == 8)
            return sp, ge, len(idx) - sp - ge
        a, b = summ(nano_depth), summ(ill_depth)
        detail[g] = {"n": len(idx), "nanopore": a, "illumina": b}
        P(f"  {g:18s} {len(idx):>6,}  {pct(a[0], len(idx))} {pct(a[1], len(idx))} {pct(a[2], len(idx))}"
          f"      {pct(b[0], len(idx))} {pct(b[1], len(idx))} {pct(b[2], len(idx))}")
    P("")

    # Helminth genera that Illumina windows lose but nanopore keeps.
    lost = Counter()
    for i, a in enumerate(amps):
        if a["group"] in ("Nematoda", "Platyhelminthes", "Acanthocephala"):
            if nano_depth[i] >= 9 and ill_depth[i] < 9:
                lost[(a["group"], (a["tax"] + [""] * 9)[7])] += 1
    P("Helminth genera with most entries resolved to species on nanopore but NOT on Illumina:")
    for (g, gen), n in lost.most_common(20):
        P(f"  {g:16s} {gen:30s} {n}")
    P("")

    # Genera of interest: how each fares, and what it collides with.
    GENERA = ["Haemonchus", "Teladorsagia", "Trichostrongylus", "Cooperia",
              "Ostertagia", "Nematodirus", "Oesophagostomum", "Chabertia",
              "Ancylostoma", "Necator", "Strongyloides", "Trichuris",
              "Ascaris", "Toxocara", "Enterobius", "Onchocerca", "Brugia",
              "Wuchereria", "Dirofilaria", "Schistosoma", "Fasciola",
              "Clonorchis", "Opisthorchis", "Taenia", "Echinococcus",
              "Hymenolepis", "Giardia", "Cryptosporidium", "Entamoeba",
              "Blastocystis", "Eimeria", "Plasmodium", "Dientamoeba"]
    keymap = {}
    for name, key_of in (("nano", lambda a: a["ins"]), ("ill", ill_key)):
        m = defaultdict(set)
        for a in amps:
            t = a["tax"] + [""] * 9
            m[key_of(a)].add(t[7])
        keymap[name] = (key_of, m)
    P("Genera of interest (entries with an informative species label):")
    P(f"  {'genus':16s} {'n':>4} {'species':>7}  {'nano sp':>7} {'nano gen':>8}  {'ill sp':>6} {'ill gen':>7}  shares a sequence with (Illumina)")
    for gname in GENERA:
        idx = [i for i, a in enumerate(amps) if (a["tax"] + [""] * 9)[7] == gname
               and informative((a["tax"] + [""] * 9)[8])]
        if not idx:
            P(f"  {gname:16s}    0  (none amplified / no species-labelled entries)")
            continue
        nsp = len({amps[i]["tax"][8] for i in idx})
        ns = sum(nano_depth[i] >= 9 for i in idx); ng = sum(nano_depth[i] >= 8 for i in idx)
        is_ = sum(ill_depth[i] >= 9 for i in idx); ig = sum(ill_depth[i] >= 8 for i in idx)
        key_of, m = keymap["ill"]
        other = Counter(g for i in idx for g in m[key_of(amps[i])] if g != gname)
        P(f"  {gname:16s} {len(idx):>4} {nsp:>7}  {pct(ns, len(idx))} {pct(ng, len(idx)):>8}  {pct(is_, len(idx))} {pct(ig, len(idx)):>7}  "
          + ", ".join(g for g, _ in sorted(other.items(), key=lambda kv: (-kv[1], kv[0]))[:6]))
    P("")

    if os.environ.get("SKIP_LOOKALIKES"):
        rows = []
    # Anchor look-alikes over unique full amplicons (primer + insert + primer).
    uniq = {}
    for a in amps:
        uniq.setdefault(a["ins"], a)
    items = [(f"u{k}", FWD + s + REV_RC.replace("N", "A").replace("Y", "C").replace("S", "C").replace("R", "A"))
             for k, s in enumerate(uniq)]
    meta = {f"u{k}": a for k, a in enumerate(uniq.values())}
    chunks = [items[i:i + 200] for i in range(0, len(items), 200)]
    if not os.environ.get("SKIP_LOOKALIKES"):
      with Pool(os.cpu_count()) as pool:
        rows = [r for part in pool.map(lookalike_job, chunks) for r in part]
    by_edit = Counter(r["edits"] for r in rows)
    seqs_hit = defaultdict(set)
    for r in rows:
        seqs_hit[r["edits"] <= 2].add(r["seq_id"])
    P(f"Anchor look-alikes (<=3 edits) in {len(items):,} unique amplicons: {len(rows):,} hits; by edits {dict(sorted(by_edit.items()))}")
    near = {r["seq_id"] for r in rows if r["edits"] <= 2}
    gc = Counter(meta[s]["group"] for s in near)
    P(f"  unique amplicons with a hit at <=2 edits (the splitter's tolerance): {len(near):,}; by group {dict(gc.most_common())}")
    hel = sorted({";".join(meta[s]["tax"][6:9]) for s in near
                  if meta[s]["group"] in ("Nematoda", "Platyhelminthes", "Acanthocephala")})
    P(f"  helminth taxa among them ({len(hel)}): " + ", ".join(hel[:40]))

    text = "\n".join(lines)
    print(text)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "summary.txt"), "w") as fh:
        fh.write(text + "\n")
    with open(os.path.join(OUT, "resolution.json"), "w") as fh:
        json.dump(detail, fh, indent=1)
    with open(os.path.join(OUT, "lengths.tsv"), "w") as fh:
        fh.write("group\tinsert_length\n")
        for a in amps:
            fh.write(f"{a['group']}\t{len(a['ins'])}\n")


if __name__ == "__main__":
    main()
