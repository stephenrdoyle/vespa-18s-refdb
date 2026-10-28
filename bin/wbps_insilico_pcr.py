#!/usr/bin/env python3
"""Stream every WBPS genome, find VESPA 18S V4 amplicons, keep only the hits.

Forward primer: found via exact 9-mer seeds from either half, verified at
<= 1 edit (pigeonhole: one edit leaves one half intact). Reverse primer:
searched at <= 2 edits only in the window downstream of a forward hit.
Both strands. Each genome is streamed, scanned and discarded.

    bin/wbps_insilico_pcr.py <genome list> <outdir> [workers]

The genome list is paths relative to the release directory, one per line:

    curl -s https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/current/CHECKSUMS \\
        | awk '/genomic.fa.gz$/ && !/masked/ {print $2}' > wbps_genomes.txt

Resumable: a genome with a .done marker in <outdir> is skipped. All 274
WBPS19 assemblies took 35 minutes on 8 workers on a laptop.
"""
import gzip, io, os, sys, time, urllib.request
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seqmatch import find_anchor_cost  # noqa: E402

BASE = "https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/current/"
FWD = "AGCAGCCGCGGTAATTCC"
REV = "TCCGTCAATTYCTTNAASTTTC"
COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
REV_RC = REV.translate(COMP)[::-1]
SEQ_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
MAX_INSERT = 2500
F_EDITS, R_EDITS = 1, 2
SEEDS = [(FWD[:9], 0), (FWD[9:], 9)]   # (seed, offset within primer)


def scan(seq):
    """Yield (f_start, f_end, f_cost, r_start, r_end, r_cost) on this strand."""
    cands = set()
    for seed, off in SEEDS:
        i = seq.find(seed)
        while i != -1:
            cands.add(max(0, i - off - 2))
            i = seq.find(seed, i + 1)
    seen = set()
    for c in sorted(cands):
        win = seq[c:c + len(FWD) + 4]
        h = find_anchor_cost(win, FWD, F_EDITS)
        if not h:
            continue
        fs, fe, fc = c + h[0], c + h[1], h[2]
        if fe in seen:
            continue
        seen.add(fe)
        down = seq[fe:fe + MAX_INSERT + len(REV_RC) + 2]
        r = find_anchor_cost(down, REV_RC, R_EDITS)
        if r:
            yield fs, fe, fc, fe + r[0], fe + r[1], r[2]


def contigs(stream):
    name, buf = None, []
    for line in stream:
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.rstrip())
    if name is not None:
        yield name, "".join(buf).upper()


def one_genome(path, outdir):
    stem = os.path.basename(path).replace(".genomic.fa.gz", "")
    done = os.path.join(outdir, stem + ".done")
    if os.path.exists(done):
        return stem, "skipped", 0, 0
    t0 = time.time()
    rows, bases = [], 0
    for attempt in range(4):
        try:
            rows, bases = [], 0
            with urllib.request.urlopen(BASE + path, timeout=120) as resp:
                text = io.TextIOWrapper(gzip.GzipFile(fileobj=resp), "ascii",
                                        errors="replace")
                for ctg, seq in contigs(text):
                    bases += len(seq)
                    L = len(seq)
                    for strand, s in (("+", seq), ("-", seq.translate(SEQ_COMP)[::-1])):
                        for fs, fe, fc, rs, re_, rc in scan(s):
                            ins = s[fe:rs]
                            # N matches anything, so a scaffold gap after a
                            # forward hit "contains" the reverse primer.
                            if len(ins) < 100 or ins.count("N") > 0.01 * len(ins) \
                                    or "N" in s[fs:fe] or "N" in s[rs:re_]:
                                continue
                            a, b = (fs, re_) if strand == "+" else (L - re_, L - fs)
                            rows.append((ctg, strand, a, b, fc, rc, s[fe:rs]))
            break
        except Exception as e:  # network hiccup: retry the whole genome
            err = e
            time.sleep(10 * (attempt + 1))
    else:
        with open(os.path.join(outdir, stem + ".failed"), "w") as fh:
            fh.write(repr(err) + "\n")
        return stem, "failed", 0, 0
    tmp = os.path.join(outdir, stem + ".tsv.partial")
    with open(tmp, "w") as fh:
        fh.write("genome\tcontig\tstrand\tstart\tend\tf_edits\tr_edits\tinsert_len\tinsert\n")
        for ctg, strand, a, b, fc, rc, ins in rows:
            fh.write(f"{stem}\t{ctg}\t{strand}\t{a}\t{b}\t{fc}\t{rc}\t{len(ins)}\t{ins}\n")
    os.replace(tmp, os.path.join(outdir, stem + ".tsv"))
    with open(done, "w") as fh:
        fh.write(f"{bases}\t{len(rows)}\t{time.time() - t0:.0f}\n")
    return stem, "ok", bases, len(rows)


def main():
    _, listing, outdir, *rest = sys.argv
    workers = int(rest[0]) if rest else 6
    os.makedirs(outdir, exist_ok=True)
    paths = [l.strip() for l in open(listing) if l.strip()]
    t0 = time.time()
    with ProcessPoolExecutor(workers) as pool:
        futs = {pool.submit(one_genome, p, outdir): p for p in paths}
        for k, f in enumerate(as_completed(futs), 1):
            stem, st, bases, n = f.result()
            print(f"[{k}/{len(paths)} {time.time() - t0:6.0f}s] {st:7s} {stem}  "
                  f"{bases / 1e6:8.1f} Mb  {n} amplicons", flush=True)


if __name__ == "__main__":
    main()
