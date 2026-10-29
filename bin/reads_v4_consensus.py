#!/usr/bin/env python3
"""V4 from short genome reads, by baiting and iterative consensus.

    bin/reads_v4_consensus.py <seed.fasta> <out prefix> <max reads> <fastq.gz URL or path> ... [--preset map-pb]

--decoy <fasta> adds sequences (host 18S, say) that baited reads are also
mapped against; a read whose best hit is a decoy is left out. For blood
or skin samples, where the host's rDNA outnumbers the parasite's.

--preset sets minimap2's mode (default sr, short reads). map-pb or map-ont
take noisy long reads - old PacBio RS subreads, nanopore - where the
primers rarely match within two edits but a consensus over many reads of
each molecule is still exact.

For species with Illumina genome reads but no rDNA in any assembly (the
2010s short-read helminth genomes). rDNA is present in tens to hundreds of
copies, so a few million reads cover V4 many times over. The seed is a
relative's 18S around V4 - both primer sites and ~250bp either side, so
reads overlapping the primers are caught. Reads are streamed (nothing is
stored but those that map), mapped to the seed with minimap2 -x sr, and
the majority base at each position (deletions and insertions included)
becomes the next reference; mapping and voting repeat until the consensus
stops changing. The V4 insert is then cut from the consensus with the
usual in-silico PCR.

Writes <prefix>.consensus.fasta, <prefix>.reads.fastq (the baited reads)
and <prefix>.report.txt: rounds, depth over the insert, and positions
where a second allele holds >= 10% (intragenomic variants, or a mixed
sample). Short-read consensus gives the majority copy only; minor rDNA
variants need long reads.
"""

from __future__ import annotations

import collections
import gzip
import io
import os
import re
import subprocess
import sys
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vespa import pcr_full  # noqa: E402

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def open_source(src):
    """A local file, or a URL. EBI's HTTPS sometimes refuses (403) under
    load: retry, then fall back to its FTP server."""
    if not src.startswith(("ftp.", "http://", "https://", "ftp://")):
        return open(src, "rb")
    import time
    path = src.split("://", 1)[-1]
    for attempt, url in enumerate(["https://" + path, "https://" + path, "ftp://" + path,
                                   "ftp://" + path]):
        try:
            return urllib.request.urlopen(url, timeout=300)
        except Exception as e:  # noqa: BLE001
            err = e
            time.sleep(10 * (attempt + 1))
    raise err


def read_fasta(path):
    seq = []
    for line in open(path):
        if not line.startswith(">"):
            seq.append(line.strip().upper())
    return "".join(seq)


def stream_reads(src, limit):
    raw = open_source(src)
    fh = io.TextIOWrapper(gzip.GzipFile(fileobj=raw), "ascii", errors="replace")
    for i in range(limit * 4):
        line = fh.readline()
        if not line:
            break
        yield line


PRESET = "sr"
DECOY = None


def bait(seed_path, sources, limit, out_fastq):
    """Stream each source's first `limit` reads through minimap2; keep mapped."""
    kept = 0
    ref_path, seed_name = seed_path, open(seed_path).readline()[1:].split()[0]
    if DECOY:
        ref_path = out_fastq + ".bait.fa"
        open(ref_path, "w").write(open(seed_path).read() + open(DECOY).read())
    with open(out_fastq, "w") as out:
        for src in sources:
            p = subprocess.Popen(["minimap2", "-ax", PRESET, "-t", "4", "--secondary=no",
                                  ref_path, "-"], stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
            import threading

            def feed():
                try:
                    for line in stream_reads(src, limit):
                        p.stdin.write(line)
                finally:
                    p.stdin.close()     # minimap2 must see EOF even on error
            t = threading.Thread(target=feed)
            t.start()
            for line in p.stdout:
                if line.startswith("@"):
                    continue
                f = line.split("\t")
                if int(f[1]) & 4 or f[2] != seed_name:
                    continue
                out.write(f"@{f[0]}\n{f[9]}\n+\n{f[10]}\n")
                kept += 1
            t.join()
            p.wait()
    return kept


def consensus(ref, fastq):
    """One round: map reads to ref, vote base by base. Returns (seq, depth, minor)."""
    rpath = fastq + ".ref.fa"
    open(rpath, "w").write(f">ref\n{ref}\n")
    sam = subprocess.run(["minimap2", "-ax", PRESET, "--secondary=no", rpath, fastq],
                         capture_output=True, text=True, check=True).stdout
    base = [collections.Counter() for _ in ref]
    ins = [collections.Counter() for _ in ref]
    for line in sam.splitlines():
        if line.startswith("@"):
            continue
        f = line.split("\t")
        flag = int(f[1])
        if flag & 4 or flag & 256 or flag & 2048:
            continue
        pos, seq = int(f[3]) - 1, f[9]
        q = 0
        for n, op in CIGAR.findall(f[5]):
            n = int(n)
            if op in "M=X":
                for k in range(n):
                    base[pos + k][seq[q + k]] += 1
                pos += n
                q += n
            elif op == "I":
                if pos > 0:
                    ins[pos - 1][seq[q:q + n]] += 1
                q += n
            elif op == "D":
                for k in range(n):
                    base[pos + k]["-"] += 1
                pos += n
            elif op == "S":
                q += n
    out, depth, minor = [], [], []
    for i, c in enumerate(base):
        d = sum(c.values())
        depth.append(d)
        if d == 0:
            out.append(ref[i])
            continue
        (b, nb), *rest = c.most_common(2) + [("", 0)]
        if rest and rest[0][1] >= 0.1 * d:
            minor.append((i, b, nb, rest[0][0], rest[0][1]))
        if b != "-":
            out.append(b)
        top = ins[i].most_common(1)
        if top and top[0][0] and top[0][1] > d / 2:
            out.append(top[0][0])
    return "".join(out), depth, minor


def main() -> int:
    global PRESET, DECOY
    args = sys.argv[1:]
    if "--decoy" in args:
        k = args.index("--decoy")
        DECOY = args[k + 1]
        del args[k:k + 2]
    if "--preset" in args:
        k = args.index("--preset")
        PRESET = args[k + 1]
        del args[k:k + 2]
    seed_path, prefix, limit, *sources = args
    limit = int(limit)
    fastq = prefix + ".reads.fastq"
    kept = bait(seed_path, sources, limit, fastq)
    ref = read_fasta(seed_path)
    log = [f"baited {kept} reads from {len(sources)} file(s), first {limit:,} reads each"]
    for rnd in range(1, 8):
        new, depth, minor = consensus(ref, fastq)
        changed = new != ref
        log.append(f"round {rnd}: length {len(new)}, {'changed' if changed else 'stable'}")
        ref = new
        if not changed:
            break
    status, ins, fc, rc, *_ = pcr_full(ref)
    if status == "amplicon":
        start = ref.find(ins)
        d = depth[start:start + len(ins)] if len(depth) == len(ref) else depth
        ds = sorted(d)
        log.append(f"V4 insert {len(ins)}bp, primer edits {fc}/{rc}; depth over insert "
                   f"min {ds[0]}, median {ds[len(ds) // 2]}")
        mi = [m for m in minor if start <= m[0] < start + len(ins)]
        log.append(f"positions in the insert with a second allele >= 10%: {len(mi)}")
        for i, b, nb, b2, n2 in mi[:20]:
            log.append(f"  insert pos {i - start + 1}: {b} {nb} / {b2} {n2}")
    else:
        log.append(f"no V4 insert in the consensus ({status})")
    open(prefix + ".consensus.fasta", "w").write(f">{os.path.basename(prefix)}\n{ref}\n")
    if status == "amplicon":
        open(prefix + ".v4.fasta", "w").write(f">{os.path.basename(prefix)}_V4\n{ins}\n")
    open(prefix + ".report.txt", "w").write("\n".join(log) + "\n")
    print("\n".join(log))
    return 0


if __name__ == "__main__":
    sys.exit(main())
