#!/usr/bin/env python3
"""The README figures, as static SVG in light and dark variants.

    bin/readme_figures.py releases/vespa18s_<version> docs/figures

Writes the four SVGs, and prints to stdout the two README <picture>
blocks with their alt text built from the same numbers - paste them over
the old ones. <outdir> is used as given in those blocks, so pass it
relative to the README (docs/figures).

1. sources_<mode>.svg  records each source contributes to, split into
                       records only that source provides and records it
                       shares with another source
2. taxa_<mode>.svg     records per major group, split into those that
                       resolve to species read whole and those that stop
                       at a higher rank

Colours are the dataviz reference palette's validated slots (categorical
blue and orange; a two-step blue ramp for resolution), stepped for each
mode's surface. Standard library only.
"""

from __future__ import annotations

import collections
import csv
import os
import sys
from xml.sax.saxutils import escape

FONT = 'system-ui, -apple-system, &quot;Segoe UI&quot;, sans-serif'
MODES = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#898781",
                  grid="#e1e0d9", axis="#c3c2b7", s1="#2a78d6", s2="#eb6834",
                  seq_dark="#256abf", seq_light="#86b6ef", ring="rgba(11,11,11,0.10)"),
    "dark": dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", muted="#898781",
                 grid="#2c2c2a", axis="#383835", s1="#3987e5", s2="#d95926",
                 seq_dark="#3987e5", seq_light="#184f95", ring="rgba(255,255,255,0.10)"),
}
SOURCE_NAMES = {"PR2_5.1.1": "PR2 v5.1.1", "GenBank": "GenBank", "WBPS19": "WormBase ParaSite 19",
                "NCBI_genomes": "NCBI protist genomes", "ParaRef": "ParaRef",
                "local": "Local (our sequences)"}


def group_of(t: dict) -> str:
    dom, sg, div, sub, cl = (t[k].split(":")[0] for k in
                             ("domain", "supergroup", "division", "subdivision", "class"))
    if ":" in t["domain"]:
        return "Organelle 16S"
    if dom in ("Bacteria", "Archaea"):
        return "Bacteria and archaea"
    if cl in ("Nematoda", "Platyhelminthes", "Acanthocephala"):
        return cl
    if cl == "Craniata":
        return "Vertebrates"
    if sub == "Metazoa":
        return "Other animals"
    if sub == "Fungi":
        return "Fungi"
    if sub == "Apicomplexa":
        return "Apicomplexa"
    if div == "Alveolata":
        return "Other alveolates"
    if div == "Stramenopiles":
        return "Stramenopiles"
    if div == "Metamonada":
        return "Metamonada"
    if div == "Discoba":
        return "Discoba"
    if sg == "Amoebozoa":
        return "Amoebozoa"
    if div == "Rhizaria":
        return "Rhizaria"
    if sg == "Archaeplastida":
        return "Plants and algae"
    return "Other eukaryotes"


def nice_ticks(vmax: float):
    for step in (1, 2, 2.5, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 2500,
                 5000, 10000, 20000, 25000, 50000):
        if vmax / step <= 5:
            n = int(vmax // step) + 1
            return [step * i for i in range(n + 1)], step * n
    return [0, vmax], vmax


def fmt(n: int) -> str:
    return f"{n:,}"


def tick_label(v: float) -> str:
    return f"{v / 1000:g}k" if v >= 1000 else f"{v:g}"


def hbar_svg(title, subtitle, rows, seg_names, seg_colors, mode, note=None,
             emphasise=(), value_label=None, right=92):
    """rows: [(label, [seg values])]. Stacked horizontal bars, one axis.
    value_label(vals) -> the text at the bar end (default: the total)."""
    c = MODES[mode]
    W, left = 760, 196
    top, row_h, bar_h = 104, 28, 16
    plot_w = W - left - right
    H = top + row_h * len(rows) + 46 + (22 if note else 0)
    ticks, vmax = nice_ticks(max(sum(v) for _, v in rows))
    x = lambda v: left + plot_w * v / vmax  # noqa: E731
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
         f'viewBox="0 0 {W} {H}" font-family="{FONT}" role="img" '
         f'aria-label="{escape(title)}">',
         f'<rect width="{W}" height="{H}" rx="8" fill="{c["surface"]}" stroke="{c["ring"]}"/>',
         f'<text x="24" y="34" font-size="17" font-weight="600" fill="{c["ink"]}">{escape(title)}</text>',
         f'<text x="24" y="56" font-size="13" fill="{c["ink2"]}">{escape(subtitle)}</text>']
    # legend, one row, above the plot
    lx = left
    for name, col in zip(seg_names, seg_colors):
        o.append(f'<rect x="{lx}" y="72" width="12" height="12" rx="3" fill="{c[col]}"/>')
        o.append(f'<text x="{lx + 18}" y="82.5" font-size="12.5" fill="{c["ink2"]}">{escape(name)}</text>')
        lx += 18 + 7.2 * len(name) + 24
    # recessive grid and axis
    y0, y1 = top - 6, top + row_h * len(rows)
    for t in ticks:
        o.append(f'<line x1="{x(t):.1f}" y1="{y0}" x2="{x(t):.1f}" y2="{y1}" '
                 f'stroke="{c["grid" if t else "axis"]}" stroke-width="1"/>')
        o.append(f'<text x="{x(t):.1f}" y="{y1 + 18}" font-size="11.5" fill="{c["muted"]}" '
                 f'text-anchor="middle" font-variant-numeric="tabular-nums">{tick_label(t)}</text>')
    for k, (label, vals) in enumerate(rows):
        yc = top + row_h * k + row_h / 2
        weight = "600" if label in emphasise else "400"
        o.append(f'<text x="{left - 12}" y="{yc + 4.5}" font-size="13" font-weight="{weight}" '
                 f'fill="{c["ink"]}" text-anchor="end">{escape(label)}</text>')
        segs = [(v, col) for v, col in zip(vals, seg_colors) if v > 0]
        xa = x(0)
        for j, (v, col) in enumerate(segs):
            x0 = xa
            # Every non-zero segment stays visible - a source's handful of
            # records only it provides must not vanish beside PR2's 89,000.
            # The label carries the exact values.
            w = max(x(v) - x(0), 4)
            xa = x0 + w
            last = j == len(segs) - 1
            gap = 0 if last else 2  # surface gap between stacked segments
            ww = max(w - gap, 1.5)
            if last:  # rounded data end, square at the baseline / joint
                r = min(4, ww / 2, bar_h / 2)
                d = (f"M{x0:.1f},{yc - bar_h / 2:.1f} h{ww - r:.1f} "
                     f"a{r},{r} 0 0 1 {r},{r} v{bar_h - 2 * r:.1f} "
                     f"a{r},{r} 0 0 1 -{r},{r} h-{ww - r:.1f} z")
                o.append(f'<path d="{d}" fill="{c[col]}"/>')
            else:
                o.append(f'<rect x="{x0:.1f}" y="{yc - bar_h / 2:.1f}" width="{ww:.1f}" '
                         f'height="{bar_h}" fill="{c[col]}"/>')
        text = value_label(vals) if value_label else fmt(sum(vals))
        o.append(f'<text x="{xa + 8:.1f}" y="{yc + 4.5}" font-size="12.5" '
                 f'fill="{c["ink2"]}" font-variant-numeric="tabular-nums">{escape(text)}</text>')
    if note:
        o.append(f'<text x="24" y="{H - 16}" font-size="12" fill="{c["muted"]}">{escape(note)}</text>')
    o.append("</svg>")
    return "\n".join(o)


def main() -> int:
    rel, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    version = os.path.basename(rel.rstrip("/")).replace("vespa18s_", "")

    sources = collections.defaultdict(set)
    with open(os.path.join(rel, "members.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["role"] == "kept":
                sources[r["id"]].add(r["source"])
    per, only = collections.Counter(), collections.Counter()
    for s in sources.values():
        for x in s:
            per[x] += 1
        if len(s) == 1:
            only[next(iter(s))] += 1
    src_rows = [(SOURCE_NAMES.get(s, s), [only[s], per[s] - only[s]])
                for s, _ in per.most_common()]

    with open(os.path.join(rel, "taxonomy.tsv")) as fh:
        tax = list(csv.DictReader(fh, delimiter="\t"))
    n, sp = collections.Counter(), collections.Counter()
    for t in tax:
        g = group_of(t)
        n[g] += 1
        sp[g] += t["deepest_rank"] == "species"
    tax_rows = [(g, [sp[g], n[g] - sp[g]]) for g, _ in n.most_common()]

    for mode in MODES:
        with open(os.path.join(out, f"sources_{mode}.svg"), "w") as fh:
            fh.write(hbar_svg(
                f"Where the {len(sources):,} reference sequences come from",
                f"vespa18s {version}: records each source contributes to",
                src_rows, ["Only from this source", "Shared with another source"],
                ["s1", "s2"], mode,
                note="A record is one unique V4 sequence; it can be backed by several sources.",
                value_label=lambda v: f"{fmt(sum(v))}  ·  {fmt(v[0])} only here", right=176))
        with open(os.path.join(out, f"taxa_{mode}.svg"), "w") as fh:
            fh.write(hbar_svg(
                "Taxonomic spread of the reference sequences",
                f"vespa18s {version}: records per group, and how many resolve to species",
                tax_rows, ["Resolves to species", "Stops at genus or above"],
                ["seq_dark", "seq_light"], mode,
                note="Bold: the parasite groups the database was built for. Resolution read "
                     "whole (nanopore); through 2x150 Illumina windows it is lower.",
                emphasise=("Nematoda", "Platyhelminthes", "Acanthocephala", "Apicomplexa")))
    print(f"wrote {out}/sources_{{light,dark}}.svg and taxa_{{light,dark}}.svg",
          file=sys.stderr)

    # The README blocks, alt text included, from the same numbers - so the
    # alt text cannot drift from the figures. Paste them over the old ones.
    src_alt = (f"Bar chart: where the {len(sources):,} reference sequences in {version} "
               "come from. Records each source backs, and how many only that source "
               "provides: " + "; ".join(f"{name} {fmt(sum(v))} ({fmt(v[0])} only here)"
                                        for name, v in src_rows) + ".")
    tax_alt = (f"Bar chart: records per major taxonomic group in {version}, and how many "
               "resolve to species read whole: " +
               "; ".join(f"{g} {fmt(sum(v))} ({fmt(v[0])} to species)" for g, v in tax_rows)
               + ".")
    for stem, alt in (("sources", src_alt), ("taxa", tax_alt)):
        print(f'<picture>\n'
              f'  <source media="(prefers-color-scheme: dark)" srcset="{out}/{stem}_dark.svg">\n'
              f'  <img alt="{escape(alt, {chr(34): "&quot;"})}" src="{out}/{stem}_light.svg" '
              f'width="760">\n'
              f'</picture>\n')
    return 0


if __name__ == "__main__":
    sys.exit(main())
