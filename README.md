# vespa-18s-refdb

Builds the 18S V4 reference database for VESPA-primer metabarcoding:
everything the primers amplify, trimmed to the amplicon, with taxonomy,
resolution tables and the adapter look-alike mask. It is consumed by
[nf-helminth-amplicon](../nf-helminth-amplicon) in its
`--reference_database` mode, and kept separate from it because building
a reference is its own job - its own inputs, its own release cycle - and
its output is a versioned artefact, not a pipeline step.

**Status: exploratory scripts and first measurements. No release has
been built yet.** See [docs/DESIGN.md](docs/DESIGN.md) for the design,
what PR2 and the WormBase ParaSite genomes look like through these
primers, and the release format the pipeline expects.

Primers, as ordered (locus-specific parts; both carry Nextera tails):

- `AMPSEQ_VESPA18S_29F` `AGCAGCCGCGGTAATTCC`
- `AMPSEQ_VESPA18S_21b8_I_R_Nmod` `TCCGTCAATTYCTTNAASTTTC`

## Layout

| path | |
|---|---|
| `bin/pr2_insilico_pcr.py` | In-silico PCR over PR2; insert lengths; resolution over the whole insert (nanopore) and the 2x150 read windows (Illumina); adapter-anchor look-alikes |
| `bin/wbps_insilico_pcr.py` | Streams every WormBase ParaSite genome from the EBI FTP, finds V4 amplicons on both strands, keeps only the hits |
| `bin/compare_wbps_pr2.py` | Per genome: is the species in PR2, does its V4 match PR2 and under whose label, nearest PR2 taxon if not, possible contaminants |
| `bin/seqmatch.py` | Approximate matching and the look-alike scan, copied from nf-helminth-amplicon's `bin/ampseqlib.py`. The look-alike scan must stay identical to the pipeline's, because the pipeline reads the mask built here |
| `results/` | Reports and data from the runs the design note quotes |
| `docs/DESIGN.md` | Design, measurements, release format |

Python 3.9+, standard library only.

## Reproducing

```bash
# PR2 v5.1.0, DADA2 format, from https://github.com/pr2database/pr2database
python3 bin/pr2_insilico_pcr.py pr2_version_5.1.0_SSU_dada2.fasta.gz pr2_out

curl -s https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/current/CHECKSUMS \
    | awk '/genomic.fa.gz$/ && !/masked/ {print $2}' > wbps_genomes.txt
python3 bin/wbps_insilico_pcr.py wbps_genomes.txt wbps_out 8

python3 bin/compare_wbps_pr2.py pr2_out/pcr_results.json wbps_out wbps_vs_pr2
```

On a laptop: PR2 in-silico PCR about 10 minutes (cached as
`pcr_results.json`, about 120MB, and reused), the look-alike scan a few
minutes, all 274 WBPS19 genomes 35 minutes on 8 workers (about 20GB
streamed, nothing kept), the comparison about 10 minutes.

## results/

| file | contents |
|---|---|
| `pr2_v5.1.0_summary.txt` | amplification, insert lengths, resolution, genera of interest, anchor look-alikes |
| `wbps19_vs_pr2.txt` | genome-versus-PR2 report |
| `wbps19_scan.tsv` | one row per WBPS19 assembly: bases scanned, amplicons, seconds (one failed download listed) |
| `wbps19_v4_amplicons.tsv.gz` | every V4 amplicon found in WBPS19: genome, contig, strand, coordinates, primer edits, insert. Raw material for the genome layer **before** the filtering rules in the design note - it still contains host, vector and other contaminant rDNA |

## Caveats

- Resolution is a best case: "every entry with an identical sequence
  shares this genus/species". Sequencing error and missing references
  both make the real figure worse.
- Primer matching allows <= 2 edits in PR2; in genomes, <= 1 edit in the
  forward primer (seeded search) and <= 2 in the reverse.
- A genome hit is discarded if a primer site contains N or the insert is
  more than 1% N: N in a scaffold gap matches a primer for free.
- "Nearest" in the comparison is by shared 12-mers, not an alignment.
