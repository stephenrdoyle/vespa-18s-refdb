# vespa-18s-refdb: design

**Status: proposal and first measurements. No release has been built yet.**

This repository builds the reference database that
[nf-helminth-amplicon](../nf-helminth-amplicon) uses in its
`--reference_database` mode: an 18S V4 database covering everything the
VESPA primers amplify, trimmed to the amplicon, with taxonomy, resolution
tables and the adapter look-alike mask the pipeline needs. How the
pipeline *uses* it - mapping, LCA, the Illumina and nanopore routes,
reporting - is designed in the pipeline's
`docs/REFERENCE_DATABASE.md`. The two meet at the release format in §5.

VESPA: Owens et al. 2024, *Nat Commun* 15:402,
doi:10.1038/s41467-023-44521-3.

---

## 1. Fixed inputs

| | |
|---|---|
| Primers | `AMPSEQ_VESPA18S_29F` `AGCAGCCGCGGTAATTCC` / `AMPSEQ_VESPA18S_21b8_I_R_Nmod` `TCCGTCAATTYCTTNAASTTTC` (locus-specific parts; both carry the standard Nextera tails). Our reverse primer differs from VESPA's published 21b8R (`TCAATTYCTTIAASTTTC`, Supplementary Table 1): a 5′ `TCCG` extension that matches the template, and N in place of deoxyinosine |
| Backbone | PR2 v5.1.0 SSU, DADA2 format: 9 ranks in the header, `Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species;` |
| Read geometry | Illumina 2x150 (fixed for now); nanopore reads the whole amplicon |

## 2. Preparation

Keyed by PR2 version and primer pair; built once and reused across runs.

1. **In-silico PCR.** Find 29F and the reverse complement of the
   reverse primer in every PR2 entry, allowing a small edit budget.
   Entries that do not span both primer sites are set aside and counted.
2. **Trim** to the insert between the primers.
3. **Collapse** identical inserts. Their taxonomy is merged to the ranks
   they share (lowest common ancestor). PR2 placeholder ranks (`_X`,
   `_XX`, `_sp.`) count as unresolved, not as taxa.
4. **Write**, alongside the trimmed database:
   - the **insert length distribution**, by lineage;
   - a **nanopore resolution table** over the whole insert;
   - an **Illumina resolution table** over the two read windows only -
     at 2x150, the first 132 and last 128 bases of the insert, after
     the 18 and 22 base primers - with entries collapsed when both
     windows match.
     Read length is a parameter, so moving to longer reads is a
     recomputation, not a redesign;
   - the **anchor look-alike mask** over every unique insert.

This answers most of the pipeline's open questions with numbers before
any read is assigned.

## 3. What PR2 v5.1.0 looks like through these primers

A first pass (`bin/pr2_insilico_pcr.py`; report in
`results/pr2_v5.1.0_summary.txt`): exact primer
match, falling back to <= 2 edits; resolution measured as "every entry
with an identical sequence shares this genus / species" - a best case,
since it assumes the sequence is read without error.

**Amplification.** 151,931 of 240,200 entries (63%) span both primer
sites. Most of the rest are partial sequences that stop short of one
site, not primer failures. 18% of amplicons needed at least one primer
edit. Nematoda 71%, Platyhelminthes 84%, Acanthocephala 86%, host
(Craniata) 85%. Bacteria and plastid 16S largely do not amplify, and
mitochondrial entries essentially never do.

**Insert length** (between the primers): median 546bp, 5th-95th
percentile mostly ~520-710; Nematoda 524/544/551, Platyhelminthes
544/640/755, Metamonada (Giardia) 403/425/799. At 2x150 the reads cover
260 insert bases, so **99.9% of inserts leave an unread gap** in the
middle - the Illumina route cannot merge, confirmed.

**Resolution** (% of species-labelled amplified entries):

| group | n | nanopore species / genus | Illumina 2x150 species / genus |
|---|---:|---|---|
| Nematoda | 2,411 | 78% / 95% | 65% / 91% |
| Platyhelminthes | 1,719 | 89% / 96% | 78% / 91% |
| Acanthocephala | 94 | 87% / 100% | 84% / 100% |
| Metamonada | 862 | 94% / 100% | 91% / 100% |
| *Entamoeba* | 353 | 96% / 100% | 96% / 100% |
| *Blastocystis* | 217 | 43% / 100% | 26% / 100% |
| Craniata (host) | 1,609 | 49% / 54% | 36% / 40% |

The Illumina windows cost about 10-15 points of species resolution for
helminths; *Schistosoma* is the largest single loss (35% species on
Illumina vs 84% on nanopore).

Genera of interest:

| genus | finding |
|---|---|
| *Trichuris*, *Strongyloides*, *Taenia*, *Echinococcus*, *Giardia*, *Entamoeba* | resolve to genus on both platforms, mostly to species |
| *Haemonchus* | genus yes, species no: the three species share a V4 sequence |
| *Ancylostoma* | genus on both; species on nanopore only |
| *Ascaris* | genus on both; species partly on nanopore (40%), not on Illumina |
| *Trichostrongylus* / *Ostertagia* | share a sequence over the Illumina windows |
| *Clonorchis* / *Opisthorchis* / *Metorchis* | not separable even at genus on Illumina |
| *Toxocara* | collides with *Toxascaris*, *Baylisascaris*, *Porrocaecum* |
| ***Teladorsagia*, *Cooperia*** | **absent from PR2** |
| *Oesophagostomum* | one partial entry, no amplicon |
| *Cryptosporidium* | 373 entries but only 1 spans the amplicon: most stop ~436bp after the forward primer. The one full-length entry matches both primers exactly - a database gap, not a primer failure |
| *Plasmodium* | 3 entries in PR2 |

**PR2 has labelling errors that become false calls.** A *Schistosoma*
entry carries a *Brugia malayi* V4 sequence, a *Haemonchus* entry is
identical to *Necator americanus*, and a *Brugia* entry is identical to
*Caenorhabditis*. Collapsing identical sequences turns each into an
ambiguity (or, if only one entry carries the sequence, a confident wrong
answer). Database preparation should therefore **flag identical
sequences whose labels disagree at a high rank** (different phyla or
classes cannot share a V4 sequence) and report them for curation, not
silently merge them.

**PR2 alone does not cover veterinary trichostrongylids.** For ruminant
work it would need supplementing, e.g. with curated 18S from NCBI for
*Teladorsagia*, *Cooperia*, *Oesophagostomum* and others. The database must
therefore accept local additions (§4).

**Anchor look-alikes.** 953 unique amplicons carry an adapter anchor
within 2 edits - the splitter's tolerance - including 32 platyhelminth
(mostly cestode) and 18 nematode amplicons; 2 are within 1 edit. Without
a database-wide mask every read of those taxa would be cut or
quarantined on nanopore. At 3 edits there are tens of thousands of
chance hits, which is the expected background for a 14-15bp pattern.

## 4. Which database: a layered build, not a single download

No single public database is adequate, so the reference is **built**,
with PR2 as the backbone and curated layers on top.

Alternatives considered:

- **SILVA NR99** is broader, but clustered at 99% identity (close
  species merged before we start) and its eukaryote taxonomy is
  inconsistent below class. Useful to fill gaps, not as the backbone.
  The DADA2 `silva_*_assignSpecies.fa.gz` files are bacteria and archaea
  only (checked: 1.47M entries, no eukaryotes) and are of no use here.
- **NCBI GenBank** is the most complete but uneven. 18S records of 1-3kb
  (Oct 2026): *Teladorsagia* 1, *Cooperia* 0, *Trichostrongylus* 3,
  *Haemonchus* 13. Allowing 400bp+ raises this (*Cooperia* 88), but many
  of those are likely 18S 3′ end plus ITS1, which does not cover V4.

The build:

1. **PR2 backbone**, trimmed to the amplicon (§2).
2. **Genome assemblies.** The same in-silico PCR on WormBase ParaSite
   helminth assemblies, and VEuPathDB for protists such as
   *Cryptosporidium*. Authoritative species labels, and genomes exist
   for exactly the taxa PR2 lacks. Risk: rDNA arrays are sometimes
   collapsed or missing in assemblies, so coverage is checked per
   species.
3. **Targeted GenBank additions** for named gap genera, kept only if
   they span the amplicon, with ranks taken from NCBI taxonomy IDs and
   mapped onto PR2's nine.
4. **Local sequences**: V4 from vouchered specimens and lab strains.
   The most trusted layer, and the one that grows.
5. **Curation across layers**: provenance per sequence, identical
   sequences collapsed, label conflicts at high rank flagged for review
   (§3), and a versioned build recorded in every result.

### 4.1 What the WBPS genomes give (WBPS19, all assemblies)

Every genomic FASTA in WormBase ParaSite release 19 was streamed and
scanned with the same primers: forward <= 1 edit, reverse <= 2 edits,
both strands; hits with N in a primer site or > 1% N in the insert
discarded, since N in a scaffold gap matches a primer for free. Scripts
are `bin/wbps_insilico_pcr.py` and `bin/compare_wbps_pr2.py`; every
hit is in `results/wbps19_v4_amplicons.tsv.gz`.

- 273 of 274 assemblies scanned (one download truncated at source).
  **185 (68%) yield at least one V4 amplicon**; 145 of 208 species.
- **42 species gain a V4 sequence PR2 does not have.** Most are
  free-living nematodes; the parasites among them include *Brugia
  pahangi*, *B. timori*, *Cylicocyclus nassatus*, *Cylicostephanus
  goldi*, *Parascaris univalens*, *Parelaphostrongylus tenuis*,
  *Paragonimus heterotremus*, *Echinococcus oligarthrus*, *Schistosoma
  guineensis* and several *Trichinella*.
- **The priority gaps mostly stay open.** No clean amplicon from
  *Teladorsagia* (the 18S is broken by scaffold gaps), *Oesophagostomum*,
  any of four *Ancylostoma* assemblies, any of five *Onchocerca*,
  *Taenia solium* or *Ascaris lumbricoides*. *Cooperia* has no WBPS
  genome. rDNA arrays are routinely collapsed or missing in older
  assemblies, as feared. For strongylids, GenBank and our own Sanger
  sequences (layers 3 and 4) are the route, not genomes.
- **Genomes adjudicate PR2's mislabels.** The *Necator americanus*
  genome's V4 is identical to PR2's "*Haemonchus* sp." entry, and the
  *Brugia malayi* and *B. timori* genomes are identical to PR2's
  "*Schistosoma mansoni*" entry - both PR2 labels are wrong.
- **Genomes carry other organisms' rDNA.** *Protopolystoma xenopodis*
  contains *Xenopus* (its host), *Wuchereria bancrofti* contains *Culex*
  (its vector), and several assemblies carry fungal or amoebal 18S.
  Seven *Trichinella* assemblies and *Meloidogyne hapla* each carry a
  553/545bp sequence identical to PR2 entries labelled as bacteria;
  at that length it is eukaryotic 18S, so those PR2 labels are probably
  wrong too.
- **Copies and variants.** Chromosome-level assemblies hold the array:
  *Fasciola hepatica* PRJEB58756 has 516 copies in 28 variants,
  *Schistosoma curassoni* 166 in 67. Draft assemblies add error: the
  older *F. hepatica* assembly's variants differ substantially from both
  PR2 and the newer assembly, which matches PR2 exactly.
- **Long V4 in cestodes.** *Taenia* inserts are 941-966bp and
  *Echinococcus* 828-889bp, well above the 546bp median.
- **Shared sequences are biology, not error,** in several places:
  *Haemonchus contortus*/*placei*/*similis*, *Trichinella* species, and
  the *Schistosoma haematobium* group (*curassoni*, *intercalatum*,
  *mattheei*, *guineensis*) are identical over V4.

### 4.2 Rules for the genome layer

1. Keep a variant only if its nearest database neighbour is in the same
   phylum as the genome's species; send the rest to review as possible
   contaminants or PR2 mislabels.
2. Prefer the best assembly per species, and variants seen in more than
   one copy, over every variant of every draft.
3. Where a genome and PR2 disagree on a label for the same sequence,
   the genome's label wins if the variant is the dominant one in a
   good assembly.

## 5. Release format

A release is one directory, named for its build version, holding:

| file | contents |
|---|---|
| `refdb.fasta` | one record per unique insert, primers trimmed; the header is the record id only |
| `taxonomy.tsv` | id, the nine PR2 ranks (merged to the ranks all members share; placeholder ranks left blank), the source layers contributing, and the number of member sequences |
| `members.tsv` | every source sequence collapsed into each id: accession or genome and contig, layer, original label |
| `resolution_nanopore.tsv` | per id, the deepest rank at which it is distinguishable over the whole insert |
| `resolution_illumina_<R1>x<R2>.tsv` | the same over the two read windows for that read length |
| `anchor_lookalikes.tsv` | adapter-anchor near-copies, in exactly the columns the pipeline's pre-flight writes (`seq_id anchor strand position edits context`), so the pipeline reads it unchanged |
| `conflicts.tsv` | identical sequences whose labels disagree at a high rank, and what the build did with each |
| `manifest.json` | build version, primers, PR2 version, WBPS release, version of every other layer, build date and code commit |

The pipeline records `manifest.json`'s build version in every result.

## 6. Order of work

1. Turn the exploratory scripts into a build: PR2 layer, collapse and
   merge, conflict flagging, release files (§5).
2. Genome layer from `results/wbps19_v4_amplicons.tsv.gz`, applying
   §4.2.
3. GenBank layer for a priority gap list; measure how many gaps it
   closes.
4. VEuPathDB genomes for protists (*Cryptosporidium* first).
5. Layer 4 format, and the first local sequences.

## 7. Open questions

- The priority gap list: ruminant trichostrongylids, human
  soil-transmitted helminths, or both. It decides which GenBank queries
  come first.
- File format for local additions (layer 4), and how each layer is
  versioned.
- Whether to filter to lineages of interest, or keep all eukaryotes so
  off-target reads are assigned rather than forced onto the nearest
  parasite. (Current lean: keep everything; the pipeline reports
  off-target fractions from it.)
- Retry the one WBPS19 assembly whose download failed
  (*Schistosoma intercalatum* TD2).
