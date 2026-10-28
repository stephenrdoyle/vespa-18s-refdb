# vespa-18s-refdb: design

**Status: proposal and first measurements. No release has been built yet.**

This repository builds the reference database that
[amphid](https://github.com/stephenrdoyle/nf-amphid-pipeline) uses in its
`--reference_database` mode: an 18S V4 database covering everything the
VESPA primers amplify, trimmed to the amplicon, with taxonomy, resolution
tables and the adapter look-alike mask the pipeline needs. How the
pipeline *uses* it - mapping, LCA, the Illumina and nanopore routes,
reporting - is designed in the pipeline's
`docs/REFERENCE_DATABASE.md`. The two meet at the release format in §6.

VESPA: Owens et al. 2024, *Nat Commun* 15:402,
doi:10.1038/s41467-023-44521-3.

---

## 1. Fixed inputs

| | |
|---|---|
| Primers | `AMPSEQ_VESPA18S_29F` `AGCAGCCGCGGTAATTCC` / `AMPSEQ_VESPA18S_21b8_I_R_Nmod` `TCCGTCAATTYCTTNAASTTTC` (locus-specific parts; both carry the standard Nextera tails). Our reverse primer differs from VESPA's published 21b8R (`TCAATTYCTTIAASTTTC`, Supplementary Table 1): a 5′ `TCCG` extension that matches the template, and N in place of deoxyinosine |
| Backbone | PR2 v5.1.1 SSU, `taxo_long` format (accession, gene, organelle, strain, then the nine ranks `Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species`). The exploratory measurements in §3 were made on v5.1.0 in DADA2 format, which has no accessions; the build uses v5.1.1, whose in-silico PCR gives 151,932 amplicons (one more than v5.1.0) and the same 89,738 unique inserts |
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

As implemented in `bin/wbps_layer.py`:

1. **Support.** Per species, pooled across its assemblies, a variant
   enters only if it is seen in at least two copies. A species with no
   such variant contributes its most common one. Of 1,034 variants
   seen, 798 single-copy ones are dropped - mostly draft-assembly noise,
   and with them most stray contaminant copies.
2. **Phylum.** A variant enters only if its nearest PR2 sequence is in
   the species' own PR2 class. Four failed, all correctly: the only V4
   in *Cylicostephanus goldi* is a chytrid, in *Protopolystoma* its host
   *Xenopus*; a *Trichinella pseudospiralis* assembly carries vertebrate
   18S in two copies; a *Meloidogyne graminicola* contig carries
   *Tetrahymena*. They go to `results/curation/wbps_review.tsv`.
3. **Ranks.** From PR2 where PR2 knows the genus (135 species);
   otherwise from the NCBI lineage, walked up to the first name PR2
   has, with PR2-style placeholders below (10 species; strain-coded
   names such as *Auanema* JU1783 are looked up by genus).

232 genome variants from 145 species enter the build.

### 4.3 GenBank layer

`bin/genbank_fetch.py` and `bin/genbank_layer.py`. 52,121 nuccore
records for Nematoda, Platyhelminthes and Acanthocephala (18S or "small
subunit ribosomal", 500-4000bp, not environmental or uncultured):

| | |
|---|---|
| already a PR2 entry (by accession) | 6,230 |
| no V4 (partial 18S, ITS) | 37,794 |
| sent to review as another class | 49 |
| **entered** | **~7,900** |

That roughly doubles helminth V4 coverage beyond PR2: Acanthocephala
alone goes from 107 PR2 amplicons to about 400. Ranks come from the NCBI
taxid. Two traps, both fixed:

- **Genus names repeat across kingdoms** - *Setaria* is a nematode and
  a grass, *Labrys* a nematode and a bacterium. A PR2 genus is used only
  if PR2's class for it appears in the record's NCBI lineage
  (`vespa.PR2Index.ranks_for`); "Metazoa" alone is not enough.
- **The first word of a name is not always a genus.** "Nematoda sp."
  and "Dicrocoeliidae sp." had become genera. The genus is now NCBI's,
  from the ranked lineage; a record with none gets an empty genus.

The review list is host and vector rDNA, as expected: trematode records
nearest molluscs (intermediate hosts), *Wuchereria* nearest mosquitoes.

**No public V4 exists for *Teladorsagia*, *Cooperia* or *Bunostomum*.**
Their GenBank records are ITS or 18S fragments outside V4.
*Onchocerca* has about 30 partial records covering V4's 3′ end, not
used for now.

### 4.4 ParaRef

ParaRef (Niemann et al. 2025; zenodo.org/records/13744644) is 645
decontaminated parasite assemblies: 221 fungi, 184 nematodes, 65
apicomplexans, 61 platyhelminths, 37 kinetoplastids and more. The
21.8GB archive was streamed through the same scan
(`bin/pararef_insilico_pcr.py`, 77 minutes; 644 genomes - the
archive lacks *Elaeophora elaphi*). Hits: `results/pararef_v4_amplicons.tsv.gz`.

**ParaRef's decontamination strips rDNA.** Of 174 species in both
ParaRef and WBPS, 82 have V4 in WBPS and none in ParaRef - every
*Schistosoma*, *Trichuris*, *Trichinella* and *Echinococcus*, and
*C. elegans* itself; *Fasciola hepatica* falls from 519 copies to 1.
Only 4 species gain. For the gap species the assemblies are the same
(same BioProject), so this is masking, not different assemblies.
Conterminator flags sequence shared across distant taxa, and 18S is
exactly that. So a missing V4 in a ParaRef genome means nothing, and
ParaRef is used only as a lower-trust genome layer. (It matters beyond
this project: metagenomic reads from rDNA will find nothing to map to
in ParaRef.)

Even so, 63 ParaRef variants are in no other source - helminths such as
*Dipylidium caninum*, *Taenia pisiformis*, *Mansonella ozzardi*,
*Ascaris lumbricoides*, *Parascaris univalens*, *Heligmosomoides
polygyrus*; protists such as *Plasmodium vivax* and *P. malariae*,
*Giardia muris*, *Entamoeba dispar*, *Naegleria fowleri*,
*Histomonas*, *Trypanosoma congolense*. Still nothing for
*Cryptosporidium*.

### 4.5 Local layer and *Teladorsagia*

`local/sequences.tsv` is hand-curated and the most trusted layer
(`bin/local_layer.py`). Its first entry is the *Teladorsagia
circumcincta* V4 from the tci2_wsi3.0 assembly - chr_1, two identical
copies, both primers exact - found with `bin/genome_insilico_pcr.py`.

It is **identical to PR2's *Trichostrongylus colubriformis*
(AJ920350)**, one base from *Ostertagia*, and two from all three
*Haemonchus* (one difference in each 2x150 read window). One
*Ostertagia ostertagi* record differs only in the unread middle, and two
*O. ostertagi* records differ from each other. So **V4 resolves the
ruminant trichostrongylids to family and seldom further**; telling
*Teladorsagia* from *Trichostrongylus* or *Ostertagia* needs ITS2. The
merged record stops at the shared rank.

### 4.6 Protists: *Cryptosporidium*

**PR2 splits apicomplexan genera into numbered clades** -
*Cryptosporidium1* (*C. parvum*, *C. hominis*), *Cryptosporidium13*
(*C. andersoni*), and likewise *Theileria*, *Babesia*, *Hepatozoon*,
*Eimeria*, *Sarcocystis*, *Plasmodium*, *Toxoplasma*: about 9,000 PR2
entries. They are now read as their base genus (`vespa.normalise_ranks`)
in the merge and in rank lookup; `members.tsv` keeps PR2's own label.
Before this, a genome's *Cryptosporidium* and PR2's *Cryptosporidium1*
were two genera and their shared records stopped at family - and the
v0.2.0 coverage table (§6.2) undercounts every one of these genera: its
"*Cryptosporidium* 1" and "*Toxoplasma* 0" are wrong.

**PR2 and NCBI name protist classes differently** (Coccidiomorphea /
Gregarinomorphea against Conoidasida), so the lineage check rejected
PR2's ranks for every *Cryptosporidium* genome - and the placeholder
class that replaced them also switched off the other-class check. The
sequence is now the witness: a PR2 genus is accepted when its class is
the class of the PR2 insert nearest the taxon's own V4.

**Sources.** Most PR2 and GenBank *Cryptosporidium* records are the
standard nested-PCR genotyping fragment: it holds VESPA's forward primer
and runs about 436 of the 513 insert bases, so it never spans the
amplicon. Full-length V4 comes from genomes:
`bin/ncbi_genomes_insilico_pcr.py` streamed all 102 NCBI assemblies (17
species; ~2 minutes), and the genome-layer rules
(`bin/pararef_layer.py ... NCBI_genomes`) keep 27 variants. GenBank adds
57 full-length records PR2 lacks (`bin/genbank_layer.py`), mostly
*C. parvum* and *C. baileyi*; 7,115 are genotyping fragments.

With PR2, all 86 *Cryptosporidium* records reach genus, and the named
species - *parvum*, *hominis*, *meleagridis*, *muris*, *andersoni*,
*bovis*, *ryanae*, *xiaoi*, *baileyi*, *serpentis*, *felis*, *canis*,
*ubiquitum*, *viatorum*, *proliferans*, *fragile* - mostly reach
species, whole or through the Illumina windows.

Open: the genotyping fragments carry far more species and genotypes
than the genomes. Admitting partial references (here: all of the R1
window, a third of the R2 window, 85% of the insert) is the same
question as *Onchocerca*'s partial records (§7).

## 5. Curation: collapsing the layers

As implemented in `bin/collapse.py` (its docstring is the reference).
Layers are given in increasing order of trust: **GenBank, PR2,
ParaRef, WBPS, local**. GenBank is below PR2 because PR2 *is* GenBank,
curated - a raw GenBank label must not overrule PR2's. ParaRef is below
WBPS because its decontaminated genomes hold less rDNA (§4.4).
Inserts shorter than 100bp are set aside: primer matches that touch or
overlap "amplify" nothing (33 GenBank records, 6 PR2 entries).

- **Merge.** Identical inserts become one record with a stable id
  derived from the sequence. At each rank the record takes the single
  informative name its members agree on; placeholders do not vote; two
  different names stop it at the rank above.
- **Conflicts at class or above** (identical V4 in two PR2 classes - a
  nematode and a flatworm, a plant and a bacterium) are resolved by the
  most trusted layer if it is consistent, else by a two-thirds
  majority of at least two, else left **unresolved**: the record is kept
  at the ranks all members share and listed for review. Kept, not
  dropped - a sequence missing from the database still sends reads,
  which would then be pushed onto its nearest neighbour.
- **Below class, a trusted layer overrules narrowly.** Species never
  overrule: a genome shows a sequence is *in* a species, not *only* in
  it. At order, family and genus, a conflicting less-trusted label is
  overruled only when a single accession carries it **and** it either
  was never identified to species (PR2's "*Haemonchus* sp." identical
  to *Necator*) or is outvoted at least three to one (a *Brugia*
  genome-project contig among nine *C. elegans* records). A named
  species standing alone is believed: the *Teladorsagia* genome's V4 is
  identical to PR2's single *Trichostrongylus colubriformis* (§4.5), and
  erasing that record on one genome's word would be wrong.
- Overruled members stay in `members.tsv`, marked.

On PR2 v5.1.1 plus WBPS19: 89,872 unique inserts; 103 class-level
conflicts resolved by majority, 8 by a genome, 39 unresolved. The
helminth cases come out as they should:

| record | outcome |
|---|---|
| *Brugia malayi* V4, one PR2 member labelled *Schistosoma mansoni* (a contig from the *S. mansoni* genome project) | majority: *Brugia*; *Wuchereria* and *B. timori* genomes share it, so the record stops above genus |
| *Necator americanus* V4, one PR2 member labelled "*Haemonchus* sp." | genome overrules the lone label: *Necator americanus* |
| *C. elegans* V4, one PR2 member from the *B. malayi* genome project | majority plus genome: *C. elegans* |
| *Tetrahymena* V4, one PR2 member from a *M. graminicola* genome contig | majority: *Tetrahymena* |
| *Haemonchus contortus*/*placei*/*similis* | genus *Haemonchus*: real sharing |
| *Toxocara canis*/*cati*, one PR2 "*Porrocaecum*" | lone label overruled; genus *Toxocara* |
| *Pristionchus*, one PR2 "*Diplogaster lethieri*" | an old synonym of *P. lheritieri*; genus *Pristionchus* |
| *Clonorchis*/*Opisthorchis* | several accessions each: real sharing; family |

Every decision is in `results/curation/`: `conflicts.tsv`,
`overruled.tsv` and `wbps_review.tsv`.

## 6. Release format

A release is one directory, named for its build version, holding:

| file | contents |
|---|---|
| `refdb.fasta` | one record per unique insert, primers trimmed; the header is the record id only |
| `taxonomy.tsv` | id, the nine PR2 ranks (merged to the ranks all members share; placeholder ranks left blank), the source layers contributing, and the number of member sequences |
| `members.tsv` | every source sequence collapsed into each id: accession or genome and contig, layer, original label, role (`kept`, `overruled`, `dropped`) and, for dropped members, the reason (`note`, from v0.5.0) |
| `resolution_nanopore.tsv` | per id, the deepest rank at which it is distinguishable over the whole insert |
| `resolution_illumina_<R1>x<R2>.tsv` | the same over the two read windows for that read length; one per geometry, listed in the manifest's `illumina_geometries` (2x150 and 2x300 from v0.5.0) |
| `anchor_lookalikes.tsv` | adapter-anchor near-copies, in exactly the columns the pipeline's pre-flight writes (`seq_id anchor strand position edits context`), so the pipeline reads it unchanged |
| `conflicts.tsv` | identical sequences whose labels disagree at a high rank, and what the build did with each |
| `primer_fit.tsv` | per id, how well the primers fit its best-matching template: forward and reverse edits, mismatch positions counted from each primer's 3′ end, the template at each site, and a class - `exact`, `away_3prime`, `near_3prime` (a mismatch in the last 5 bases), `indel`, `edits_only` (a genome-only record: counts known, sites not), `unknown` (from v0.4.0) |
| `manifest.json` | build version, primers, PR2 version, WBPS release, version of every other layer, build date and code commit |

The pipeline records `manifest.json`'s build version in every result.

### 6.1 v0.1.0

PR2 v5.1.1 plus WBPS19, built by commit `930d392`;
`releases/vespa18s_v0.1.0/summary.txt` has the full counts and
`manifest.json` the checksums. The data (19.4MB tarball) are not in git.

- 89,872 records; 69,419 distinct Illumina 150x150 window groups.
- Deepest rank reached, all records: species 49% / genus 25% read
  whole (nanopore); species 46% / genus 25% through the Illumina
  windows.
- Helminth records (4,678): species 67% / genus 27% nanopore; species
  63% / genus 30% Illumina. (Per record, not per PR2 entry, so lower
  than the per-entry figures in §3.)
- Checked against this format: checksums, every id in every table,
  look-alike columns, and the pipeline's own `LookalikeMask` loads the
  mask.

**The mask is valid but the pipeline cannot use it yet.** 43,213
look-alikes (42,244 of them at 3 edits, the chance background for a
14-15bp anchor) load as 29,332 contexts, and the pipeline's
`LookalikeMask.covers()` tries every context against every anchor hit:
measured at 3 seconds per hit, about 350 CPU-days for a 5M-read run.
That is the pipeline's to fix - index the contexts, for example by
k-mer, so a hit is checked only against contexts that share sequence
with it. The release format stays as it is.

**Update:** fixed on the pipeline's `lookalike-mask-index` branch
(commit `e9c5245`, not yet merged): a seed-piece index plus the q-gram
lemma, both exact bounds, so the answers are unchanged. On the v0.4.0
mask a real junction hit takes 7.8 ms instead of 389 ms - about 22
CPU-hours for a 5M-read run instead of about 740.

**Erratum (found while cutting v0.2.0).** v0.1.0's
`resolution_illumina_150x150.tsv` overstates the rank reached for
**8,775 of 89,872 records (227 helminth)**. Records sharing read windows
were merged with placeholders abstaining, which is right for the labels
of one sequence but wrong for different sequences a read cannot tell
apart. `bin/release.py` now takes a strict common ancestor and refuses
to release if any record resolves deeper through the windows than
whole. v0.1.0 is left as tagged; use v0.2.0.

### 6.2 v0.2.0

GenBank + PR2 v5.1.1 + ParaRef + WBPS19 + local, in that order of trust;
built by commit `af46077`. `releases/vespa18s_v0.2.0/summary.txt` and
`manifest.json` have the details; data in the tarball.

- **94,028 records** (v0.1.0: 89,872); 72,727 Illumina window groups.
- **Helminth records 8,808** (v0.1.0: 4,678), almost all the gain from
  GenBank. Species reached: 64% nanopore, 55% Illumina; genus or better
  94% / 89%.
- Conflicts: 103 resolved by majority, 20 by trust, 39 unresolved.
- Format checks pass, and no record resolves deeper through the Illumina
  windows than whole.

Genera of interest (`results/v0.2.0_genera_of_interest.tsv`; records,
then how many reach species, nanopore/Illumina):

| genus | records | species N/I | notes |
|---|---:|---|---|
| *Haemonchus* | 4 | 2/1 | three species share one V4 |
| *Teladorsagia* | 1 | 0/0 | local only; shares its record with *Trichostrongylus* |
| *Ostertagia* | 4 | 3/2 | |
| *Cooperia*, *Bunostomum* | 0 | | no V4 in any public source |
| *Oesophagostomum*, *Nematodirus*, *Chabertia* | 1 each | 1/1 | one record each |
| *Ancylostoma* / *Necator* | 4 / 1 | 2/1, 1/1 | |
| *Strongyloides* | 25 | 18/16 | |
| *Trichuris* | 29 | 23/21 | |
| *Ascaris* / *Toxocara* | 4 / 6 | 3/1, 4/2 | |
| *Onchocerca* | 1 | 0/0 | 30 partial GenBank records unused |
| *Brugia* / *Wuchereria* / *Mansonella* / *Dirofilaria* | 8 / 3 / 2 / 7 | 6/4, 2/2, 2/2, 7/7 | |
| *Schistosoma* | 87 | 53/31 | the largest Illumina loss |
| *Fasciola* / *Taenia* / *Echinococcus* | 10 / 42 / 18 | 8/4, 40/33, 18/15 | |
| *Moniezia* / *Dipylidium* / *Hymenolepis* | 10 / 3 / 11 | 8/8, 3/3, 10/8 | |
| *Giardia* / *Entamoeba* / *Blastocystis* | 20 / 146 / 294 | 15/14, 98/94, 89/55 | |
| *Cryptosporidium* | 1 | 0/0 | still effectively absent |
| *Toxoplasma* | 0 | | no amplicon-spanning entry |

### 6.3 v0.3.0

Seven layers, in increasing trust: GenBank helminths and GenBank
protists, PR2 v5.1.1, ParaRef, NCBI protist genomes, WBPS19, local.
Built by commit `dcf9249`.

- **97,559 records**; helminth records unchanged at 8,808.
- **Protists**: the NCBI genome layer (20 genera, three best assemblies
  per species, 260 variants) and the GenBank protist layer (60 genera
  of gut, blood, tissue and zoo-relevant protists and protist-like fungi;
  8,126 records PR2 lacked, collapsing to about 3,500 new sequences,
  because many *Theileria*, *Babesia* and *Sarcocystis* records are
  identical).
- PR2's numbered clade genera (§4.6) are now counted with their base
  genus, which corrects the v0.2.0 table: *Cryptosporidium* has 102
  records, not 1; *Toxoplasma* 48, not 0.
- 9 GenBank protist records went to review: apicoplast 16S labelled as
  *Neospora*, and fungal or plant sequences under parasite names.

The full table is `results/v0.3.0_genera_of_interest.tsv`. Selected
rows (records; reaching species whole / through the Illumina windows):

| genus | records | species N/I |
|---|---:|---|
| *Cryptosporidium* | 102 | 82/70 |
| *Giardia* | 40 | 32/31 |
| *Entamoeba* | 217 | 159/157 |
| *Blastocystis* | 638 | 160/113 (subtypes, mostly not named species) |
| *Balantioides* | 18 | 14/14 |
| *Eimeria* | 534 | 401/334 |
| *Toxoplasma* | 48 | 41/30 |
| *Sarcocystis* | 801 | 651/511 |
| *Plasmodium* | 447 | 354/337 |
| *Babesia* / *Theileria* | 778 / 945 | 498/432, 684/504 |
| *Leucocytozoon* / *Haemoproteus* | 145 / 7 | 40/40, 7/7 |
| *Trypanosoma* / *Leishmania* | 565 / 20 | 256/181, 14/7 |

Still absent: *Cooperia*, *Bunostomum*, *Hammondia* (probably identical
to *Toxoplasma* over V4, and in PR2 if at all), and among zoo taxa
*Murshidia*, *Prosthenorchis*, *Rhabdias*, *Anoplocephala* - none of
which has 18S in GenBank.

### 6.4 v0.4.0

The same 97,559 records as v0.3.0 (`refdb.fasta` and `taxonomy.tsv`
are byte-identical), plus `primer_fit.tsv`. The in-silico PCR now keeps
the template at each primer site (`vespa.pcr_full`), layers carry them
(`f_site`, `r_site`), and the genome layers keep the edit counts their
scans recorded. Built by commit `9110efd`.

Every record is an *in-silico* target - both primer sites present within
two edits each - but not every one is equally likely to amplify:

| | exact | mismatch, not near 3′ | mismatch in last 5 bases | indel / sites unknown |
|---|---:|---:|---:|---:|
| all records | 78% | 8% | 14% | 0.4% |
| Nematoda | 91% | 3% | 6% | <1% |
| Platyhelminthes | 32% | 1% | **65%** | 2% |
| Acanthocephala | 87% | 3% | 10% | 0% |

- **Flatworms carry a systematic forward-primer mismatch.** Almost all
  the near-3′ cases (2,019 records) are one base, 4 from 29F's 3′ end:
  template `...GGTAACTCC` against primer `...GGTAATTCC`. In the duplex it
  is a T·G pair, the mildest mismatch and generally tolerated, so
  flatworms should amplify - but possibly less efficiently than
  nematodes in the same tube, so under-represented. Untested: a mock
  with a flatworm carrying this variant (most *Schistosoma* and
  *Fasciola*) beside a nematode would settle it.
- **187 nematode records mismatch the reverse primer's last base** - the
  mismatch most likely to block extension. Those organisms may be
  detected poorly or not at all.
- The pipeline can use `primer_fit.tsv` to say, beside a call or a
  non-detection, that amplification was predicted to be weak.

### 6.5 v0.5.0

Curation from the first use of a release in the pipeline
(`docs/feedback/v0.4.0_feedback.md`, response in
`docs/feedback/v0.4.0_response.md`), the *D. immitis* WSI 2.2 genome's
three V4 copies as local sequences, and a 2x300 resolution table.

- 97,537 records: 24 removed (eukaryote 18S carrying only a bacterial
  label), 2 added (*D. immitis*), 21 relabelled.
- New rules in `collapse.py`: prokaryote label on a eukaryote-length
  insert dropped; organelle label beside a nuclear one dropped; directly
  sequenced isolate over WGS/TSA contig; host-named species to Genus_sp.;
  and `curation/overrides.tsv` for judgements no rule can make.
- `members.tsv` gains a `note` column; dropped members are listed with
  role `dropped` and the reason.
- `resolution_illumina_300x300.tsv` (250/200-base windows) beside
  `resolution_illumina_150x150.tsv`; `manifest.json` lists them under
  `illumina_geometries`. Species-level records: 38,479 at 2x150, 46,566
  at 2x300, 48,424 read whole.

### 6.6 v0.6.0

Long-read helminth genomes, mostly from the Darwin Tree of Life, scanned
for V4 and built into the NCBI genome layer under the genome rules (two
copies, own class). Built by commit `da7d3a4`.

- 97,708 records (171 added, none removed, 16 relabelled); helminth
  records 8,810 -> 8,981.
- `curation/longread_genomes.txt`: 16 assemblies for the trichostrongylid
  and other gaps. *Cooperia oncophora*, *Oesophagostomum venulosum*,
  *Nematodirus battus*, *Anoplocephala perfoliata* and *Necator
  americanus* now resolve to species read whole and through both
  Illumina geometries; *Trichostrongylus axei* mostly.
- `curation/tol_genomes.txt`: 65 more - 13 vertebrate parasites (new
  genera *Heligmosomum*, *Litomosoides* and an Anisakidae; *Trichuris
  ovis*, *Opisthorchis viverrini*, *Dictyocaulus viviparus*,
  *Spirocerca*, *Thelazia*, fish tapeworms) and 52 free-living or
  invertebrate-parasite species, as references for faecal colonisers.
  No free-living record shares a sequence with a parasite.
- Independent genomes settle two limits of V4: *Teladorsagia*,
  *Ostertagia ostertagi* and *Trichostrongylus colubriformis* share one
  sequence (family), and so do *Ancylostoma caninum* and *A. duodenale*
  (genus). *Trichuris ovis* is identical to *T. discolor*.
- The *D. immitis* Tree of Life assembly carries the V4 of AB973230/1,
  the record VESPA's human reads matched - a genuine *D. immitis*.
- Two *Sphaerolaimus hirsutus* variants sent to review
  (`results/curation/tol_genome_review.tsv`): weakly nearest a fungus,
  more likely a divergent nematode V4 than a contaminant.

Still no V4 anywhere: *Bunostomum*, *Chabertia*, *Murshidia*,
*Prosthenorchis*, *Rhabdias* - none has an assembly on NCBI.

## 7. Order of work

1. ~~PR2 layer, genome layer, collapse with conflict handling~~ (done:
   §4.2, §5).
2. ~~The rest of a release (§6)~~ (done: `bin/release.py`; first
   release **v0.1.0**, §6.1).
3. ~~GenBank layer~~ (helminths done, §4.3); ~~ParaRef~~ (§4.4);
   ~~local layer~~ (§4.5); **v0.2.0** (§6.2).
4. ~~Protists~~ (done: NCBI genomes and GenBank, §4.6, §6.3).
5. Local sequences for *Cooperia*, *Bunostomum* and other trichostrongylid
   gaps (Sanger, or our own assemblies via `bin/genome_insilico_pcr.py`).
6. Decide whether partial GenBank records (V4's 3′ end only, as for
   *Onchocerca*) can enter as partial references.
7. ~~Make the GenBank other-class check use all cores~~ (done).
8. Zoo gaps with no 18S anywhere (elephant and rhino strongyles,
   *Prosthenorchis*, *Rhabdias*, anoplocephalids): genomes or Sanger.

## 8. Open questions

- The priority gap list: ruminant trichostrongylids, human
  soil-transmitted helminths, or both. It decides which GenBank queries
  come first.
- File format for local additions (layer 4), and how each layer is
  versioned.
- Whether to filter to lineages of interest, or keep all eukaryotes so
  off-target reads are assigned rather than forced onto the nearest
  parasite. (Current lean: keep everything; the pipeline reports
  off-target fractions from it.)
