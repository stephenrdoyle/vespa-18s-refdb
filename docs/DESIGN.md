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

### 4.7 *Blastocystis* subtypes

*Blastocystis* is reported by subtype (ST1, ST2 ...), which are not
species, and PR2 and GenBank name few: v0.6.0's 639 *Blastocystis*
records were "*Blastocystis* sp." (463) or the host-based
"*B. hominis*" (153), so even an exact match named no subtype.

- **References.** `curation/blastocystis_subtypes.tsv`: published
  full-length references for every subtype, ST1-ST50, with source and
  DOI - the backbone is Santin et al. 2024 (*J Eukaryot Microbiol*,
  Table 2), plus each describing paper. 45 subtypes are accepted (ST1-17,
  ST21, ST23-49); ST18-20 and ST22 (Zhao 2017, possibly chimeric) and an
  unpublished dugong "ST45" are invalid; ST50 (unpublished, 2025) is
  disputed. ST10 was split in 2024 (old ST10c-f are ST42-44), so older
  "ST10" titles cannot be trusted. Checked against Gentekaki et al. 2026
  (*Trends Parasitol* 42:322, doi:10.1016/j.pt.2026.03.006): its Table 1
  lists ST27-ST49 with the same representative accessions as ours and
  invalidates none (it says so). It argues against the 2% criterion that
  made ST42-ST44 (from ST10) and ST49, and would treat such finer
  structure as allelic groups within established subtypes. Those four
  stay named here - V4 separates all of them - with the dispute noted in
  the table; folding them into their parent subtypes is a change to
  `blastocystis_v4_subtypes.tsv` alone. Its criterion (ii) - fragments
  may assign isolates to known subtypes but not name new ones - is what
  this section does with V4.
- **Assignment.** `bin/blastocystis_subtypes.py` trims the references
  to V4 (82 inserts) and compares every *Blastocystis* insert we hold
  with each, by exact edit distance. An insert takes the nearest
  subtype when it is within 8% of the insert's length and every other
  subtype is at least 10 edits further. V4 varies a lot inside some
  subtypes - GenBank records titled ST3 run to 7.2% from the ST3
  references - but subtypes stay apart (those ST3 records are >= 11%
  from any other), hence a wide distance with a firm margin. Against
  GenBank titles in the current numbering it contradicts none; the
  titles it overrules are old ST10/ST12/ST14 numbers since renamed
  (ST42/44, ST26, ST24/25). Output: `curation/blastocystis_v4_subtypes.tsv`
  (record id -> subtype), which `bin/collapse.py` applies - every
  *Blastocystis* member of the record takes `Blastocystis_STn`; members
  of other genera keep their labels, so the conflict rules still judge.
  On records with no subtype, "*B. hominis*" becomes "*Blastocystis* sp.".
- **The references as records.** Some never reach a GenBank title
  search, so they are a layer of their own (`genbank_layer.py` on the
  script's `references_genbank.fasta.gz`): every subtype has at least
  its reference.
- **Result** (trial build on the v0.6.0 layers): 653 *Blastocystis*
  records, 581 named to subtype, all 45 subtypes present, and every
  subtyped record keeps its subtype read whole and through both Illumina
  geometries. 72 stay at genus: the reptile, amphibian and insect
  lineages (*B. lapemi*, *B. cycluri* ...) and possible undescribed
  subtypes, > 8% from every reference, and a few between close subtypes
  (ST14/24/25, ST10/44).
- **PubMLST and all the genomes** (after v0.7.0). PubMLST's curated
  `18S_rRNA_full_length` alleles (125, each typed; 87 span V4; public
  API, data to 2024-12-31) agree with every subtype this section had
  given (82 of 82) and are now extra references
  (`curation/blastocystis_pubmlst_18S.tsv`, `--extra`): they settle the
  ST14/ST24 cases left between the two (10 more ST24), 2 ST17 and 1
  ST1; none of the earlier assignments moved. The genome scan had
  capped assemblies per species, and most are filed as "*Blastocystis*
  sp.", so 9 of the 100 on NCBI had been scanned; all are now. 27 carry
  V4; the two-copy rule keeps 19 variants (was 6), 5 of them new - ST1
  variants and the RusS/HerGr lineage, 15% from every subtype. The
  subtype script now also types the genome layer's rows. Trial build:
  658 *Blastocystis* records, 596 named to subtype (was 653 and 581).
  They add little for VESPA: of ~100,000 reads at 85-99% to a
  *Blastocystis* record, 310 come closer.
- **VESPA's candidates.** The 26 *Blastocystis* candidates from the
  v0.4.0 feedback, compared over their two read windows, are all clearly
  one subtype - ST1 (10), ST2 (9), ST3 (7) - each 0.7-3.6% from a
  subtyped record and >= 13 edits nearer it than any other subtype. They
  were unnamed subtype variants, not missing taxa.

### 4.8 Genera found missing against the VESPA samples

The VESPA samples' unmatched clusters (v0.7.0, 2x301, identity < 97%)
pointed at protist genera never fetched from GenBank: gut trichomonads
and diplomonads, archamoebae and trichostome ciliates. Fetched
(`work/genbank_targeted_layer.tsv.gz`, a GenBank layer like the protist
one): 724 records in 15 genera, 451 already in PR2, 115 not spanning V4,
157 entered - **61 V4 sequences new to the release**: *Trichomitus
batrachorum* (26), *Tetratrichomonas* (11, among them *T. buttreyi* and
*T. gallinarum*), *Hypotrichomonas* (3), *Trepomonas* (3), *Trimitus*
(2), *Mastigamoeba* (4), *Rhizomastix libera* (3), *Buxtonella
sulcata* (5), *Cycloposthium* (3), *Dasytricha* (1).

They do not close VESPA's gaps: no unmatched cluster comes more than 1%
closer (the 83% "*Hexamita*", 87% "*Rhizomastix*" and 88%
"*Lacusteria*" clusters stay where they were), so those organisms are
not in GenBank under these names with V4-spanning records. Nematoda and
Platyhelminthes were already fetched whole; *Molineus*, a candidate for
the "*Viannaia*" cluster, has only a 36bp partial over V4.

NCBI BLAST of 15 representative unmatched clusters (largest per gap,
`results/vespa_gaps/`; core_nt, unfiltered, 2026-10-29):
- **recovered** (added to the same layer, `VESPA_blast_hits.fasta.gz`):
  a primate trichostome 98.0% to "*Buxtonella*-like sp." OM423619 (was
  93% to *Neobalantidium*); *Paramicrosporidium* sp. MW898450 (90 ->
  93%); "Hexamitinae sp." OR291217 (83 -> 93%).
- **to review** (the layer's class check sent them there): two
  "uncultured fungus" clones, LC365266/61, 97.1% to a 5,148-read
  cluster PR2 places nearest a colpodellid - NCBI says fungus, the
  sequence's nearest PR2 class is Coccidiomorphea; and "Parabasalia sp."
  OK584312 (88 -> 91%), held only because its NCBI class is a
  placeholder.
- **nothing closer anywhere:** *Entamoeba* and *Iodamoeba* (best hits
  are records we hold), an unknown alveolate (82% at best).
- **not gaps:** the "*Viannaia*" clusters are *Pichia* yeast and
  *Blastocystis* fragments, and the low-identity "*Blastocystis*",
  "*Theileria*", "*Endolimax*" and "*Rhizomastix*" ones yeast or one-window
  artefacts - labels like "unresolved below domain (94% to X)" come from
  windows that disagree, not from a missing organism.

### 4.9 A watch list of important helminths

`bin/watchlist_audit.py` reports, for the 122 helminths in
`curation/helminth_watchlist.tsv` (livestock, horse, human, companion
animal, zoo and wildlife), whether each is present, how deep it resolves
(whole, 2x300, 2x150) and its species margin. Rerun it when genomes
arrive. `results/audits/`: v0.8.0 and v0.9.0.

| | species | only in shared records | missing |
|---|---:|---:|---:|
| v0.8.0 | 80 | 10 | 32 |
| v0.9.0 | 89 | 11 | 22 |

What the gaps are:

- **Shared V4, not missing:** *Baylisascaris* (order), *Teladorsagia*
  and *T. colubriformis* (family), the rumen flukes (15 members from
  about 10 paramphistome genera, *Calicophoron* among them), the
  cyathostomins (§4.10), *Ancylostoma ceylanicum* / *caninum* /
  *duodenale*, *Mansonella perstans* / *ozzardi*, *Onchocerca volvulus* /
  *ochengi* / *lupi* / *cervicalis*. A call there names the genus or
  above, which is what V4 allows.
- **Fragile, 1 edit:** *H. contortus*, *Ostertagia*, *Chabertia*,
  *Cylicocyclus*, *A. lumbricoides*, *A. duodenale*, *Loa*,
  *Dicrocoelium*, *O. stephanostomum*, and now *O. volvulus* and *D.
  repens* (§4.10).
- **Still missing (22), no assembly and no genome reads in ENA:**
  *Cooperia punctata* / *curticei*, *Oesophagostomum dentatum* /
  *radiatum* / *columbianum* / *bifurcum*, *Nematodirus filicollis* /
  *spathiger*, *Trichostrongylus vitrinus*, *Bunostomum phlebotomum*,
  *Strongylus edentatus*, *Anoplocephala magna*, *Strongyloides
  westeri*, *A. braziliense* / *tubaeforme*, *Capillaria aerophila* /
  *obsignata*, *Ternidens*, *Pterygodermatites*, *Molineus*,
  *Enterobius anthropopitheci*, *Prosthenorchis*. Several have partial
  18S only (§8).

### 4.10 V4 from raw reads, where assemblies have none

Most 2010s short-read helminth genomes dropped their rDNA (an array of
identical repeats is what short reads cannot assemble): the WBPS genomes
of *A. ceylanicum*, *A. duodenale*, *O. dentatum*, *S. vulgaris*,
*H. placei*, *P. equorum* and every *Onchocerca* yield no V4. The raw
reads still hold it, at tens to hundreds of copies per genome, and ENA
holds the reads - including Darwin Tree of Life projects not yet
assembled. Three routes, all streaming only the start of a run:

- **Long accurate reads (PacBio HiFi):** `bin/reads_insilico_pcr.py`
  runs the usual scan on each read (100,000 reads give hundreds to
  thousands of V4 copies). `bin/read_variants.py` keeps a variant at
  >= 2% of a run's V4 reads and >= 5 reads (HiFi errors recur at
  homopolymers) and writes genome-layer input, so the genome rules and
  class check build `work/reads_layer.tsv.gz` (source `ENA_reads`;
  runs in `curation/longread_runs.tsv`). The class check sent the two
  contaminants to review: a *Bullera* yeast with *Bunostomum*, an
  unnamed *Entamoeba*-like amoeba with *Coronocyclus*.
- **Short reads, or noisy long reads (old PacBio subreads):**
  `bin/reads_v4_consensus.py` maps reads to a relative's 18S around V4
  (primers plus 250bp flanks), takes the majority base position by
  position and repeats until stable; `--decoy` drops reads nearer host
  18S (blood, skin nodules), `--preset map-pb` takes noisy long reads.
  Checked on *A. duodenale*: it recovers the PacBio-confirmed V4
  exactly, and corrects a seed with 8 planted errors in one round; from
  *O. cervicalis* as seed it recovers *D. repens* from RS II subreads.
  It gives the majority copy only: minor alleles from short reads were
  noise (absent from long reads of the same species), except where two
  independent isolates agree (*H. placei*, *O. volvulus* below).
  Results go in as local sequences, evidence in each row;
  `results/reads/` holds every run's report.
- **Assemblies:** further long-read assemblies scanned as before
  (`curation/longread_genomes.txt`), and two pre-release Tree of Life
  assemblies as their own layer (`work/tol_prerelease_layer.tsv.gz`,
  source `ToL_prerelease`: *Cyathostoma americana* nxCyaAmer1.1,
  *Crenosoma striatum* nxCreStri23; two *Angiostrongylus* ones held
  back).

What it found:

- **New species:** *Bunostomum trigonocephalum*, *Cyathostoma
  americana*, *Strongylus vulgaris*, and - by minor variants of their
  own - *Cylicostephanus goldi* / *longibursatus*, *Cylicocyclus
  insigne*, *Coronocyclus labiatus*, *Calicophoron daubneyi*.
- **Cyathostomins share one main V4.** Six species (*Cyathostomum*,
  *Cylicostephanus* x2, *Cylicocyclus* x2, *Coronocyclus*) carry
  `V4_c474d5b76665`, which v0.8.0 named *Cylicocyclus* - a genus call
  that would have misnamed every other cyathostomin. It now stops at
  Strongylida.
- **Filarial worms sit within a few bases of each other** (*Onchocerca*
  / *Loa* 1 edit; *A. viteae* / *D. repens* 1-2; *Loa* / *D. repens* 3).
  Confirmed from independent data: *Loa* (long-read assembly
  GCA_053739935.1), *D. repens* (RS II subreads), *A. viteae*
  (GCA_046563165.1), *Onchocerca* (*O. lupi* HiFi; *O. volvulus* from
  Cameroon, Ghana and Mali reads; *O. ochengi*, identical to *O.
  volvulus*). The base separating *Onchocerca* from *Loa* (insert 347)
  is real, but lies between the 2x150 read windows, so at 2x150 the
  three genera cannot be told apart. *O. volvulus* is polymorphic at
  insert 123 (Ghana: A 64%, C 36%); *D. repens* has A there, so the A
  variant went in as an *O. volvulus* record - without it those reads
  tied with *D. repens*.
- ***Mansonella*:** four read sets (two *M. perstans*, two *M.
  ozzardi*) agree on one shared main V4, so a call stops at genus. The
  v0.8.0 *M. ozzardi* record was a divergent minor copy (11% from every
  filarid, about 2% of reads) - kept, as what it is; the *M. perstans*
  GenBank record, 3 bases off, is dropped (`overrides.tsv`).
- ***H. placei*** carries two V4 variants about 50:50 at insert 138 in
  two isolates (2012, 2023). The A variant is the V4 it shares with *H.
  contortus*; the G variant is in no other record (the *H. contortus*
  genome has A in 12 of 13 copies), so it names *H. placei*: reads from
  a cattle *H. placei* split between *Haemonchus* and *H. placei*.
- ***A. ceylanicum*** (two isolates) has the common *Ancylostoma* V4 -
  not missing, but not separable from *A. caninum* / *duodenale*.
- **A collapse fix:** the narrow trust rule (a lone lower-trust label
  outvoted 3:1 below class is overruled) had overruled *Cyathostomum*'s
  and *Coronocyclus*'s reads. Members from genome and read sources
  (`GENOME_SOURCES` in `bin/collapse.py`) are now never outvoted there;
  they still face the class check and the screens.

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

### 5.1 The neighbour check

The conflict rules see only identical inserts. A mislabelled sequence a
few edits from its true relatives is a record of its own, and names
reads wrongly as confidently as any other. `bin/neighbour_check.py`
finds every pair of records >= 99% identical (exact edit distance over
V4; minimap2 proposes the pairs) whose labels disagree at class or
above - the check the v0.4.0 feedback asked for (§1a there).

On v0.6.0: 2,896 pairs, 2,410 records (`results/neighbours/`).

- **Class-level disagreements (1,771 pairs) are mostly biology.**
  Ascomycota/Basidiomycota (710) and the green-algal classes (484):
  V4 is conserved enough there that different classes sit within five
  edits. Left to the conflict rules.
- **Across supergroups (950 pairs) one side is almost always wrong.**
  A record outnumbered 3:1 in members by such neighbours flags 39
  records - but the rule only finds suspects: it pointed the wrong way
  where one contaminated study supplied several copies (*Betula*,
  *Alnus*, outnumbered by birch "*Eimeria*"), and some are dinoflagellates
  whose V4 really is their diatom or haptophyte endosymbiont's
  (*Durinskia*, *Unruhdinium*, *Karlodinium*).
- **Settled with NCBI BLAST of the whole source record** (15 queries,
  `results/neighbours/suspects_blast.xml.gz`), then
  `curation/overrides.tsv`, each line with its evidence:
  - dropped: the whole "*Eimeria* sp. T72184" study, OP407714-OP407729
    (PR2: *E. lateralis*, *callospermophili*, *monacis*), which is yeast,
    hazel/birch and agaric mushrooms with no *Eimeria* hit at all - reads
    of diet would have been named *Eimeria*; misidentified strains
    ("*Klebsormidium subtilissimum*", "*Nannochloropsis oculata*" x2),
    a bacterium labelled as a plant, an alveolate labelled
    *Dunaliella*; and WGS contigs or clones carrying another organism's
    18S (*Panicum hallii* x3, *S. cerevisiae*, *Rhizophagus*,
    *Phytophthora ramorum*, *Physarum*, *C. elegans* AC006763, an *Apis*
    patent sequence);
  - **reversed:** v0.5.0 kept *Morchella* and dropped the *Vigna*
    contig on `V4_54f09b637b29`; the *Morchella* scaffold is *Vigna*
    18S, 100% over 1,398bp, so now the *Morchella* member is dropped;
  - **relabelled** (a new override action - the sequence is right, the
    name is not): HM017144, filed by PR2 as a red alga, is
    *Angulamoeba microcystivorans* as its own title says; MG452736 is
    *Uronema* the ciliate, not the green alga; OP184778 is *Lunella*
    the diatom, not the sea snail; JN090888 is *Closteriopsis*, not
    *Alphamonas*.
- **Effect** (trial build): 31 more members dropped, 4 relabelled, 29
  records gone (they had no other member). *Betula* and *Alnus* keep
  their records and lose their false *Eimeria* neighbours.

Rerun on v0.7.0 (`results/neighbours/disagreements_v0.7.0.tsv`):
domain/supergroup pairs fell from 956 to 48 - 768 of them had involved
the 29 removed records (one contaminant can have dozens of neighbours),
141 more now agree. Of the 155 pairs at subdivision or above, 11 touch
parasite or pathogen groups; BLAST of the whole record
(`suspects_v0.7.0_blast.xml.gz`) settled five more overrides for the
next release: two more *Angulamoeba* strains PR2 files as a red alga
(HM017143, KP864083) and a cryptophyte it files as a ciliate (PP101464,
99.8% to *Chroomonas*) relabelled; a "*Trachelomonas* chloroplast" that
is betaproteobacterial 16S and two *Thalassiosira* transcriptome
contigs carrying *Bodo* and *Rhynchomonas* 18S dropped. Left alone:
*Pyrsonympha* against *Laeohelix* (both oxymonads; one has only PR2's
placeholder ranks), and two environmental clones.

Still open: the remaining cross-supergroup pairs outside the parasite
groups - mostly single records against single records.

## 6. Release format

A release is one directory, named for its build version, holding:

| file | contents |
|---|---|
| `refdb.fasta` | one record per unique insert, primers trimmed; the header is the record id only |
| `taxonomy.tsv` | id, the nine PR2 ranks (merged to the ranks all members share; placeholder ranks left blank), the source layers contributing, and the number of member sequences |
| `members.tsv` | every source sequence collapsed into each id: accession or genome and contig, layer, original label, role (`kept`, `overruled`, `dropped`) and, for dropped members, the reason (`note`, from v0.5.0) |
| `resolution_nanopore.tsv` | per id, the deepest rank at which it is distinguishable over the whole insert; for species-level ids, `margin`: edits to the nearest record of another species (blank: none near; from v0.8.0, `bin/species_margins.py`) |
| `resolution_illumina_<R1>x<R2>.tsv` | the same over the two read windows for that read length; one per geometry, listed in the manifest's `illumina_geometries` (2x150 and 2x300 from v0.5.0); `margin` over the windows |
| `rank_names.tsv` | what each positional taxonomy column holds, by clade - in animals and most fungi `class` is the phylum (from v0.8.0, `curation/rank_names.tsv`; v0.4.0 feedback §5g) |
| `anchor_lookalikes.tsv` | adapter-anchor near-copies, in exactly the columns the pipeline's pre-flight writes (`seq_id anchor strand position edits context`), so the pipeline reads it unchanged |
| `conflicts.tsv` | identical sequences whose labels disagree at a high rank, and what the build did with each |
| `primer_fit.tsv` | per id, how well the primers fit its best-matching template: forward and reverse edits, mismatch positions counted from each primer's 3′ end, the template at each site, and a class - `exact`, `away_3prime`, `near_3prime` (a mismatch in the last 5 bases), `indel`, `edits_only` (a genome-only record: counts known, sites not), `unknown` (from v0.4.0) |
| `manifest.json` | build version, primers, PR2 version, WBPS release, version of every other layer, build date and code commit; counts, among them primer fit by domain and species margins (from v0.8.0) |
| `summary.txt` | the same for a person; from v0.8.0 also primer fit by domain, species margins, and the clades V4 splits least |

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
- *Trichuris incognita* (the STPHxBioinformatics scaffolded assembly on
  GitHub) carries three V4 copies on one scaffold, both primers exact,
  7-20 edits apart; all three are local sequences. Each is nearest
  *T. suis* (4, 11, 17 edits), then *T. trichiura* (10-23), as expected
  for *T. suis*' sister species. A fourth hit, 1467bp with a ~900bp
  insertion, is left out as a likely misassembly. Whether the spread
  between copies is real rDNA variation or assembly error is not
  resolved; reads or a second assembly would settle it.
- Two *Sphaerolaimus hirsutus* variants sent to review
  (`results/curation/tol_genome_review.tsv`, one contig of
  GCA_964036195.1, 15 edits apart) are **not** the nematode's 18S. They
  are about 170 edits (~69% identity) from all four *S. hirsutus*
  GenBank/PR2 records and from every other *Sphaerolaimus* record, and
  nearest nothing closely in the release (unnamed environmental fungi,
  ~89% over ~380bp). Most likely another eukaryote carried with the
  specimen - gut content or an epibiont. The assembly holds no copy of
  the worm's own V4. NCBI BLAST (core_nt, 2026-10-28) finds nothing
  named: best hits are uncultured fungi and marine eukaryotes at 86% over
  the whole insert, 89-91% over ~380bp - an early-diverging fungal or
  unknown lineage with no close reference anywhere. Kept out;
  *Sphaerolaimus* stays named by its GenBank/PR2 records.

Still no V4 anywhere: *Bunostomum*, *Chabertia*, *Murshidia*,
*Prosthenorchis*, *Rhabdias* - none has an assembly on NCBI.

### 6.7 v0.7.0

Curation release: no new genomes or GenBank fetches. Built by commit
`ef97680`, with one more layer, `blastocystis_ref_layer.tsv.gz`.

- 97,696 records (17 added, 29 removed); helminth records 8,981 ->
  8,983.
- ***Blastocystis* subtypes** (§4.7): 581 of 653 *Blastocystis* records
  named `Blastocystis_STn`, all 45 accepted subtypes (ST1-ST17, ST21,
  ST23-ST49) present, each kept read whole and through both Illumina
  geometries. 14 of the added records are subtype references no title
  search had found. "*B. hominis*" is gone.
- **The neighbour check** (§5.1): 31 members dropped and 4 relabelled
  in `curation/overrides.tsv`, each with BLAST or neighbour evidence;
  the 29 removed records had no other member. Among them the whole
  "*Eimeria* sp. T72184" study (yeast, hazel/birch and agaric mushroom
  18S), so diet reads are no longer named *Eimeria*. `V4_54f09b637b29`
  is now *Vigna radiata*, reversing v0.5.0. `members.tsv`'s `note`
  column now also says "relabelled: <reason>".
- ***Trichuris incognita***: three local records from the
  STPHxBioinformatics assembly, each nearest *T. suis* (§6.6 notes).
- *Schistosoma intercalatum* TD2, scanned from NCBI because the WBPS
  copy is truncated, supports the shared *S. haematobium*-group record;
  no new record.

### 6.8 v0.8.0

Curation from the v0.7.0 feedback and the VESPA gaps, and three
additions to the release format. Built by commit `e8d6145`, with one
more layer, `genbank_targeted_layer.tsv.gz`.

- 97,739 records (46 added, 3 removed, 31 relabelled or resolved
  differently); helminth records 8,983; species-level 49,040.
- **Curation:** five more neighbour-check overrides (§5.1: two
  *Angulamoeba* strains and a cryptophyte relabelled; a "*Trachelomonas*
  chloroplast" that is bacterial 16S and two diatom transcriptome
  contigs carrying kinetoplastid 18S dropped); *Balantidium coli*
  harmonised under *Balantioides* (39 relabels), so 4 records that
  stopped above genus now name it.
- **New sequences (§4.8):** a targeted GenBank layer - trichomonads,
  diplomonads, archamoebae and trichostome ciliates never fetched, and
  six accessions BLAST found for unmatched VESPA clusters (among them a
  primate *Buxtonella*-like ciliate, now 98%).
- ***Blastocystis* (§4.7):** all 100 NCBI genomes scanned (19 genome
  variants, was 6) and PubMLST's 87 typed alleles as extra subtype
  references - 658 records, 596 named to subtype (was 653 and 581).
- **Release format (§6):** a `margin` column in each resolution table
  (`bin/species_margins.py`; ~30% of species calls rest on 1-2 edits);
  `rank_names.tsv`, what PR2's positional columns hold by clade; and in
  `summary.txt`, primer fit by domain, the margins, and the clades V4
  splits least.

### 6.9 v0.9.0

V4 from raw reads and new assemblies (§4.10), and the watch list
(§4.9). Built by commit `a6ac2b9`, with two more layers,
`reads_layer.tsv.gz` (source `ENA_reads`) and
`tol_prerelease_layer.tsv.gz` (`ToL_prerelease`).

- 97,784 records (46 added, 1 removed); helminth records 8,983 ->
  9,028; species-level 49,040 -> 49,083.
- Watch list: 89 of 122 helminths to species (v0.8.0: 80), 11 only in
  shared records, 22 missing (32).
- New species: *Bunostomum trigonocephalum*, *Cyathostoma americana*,
  *Strongylus vulgaris*; *Cylicostephanus goldi* / *longibursatus*,
  *Cylicocyclus insigne*, *Coronocyclus labiatus*, *Calicophoron
  daubneyi* and *Haemonchus placei* by variants of their own.
- Corrected: the cyathostomin main V4 (`V4_c474d5b76665`, six species)
  no longer named *Cylicocyclus*; *Mansonella*'s main V4, shared by
  *M. perstans* and *M. ozzardi*, added from four read sets, and the
  *M. perstans* GenBank record 3 bases off it dropped; an *O. volvulus*
  variant added so its reads no longer tie with *D. repens*.
- Confirmed from independent data: *Loa*, *D. repens*, *A. viteae*,
  *Onchocerca* (with *O. ochengi* and *O. lupi* sharing its V4), *A.
  ceylanicum* (sharing *Ancylostoma*'s), *Heligmosomoides*, *D. immitis*.
- `bin/collapse.py`: genome and read members are never outvoted by the
  narrow trust rule.

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
- **Partial references.** A record must contain both primer sites, so
  species known only from partial 18S are absent: about 8,200 named
  species in PR2 and our GenBank downloads (v0.7.0), 2,900 of them
  helminths (`results/partials/partial_only_species.tsv`, best record per
  helminth, apicomplexan and protist species). By their best record,
  helminths split roughly as: one primer plus 350-699bp of insert (559),
  one primer plus one Illumina window's worth (128-349bp, 570), less
  (122), no primer site at all (1,634 - outside V4 or a short internal
  stretch; not yet aligned), and full-length sequences whose other
  primer site is > 2 edits off (69, mostly land planarians - arguably
  rightly out).
  - *Near-complete tier* (measured, `results/partials/measure_partials.py`
    -> `near_complete_partials.tsv`): aligning the 648 "350-699bp"
    candidates to full V4 of their genus or family, only **49 species
    are complete or <= 10bp short** (66 at <= 30bp); most nematode
    partials stop > 100bp before the reverse primer. The 49 include
    *Acanthocheilonema reconditum*, *Mansonella atelensis*, *Gurltia
    paralysans*, *Angiostrongylus mackerrasae*, *Setaria tundra*, the pig
    *Eimeria* (*E. suis*, *debliecki*, *neodebliecki*, *spinosa*,
    *perminuta*), *Cryptosporidium scrofarum* and *C. equi*, *Cyclospora
    macacae*, *Moniliformis tupaia*. Proposed handling, not built: admit
    them trimmed to the insert, flagged in a `partial.tsv` (missing end,
    bp), primer fit "unknown" at the missing site, joining an existing
    record where contained in it. Before building: GenBank candidates
    with no full record in their genus were measured against a distant
    family member (no lineage lookup) and may be misjudged.
  - *One-window tier* (~570 helminth species, plus the truncated
    "350-699bp" ones): useful to Illumina reads only, in the window they
    cover - needs partial-record support in the pipeline (amphid).
- Whether to filter to lineages of interest, or keep all eukaryotes so
  off-target reads are assigned rather than forced onto the nearest
  parasite. (Current lean: keep everything; the pipeline reports
  off-target fractions from it.)
