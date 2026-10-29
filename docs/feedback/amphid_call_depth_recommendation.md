# Recommendation to amphid: cap each call's depth by its distance to the best hit

From vespa-18s-refdb, 2026-10-29, on release v0.7.0. Evidence:
`bin/simulate_leave_one_out.py`, results in `results/leave_one_out/`.

## The problem

The database-mode assignment names the deepest rank its best hits share,
at any distance. When a read's species is in the database that is
right. When it is not - and about 8,200 named species have no V4 record
(`docs/DESIGN.md` §8) - the read borrows its nearest relative's name.

Leave-one-out, each of about 4,000 helminth and apicomplexan species
removed in turn, its genus still in the database:

| reads | named as a sister species | named to a species of another genus | stops at the right genus (best possible) |
|---|---:|---:|---:|
| exact, whole insert | 35-46% | 3-15% | 29-43% |
| 2x150 Illumina (0.5% substitutions) | 27-41% | 3-15% | 31-48% |
| nanopore (2% errors, tolerance 2) | 16-30% | 1-11% | 36-55% |

Ranges are across Nematoda, Platyhelminthes, Acanthocephala and
Apicomplexa. Where the species was its genus's only one, nearly half
of its reads are still named to some species or genus, all wrongly.

## The fix: a distance cap on depth

A wrong call sits further from its best hit than a right one. With the
same reads scored against the full release (species present,
`leave_one_out_present.tsv`) and with their species removed:

| reads | species-level only if within | correct species calls kept | wrong species calls turned into genus or above |
|---|---|---:|---:|
| 2x150 Illumina | 1.0% | 86% | 83% |
| | **1.5%** | **96%** | **72%** |
| | 2.0% | 99.8% | 56% |
| nanopore, 2% errors | 2.0% | 55% | 92% |
| | **3.0%** | **96%** | **65%** |
| | 4.0% | 99.9% | 45% |

| reads | genus-level only if within | correct genus-or-better calls kept | wrong genus calls demoted |
|---|---|---:|---:|
| 2x150 Illumina | **5%** | 100% | 39% |
| nanopore, 2% errors | **5%** | 100% | 45% |

Suggested starting values, as parameters: **species within 1.5%
(Illumina) or 3% (nanopore) of the best hit; genus within 5%;** beyond
that, the shared rank above genus. Distance = edits / aligned read
length, as the simulation measured it (Illumina: both windows summed
over 260 bases).

## Caveats

- **"Kept" is optimistic.** Present-species reads came from a record
  that is itself in the database. A real read from a variant of a known
  species that the database lacks sits a little further away - within-
  species V4 variation is typically under 1% - so some correct calls
  near the cap would be demoted to genus. That is the safe direction.
- **Some wrong calls no cap can catch.** About 20% of wrong species
  calls from exact reads are within 0.5% of a sister: V4 does not
  separate those species, and only adding the missing species helps.
- **Calibrated on helminths and Apicomplexa.** Groups where V4 varies
  more or less within a genus (fungi, green algae) may want different
  values; the scripts rerun on any release.
- **Nanopore error rates differ by chemistry and basecaller.** The 3%
  assumes about 2% read error; the species cap should track the run's
  error rate plus about 1%.

## Reporting

A capped call should say so - for example "genus *Trichuris*, best hit
*T. suis* at 2.4%" - so a reader can tell a call limited by V4 from one
limited by distance. The *Trichinella nelsoni* calls in the VESPA
primate samples (v0.4.0 feedback §5d) are the test case: their distance
to the *T. nelsoni* records decides whether they are a finding or a
missing capillariid.
