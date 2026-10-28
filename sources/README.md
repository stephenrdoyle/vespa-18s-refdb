# sources/

Downloaded source databases. Nothing here is committed except this file;
each release's `manifest.json` records exactly which versions it used.

| file | from |
|---|---|
| `pr2_version_5.1.1_SSU_taxo_long.fasta.gz` | https://github.com/pr2database/pr2database/releases/tag/v5.1.1 - headers carry accession, gene, organelle, strain and the nine ranks |
| `pr2_version_5.1.1_chimera.xlsx` | same release. PR2's list of sequences removed as chimeras; checked, none of them remain in the 5.1.1 FASTA, so the build does not need it |

WormBase ParaSite genomes are streamed from
https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/ and not
stored.
