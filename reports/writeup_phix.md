# The complete genome that four methods missed and the fifth mis-named

**Detection:** PhiX174, complete genome, in the *Aedes aegypti* WGS library `SRR39197041`
**Final call:** `biological` · source `lab-organism` (Illumina sequencing control spike-in)
**What the pipeline said:** `noise — fewer than 50 mapped reads`

---

## What was found

The assembly arm host-depleted the library to 3,873,512 read pairs and assembled them
twice, with MEGAHIT and with rnaSPAdes. Two contigs came back that looked, on every
measure the pipeline applies, like a significant viral discovery:

| | rnaSPAdes | MEGAHIT |
|---|---|---|
| contig | `NODE_2_length_5459_cov_5.27` | `k141_52513` |
| length | 5,459 bp | 5,527 bp |
| topology (geNomad) | **DTR** (circular) | **DTR** (circular) |
| geNomad virus score | 0.9826 | 0.9822 |
| hallmark genes | 2 | 2 |
| geNomad taxonomy | Monodnaviria → Petitvirales → **Microviridae** | same |
| CheckV | **Complete**, 100.0% | **Complete**, 100.0% |
| CheckV contamination | 0.0% | 0.0% |

A circular 5.5 kb ssDNA phage genome, 100% complete, no host contamination, two
hallmark genes, recovered **independently by two assemblers** — which is exactly the
corroboration the pipeline design treats as its strongest evidence.

It is PhiX174. Aligned to `NC_001422.1` the contig matches over its full length in two
colinear blocks (the assembler opened the circle at a different point): 5,468 of 5,504
aligned bases identical, **99.3% nucleotide identity**, GC 44.6% against PhiX's 44.8%.
The 141 bp by which the contig exceeds the 5,386 bp reference is the terminal repeat
that told CheckV the genome was circular and complete. The "novel complete Microviridae
genome" is Illumina's sequencing control.

## Why every database method missed it

PhiX is in the reference set. It is in the Kraken2 database. It is in the decoy panel,
deliberately, as a named lab organism. All of them failed, and for one shared reason.

**Competitive mapping reported 0 reads at MAPQ ≥ 20.** `UniVec_Core` — a decoy added to
absorb vector and adapter reads — contains the entire phi-X174 genome as
`gnl|uv|J02482.1:1-5386-49`, **5,375 of 5,375 bases identical** to the PhiX decoy. Every
PhiX read therefore has two perfect placements, gets MAPQ 0, and is deleted by the
MAPQ ≥ 20 filter that the whole measurement rests on.

**Kraken2 reported 0 reads.** Of the 394 reads that align to PhiX in the BAM, Kraken2
called 235 *unclassified* and 159 to *taxid 1 — the root of the taxonomy*. None reached
taxid 10847. PlusPF-16 also carries phi-X174 twice, once as the phage and once inside
UniVec_Core as a synthetic construct, so the lowest common ancestor of a PhiX read's
k-mers is the root of the tree.

**sendsketch** returned only "synthetic construct", 109 matches, 0 unique.

These are not three separate bugs. Kraken2's LCA and minimap2's MAPQ are the same idea
applied to k-mers and to alignments: *score a read by how much better its best
placement is than its second best.* A sequence stored twice has no second-best
placement, so both methods score it at zero and discard it. Any method built on that
principle is structurally blind to a duplicated reference entry — and the more
carefully the panel is stocked with decoys, the more likely the duplication becomes.

Assembly was the only arm that saw it, because assembly never asks "where else could
this read go?" It asks "what do these reads spell?" That is also why it could not name
what it found: with the genome rebuilt, naming it means going back to the same
databases that were blind to it.

## What the numbers actually are

Read at MAPQ 0, with the two copies summed:

| | |
|---|---|
| reads | 608 (316 on the decoy + 292 on the UniVec copy) |
| read fraction | **0.0005%** of 120,788,540 reads |
| breadth | **100%** of the 5,386 bp genome |
| depth | **~14×** |

This is the clearest case in the project for why **read fraction is the weakest of the
three statistics**. At 0.0005% of the library PhiX sits in the bottom third of the
library's twelve non-host signals by read count — 6,600 times fewer reads than the
adapter read-through at 3.33% — and it is the only non-host organism in the library
whose genome is present in its entirety. Breadth is what distinguishes "a complete organism is here" from "some reads
landed here", and it is the statistic the automatic rules discarded first.

## What follows from it

1. **Check the panel for redundancy before mapping, not after.** An all-vs-all
   `minimap2` of the reference panel against itself takes minutes and would have caught
   `UniVec_Core` ⊃ `PhiX174` before a single read was aligned. This is now the one
   pre-flight check the design was missing.
2. **De-duplicate rather than drop.** PhiX must stay in the panel — it is a real decoy
   doing real work — so the phi-X174 entry should be removed from `UniVec_Core`, not the
   other way round.
3. **A reference that is a collection is not a genome.** The same mistake, in a second
   form, made the pipeline call 4,016,896 adapter reads (3.33% of the library) "noise"
   because they covered only 0.23% of UniVec's 688 kb — when they covered 100% of the
   58 bp Illumina PCR primer they actually came from.
4. **Report PhiX as a QC number, not a detection.** 0.0005% residual after the
   instrument's own PhiX filter is normal, and worth recording as such. The chicken
   RNA-seq library, by contrast, contains zero PhiX reads — the spike-in is specific to
   this DNA run.
