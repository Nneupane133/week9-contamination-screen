# The virus in the "ARV-infected" bird is part of the bird

**Detection:** EAV-HP (endogenous avian retrovirus, `NC_005947.1`), in the chicken jejunum RNA-seq library `BYJ10`
**Final call:** `biological` · source `integrated-in-host` (transcription of an endogenous retroviral element)
**What it is not:** an infection, a contaminant, or the avian orthoreovirus the sample was submitted for

---

## What was found

Among all the non-host references in the chicken library, only one is covered almost end to end:

| | MAPQ ≥ 0 | MAPQ ≥ 20 |
|---|---|---|
| reads | 3,069 | 467 |
| read fraction | 0.0036% of 84.7 M reads | — |
| breadth | **98.8%** of 4,302 bp | 95.0% |
| depth | **105.7×** | 16.1× |

By the usual criteria (near-complete breadth at >100× depth) this looks like an actively
replicating retrovirus. The sample's stated context is "ARV-infected", but the expected virus,
**avian orthoreovirus, has 0 reads at any MAPQ across all 10 segments** and 0 Kraken2
reads in all four chicken libraries. So the one complete-looking virus in the library is not
the one the bird was supposed to carry.

## Why it is the chicken, not an infection

**1. It competes with the host genome.** In `04b`, the 3,217 reads assigned to EAV-HP were
re-aligned against the chicken genome plus transcriptome. **2,544 of them (79.1%) align at
least as well to chicken.** EAV-HP is a family of endogenous retroviruses with many copies
in the *Gallus* genome, so these reads are transcripts from chicken loci. The low MAPQ
survival (15%) is the same signal seen from the other side: most reads have an equally good
second placement.

**2. It tracks the host species, not the batch.** If this were an infection, the signal
would follow the sample. If it were index hopping, it would follow the lane. If it were a
reagent, it would follow the kit. Kraken2 places retroviral reads at the same order of
magnitude in every chicken library, across both flowcells:

| Library | Flowcell / lane | Retroviridae reads | per million reads |
|---|---|---|---|
| **BYJ10** | `23C55KLT4` L5 | 8,992 | 212 |
| AYJ9 (same lane) | `23C55KLT4` L5 | 12,828 | 335 |
| SYJ7 | `237VF5LT4` L7 | 16,074 | 447 |
| SYJ10 | `237VF5LT4` L7 | 11,318 | 288 |

The same-lane partner AYJ9 is not 100–1,000× higher than BYJ10, which rules out bleed-through.
The other flowcell has the same signal, so it is not a single-batch artifact. A ~2-fold
spread between birds is what expression of an endogenous locus looks like.

**3. Its relatives fall apart the same way.** The other avian retroviruses in the panel give
the control result. **Y73 sarcoma virus** attracts 576 reads, but 98.3% of them prefer the
host, and its coverage covers only the 34% of the reference that is the transduced chicken
*YES1* gene (v-*yes*). Those reads are chicken *YES1* mRNA. **Rous sarcoma virus** (92%
host-competing, transduced v-*src*) and **avian leukosis virus** (78% host-competing) keep
0 reads at MAPQ ≥ 20. EAV-HP differs from them in degree, not in kind. It is the most
complete because the whole element, not one captured gene, is in the chicken genome.

**4. The assembly arm never saw it.** Host depletion left 1,199,736 read pairs for
assembly, and none of the 6,796 contigs aligns to EAV-HP. The reads were removed as chicken
before MEGAHIT and rnaSPAdes ran. This is the opposite of PhiX, where assembly was the only
method that found the genome.

## Why the call was hard

Without `04b`, `05_score_detections.py` bins EAV-HP as `noise` because of the low MAPQ
survival, and bins Y73 sarcoma virus as `consequential` because of its 380 high-MAPQ reads.
Both calls are wrong in opposite directions, and they fail for the same reason: the
competitive reference had the chicken *genome* but not its *transcriptome*. Spliced host
mRNA (*YES1*) and expressed endogenous proviruses had nowhere better to go. The
host-competition test is the only step that separates "a viral sequence is expressed" from
"a virus is present".

Kraken2/Bracken does not settle it either. Bracken pushes 1,223 reads down to "Avian
leukosis virus" in BYJ10, where Kraken2 itself assigned only 10. This is redistribution, not
evidence.

## What follows from it

1. **Report EAV-HP as host biology, not as a finding about this bird.** It is real
   (`biological`), but it is the chicken's own retroviral load. It is expected in every
   chicken RNA-seq library and says nothing about infection status.
2. **For RNA-seq, compete against the host transcriptome, not only the genome.** Otherwise
   spliced host genes carried by oncogenic retroviruses (v-*yes*, v-*src*) look like viral
   reads.
3. **The negative result matters more.** The library contains no detectable avian
   orthoreovirus at any depth. If BYJ10 is meant to be an ARV-infected bird, that label needs
   checking: by the time of sampling, virus was either absent from jejunum or below detection
   in 42 M read pairs.
