# Week 9 — Screening two readsets for non-target sequence

Detect, measure and attribute non-target (contaminating or unexpected) sequence
in one RNA-seq readset and one DNA WGS readset, and decide for each detection
how much it matters and where it most likely came from.

| | RNA-seq | DNA WGS |
|---|---|---|
| Sample | `BYJ10` | `SRR39197041` |
| Host | *Gallus gallus* (chicken), jejunum | *Aedes aegypti* (mosquito) |
| Library | paired-end RNA-seq, 2×150 bp | paired-end WGS, 2×151 bp |
| Reads | 42,352,610 pairs | 60,394,270 pairs |
| Flowcell / lane | `23C55KLT4` run 534 lane 5 | `H5F2GDSX3` lane 4 |
| Stated context | ARV-infected bird | field/colony mosquito, SRA |

Everything was run first on a **100,000-read-pair test subsample** of each
readset to validate the pipeline, then on the full readsets. Test outputs are
kept under `*/test/`, full outputs under `*/full/`.

---

## How to reproduce

```bash
scripts/00_fetch_kraken2_db.sh      # prebuilt Kraken2 PlusPF-16 index
scripts/00_make_conda_envs.sh       # kraken2 2.1.3 + bracken; genomad + checkv
scripts/00_fetch_host_refs.sh       # chicken + mosquito genomes and transcriptomes
scripts/00_fetch_virus_dbs.sh       # geNomad DB, CheckV DB, RefSeq viral proteins

qsub -v LEVEL=test scripts/01_sendsketch.sh          # screen A
qsub -v LEVEL=test scripts/02_kraken2_bracken.sh     # screen B
qsub               scripts/02b_kraken2_controls.sh   # co-sequenced controls
bash               scripts/03_fetch_candidates.sh    # candidate + decoy genomes
qsub -v LEVEL=test scripts/04_competitive_map.sh     # competitive mapping
python3            scripts/05_score_detections.py --sample BYJ10_RNA --level test
qsub -v LEVEL=test,SAMPLE=BYJ10_RNA scripts/06_assemble_and_virus.sh

# then repeat with LEVEL=full
```

All jobs are PBS scripts for the ASC `asax` cluster. `scripts/env.sh` loads the
Lmod modules and is safe to source under `set -e`.

---

## Design: why these tools, in this order

The pipeline is deliberately built as **three independent lines of evidence**
that fail in different ways, followed by a **measurement** step and an
**attribution** step. No detection is reported on the strength of one method.

### Step 1 — `sendsketch.sh` (BBTools): database-free hypothesis generation

`sendsketch.sh` builds a MinHash sketch of the reads locally and sends only the
sketch to JGI's server, which compares it against the complete RefSeq and nt
sketch collections. That gives a genome-wide, **uncapped** answer in seconds
with no local database at all.

It is the right instrument for *"what might be in here?"* and the wrong
instrument for *"how much?"* — a sketch hit of four shared hashes looks
identical whether it came from a real organism at low abundance or from four
spurious reads. So every hit it raises is treated strictly as a hypothesis.

Both `refseq` and `nt` are queried, because they fail differently: RefSeq is
curated and has clean taxonomy but misses unusual taxa; nt is far broader but
is full of mislabelled entries. **A hit present in both is much stronger
evidence than a hit in either alone**, and that disagreement is used explicitly
in the detection table.

### Step 2 — Kraken2 + Bracken: exact k-mer classification, per read

Kraken2 classifies *every* read against an exact k-mer index, so unlike
sendsketch it produces countable per-taxon read numbers — the quantity the four
significance bins actually need.

* **Database: PlusPF-16GB** (bacteria, archaea, viral, plasmid, human,
  UniVec_Core, protozoa, fungi). Fungi and protozoa are included deliberately:
  reagent and lab-organism contaminants are frequently fungal, and a
  bacteria-only database would simply be blind to them.
* **`--confidence 0.1`.** At the default of 0 a read is assigned if *any* k-mer
  matches, which generates a long tail of single-k-mer false positives —
  exactly the "noise" bin we are trying to exclude. Requiring 10 % of a read's
  k-mers to support the call costs a little sensitivity, which is acceptable
  because sendsketch already covers the sensitive end.
* **`--report-minimizer-data`.** This is the single most useful
  anti-false-positive flag Kraken2 has. It reports the number of *distinct*
  minimizers supporting each taxon. A genuinely present genome is hit across
  many different minimizers; a spurious taxon is hit by a handful of minimizers
  over and over (one repeat, one conserved rRNA). Nothing in the plain report
  distinguishes these two cases.
* **Bracken** then redistributes reads that Kraken2 parked at genus/family
  (because the species has close relatives in the database) down to species
  level. **See the caveat below — in this dataset Bracken is actively
  misleading, and that is itself one of the findings.**

### Step 2b — control screens on the co-sequenced samples

This is what makes source attribution falsifiable rather than a plausibility
argument. The read headers reveal the batch structure:

| Group | Samples | Relationship |
|---|---|---|
| A | **BYJ10**, AYJ9 | same flowcell `23C55KLT4`, **same lane 5** |
| B | SYJ7, SYJ10 | **different** flowcell `237VF5LT4` lane 7 |
| C | **SRR39197041**, SRR39197043 | same flowcell `H5F2GDSX3`, **same lane 4** |

Those groupings separate source categories that read fraction, depth and
breadth cannot separate on their own:

* **Bleed-through (index hopping)** can only come from a library on the *same
  lane*. If a taxon is abundant in AYJ9 and present in BYJ10 at ~0.1–1 % of
  AYJ9's level, index hopping is the leading explanation. If the taxon is
  absent from AYJ9, index hopping is *ruled out*.
* **Reagent contamination** tracks the kit, not the flowcell — it should appear
  across group A *and* group B.
* **A real biological associate** tracks the sample, not the batch.

Each source call therefore makes a prediction that the controls confirm or
refute.

### Step 3 — competitive mapping (minimap2 + samtools): the measurement

**Why minimap2 rather than bowtie2.** The combined reference is ~4 Gb (host +
human decoy + every candidate). minimap2 indexes that in minutes where
`bowtie2-build` takes many hours, and the ASC bowtie2 module is 2.2.9 (2016).
The `-ax sr` preset is a short-read mode as accurate as bowtie2 for 150 bp
pairs, and it tolerates **divergent** alignment better — which matters a great
deal here, because a candidate reference is usually a database *relative* of
the organism actually present, not the exact strain. bowtie2 end-to-end would
drop those reads and we would under-call real detections.

**Why one combined reference rather than one BAM per candidate.** Mapping to a
candidate in isolation is the commonest way to manufacture a false
contamination result: conserved rRNA, ribosomal-protein and transposon reads
from the host and from everything else have nowhere else to go, so they pile
onto the candidate's conserved loci and produce impressive-looking depth over a
few narrow windows. In a competitive reference every read goes to its single
best match, so a read is credited to a candidate only when that candidate
explains it better than the host, the human decoy and every lab-organism decoy.

**The decoy panel** (human, PhiX174, *E. coli* K-12, *S. cerevisiae*,
UniVec_Core) exists precisely to absorb reads that would otherwise be forced
onto a candidate, and it maps one-to-one onto the source categories we have to
adjudicate between (human / lab organism / reagent-vector).

**Why MAPQ ≥ 20.** MAPQ encodes how much better the best alignment is than the
second best. Requiring ≥ 20 discards reads that map about equally well to the
candidate and to something else — precisely the ambiguous reads that inflate
breadth. This filter is what makes breadth mean *"sequence unique to this
organism is present"*.

**`-I 12G`** forces a single-part minimap2 index. With a multi-part index
minimap2 maps against each part independently and computes MAPQ within a part,
which silently destroys the competitive comparison the whole design rests on.

### Step 4 — scoring (`05_score_detections.py`)

Per candidate, from `samtools coverage` over the whole reference:

* **read fraction** = reads at MAPQ ≥ 20 ÷ total reads in the library
* **depth** = mean per-base depth across the entire reference
* **breadth** = fraction of reference bases covered at ≥ 1×

plus three statistics that exist only to catch false positives:

* **`mapq_survival`** = reads at MAPQ ≥ 20 ÷ reads at MAPQ ≥ 0. A real
  organism's reads are unique to it and keep a high MAPQ; a cross-mapping
  artifact collapses when ambiguous reads are removed.
* **`gini`** = Gini coefficient of per-base depth. Breadth and mean depth cannot
  distinguish an even genome-wide 0.5× from a 25× spike over 2 % of the genome,
  yet those mean opposite things. A genuine organism gives low-to-moderate
  Gini; a conserved-locus or database-error artifact gives Gini near 1.
* **`top1pct_share`** = share of aligned bases in the most-covered 1 % of the
  reference — a blunter, more explainable version of the same idea.

### Step 5 — RNA-virus arm (`06_assemble_and_virus.sh`)

Everything above compares reads to references that already exist. A divergent
RNA virus produces reads that sendsketch and Kraken2 simply leave unclassified
— the commonest way a real viral finding is missed. Assembly reconstructs the
genome from the data itself.

* **Host depletion first**, against the host **genome and transcriptome**. With
  the `-ax sr` short-read preset a read spanning an mRNA splice junction aligns
  poorly to the genome and would survive depletion as a false "non-host" read;
  the transcriptome catches exactly those. Only pairs where *neither* mate maps
  are carried forward, which is the conservative choice.
* **MEGAHIT and rnaSPAdes, both.** They fail in different places and after
  depletion the read set is small enough that running both is cheap. MEGAHIT
  copes with the enormous dynamic range of a contamination screen; rnaSPAdes
  reconstructs full-length RNA virus genomes markedly better. **A contig
  recovered independently by both is much stronger evidence than one recovered
  by either alone.**
* **geNomad** scores contigs as chromosome / plasmid / virus from marker-gene
  content rather than nucleotide similarity, so it can flag a virus with no
  close database relative — and it reports host-vs-viral gene content along a
  contig, which is direct evidence for calling something an *integrated
  element* rather than a free virus.
* **CheckV** estimates viral-contig completeness, which is exactly what
  separates "a viral fragment is present" (trace) from "a near-complete viral
  genome is present" (biological). It also flags host sequence *within* a
  contig — a second, independent read on the same EVE question.
* **DIAMOND blastx** against RefSeq viral proteins is independent of geNomad's
  marker set; protein similarity survives far greater nucleotide divergence, so
  it is the most sensitive way to name a divergent virus.

---

## The four significance bins

Thresholds are applied in this order, so a hit can never be promoted on the
strength of numbers the credibility checks have already discredited.

| Bin | Rule | Meaning |
|---|---|---|
| **noise** | < 50 unique reads, **or** breadth < 0.5 %, **or** `mapq_survival` < 0.25, **or** (Gini > 0.95 and breadth < 5 %) | Not credible. Consistent with cross-mapping or a database artifact. |
| **trace** | passes the above, but read fraction < 0.1 % **and** depth < 1× **and** breadth < 10 % | Real, but there is not enough of it to change any downstream result. |
| **consequential** | read fraction ≥ 0.1 % **or** depth ≥ 1× **or** breadth ≥ 10 % | Enough to shift expression estimates, pull mis-assembled contigs into an assembly, or seed false variant calls. |
| **biological** | breadth ≥ 50 % at depth ≥ 1×, **or** CheckV completeness ≥ 90 % | A near-complete genome is present; this is a finding about the sample, not a QC footnote. |

## Source categories and the evidence that distinguishes them

| Source | Evidence that supports it over the alternatives |
|---|---|
| **real associate / cobiont** | Broad, even coverage (low Gini) at a depth consistent with a real population; tracks the *sample* not the batch; absent or different in co-lane controls; ecologically coherent with the host. |
| **integrated in host genome (EVE)** | Reads also align to the host assembly; depth tracks host depth at a fixed ratio; contigs carry viral sequence flanked by host sequence; present at the same level in *every* sample of that host regardless of flowcell. |
| **reagent** | Taxon is a known kit organism; appears across *different* flowcells and library batches; abundance inversely related to input biomass; low, patchy breadth. |
| **bleed-through (index hopping)** | Taxon is the *dominant* organism of another library on the **same lane**, present here at ~0.1–1 % of that level. Excluded if the co-lane sample lacks it or has it at equal level. |
| **library chimera** | Reads are internally chimeric — part host, part other — with breakpoints at ligation junctions; discordant pairs; coverage islands at the junction only. |
| **human** | Clean primate phylogeny with *Homo sapiens* on top and unique (non-shared) matches; even genome-wide coverage distinguishes a bloodmeal or real contamination from conserved-vertebrate mis-assignment. |
| **lab organism** | PhiX174, *E. coli* K-12, *S. cerevisiae*, cloning vectors — present in the decoy panel so they claim their own reads instead of contaminating a candidate's count. |
| **database error** | Reads pile onto a *narrow region* of the candidate reference (high Gini, low breadth) and that region is itself a contaminant insert in the reference; `mapq_survival` collapses once the true source is in the competitive reference. |

---

## Repository layout

```
data/                 raw readsets (not tracked) + 100k-pair test subsamples
metadata/
  datasets.txt        provenance of both readsets
  candidates.tsv      reference panel for competitive mapping (name/role/source/id/restrict)
references/
  host/               chicken + mosquito genomes and transcriptomes
  candidates/         candidate and decoy genomes
  combined/           per-readset competitive references and cached minimap2 indexes
  db/                 Kraken2 PlusPF-16, geNomad, CheckV, RefSeq viral proteins
scripts/              numbered pipeline stages, each with its rationale in the header
results/
  sendsketch/{test,full}/
  kraken2/{test,full,controls}/
  coverage/{test,full}/     BAMs, samtools coverage, per-base depth, detection tables
  assembly/{test,full}/     contigs, geNomad, CheckV, DIAMOND
reports/              detection table and the one-page write-up
logs/                 PBS job logs
```
