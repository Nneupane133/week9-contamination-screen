# Detection table (full)

`bin` is the final call; where the automatic rule in `05_score_detections.py` disagreed, its verdict is shown as `auto: <bin>` and the reason is in the evidence column.

`evenness` is observed breadth divided by the breadth a randomly placed library of that depth would give (`1 - e^-depth`). Near 1 means the reads are spread across the reference like a real population's; near 0 means they are piled on a few loci.

† measured at MAPQ 0 rather than MAPQ>=20, because the sequence is integrated in the host genome or present twice in the reference panel — in both cases a MAPQ filter throws away its real reads. MAPQ survival is left blank for those rows, where it would be 1.0 by construction.

## BYJ10_RNA — Gallus gallus (chicken)

*jejunum RNA-seq, 2x150 bp, flowcell 23C55KLT4 lane 5; 84,705,220 reads; host claims 79,175,062 (93.5%) at MAPQ>=20.*

| organism | bin | likely source | reads | read fraction | depth | breadth | evenness | MAPQ surv. | why |
|---|---|---|---:|---:|---:|---:|---:|---:|---|
| EAV-HP endogenous retrovirus † | **biological** | integrated-in-host | 3,069 | 0.00362% | 106 | 98.8% | 0.99 | - | 98.8% of the 4,302 bp genome covered at 105.7x; 79.1% of its reads align at least as well to the chicken genome+transcriptome (04b); alpharetroviral signal is equal in all four chicken libraries across two flowcells, so it tracks the host species, not the batch |
| Homo sapiens | **trace**<br>`auto: noise` | human | 17,941 | 0.0212% | 0.000699 | 0.00446% | 0.064 | 0.12 | 17,941 reads (0.021%) survive MAPQ>=20 and 69% of them are not host-competing; Kraken2 puts human at 0.08-0.10% in all six libraries in this project across two host species, two labs and three flowcells -- a constant floor, not a property of this sample |
| NEBNext adaptor (UniVec NGB00846) | **trace**<br>`auto: noise` | reagent | 251 | 0.000296% | 125 | 100% | 1 | 0.0023 | 295 reads at MAPQ>=20, of which 251 sit on a single 65 bp NEBNext adaptor at 100% breadth; adapter read-through, not an organism |
| Y73 sarcoma virus † | **noise** | database-error | 576 | 0.000680% | 16.4 | 44.2% | 0.44 | - | 98.3% of its reads prefer the host once the chicken transcriptome is in the competition (04b); coverage is confined to the 34% of the reference that is the transduced host v-yes/YES1 gene, so these are chicken transcripts, not virus |
| UniVec_Core, all other entries (vectors, other adapters) | **noise** | reagent | 44 | 0.000052% | - | - | - | - | the remainder of the UniVec_Core collection once the one entry that carries the bulk of the reads is counted separately; spread thinly over ~3,000 unrelated vector and adapter sequences |
| Rous sarcoma virus † | **noise** | database-error | 15 | 0.000018% | 0.227 | 9.60% | 0.47 | - | 0 reads at MAPQ>=20; 92% of the 77 reads it attracts are host-competing (transduced v-src) |
| Avian leukosis virus † | **noise** | integrated-in-host | 7 | 0.000008% | 0.0498 | 2.73% | 0.56 | - | 0 reads at MAPQ>=20; the 18 reads it attracts are 78% host-competing -- endogenous ALV-like background, the same signal EAV-HP captures |
| Escherichia coli K-12 | **noise** | lab-organism | 6 | 0.000007% | 0.000176 | 0.0124% | 0.71 | 0.023 | 6 reads at MAPQ>=20; decoy working as intended |
| Aeromonas hydrophila | **noise** | database-error | 2 | 0.000002% | 2.59e-05 | 0.00133% | 0.51 | 0.0046 | 2 reads at MAPQ>=20 out of 433; survival 0.005 |
| Saccharomyces cerevisiae | **noise** | lab-organism | 1 | 0.000001% | 1.15e-05 | 0.00115% | 1 | 0.0086 | 1 read at MAPQ>=20; decoy working as intended |
| Avian orthoreovirus | **noise** | absent | 0 | 0% | 0 | 0% | - | - | 0 reads at any MAPQ over all 10 segments, and 0 Kraken2 reads in all four chicken libraries -- despite the stated "ARV-infected" sample context |
| Chlamydia abortus | **noise** | database-error | 0 | 0% | 0 | 0% | - | 0 | 1 read at MAPQ 0, 0 survive |
| Chlamydia psittaci | **noise** | database-error | 0 | 0% | 0 | 0% | - | 0 | 3 reads at MAPQ 0, 0 survive |
| PhiX174 | **noise** | absent | 0 | 0% | 0 | 0% | - | - | 0 reads at any MAPQ -- this RNA run carried no PhiX spike-in, in contrast to the DNA library |

## SRR39197041_DNA — Aedes aegypti (mosquito)

*whole-genome DNA, 2x151 bp, flowcell H5F2GDSX3 lane 4; 120,788,540 reads; host claims 52,500,827 (43.5%) at MAPQ>=20.*

| organism | bin | likely source | reads | read fraction | depth | breadth | evenness | MAPQ surv. | why |
|---|---|---|---:|---:|---:|---:|---:|---:|---|
| PhiX174 (+ its duplicate inside UniVec_Core) † | **biological**<br>`auto: noise` | lab-organism | 608 | 0.000503% | 14.5 | 100% | 1 | - | Complete 5,386 bp genome at 100% breadth and ~14x once its duplicate copy inside UniVec_Core is counted back in; reassembled end-to-end with direct terminal repeats by MEGAHIT and rnaSPAdes independently and scored 100% complete by CheckV; missed entirely by Kraken2 (235/394 reads unclassified, 159 parked at root) and zeroed by competitive mapping, both for the same reason -- the sequence is in the database twice |
| Illumina PCR primer (UniVec NGB00360) | **consequential**<br>`auto: noise` | reagent | 4,016,896 | 3.33% | 999,846 | 100% | 1 | 0.58 | 4,016,896 reads, 3.33% of the library, on one 58 bp Illumina PCR primer entry at 100% breadth and ~1e6 x depth; the library was not adapter-trimmed before submission |
| Homo sapiens | **consequential**<br>`auto: noise` | human | 161,729 | 0.134% | 0.00734 | 0.473% | 0.65 | 0.62 | 161,729 reads (0.134%); observed breadth 0.473% is 65% of the Poisson expectation for 0.0073x, so the reads are spread genome-wide rather than piled on conserved loci; the same 0.08-0.10% floor appears in every library in this project, so handler/reagent -- not a bloodmeal and not index hopping |
| Leishmania major (trypanosomatid relative) | **trace**<br>`auto: noise` | real-associate | 11,537 | 0.00955% | 0.0391 | 2.13% | 0.56 | 0.58 | 11,537 reads survive MAPQ>=20 and are spread over 35 of the 36 chromosomes at a near-identical 2-3% breadth and 0.04-0.07x depth each -- the signature of a DIVERGENT relative, where only the conserved few per cent of each chromosome aligns, not of a spike artifact; Trypanosomatida sits at the same level in the co-lane control (4,073 vs 4,295), which is what a cobiont shared by the mosquito population looks like. The organism present is a trypanosomatid, most plausibly a monoxenous mosquito-infecting genus; Leishmania major is the database stand-in, not the identification |
| Acinetobacter schindleri | **trace**<br>`auto: noise` | reagent | 200 | 0.000166% | 0.00696 | 0.374% | 0.54 | 0.75 | 200 reads with 75% MAPQ survival and breadth at 54% of Poisson expectation -- credible but minute; Acinetobacter is a canonical extraction-kit organism |
| Cell fusing agent virus | **trace** | real-associate | 60 | 0.000050% | 0.572 | 6.52% | 0.15 | 0.91 | 60 reads, 91% of which survive MAPQ>=20, over 6.5% of the genome at 0.57x; CFAV is the canonical insect-specific flavivirus of Aedes aegypti, present both as free virus and as EVEs in the host genome |
| UniVec_Core, all other entries (vectors, other adapters) | **noise** | reagent | 60,571 | 0.0501% | - | - | - | - | the remainder of the UniVec_Core collection once the one entry that carries the bulk of the reads is counted separately; spread thinly over ~3,000 unrelated vector and adapter sequences |
| Lodderomyces elongisporus | **noise** | database-error | 26,496 | 0.0219% | 0.1 | 0.0219% | 0.0023 | 0.037 | 714,810 reads at MAPQ 0 collapse to 26,496 at MAPQ>=20, all inside 0.02% of the genome (2/1000 of Poisson expectation) -- rDNA and repeat cross-mapping |
| Aspergillus flavus | **noise** | database-error | 3,595 | 0.00298% | 0.00411 | 0.00660% | 0.016 | 0.17 | survival 0.17, breadth at 1.6% of expectation |
| Plasmodium relictum | **noise** | database-error | 1,109 | 0.000918% | 0.00403 | 0.0328% | 0.082 | 0.038 | survival 0.04 and breadth at 8% of expectation; Kraken2 reports it at the same level in the co-lane control (10,316 vs 8,291) |
| Saccharomyces cerevisiae | **noise** | lab-organism | 7 | 0.000006% | 3.46e-05 | 0.00197% | 0.57 | 0.01 | 7 reads at MAPQ>=20 out of 702; decoy working as intended |
| Escherichia coli K-12 | **noise** | lab-organism | 4 | 0.000003% | 8.88e-05 | 0.00444% | 0.5 | 0.57 | 4 reads at MAPQ>=20; decoy working as intended |
| Achromobacter xylosoxidans | **noise** | absent | 0 | 0% | 0 | 0% | - | 0 | 4 reads at MAPQ 0, 0 survive |
| Lactobacillus phage Semele | **noise** | absent | 0 | 0% | 0 | 0% | - | - | 0 reads at any MAPQ |
| Phasi Charoen-like virus | **noise** | absent | 0 | 0% | 0 | 0% | - | - | 0 reads on all three segments -- the commonest Aedes aegypti insect-specific virus is not present |
| Pseudomonas aeruginosa | **noise** | absent | 0 | 0% | 0 | 0% | - | 0 | 5 reads at MAPQ 0, 0 survive |
| Wolbachia pipientis | **noise** | absent | 0 | 0% | 0 | 0% | - | 0 | 0 reads at MAPQ>=20 and 4 at MAPQ 0; Kraken2 sees 0 here and 3 in the control -- this mosquito line is Wolbachia-free |

