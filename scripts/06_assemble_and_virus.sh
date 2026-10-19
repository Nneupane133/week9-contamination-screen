#!/bin/bash
#PBS -N w9_asm
#PBS -q large
#PBS -l select=1:ncpus=32:mem=180gb
#PBS -l walltime=48:00:00
#PBS -j oe
#PBS -o /home/aubnxn001/week9/logs/
# ---------------------------------------------------------------------------
# STEP 5: host depletion -> assembly -> virus calling -> completeness
#
# WHY ASSEMBLY IS REQUIRED AND READ MAPPING IS NOT ENOUGH (esp. for RNA viruses):
#   Everything up to here compares reads against references that already exist.
#   An RNA virus that is divergent from anything in RefSeq produces reads that
#   sendsketch and Kraken2 simply leave unclassified -- the commonest way a real
#   viral finding is missed. Assembling the non-host reads reconstructs the
#   genome from the data itself, so the detection no longer depends on the virus
#   already being in a database. Contigs also carry information single reads
#   cannot: genome length, ORF structure, terminal repeats, and whether viral
#   sequence sits on the same contig as host sequence (the signature of an
#   endogenous viral element rather than a replicating virus).
#
# WHY HOST-DEPLETE FIRST:
#   ~95 % of these libraries are host. Assembling that is an enormous waste of
#   memory and it buries low-coverage non-host contigs in the host graph.
#   Depletion is done against the host GENOME **and** the host TRANSCRIPTOME:
#   with the -ax sr short-read preset a read spanning an mRNA splice junction
#   aligns poorly to the genome and would survive depletion as a false
#   "non-host" read. The transcriptome catches exactly those.
#   Only pairs where NEITHER mate maps are carried forward, which is the
#   conservative choice -- a chimeric host/non-host pair stays in.
#
# WHY BOTH MEGAHIT AND rnaSPAdes:
#   They fail in different places, and after depletion the read set is small
#   enough that running both is cheap.
#     MEGAHIT   - succeeds-or-dies on very uneven coverage, low memory, handles
#                 the 0.01x-to-1000x dynamic range of a contamination screen.
#     rnaSPAdes - built for the non-uniform coverage of transcriptomes and
#                 reconstructs full-length RNA virus genomes markedly better.
#   A contig recovered by BOTH assemblers independently is much stronger
#   evidence than one recovered by either alone.
#
# WHY geNomad:
#   It scores each contig as chromosome / plasmid / virus from marker-gene
#   content rather than from nucleotide similarity to a known genome, so it can
#   flag a virus that has no close database relative. Crucially for this task,
#   geNomad also reports host-vs-viral gene content along a contig, which is the
#   direct evidence for calling something an integrated element (viral genes
#   embedded in host-gene context) instead of a free virus.
#
# WHY CheckV:
#   It estimates how complete each viral contig is against reference genomes.
#   Completeness is what separates "a viral fragment is present" (trace) from
#   "a near-complete viral genome is present" (biological), so it feeds the bin
#   assignment directly. CheckV also flags host-contamination within a contig,
#   a second independent read on the EVE question.
#
# WHY DIAMOND blastx AS WELL:
#   A protein-level search against RefSeq viral proteins is independent of
#   geNomad's marker set. Protein similarity persists across far greater
#   nucleotide divergence, so it is the most sensitive way to identify a
#   divergent virus and to put a name on a contig geNomad only scored.
# ---------------------------------------------------------------------------
set +e; source /apps/profiles/modules_asax.sh.dyn >/dev/null 2>&1
module load bbmap samtools/1.18 minimap2 seqkit megahit_voutcn spades diamond >/dev/null 2>&1
set -euo pipefail

WK=/home/aubnxn001/week9
LEVEL=${LEVEL:-test}
SAMPLE=${SAMPLE:-BYJ10_RNA}
THREADS=${THREADS:-32}
MEMGB=${MEMGB:-170}
OUT=$WK/results/assembly/$LEVEL/$SAMPLE
mkdir -p "$OUT"
export PATH="$WK/envs/w9virus/bin:$PATH"

case "$SAMPLE" in
  BYJ10_RNA)        HOSTREF=$WK/references/host/Gallus_gallus_GRCg7b.fna.gz
                    HOSTTX=$WK/references/host/Gallus_gallus_GRCg7b_rna.fna.gz
                    if [ "$LEVEL" = test ]; then R1=$WK/data/test_subsample/BYJ10_test_1.fq.gz; R2=$WK/data/test_subsample/BYJ10_test_2.fq.gz
                    else R1=$WK/data/RNAseq/BYJ10_1.fq.gz; R2=$WK/data/RNAseq/BYJ10_2.fq.gz; fi ;;
  SRR39197041_DNA)  HOSTREF=$WK/references/host/Aedes_aegypti_AaegL5.fna.gz
                    HOSTTX=$WK/references/host/Aedes_aegypti_AaegL5_rna.fna.gz
                    if [ "$LEVEL" = test ]; then R1=$WK/data/test_subsample/SRR39197041_test_1.fq.gz; R2=$WK/data/test_subsample/SRR39197041_test_2.fq.gz
                    else R1=$WK/data/DNA_WGS/SRR39197041_1.fastq.gz; R2=$WK/data/DNA_WGS/SRR39197041_2.fastq.gz; fi ;;
  *) echo "unknown SAMPLE $SAMPLE" >&2; exit 1 ;;
esac

# ---- 1. host depletion ---------------------------------------------------
NH1=$OUT/nonhost_1.fq.gz; NH2=$OUT/nonhost_2.fq.gz
if [ ! -s "$NH1" ]; then
  DEPREF=$OUT/hostdepl.fna
  if [ ! -s "$DEPREF" ]; then
    zcat "$HOSTREF" > "$DEPREF"
    [ -s "$HOSTTX" ] && zcat "$HOSTTX" >> "$DEPREF"
  fi
  echo "[$(date +%T)] host depletion ($SAMPLE)"
  # -f 12 : both mates unmapped;  -F 256 : drop secondary alignments
  minimap2 -ax sr -t "$THREADS" -I 12G --secondary=no "$DEPREF" "$R1" "$R2" \
      2> "$OUT/hostdepl.minimap2.log" \
    | samtools view -@ 4 -b -f 12 -F 256 - \
    | samtools sort -@ 8 -m 2G -n -o "$OUT/nonhost.bam" -
  samtools fastq -@ 8 -1 "$NH1" -2 "$NH2" -0 /dev/null -s /dev/null -n "$OUT/nonhost.bam"
  rm -f "$OUT/nonhost.bam" "$DEPREF"
fi
seqkit stats "$NH1" "$NH2" | tee "$OUT/nonhost.stats.txt"

# ---- 2. assembly ---------------------------------------------------------
if [ ! -s "$OUT/megahit/final.contigs.fa" ]; then
  echo "[$(date +%T)] MEGAHIT"
  rm -rf "$OUT/megahit"
  megahit -1 "$NH1" -2 "$NH2" -o "$OUT/megahit" -t "$THREADS" \
          --min-contig-len 300 -m $((MEMGB*1000000000)) > "$OUT/megahit.log" 2>&1 \
    || echo "  WARN: MEGAHIT failed"
fi
if [ ! -s "$OUT/rnaspades/transcripts.fasta" ]; then
  echo "[$(date +%T)] rnaSPAdes"
  rnaspades.py -1 "$NH1" -2 "$NH2" -o "$OUT/rnaspades" \
          -t "$THREADS" -m "$MEMGB" > "$OUT/rnaspades.log" 2>&1 \
    || echo "  WARN: rnaSPAdes failed"
fi

# Pool both assemblies, tagging each contig with its assembler so the two can
# be told apart (and cross-confirmed) downstream.
CONTIGS=$OUT/contigs.pooled.fa
: > "$CONTIGS"
[ -s "$OUT/megahit/final.contigs.fa" ] && \
  seqkit replace -p '^' -r 'megahit_' "$OUT/megahit/final.contigs.fa" >> "$CONTIGS"
[ -s "$OUT/rnaspades/transcripts.fasta" ] && \
  seqkit seq -m 300 "$OUT/rnaspades/transcripts.fasta" | seqkit replace -p '^' -r 'rnaspades_' >> "$CONTIGS"
seqkit stats "$CONTIGS" | tee "$OUT/contigs.stats.txt"

# ---- 3. virus calling ----------------------------------------------------
echo "[$(date +%T)] geNomad"
genomad end-to-end --cleanup --threads "$THREADS" \
        "$CONTIGS" "$OUT/genomad" "$WK/references/db/genomad_db" \
        > "$OUT/genomad.log" 2>&1 || echo "  WARN: geNomad failed"

echo "[$(date +%T)] CheckV"
VIR=$OUT/genomad/contigs.pooled_summary/contigs.pooled_virus.fna
if [ -s "$VIR" ]; then
  checkv end_to_end "$VIR" "$OUT/checkv" -t "$THREADS" \
         -d "$WK/references/db/checkv/checkv-db-v1.5" \
         > "$OUT/checkv.log" 2>&1 || echo "  WARN: CheckV failed"
else
  echo "  no geNomad virus contigs; running CheckV on all contigs instead"
  checkv end_to_end "$CONTIGS" "$OUT/checkv" -t "$THREADS" \
         -d "$WK/references/db/checkv/checkv-db-v1.5" \
         > "$OUT/checkv.log" 2>&1 || echo "  WARN: CheckV failed"
fi

echo "[$(date +%T)] DIAMOND blastx vs RefSeq viral proteins"
diamond blastx -q "$CONTIGS" -d "$WK/references/db/refseq_viral_protein.dmnd" \
        -o "$OUT/contigs.viralprot.tsv" -p "$THREADS" \
        --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
        --max-target-seqs 5 --evalue 1e-5 --sensitive \
        > "$OUT/diamond.log" 2>&1 || echo "  WARN: DIAMOND failed"

echo "[$(date +%T)] $SAMPLE assembly/virus arm complete"
