#!/bin/bash
#PBS -N w9_evetest
#PBS -q large
#PBS -l select=1:ncpus=16:mem=64gb
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -o /home/aubnxn001/week9/logs/
# ---------------------------------------------------------------------------
# STEP 3c: the host-competition (EVE) test.
#
# WHY THIS STEP EXISTS:
#   In the competitive BAM a candidate's reads can carry low MAPQ for two
#   completely different reasons, and the significance bins get the answer
#   backwards if we cannot tell them apart:
#
#     (a) the reads align about equally well to some UNRELATED reference
#         -> cross-mapping artifact -> genuinely noise; or
#     (b) the reads align about equally well to the HOST GENOME, because the
#         sequence is physically integrated in the host (an endogenous viral
#         element, an endosymbiont-derived insertion, a nuclear mitochondrial
#         insert) -> not an artifact at all, but one of the most interesting
#         findings available, and precisely the "integrated in host genome"
#         source category.
#
#   MAPQ alone cannot distinguish them, because MAPQ only says "a near-equal
#   alignment exists somewhere", never where. This step asks where.
#
# HOW:
#   1. Pull every read that aligns to a candidate in the competitive BAM, at
#      ANY MAPQ (the low-MAPQ ones are exactly the ones in question).
#   2. Re-map just those reads to the SAME combined reference, but with
#      secondary alignments switched ON (-N 10 --secondary=yes -p 0.5), so
#      every near-equal placement is recorded rather than silently discarded.
#   3. For each read, compare its best alignment score (AS) on the candidate
#      with its best AS on the host. Reads scoring >= as well on the host are
#      "host-competing".
#
#   A candidate whose reads are overwhelmingly host-competing is integrated in
#   (or indistinguishable from) the host genome. A candidate whose reads are
#   NOT host-competing is genuinely exogenous, and its low MAPQ (if any) must
#   be explained by something else.
#
#   The re-mapping reference deliberately ADDS THE HOST TRANSCRIPTOME, which the
#   main competitive reference does not contain. Reason: several viral reference
#   genomes carry a transduced HOST gene (Y73 sarcoma virus carries the host
#   src oncogene; many acutely transforming retroviruses are built this way).
#   That transduced copy is an intronless cDNA, so an RNA-seq read spanning an
#   exon-exon junction aligns BETTER to the "virus" than to the host genomic
#   locus, and the host genome alone loses the competition. Adding the host
#   transcriptome restores it, and reads that are really host mRNA are then
#   correctly scored as host-competing instead of being reported as a virus.
#   The transcriptome is used ONLY here, not in the main measurement, because
#   genome+transcriptome are redundant with each other and would split the
#   host's own reads and destroy its MAPQ.
# ---------------------------------------------------------------------------
set +e; source /apps/profiles/modules_asax.sh.dyn >/dev/null 2>&1
module load samtools/1.18 minimap2 >/dev/null 2>&1; set -euo pipefail

WK=/home/aubnxn001/week9
LEVEL=${LEVEL:-test}
THREADS=${THREADS:-16}
COV=$WK/results/coverage/$LEVEL
OUT=$WK/results/host_competition/$LEVEL
mkdir -p "$OUT"

run_one () {
  local SAMP=$1 HOSTPFX=$2 HOSTTX=$3
  local BAM=$COV/${SAMP}.sorted.bam
  [ -s "$BAM" ] || { echo "no BAM for $SAMP"; return; }

  # Competition reference = the competitive reference + the host transcriptome,
  # with transcripts carried under the host's own prefix so they count as host.
  local CREF=$WK/references/combined/${SAMP}.competition.fna
  local MMI=$WK/references/combined/${SAMP}.competition.mmi
  if [ ! -s "$MMI" ]; then
    if [ ! -s "$CREF" ]; then
      cp "$WK/references/combined/${SAMP}.fna" "$CREF"
      if [ -s "$HOSTTX" ]; then
        zcat "$HOSTTX" | awk -v n="$HOSTPFX" '/^>/{sub(/^>/,"");print ">" n "::TX_" $0; next}{print}' >> "$CREF"
      fi
    fi
    echo "[$(date +%T)] building competition index for $SAMP (genome + transcriptome + candidates)"
    minimap2 -x sr -t "$THREADS" -I 12G -d "$MMI" "$CREF" 2>> "$OUT/${SAMP}.index.log"
  fi

  echo "[$(date +%T)] $SAMP: collecting reads that touch any non-host reference"
  # Reference sequences belonging to candidates/decoys (everything but the host).
  samtools idxstats "$BAM" | awk -F'\t' -v h="$HOSTPFX" '$1!="*" && $1 !~ "^"h"::" && $3>0 {print $1}' \
    > "$OUT/${SAMP}.nonhost_refs.txt"
  [ -s "$OUT/${SAMP}.nonhost_refs.txt" ] || { echo "  none"; return; }

  # Collect those reads (both mates) by name, then pull them out of the BAM.
  samtools view -@ 4 "$BAM" $(tr '\n' ' ' < "$OUT/${SAMP}.nonhost_refs.txt") \
    | cut -f1 | sort -u > "$OUT/${SAMP}.readnames.txt"
  echo "  $(wc -l < "$OUT/${SAMP}.readnames.txt") read names"

  samtools view -@ 4 -N "$OUT/${SAMP}.readnames.txt" -b "$BAM" \
    | samtools sort -@ 4 -n -o "$OUT/${SAMP}.subset.bam" -
  samtools fastq -@ 4 -1 "$OUT/${SAMP}.sub_1.fq" -2 "$OUT/${SAMP}.sub_2.fq" \
                 -0 /dev/null -s "$OUT/${SAMP}.sub_s.fq" -n "$OUT/${SAMP}.subset.bam" 2>/dev/null

  echo "[$(date +%T)] $SAMP: re-mapping that subset with secondary alignments ON"
  # -N 10 keeps up to 10 placements, -p 0.5 reports any whose score is within
  # 50% of the best. Without these the host placement we need to see is dropped.
  minimap2 -ax sr -t "$THREADS" -I 12G -N 10 -p 0.5 --secondary=yes \
      "$MMI" "$OUT/${SAMP}.sub_1.fq" "$OUT/${SAMP}.sub_2.fq" \
      2> "$OUT/${SAMP}.remap.log" \
    | samtools view -@ 4 -b -o "$OUT/${SAMP}.remap.bam" -

  echo "[$(date +%T)] $SAMP: scoring host competition"
  # For every alignment record: read name, flag(mate), reference source, AS.
  samtools view "$OUT/${SAMP}.remap.bam" \
    | awk -F'\t' -v h="$HOSTPFX" '
        $3=="*" { next }                       # unmapped record, nothing to compare
        {
          as=0; for(i=12;i<=NF;i++) if ($i ~ /^AS:i:/) { split($i,a,":"); as=a[3] }
          split($3, r, "::"); src=r[1]
          if (src=="") next
          mate = and($2,64) ? 1 : 2
          key = $1 "/" mate
          if (src==h) { if (as>hostbest[key]) hostbest[key]=as }
          else        { if (as>candbest[key]) { candbest[key]=as; candsrc[key]=src } }
        }
        END{
          for (k in candbest) {
            s=candsrc[k]; tot[s]++
            if (k in hostbest && hostbest[k] >= candbest[k]) hostwin[s]++
          }
          for (s in tot) printf "%s\t%d\t%d\t%.4f\n", s, tot[s], hostwin[s]+0, (hostwin[s]+0)/tot[s]
        }' | sort -k2,2 -rn \
    | cat <(printf "source\treads\thost_competing\thost_competing_frac\n") - \
    > "$OUT/${SAMP}.host_competition.tsv"
  cat "$OUT/${SAMP}.host_competition.tsv"
  rm -f "$OUT/${SAMP}.sub_"*.fq "$OUT/${SAMP}.subset.bam"
}

run_one BYJ10_RNA       Gallus_gallus_GRCg7b  "$WK/references/host/Gallus_gallus_GRCg7b_rna.fna.gz"
run_one SRR39197041_DNA Aedes_aegypti_AaegL5  "$WK/references/host/Aedes_aegypti_AaegL5_rna.fna.gz"
echo "[$(date +%T)] host-competition test complete"
