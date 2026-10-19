#!/bin/bash
#PBS -N w9_map
#PBS -q large
#PBS -l select=1:ncpus=32:mem=120gb
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -o /home/aubnxn001/week9/logs/
# ---------------------------------------------------------------------------
# STEP 3b: COMPETITIVE mapping + per-candidate depth / breadth / read fraction
#
# WHY minimap2 AND NOT bowtie2:
#   The combined reference here is several Gb (host + human decoy + candidates).
#   minimap2 indexes that in minutes where bowtie2-build takes many hours, and
#   the ASC bowtie2 module is 2.2.9 (2016). minimap2's -ax sr preset is a
#   short-read mode that is as accurate as bowtie2 for 150 bp pairs, and it
#   tolerates DIVERGENT alignment better -- which matters a lot here, because
#   our "candidate" reference is usually a database relative of the organism
#   actually present, not the exact strain. bowtie2 end-to-end would simply
#   drop those reads and we would under-call a real detection.
#
# WHY ONE COMBINED REFERENCE (competitive) RATHER THAN ONE BAM PER CANDIDATE:
#   Mapping to a candidate in isolation is the single commonest way to
#   manufacture a false contamination result: conserved rRNA / ribosomal
#   protein / transposon reads from the host and from everything else in the
#   library have nowhere else to go, so they pile onto the candidate's
#   conserved loci and produce high apparent depth over a few narrow windows.
#   In a competitive reference every read is assigned to its single best match,
#   so a read is only credited to a candidate when that candidate explains it
#   better than the host, the human decoy and every lab-organism decoy do.
#
# WHY MAPQ >= 20:
#   MAPQ encodes "how much better is the best alignment than the second best".
#   Requiring >= 20 discards reads that map about equally well to the candidate
#   and to something else in the reference -- precisely the ambiguous reads that
#   would otherwise inflate breadth. This is the filter that makes the breadth
#   number mean "sequence unique to this organism is present".
#
# WHY samtools coverage FOR THE MEASUREMENT:
#   It reports, per reference sequence and in one pass, exactly the three
#   quantities the assignment asks for: numreads (-> read fraction), covbases /
#   length (-> breadth) and meandepth (-> depth). Computing them from samtools
#   depth would need a full per-base dump of a multi-Gb reference.
#
# WHY WE ALSO DUMP PER-BASE DEPTH FOR THE SMALL REFERENCES:
#   Breadth and mean depth cannot tell an even, genome-wide 0.5x from a spike
#   over 2 % of the genome at 25x -- but those two mean completely different
#   things (real organism vs conserved-region cross-mapping / database error).
#   So for the non-host references we keep per-base depth and compute an
#   evenness statistic in step 5.
# ---------------------------------------------------------------------------
set +e; source /apps/profiles/modules_asax.sh.dyn >/dev/null 2>&1
module load bbmap samtools/1.18 minimap2 seqkit >/dev/null 2>&1; set -euo pipefail

WK=/home/aubnxn001/week9
LEVEL=${LEVEL:-test}
SAMPLE=${SAMPLE:-}            # BYJ10_RNA | SRR39197041_DNA  (both if empty)
THREADS=${THREADS:-32}
REFDIR=$WK/references
OUT=$WK/results/coverage/$LEVEL
mkdir -p "$OUT" "$REFDIR/combined"

# --- which reference pieces belong to which readset -----------------------
# Column 2 of the manifest marks role; column 5 (optional) restricts a piece to
# one readset, so the chicken host is not in the mosquito reference and vice versa.
build_combined () {
  local SAMP=$1 COMB=$REFDIR/combined/${SAMP}.fna
  if [ -s "$COMB" ]; then echo "  combined reference already built: $COMB"; return; fi
  echo "[$(date +%T)] building combined reference for $SAMP"
  : > "$COMB"
  while IFS=$'\t' read -r name role taxid source id restrict; do
    case "$name" in ''|'#'*) continue;; esac
    # restrict: which readset this reference piece belongs in. Written as an
    # explicit if, not a && chain -- under set -e a failing && chain would end
    # the script rather than skip the row.
    if [ -n "${restrict:-}" ] && [ "$restrict" != "both" ] && [ "$restrict" != "$SAMP" ]; then
      continue
    fi
    local f
    if   [ -s "$REFDIR/host/${name}.fna.gz" ];       then f=$REFDIR/host/${name}.fna.gz
    elif [ -s "$REFDIR/candidates/${name}.fna.gz" ]; then f=$REFDIR/candidates/${name}.fna.gz
    else echo "  WARN: no fasta for $name"; continue; fi
    # Prefix every header with the source name so alignments can be attributed
    # back to an organism regardless of how many contigs that organism has.
    zcat "$f" | awk -v n="$name" '/^>/{sub(/^>/,""); print ">" n "::" $0; next} {print}' >> "$COMB"
    echo "  + $name ($role)"
  done < "$WK/metadata/candidates.tsv"
  seqkit stats "$COMB" | sed 's/^/  /'
}

map_one () {
  local SAMP=$1 R1=$2 R2=$3
  local COMB=$REFDIR/combined/${SAMP}.fna
  local BAM=$OUT/${SAMP}.sorted.bam
  build_combined "$SAMP"

  # Build the minimap2 index ONCE and cache it: the combined reference is ~4 Gb
  # and the test-level and full-level runs would otherwise rebuild it twice.
  # -I 12G forces a SINGLE-PART index. With a multi-part index minimap2 maps
  # against each part independently and MAPQ is computed within a part, which
  # silently destroys the competitive comparison this whole design rests on.
  local MMI=$REFDIR/combined/${SAMP}.sr.mmi
  if [ ! -s "$MMI" ]; then
    echo "[$(date +%T)] building minimap2 index for $SAMP (single part, -I 12G)"
    minimap2 -x sr -t "$THREADS" -I 12G -d "$MMI" "$COMB" 2>> "$OUT/${SAMP}.index.log"
  fi

  if [ ! -s "$BAM" ]; then
    echo "[$(date +%T)] minimap2 $SAMP -> $(basename "$MMI")"
    minimap2 -ax sr -t "$THREADS" -I 12G --secondary=no "$MMI" "$R1" "$R2" \
      2> "$OUT/${SAMP}.minimap2.log" \
      | samtools sort -@ 8 -m 2G -o "$BAM" -
    samtools index -@ 8 "$BAM"
  fi

  echo "[$(date +%T)] samtools stats/coverage $SAMP"
  samtools flagstat -@ 8 "$BAM" > "$OUT/${SAMP}.flagstat.txt"
  samtools idxstats "$BAM"      > "$OUT/${SAMP}.idxstats.tsv"
  # Unfiltered and MAPQ>=20 coverage, so the report can show how much of a
  # detection survives the uniqueness filter -- that difference is itself
  # diagnostic (a database-error hit loses almost everything at MAPQ>=20).
  samtools coverage               "$BAM" > "$OUT/${SAMP}.coverage.mq0.tsv"
  samtools coverage --min-MQ 20   "$BAM" > "$OUT/${SAMP}.coverage.mq20.tsv"

  # Per-base depth for the non-host references only (host is multi-Gb).
  echo "[$(date +%T)] per-base depth for non-host references, $SAMP"
  # idxstats columns: ref, length, mapped, unmapped. Keep every non-host
  # reference sequence (regardless of whether it got reads) so that breadth is
  # computed over the FULL candidate genome, not only its covered contigs.
  awk -F'\t' '$1!="*" && $1 !~ /^(Gallus_gallus_GRCg7b|Aedes_aegypti_AaegL5|Homo_sapiens)::/ {print $1"\t0\t"$2}' \
      "$OUT/${SAMP}.idxstats.tsv" > "$OUT/${SAMP}.small_refs.bed"
  if [ -s "$OUT/${SAMP}.small_refs.bed" ]; then
    samtools depth -a -Q 20 -@ 8 -b "$OUT/${SAMP}.small_refs.bed" "$BAM" \
      | gzip > "$OUT/${SAMP}.depth.nonhost.tsv.gz"
  fi
  echo "[$(date +%T)] $SAMP done"
}

if [ "$LEVEL" = "test" ]; then
  D=$WK/data/test_subsample
  RNA1=$D/BYJ10_test_1.fq.gz; RNA2=$D/BYJ10_test_2.fq.gz
  DNA1=$D/SRR39197041_test_1.fq.gz; DNA2=$D/SRR39197041_test_2.fq.gz
else
  RNA1=$WK/data/RNAseq/BYJ10_1.fq.gz; RNA2=$WK/data/RNAseq/BYJ10_2.fq.gz
  DNA1=$WK/data/DNA_WGS/SRR39197041_1.fastq.gz; DNA2=$WK/data/DNA_WGS/SRR39197041_2.fastq.gz
fi

if [ -z "$SAMPLE" ] || [ "$SAMPLE" = BYJ10_RNA ]; then
  map_one BYJ10_RNA "$RNA1" "$RNA2"
fi
if [ -z "$SAMPLE" ] || [ "$SAMPLE" = SRR39197041_DNA ]; then
  map_one SRR39197041_DNA "$DNA1" "$DNA2"
fi
echo "[$(date +%T)] competitive mapping $LEVEL complete"
