#!/bin/bash
#PBS -N w9_kraken2
#PBS -q large
#PBS -l select=1:ncpus=16:mem=48gb
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -o /home/aubnxn001/week9/logs/
# ---------------------------------------------------------------------------
# STEP 2 (screen B): Kraken2 + Bracken -- exact k-mer classification, per read
#
# WHY Kraken2 AS A SECOND, INDEPENDENT SCREEN:
#   sendsketch answers "is this genome present?" from a MinHash sketch. It
#   cannot give a read fraction, and a sketch hit of 4 shared hashes is
#   indistinguishable from a 4-read artifact. Kraken2 classifies EVERY read
#   against an exact k-mer index, so it yields countable per-taxon read numbers
#   -- the quantity the bins (noise / trace / consequential / biological)
#   actually need. The two methods also fail differently: sendsketch queries
#   the complete uncapped RefSeq/nt sketch set (broad, no local DB), Kraken2
#   queries a 16 GB capped index (narrower, but quantitative). Agreement
#   between them is real evidence; disagreement flags a candidate as shaky.
#
# WHY --confidence 0.1:
#   At the default 0 a read is assigned if ANY k-mer matches, which produces a
#   long tail of single-k-mer false positives -- exactly the "noise" bin we are
#   trying to keep out. 0.1 requires 10 % of a read's k-mers to support the
#   assignment. It costs a little sensitivity, which is acceptable because
#   sendsketch is already covering the sensitive end.
#
# WHY --report-minimizer-data (the key anti-false-positive trick):
#   It adds the count of DISTINCT minimizers supporting each taxon. A genuinely
#   present genome is hit across many different minimizers; a spurious taxon is
#   hit by a handful of minimizers over and over (one repeat, one conserved
#   rRNA). The ratio taxon_reads / distinct_minimizers separates these two
#   cases, and nothing in the plain report can.
#
# WHY Bracken afterwards:
#   Kraken2 assigns each read to the LCA of its k-mer hits, so reads from a
#   species with close DB relatives pile up at genus/family and the species row
#   under-counts. Bracken redistributes them to species using the DB's own
#   k-mer distribution, giving the honest species-level read fraction.
# ---------------------------------------------------------------------------
set +e; source /apps/profiles/modules_asax.sh.dyn >/dev/null 2>&1; set -e

WK=/home/aubnxn001/week9
LEVEL=${LEVEL:-${1:-test}}          # test | full
THREADS=${THREADS:-16}
DB=$WK/references/db/k2_pluspf_16gb
OUT=$WK/results/kraken2/$LEVEL
mkdir -p "$OUT"

export PATH="$WK/envs/w9tax/bin:$PATH"   # kraken2 2.1.3 + bracken 3.0

if [ "$LEVEL" = "test" ]; then
  S1=$WK/data/test_subsample; RNA1=$S1/BYJ10_test_1.fq.gz; RNA2=$S1/BYJ10_test_2.fq.gz
  DNA1=$S1/SRR39197041_test_1.fq.gz;  DNA2=$S1/SRR39197041_test_2.fq.gz
else
  RNA1=$WK/data/RNAseq/BYJ10_1.fq.gz;         RNA2=$WK/data/RNAseq/BYJ10_2.fq.gz
  DNA1=$WK/data/DNA_WGS/SRR39197041_1.fastq.gz; DNA2=$WK/data/DNA_WGS/SRR39197041_2.fastq.gz
fi

run_one () {
  local NAME=$1 R1=$2 R2=$3
  echo "[$(date +%T)] kraken2 $NAME ($LEVEL)"
  kraken2 --db "$DB" --threads "$THREADS" --paired \
          --confidence 0.1 --report-minimizer-data \
          --report "$OUT/${NAME}.minimizer.report" \
          --output "$OUT/${NAME}.kraken2.out" \
          "$R1" "$R2" 2> "$OUT/${NAME}.kraken2.log"

  # Strip the two minimizer columns to get a report Bracken can parse.
  cut -f1,2,3,6,7,8 "$OUT/${NAME}.minimizer.report" > "$OUT/${NAME}.report"

  # Bracken at species (S) level. -r 150 matches the 150 bp read length of both
  # libraries; -t 10 requires >=10 reads before a species is re-estimated at
  # all, which suppresses the single-read tail.
  echo "[$(date +%T)] bracken $NAME"
  bracken -d "$DB" -i "$OUT/${NAME}.report" \
          -o "$OUT/${NAME}.bracken.species.tsv" \
          -w "$OUT/${NAME}.bracken.report" \
          -r 150 -l S -t 10 > "$OUT/${NAME}.bracken.log" 2>&1 || \
          echo "  WARN: bracken failed for $NAME (see ${NAME}.bracken.log)"

  gzip -f "$OUT/${NAME}.kraken2.out" &
}

run_one BYJ10_RNA       "$RNA1" "$RNA2"
run_one SRR39197041_DNA "$DNA1" "$DNA2"
wait
echo "[$(date +%T)] kraken2/bracken $LEVEL complete"
