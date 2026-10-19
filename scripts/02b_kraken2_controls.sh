#!/bin/bash
#PBS -N w9_k2_ctrl
#PBS -q large
#PBS -l select=1:ncpus=16:mem=48gb
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -o /home/aubnxn001/week9/logs/
# ---------------------------------------------------------------------------
# STEP 2b: screen the CO-SEQUENCED samples as controls.
#
# WHY THIS EXISTS -- it is what makes source attribution falsifiable:
#   Read headers show the batch structure of this data set.
#     BYJ10  + AYJ9          : flowcell 23C55KLT4, run 534, LANE 5  (same lane)
#     SYJ7   + SYJ10         : flowcell 237VF5LT4, run 167, lane 7  (different flowcell)
#     SRR39197041 + ...043   : flowcell H5F2GDSX3, lane 4           (same lane)
#
#   Those three groupings separate the candidate source categories that read
#   fraction / depth / breadth alone CANNOT separate:
#
#   * bleed-through (index hopping) can only come from a library on the SAME
#     LANE. If a taxon is abundant in AYJ9 and present in BYJ10 at roughly
#     0.1-1 % of AYJ9's level, index hopping is the leading explanation. If the
#     taxon is absent from AYJ9, index hopping is ruled out -- which is a much
#     stronger statement than "it looks low abundance".
#
#   * reagent contamination tracks the KIT, not the flowcell. A taxon that
#     appears in BYJ10, AYJ9, SYJ7 and SYJ10 -- across two different flowcells
#     -- but is unrelated to chicken biology, behaves like a reagent organism.
#     One confined to a single flowcell does not.
#
#   * a real biological associate should track the SAMPLE, not the batch: present
#     in the samples whose tissue/host justifies it and absent from the others,
#     at a level unrelated to which lane it was sequenced on.
#
#   Without these controls every source call would rest on plausibility
#   argument alone. With them, each call makes a prediction that the control
#   samples either confirm or refute.
# ---------------------------------------------------------------------------
set +e; source /apps/profiles/modules_asax.sh.dyn >/dev/null 2>&1; set -e
WK=/home/aubnxn001/week9
THREADS=${THREADS:-16}
DB=$WK/references/db/k2_pluspf_16gb
OUT=$WK/results/kraken2/controls
mkdir -p "$OUT"
export PATH="$WK/envs/w9tax/bin:$PATH"

declare -A R1=(
  [AYJ9_RNA]=$WK/data/RNAseq/AYJ9_1.fq.gz
  [SYJ7_RNA]=$WK/data/RNAseq/SYJ7_1.fq.gz
  [SYJ10_RNA]=$WK/data/RNAseq/SYJ10_1.fq.gz
  [SRR39197043_DNA]=$WK/data/DNA_WGS/SRR39197043_1.fastq.gz
)
declare -A R2=(
  [AYJ9_RNA]=$WK/data/RNAseq/AYJ9_2.fq.gz
  [SYJ7_RNA]=$WK/data/RNAseq/SYJ7_2.fq.gz
  [SYJ10_RNA]=$WK/data/RNAseq/SYJ10_2.fq.gz
  [SRR39197043_DNA]=$WK/data/DNA_WGS/SRR39197043_2.fastq.gz
)

for S in AYJ9_RNA SYJ7_RNA SYJ10_RNA SRR39197043_DNA; do
  [ -s "${R1[$S]}" ] || { echo "skip $S (missing fastq)"; continue; }
  [ -s "$OUT/${S}.report" ] && { echo "have $S"; continue; }
  echo "[$(date +%T)] kraken2 $S"
  kraken2 --db "$DB" --threads "$THREADS" --paired \
          --confidence 0.1 --report-minimizer-data \
          --report "$OUT/${S}.minimizer.report" --output /dev/null \
          "${R1[$S]}" "${R2[$S]}" 2> "$OUT/${S}.kraken2.log"
  cut -f1,2,3,6,7,8 "$OUT/${S}.minimizer.report" > "$OUT/${S}.report"
  bracken -d "$DB" -i "$OUT/${S}.report" -o "$OUT/${S}.bracken.species.tsv" \
          -w "$OUT/${S}.bracken.report" -r 150 -l S -t 10 \
          > "$OUT/${S}.bracken.log" 2>&1 || echo "  WARN bracken $S"
done
echo "[$(date +%T)] control screens complete"
