#!/bin/bash
# ---------------------------------------------------------------------------
# STEP 1 (screen A): BBTools sendsketch.sh -- fast, database-free taxonomic screen
#
# WHY sendsketch FIRST:
#   It builds a MinHash sketch of the reads locally and sends only that sketch
#   (a few hundred kB) to JGI's server, which compares it against the full
#   RefSeq / nt sketch collections. So we get a genome-wide, *uncapped* answer
#   in seconds with zero local database. That makes it the right instrument for
#   HYPOTHESIS GENERATION -- "what organisms might be in here at all?"
#
# WHY IT IS NOT ENOUGH ON ITS OWN:
#   A sketch comparison reports presence/ANI, not read counts, and it cannot
#   tell a real 0.1 % cobiont from a 2-read database artifact. That is why
#   every hit it raises is carried forward to Kraken2/Bracken (step 2) for a
#   second independent opinion and then to alignment (step 3) for the actual
#   read fraction / depth / breadth measurements.
#
#   Two databases are queried because they fail differently:
#     refseq -- curated, clean taxonomy, but misses unusual/unsubmitted taxa
#     nt     -- vastly broader (catches odd viruses, EVEs, vectors), but is
#               full of mislabelled entries, so an nt-only hit is weak evidence
#   A hit present in BOTH is much stronger than a hit in either alone.
# ---------------------------------------------------------------------------
set -euo pipefail
source /home/aubnxn001/week9/scripts/env.sh
LEVEL=${1:-test}           # test | full
OUT=$WK/results/sendsketch/$LEVEL
mkdir -p "$OUT"

if [ "$LEVEL" = "test" ]; then
  declare -A R1=( [BYJ10_RNA]=$WK/data/test_subsample/BYJ10_test_1.fq.gz \
                  [SRR39197041_DNA]=$WK/data/test_subsample/SRR39197041_test_1.fq.gz )
  declare -A R2=( [BYJ10_RNA]=$WK/data/test_subsample/BYJ10_test_2.fq.gz \
                  [SRR39197041_DNA]=$WK/data/test_subsample/SRR39197041_test_2.fq.gz )
else
  declare -A R1=( [BYJ10_RNA]=$WK/data/RNAseq/BYJ10_1.fq.gz \
                  [SRR39197041_DNA]=$WK/data/DNA_WGS/SRR39197041_1.fastq.gz )
  declare -A R2=( [BYJ10_RNA]=$WK/data/RNAseq/BYJ10_2.fq.gz \
                  [SRR39197041_DNA]=$WK/data/DNA_WGS/SRR39197041_2.fastq.gz )
fi

for S in "${!R1[@]}"; do
  for DB in refseq nt; do
    echo "[$(date +%T)] sendsketch $S vs $DB"
    # Both mates are streamed into one sketch via stdin: sendsketch.sh ignores
    # in2= / comma-lists and would otherwise silently sketch R1 only.
    # records=100 keeps a deep tail so low-abundance candidates stay visible;
    # printall disables the default display cutoff.
    cat "${R1[$S]}" "${R2[$S]}" \
      | sendsketch.sh in=stdin.fq.gz out="$OUT/${S}.${DB}.tsv" \
          name="$S" $DB records=100 printall ow=t \
          > "$OUT/${S}.${DB}.stderr" 2>&1 || echo "  WARN: $S/$DB failed"
  done
done
echo "[$(date +%T)] sendsketch $LEVEL complete"
