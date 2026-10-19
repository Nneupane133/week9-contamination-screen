#!/bin/bash
# Download the prebuilt Kraken2 PlusPF-16GB index (Langmead lab, genome-idx S3).
# WHY a prebuilt index: building Standard/PlusPF from scratch needs ~500 GB of
# downloaded genomes and many hours; the prebuilt index is the identical product
# and ships the Bracken *.kmer_distrib files we need for abundance re-estimation.
# WHY PlusPF-16: bacteria+archaea+viral+plasmid+human+UniVec_Core+protozoa+fungi.
# Fungi/protozoa matter here because reagent and lab-organism contaminants are
# frequently fungal, and the 16 GB cap keeps it runnable on an ordinary node.
set -euo pipefail
DB=/home/aubnxn001/week9/references/db/k2_pluspf_16gb
BUILD=20250402
mkdir -p "$DB"
cd "$DB"
URL="https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16gb_${BUILD}.tar.gz"
echo "[$(date)] downloading $URL"
curl -fL --no-progress-meter -C - -o "k2_pluspf_16gb_${BUILD}.tar.gz" "$URL"
echo "[$(date)] extracting"
tar -xzf "k2_pluspf_16gb_${BUILD}.tar.gz"
rm -f "k2_pluspf_16gb_${BUILD}.tar.gz"
echo "[$(date)] done"; ls -la "$DB"
