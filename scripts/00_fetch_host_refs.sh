#!/bin/bash
# Fetch the two HOST reference genomes.
#
# WHY we need the host genome at all, when the question is about non-target DNA:
#   1. Denominator. "Read fraction" is only meaningful next to how much of the
#      library is host; a 0.1 % hit means different things at 95 % vs 40 % host.
#   2. EVE discrimination. The single hardest call in this assignment is
#      "free virus / live organism" vs "sequence integrated in the host genome".
#      The decisive test is whether the reads that hit a viral or bacterial
#      reference ALSO align, at equal or better identity, to the host assembly.
#      That test is impossible without the host genome.
#   3. Host depletion before assembly, so MEGAHIT/rnaSPAdes spend their memory
#      on the ~1 % that is interesting instead of the ~95 % that is chicken.
#
# Assemblies chosen:
#   GRCg7b (GCF_016699485.2) -- current RefSeq reference for Gallus gallus.
#   AaegL5.0 (GCF_002204515.2) -- current RefSeq reference for Aedes aegypti.
set -euo pipefail
REF=/home/aubnxn001/week9/references/host
mkdir -p "$REF"
B=https://ftp.ncbi.nlm.nih.gov/genomes/all
declare -A G=(
 [Gallus_gallus_GRCg7b]="$B/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz"
 [Aedes_aegypti_AaegL5]="$B/GCF/002/204/515/GCF_002204515.2_AaegL5.0/GCF_002204515.2_AaegL5.0_genomic.fna.gz"
)
for n in "${!G[@]}"; do
  echo "[$(date +%T)] $n"
  curl -fL --no-progress-meter -C - -o "$REF/$n.fna.gz" "${G[$n]}"
  ls -la "$REF/$n.fna.gz"
done
echo "[$(date +%T)] host refs done"
