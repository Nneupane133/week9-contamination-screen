#!/bin/bash
# Databases for the RNA-virus arm (step 6).
#   geNomad DB (~1.4 GB) : marker-gene profiles that let geNomad score a contig
#                          as virus / plasmid / chromosome. This is what lets us
#                          separate a genuine replicating virus from an
#                          endogenous viral element carried in the host genome.
#   CheckV  DB (~5 GB)   : reference viral genomes + HMMs used to estimate how
#                          COMPLETE a recovered viral contig is. Completeness is
#                          the difference between "a viral fragment is present"
#                          (trace) and "a near-full viral genome is present"
#                          (biological) -- i.e. it drives the bin directly.
#   RefSeq viral protein : DIAMOND blastx target, an independent protein-level
#                          opinion that does not depend on geNomad's markers.
set -euo pipefail
WK=/home/aubnxn001/week9
export PATH="$WK/envs/w9virus/bin:$PATH"
mkdir -p "$WK/references/db"

echo "[$(date)] genomad database"
genomad download-database "$WK/references/db/"

echo "[$(date)] checkv database"
checkv download_database "$WK/references/db/checkv"

echo "[$(date)] RefSeq viral proteins -> DIAMOND"
cd "$WK/references/db"
curl -fL --no-progress-meter -o viral.protein.faa.gz \
  https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz
curl -fL --no-progress-meter -o viral.2.protein.faa.gz \
  https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.protein.faa.gz || true
cat viral.protein.faa.gz viral.2.protein.faa.gz 2>/dev/null > viral.all.protein.faa.gz || \
  mv viral.protein.faa.gz viral.all.protein.faa.gz
set +e; source /apps/profiles/modules_asax.sh.dyn >/dev/null 2>&1
module load diamond >/dev/null 2>&1; set -e
diamond makedb --in viral.all.protein.faa.gz -d refseq_viral_protein --threads 8
echo "[$(date)] virus DBs done"
ls -la "$WK/references/db"
