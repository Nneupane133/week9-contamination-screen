#!/bin/bash
# Build conda envs for tools ASC does not provide.
# Envs live under $HOME because the system anaconda prefix is read-only.
#
# WHY a modern kraken2 here even though ASC ships one: the ASC module is
# 2.0.9-beta (2019). The prebuilt PlusPF index is built with 2.1.x and Bracken
# expects 2.1.x behaviour, so we pin a matching pair rather than risk a silent
# format mismatch.
# WHY bracken at all: kraken2 assigns each read to the LCA of its k-mer hits,
# so reads from a species with close relatives in the DB pile up at genus or
# family level and the species row under-counts. Bracken redistributes those
# reads to species using the DB's k-mer distribution, which is what an honest
# "read fraction" needs.
set -euo pipefail
export CONDA_PKGS_DIRS=/home/aubnxn001/.conda/pkgs
ENVROOT=/home/aubnxn001/week9/envs
mkdir -p "$ENVROOT"
CH="-c conda-forge -c bioconda --override-channels -y"
mamba create -p "$ENVROOT/w9tax"   $CH kraken2=2.1.3 bracken=3.0 seqtk=1.4
echo "=== w9tax done $(date) ==="
mamba create -p "$ENVROOT/w9virus" $CH genomad=1.8.0 checkv=1.0.3
echo "=== w9virus done $(date) ==="
