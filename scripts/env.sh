#!/bin/bash
# Shared environment loader (ASC asax / Lmod).
# NOTE: /apps/profiles/modules_asax.sh.dyn returns a non-zero status internally,
# so errexit has to be suspended across the source or every pipeline script dies
# on its first line.
__w9_errexit=0; case $- in *e*) __w9_errexit=1;; esac
set +e
source /apps/profiles/modules_asax.sh.dyn >/dev/null 2>&1
module load bbmap kraken2 samtools/1.18 minimap2 seqkit megahit_voutcn spades diamond blast+ >/dev/null 2>&1
[ "$__w9_errexit" = 1 ] && set -e
export WK=/home/aubnxn001/week9
true
