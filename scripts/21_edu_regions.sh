#!/usr/bin/env bash
# EdU-seq -> regiões replicadas no pulso (Bruno et al. 2024, GSE220088/90/91)
# Uso: bash scripts/21_edu_regions.sh EdUseq     (nome da amostra já processada pelo 20_chorseq.sh)
set -euo pipefail
name=${1:?"nome da amostra EdU-seq"}
THREADS=${THREADS:-8}
mkdir -p tracks results/edu
# [DESCRITO] bamCoverage 3.5.1, bin 1 bp, CPM
bamCoverage -p "$THREADS" -b bam/$name.hg38.bam --binSize 1 --normalizeUsing CPM \
  --outFileFormat bedgraph -o tracks/$name.edu.cpm.bedGraph
# [DESCRITO] MACS2 bdgbroadcall 2.2.7.1 --cutoff 0.85 --min-length 3000 --max-gap 20000
macs2 bdgbroadcall -i tracks/$name.edu.cpm.bedGraph \
  -c 0.85 -l 3000 -g 20000 -o results/edu/$name.replicated.bed
# bdgbroadcall escreve um BED12 com regiões fortes/fracas; guarde as 3 primeiras colunas
cut -f1-3 results/edu/$name.replicated.bed | LC_ALL=C sort -k1,1 -k2,2n | bedtools merge \
  > results/edu/$name.replicated.bed3
awk '{s+=$3-$2} END{printf "Regiões: %d  |  total replicado: %.1f Mb\n", NR, s/1e6}' results/edu/$name.replicated.bed3
