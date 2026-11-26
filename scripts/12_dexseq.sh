#!/usr/bin/env bash
# Uso diferencial de éxons (DEXSeq) - preparação da anotação e contagem por éxon.
# Bruno et al. 2024 relatam uso diferencial de éxons/íntrons com DEXSeq
# (fonte secundária; confirme os parâmetros nos STAR Methods).
set -euo pipefail
THREADS=${THREADS:-8}
GTF=ref/Homo_sapiens.GRCh38.108.gtf
PYDIR=$(Rscript -e 'cat(system.file("python_scripts", package="DEXSeq"))')
mkdir -p counts/dexseq
[[ -s ref/dexseq.gff ]] || python "$PYDIR/dexseq_prepare_annotation.py" "$GTF" ref/dexseq.gff
ls bam/*.bam | xargs -P "$THREADS" -I{} sh -c '
  id=$(basename {} .bam)
  python '"$PYDIR"'/dexseq_count.py -p yes -r pos -s reverse -f bam \
     ref/dexseq.gff {} counts/dexseq/$id.txt'
echo "Contagens por éxon em counts/dexseq/"
