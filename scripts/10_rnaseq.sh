#!/usr/bin/env bash
# RNA-seq / ChrRNA-seq paired-end, fita-específica (dUTP, RF)
# FASTQ -> QC -> HISAT2 -> BAM -> featureCounts + bigWig por fita
# Uso: bash scripts/10_rnaseq.sh config/samples_rnaseq.tsv
set -euo pipefail
SHEET=${1:?"samples.tsv"}
THREADS=${THREADS:-8}
IDX=ref/hisat2/GRCh38
GTF=ref/Homo_sapiens.GRCh38.108.gtf
mkdir -p qc/fastqc qc/hisat2 bam counts tracks

for name in $(grep -v '^#' "$SHEET" | cut -f1); do
  # cada SRR (réplica/lane) é alinhado separadamente
  for r1 in raw/${name}__*_1.fastq.gz; do
    r2=${r1/_1.fastq.gz/_2.fastq.gz}
    id=$(basename "$r1" _1.fastq.gz)                    # <nome>__<SRR>
    [[ -s bam/$id.bam ]] && continue
    fastqc -t "$THREADS" -o qc/fastqc "$r1" "$r2"        # [SUGERIDO]
    # [DESCRITO] HISAT2 2.2.1, padrão + --rna-strandness RF --fr
    hisat2 -p "$THREADS" -x "$IDX" --rna-strandness RF --fr \
           -1 "$r1" -2 "$r2" --summary-file qc/hisat2/$id.txt \
      | samtools sort -@ "$THREADS" -o bam/$id.bam -
    samtools index bam/$id.bam
    samtools flagstat bam/$id.bam > qc/hisat2/$id.flagstat  # [SUGERIDO]
  done
done

# --- Contagem por gene ----------------------------------------------------
# [DESCRITO] featureCounts + Ensembl 108 ; [INFERIDO] -s 2 (RF = reversa) e pares
featureCounts -T "$THREADS" -p --countReadPairs -s 2 -t exon -g gene_id \
  -a "$GTF" -o counts/featureCounts_genes.txt bam/*.bam

# Checagem de fita [SUGERIDO]: -s 2 deve atribuir muito mais leituras que -s 1
featureCounts -T "$THREADS" -p --countReadPairs -s 1 -t exon -g gene_id \
  -a "$GTF" -o counts/strandcheck_s1.txt $(ls bam/*.bam | head -1) >/dev/null 2>&1 || true

# --- Cobertura por fita, 1 bp, CPM ---------------------------------------
# [DESCRITO] bamCoverage 3.5.1, binSize 1, CPM, fita-específico
for bam in bam/*.bam; do
  id=$(basename "$bam" .bam)
  for s in forward reverse; do
    [[ -s tracks/$id.$s.bw ]] && continue
    bamCoverage -p "$THREADS" -b "$bam" -o tracks/$id.$s.bw \
      --binSize 1 --normalizeUsing CPM --filterRNAstrand $s
  done
done

multiqc -f -o qc qc counts                                # [SUGERIDO]
echo "RNA-seq concluído: counts/featureCounts_genes.txt e tracks/*.bw"
