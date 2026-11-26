#!/usr/bin/env bash
# QC de ChIP/ChOR-seq e perfis médios sobre genes replicados no pulso de EdU [SUGERIDO]
# Uso: bash scripts/22_metaplots_qc.sh EdUseq "tot_T0 tot_T1 tot_T2" rrpm
set -euo pipefail
edu=${1:?"amostra EdU-seq"}; samples=${2:?"lista de amostras"}; norm=${3:-cpm}
THREADS=${THREADS:-8}
GENES=ref/genes.hg38.chr.bed
mkdir -p results/qc results/metaplots

bams=$(for s in $samples; do echo bam/$s.hg38.bam; done)
bws=$(for s in $samples; do echo tracks/$s.$norm.bw; done)

# --- QC -------------------------------------------------------------------
plotFingerprint -p "$THREADS" -b $bams --labels $samples -o results/qc/fingerprint.pdf \
  --outRawCounts results/qc/fingerprint.tsv
multiBamSummary bins -p "$THREADS" -b $bams --labels $samples --binSize 1000 -o results/qc/bins.npz
plotCorrelation -in results/qc/bins.npz --corMethod pearson --whatToPlot heatmap \
  --skipZeros --removeOutliers -o results/qc/correlation.pdf

# --- genes de proteína inteiramente dentro das regiões replicadas [INFERIDO do desenho] ---
bedtools intersect -u -f 1.0 -a "$GENES" -b results/edu/$edu.replicated.bed3 \
  | awk '$3-$2 > 5000' > results/metaplots/genes_replicated.bed
echo "Genes replicados (>5 kb): $(wc -l < results/metaplots/genes_replicated.bed)"

# --- perfil ao longo do gene (TSS -> TES) ------------------------------------
computeMatrix scale-regions -p "$THREADS" -S $bws --samplesLabel $samples \
  -R results/metaplots/genes_replicated.bed -b 2000 -a 2000 --regionBodyLength 10000 \
  --skipZeros -o results/metaplots/genebody.$norm.mat.gz
plotProfile -m results/metaplots/genebody.$norm.mat.gz --perGroup \
  -o results/metaplots/genebody.$norm.profile.pdf
plotHeatmap -m results/metaplots/genebody.$norm.mat.gz \
  -o results/metaplots/genebody.$norm.heatmap.pdf

# --- centrado no TSS (pausa promotor-proximal) -----------------------------
computeMatrix reference-point --referencePoint TSS -p "$THREADS" -S $bws --samplesLabel $samples \
  -R results/metaplots/genes_replicated.bed -b 1000 -a 3000 --skipZeros \
  -o results/metaplots/tss.$norm.mat.gz
plotProfile -m results/metaplots/tss.$norm.mat.gz --perGroup -o results/metaplots/tss.$norm.profile.pdf
echo "Pronto: results/qc e results/metaplots"
