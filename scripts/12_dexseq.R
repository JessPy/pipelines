#!/usr/bin/env Rscript
# DEXSeq: uso diferencial de éxons entre condições (precisa de >=2 réplicas por condição)
# Uso: Rscript scripts/12_dexseq.R config/dexseq_design.tsv results/dexseq
#   design.tsv: sample<TAB>condition<TAB>rep   (sample = nome do arquivo em counts/dexseq sem .txt)
suppressPackageStartupMessages(library(DEXSeq))
args <- commandArgs(TRUE)
design <- read.delim(args[1]); outdir <- ifelse(length(args) > 1, args[2], "results/dexseq")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
design$condition <- factor(design$condition); design$rep <- factor(design$rep)
files <- file.path("counts/dexseq", paste0(design$sample, ".txt"))
dxd <- DEXSeqDataSetFromHTSeq(files, sampleData = design,
         design = ~ sample + exon + condition:exon, flattenedfile = "ref/dexseq.gff")
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, BPPARAM = BiocParallel::MulticoreParam(4))
dxd <- testForDEU(dxd)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition")
res <- DEXSeqResults(dxd)
tab <- as.data.frame(res)[, c("groupID", "featureID", "exonBaseMean", "dispersion", "stat", "pvalue", "padj")]
write.csv(tab[order(tab$padj), ], file.path(outdir, "DEU_results.csv"), row.names = FALSE)
cat("Éxons com padj < 0.05:", sum(res$padj < 0.05, na.rm = TRUE), "\n")
DEXSeqHTML(res, path = file.path(outdir, "html"), FDR = 0.05)
