rm(list = ls())

#
# -1. install libraries
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")

#
# 0. load libraries
#
library(DESeq2)
library(tximport)
library(biomaRt)
library(BiocParallel)
library(crayon) 
library(ggplot2)

#
# 0. user-defined variables
#
setwd("~/scratch/")
kallisto_dir = "~/RNAseq/kallisto.100"
results_dir = '~/RNAseq/deseq2'

#
# 1. generate gene to transcript mapping
#
mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                        dataset="hsapiens_gene_ensembl",
                        host = 'https://www.ensembl.org',
                        verbose = TRUE)
working_attributes = c('ensembl_transcript_id', 
                       'ensembl_gene_id', 
                       'external_gene_name',
                       'gene_biotype',
                       'description')
t2g = biomaRt::getBM(attributes=working_attributes, 
                     mart=mart,
                     verbose=TRUE)
dim(t2g)
View(t2g)

#
# 2. define metadata
#
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
paths = file.path(dirnames, 'abundance.h5')
labels = sapply(strsplit(paths, split='/',fixed=TRUE), function(x) (x[9]))
condition = c(rep('17-1.8', 3), rep('17-4.2', 3), rep('KO', 3), rep('Ctrl', 3))

metadata = data.frame(labels)
metadata$condition = condition
metadata$path = paths
View(metadata)

#
# 3. contrasts
#
threshold = 10
effect_size_threshold = log2(2)

#
# 3.1. contrast Ctrl vs KO
#
rule = (metadata$condition == 'KO') | (metadata$condition == 'Ctrl')
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~condition) 
dds$time = relevel(dds$condition, ref="KO")

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('contrast Ctrl vs KO:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, 
            file=paste(results_dir, '/KO_VS_Ctrl.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, 
                       '/KO_VS_Ctrl.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('condition')) + ggtitle('effect Ctrl vs KO')
#ggsave(file.path(results_dir, 'KO_VS_Ctrl.png'))

#
# 3.2. contrasts 17-1.8 vs KO
#
rule = (metadata$condition == '17-1.8') | (metadata$condition == 'KO')
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~condition) 
dds$time = relevel(dds$condition, ref="KO")

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('contrast 17-1.8 vs KO:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, 
            file=paste(results_dir, '/KO_VS_17-1.8.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, 
                       '/KO_VS_17-1.8.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('condition')) + ggtitle('effect 17-1.8 vs KO')
#ggsave(file.path(results_dir, '17-1.8_VS_KO.png'))

#
# 3.3. contrasts 17-4.2 vs KO
#
rule = (metadata$condition == '17-4.2') | (metadata$condition == 'KO')
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~condition) 
dds$time = relevel(dds$condition, ref="KO")

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('contrast 17-4.2 vs KO:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, 
            file=paste(results_dir, '/KO_VS_17-4.2.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, 
                       '/KO_VS_17-4.2.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('condition')) + ggtitle('effect 17-4.2 vs KO')
#ggsave(file.path(results_dir, 'KO_VS_17-4.2.png'))

#
# 3.3. contrasts 17-1.8 vs 17-4.2
#
rule = (metadata$condition == '17-1.8') | (metadata$condition == '17-4.2')
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~condition) 
dds$time = relevel(dds$condition, ref="17-4.2")

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('contrast 17-1.8 vs 17-4.2:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, 
            file=paste(results_dir, '/17-4.2_VS_17-1.8.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, 
                       '/17-4.2_VS_17-1.8.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('condition')) + ggtitle('effect 17-1.8 vs 17-4.2')
#ggsave(file.path(results_dir, '17-4.2_VS_17-1.8.png'))

