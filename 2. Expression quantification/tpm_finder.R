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
kallisto_dir = "~/RNAseq/Kallisto"
results_dir = '~/RNAseq/Kallisto'

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
path = file.path(dirnames, 'abundance.h5')
labels = sapply(strsplit(path, split='/',fixed=TRUE), function(x) (x[9]))
condition = c(rep('17-1.8', 3), rep('17-4.2', 3), rep('KO', 3), rep('Ctrl', 3))

metadata = data.frame(labels)
metadata$condition = condition
metadata$path = path
View(metadata)

#
# 2.1. read files
#
txi = tximport(metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

#
# 2.2. find abundance
#
tpm = txi$abundance
colnames(tpm) = metadata$condition
dim(tpm)
View(tpm)

#
# 2.3. store
#
store = paste(results_dir, '/DESeq2_TPM_values.tsv', sep='')
write.table(tpm, file=store, quote=FALSE, sep='\t', col.names=NA)

store = paste(results_dir, '/annotation.tsv', sep='')
write.table(t2g, file=store, quote=FALSE, sep='\t', col.names=NA)
#
# 3. contrasts
#
threshold = 10
effect_size_threshold = 1

#
# 3.1. contrast effect of ATG7 expression in Control VS. 17-1.8
#
rule = (metadata$condition == 'Ctrl') | (metadata$condition == '17-high')
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~condition) 
dds$time = relevel(dds$condition, ref="Ctrl")

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('Differential expression between Ctrl and 17-1.8:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, 
            file=paste(results_dir, '/Ctrl_VS_17-1.8.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, 
                       '/Ctrl_VS_17-1.8.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('condition')) + ggtitle('Ctrl VS 17-1.8')
ggsave(file.path(results_dir, 'Ctrl VS 17-1.8.png'))

# 3.2. contrast effect of ATG7 expression in Control VS. 17-4.2
#
rule = (metadata$condition == 'Ctrl') | (metadata$condition == '17-low')
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~condition) 
dds$time = relevel(dds$condition, ref="Ctrl")

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('Differential expression between Ctrl and 17-4.2:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, 
            file=paste(results_dir, '/Ctrl_VS_17-4.2.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, 
                       '/Ctrl_VS_17-4.2.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('condition')) + ggtitle('Ctrl VS 17-4.2')
ggsave(file.path(results_dir, 'Ctrl VS 17-4.2.png'))

# 3.3. contrast effect of ATG7 expression in Control VS. KO
#
rule = (metadata$condition == 'Ctrl') | (metadata$condition == 'KO')
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~condition) 
dds$time = relevel(dds$condition, ref="Ctrl")

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('Differential expression between Ctrl and KO:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, 
            file=paste(results_dir, '/Ctrl_VS_KO.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, 
                       '/Ctrl_VS_KO.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('condition')) + ggtitle('Ctrl VS KO')
ggsave(file.path(results_dir, 'Ctrl VS KO.png'))

# 3.4. contrast effect of ATG7(2) expression in KO VS. 17-1.8
#
rule = (metadata$condition == 'KO') | (metadata$condition == '17-high')
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
cat(blue(paste('Differential expression between KO and 17-1.8:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, 
            file=paste(results_dir, '/KO_VS_17-1.8.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, 
                       '/KO_VS_17-1.8.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('condition')) + ggtitle('KO VS 17-1.8')
ggsave(file.path(results_dir, 'KO VS 17-1.8.png'))

# 3.5. contrast effect of ATG7(2) expression in KO VS. 17-4.2
#
rule = (metadata$condition == 'KO') | (metadata$condition == '17-low')
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
cat(blue(paste('Differential expression between KO and 17-4.2:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, 
            file=paste(results_dir, '/KO_VS_17-4.2.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, 
                       '/KO_VS_17-4.2.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('condition')) + ggtitle('KO VS 17-4.2')
ggsave(file.path(results_dir, 'KO VS 17-4.2.png'))

# 3.6. contrast effect of ATG7(2) expression in 17-4.2 VS. 17-1.8
#
rule = (metadata$condition == '17-low') | (metadata$condition == '17-high')
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~condition) 
dds$time = relevel(dds$condition, ref="17-low")

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('Differential expression between 17-4.2 and 17-1.8:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, 
            file=paste(results_dir, '/17-4.2_VS_17-1.8.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, 
                       '/17-4.2_VS_17-1.8.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('condition')) + ggtitle('17-4.2 VS 17-1.8')
ggsave(file.path(results_dir, '17-4.2 VS 17-1.8.png'))