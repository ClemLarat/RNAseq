#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("ReactomePA")
#BiocManager::install("tictoc")
#BiocManager::install("wesanderson")

#
# 0. load libraries
#
library(crayon)
library(clusterProfiler)
library(enrichplot)
library(tictoc)
library(viridis)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)

#
# 1. read files and generate lists of genes
#
filename = 'KO_VS_17-1.8.tsv'
df = read.table(filename, sep='\t', header=TRUE)
df_up = df[df$log2FoldChange > 0, ] 
df_down = df[df$log2FoldChange < 0, ] 

ensemblIDs = row.names(df_up)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_one_up = convertedIDs$ENTREZID
length(list_one_up)

ensemblIDs = row.names(df_down)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_one_down = convertedIDs$ENTREZID
length(list_one_down)

filename = 'KO_VS_17-4.2.tsv'
df = read.table(filename, sep='\t', header=TRUE)
df_up = df[df$log2FoldChange > 0, ] 
df_down = df[df$log2FoldChange < 0, ] 

ensemblIDs = row.names(df_up)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_two_up = convertedIDs$ENTREZID
length(list_two_up)

ensemblIDs = row.names(df_down)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_two_down = convertedIDs$ENTREZID
length(list_two_down)

filename = 'KO_VS_Ctrl.tsv'
df = read.table(filename, sep='\t', header=TRUE)
df_up = df[df$log2FoldChange > 0, ] 
df_down = df[df$log2FoldChange < 0, ] 

ensemblIDs = row.names(df_up)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_three_up = convertedIDs$ENTREZID
length(list_three_up)

ensemblIDs = row.names(df_down)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_three_down = convertedIDs$ENTREZID
length(list_three_down)

filename = '17-4.2_VS_17-1.8.tsv'
df = read.table(filename, sep='\t', header=TRUE)
df_up = df[df$log2FoldChange > 0, ] 
df_down = df[df$log2FoldChange < 0, ] 

ensemblIDs = row.names(df_up)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_four_up = convertedIDs$ENTREZID
length(list_four_up)

ensemblIDs = row.names(df_down)
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_four_down = convertedIDs$ENTREZID
length(list_four_down)

geneLists = list('17-1.8 vs KO down-regulation'=list_one_up, 
                '17-1.8 vs KO up-regulation'=list_one_down,
                '17-4.2 vs KO down-regulation'=list_two_up, 
                '17-4.2 vs KO up-regulation'=list_two_down,
                'Ctrl vs KO down-regulation'=list_three_up,
                'Ctrl vs KO up-regulation'=list_three_down,
                '17-1.8 vs 17-4.2 down-regulation'=list_four_up,
                '17-1.8 vs 17-4.2 up-regulation'=list_four_down)

#
# 2. run the analysis on pathway enrichment
#
tic()
ck = compareCluster(geneLists, 
                    fun="enrichPathway", 
                    pvalueCutoff=0.05)
toc()
# this step takes a long time

#
# 3. plot the enrichment results
#
p1 = dotplot(ck, size='count', showCategory=10, font.size=3)

my_log_breaks = round(log10(0.05)):round(log10(min(ck@compareClusterResult$p.adjust)))
my_breaks = 10**my_log_breaks
p2 = p1 +  scale_fill_gradientn(colours = pal1, trans="log", breaks = my_breaks)
print(p2)

#
# 4. store the enrichment results in a table
#
storage_file = 'clusterProfiler_enrichments_RP.tsv'
write.table(ck@compareClusterResult, storage_file, quote=FALSE, sep='\t')
