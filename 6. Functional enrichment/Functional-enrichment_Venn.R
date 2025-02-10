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
filename = 'ATG7(2)-high_no-ATG7(1).tsv'
df = read.csv(filename, sep='\t', header=TRUE)

ensemblIDs = df$ENSEMBL
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_one = convertedIDs$ENTREZID
length(list_one)

filename = 'ATG7(2)-low_no-ATG7(1).tsv'
df = read.csv(filename, sep='\t', header=TRUE)

ensemblIDs = df$ENSEMBL
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_two = convertedIDs$ENTREZID
length(list_two)

filename = 'ATG7(2)-dose-independent_no-ATG7(1).tsv'
df = read.csv(filename, sep='\t', header=TRUE)

ensemblIDs = df$ENSEMBL
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_three = convertedIDs$ENTREZID
length(list_three)

filename = 'ATG7(1).tsv'
df = read.csv(filename, sep='\t', header=TRUE)

ensemblIDs = df$ENSEMBL
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_four = convertedIDs$ENTREZID
length(list_four)

filename = 'ATG7(2)-dose-independent_ATG7(1).tsv'
df = read.csv(filename, sep='\t', header=TRUE)

ensemblIDs = df$ENSEMBL
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_five = convertedIDs$ENTREZID
length(list_five)

filename = 'ATG7(2)-low_ATG7(1).tsv'
df = read.csv(filename, sep='\t', header=TRUE)

ensemblIDs = df$ENSEMBL
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_six = convertedIDs$ENTREZID
length(list_six)

filename = 'ATG7(2)-high_ATG7(1).tsv'
df = read.csv(filename, sep='\t', header=TRUE)

ensemblIDs = df$ENSEMBL
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_seven = convertedIDs$ENTREZID
length(list_seven)

filename = 'ATG7(2)-dose-independent_overall.tsv'
df = read.csv(filename, sep='\t', header=TRUE)

ensemblIDs = df$ENSEMBL
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_eight = convertedIDs$ENTREZID
length(list_eight)

filename = 'ATG7(2)-high_vs_KO-17-4.2.tsv'
df = read.csv(filename, sep='\t', header=TRUE)

ensemblIDs = df$ENSEMBL
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_nine = convertedIDs$ENTREZID
length(list_nine)

filename = 'ATG7(2)-low_vs_KO-17-1.8.tsv'
df = read.csv(filename, sep='\t', header=TRUE)

ensemblIDs = df$ENSEMBL
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_ten = convertedIDs$ENTREZID
length(list_ten)

filename = 'ATG7(2)-dose-gradient.tsv'
df = read.csv(filename, sep='\t', header=TRUE)

ensemblIDs = df$ENSEMBL
convertedIDs = bitr(ensemblIDs, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb='org.Hs.eg.db')
list_eleven = convertedIDs$ENTREZID
length(list_eleven)

geneLists = list('ATG7(2)high no ATG7(1)'=list_one, 
                'ATG7(2)low no ATG7(1)'=list_two,
                'ATG7(2) no ATG7(1)'=list_three,
                'ATG7(1)'=list_four,
                'ATG7(2) ATG7(1)'=list_five,
                'ATG7(2)low ATG7(1)'=list_six,
                'ATG7(2)high ATG7(1)'=list_seven,
                'ATG7(2) dose independent'=list_eight,
                'ATG7(2)high vs KO and ATG7(2)low'=list_nine,
                'ATG7(2)low vs KO and ATG7(2)high'=list_ten,
                'ATG7(2) dose gradient' =list_eleven)

#
# 2. run the analysis on pathway enrichment
#
tic()
ck = compareCluster(geneLists, 
                    fun="enrichPathway", 
                    pvalueCutoff=0.05)
toc()

# this long step takes a long time

#
# 3. plot the enrichment results
#
p1 = dotplot(ck, size='count', showCategory=10, font.size=4)

my_log_breaks = round(log10(0.05)):round(log10(min(ck@compareClusterResult$p.adjust)))
my_breaks = 10**my_log_breaks
p2 = p1 +  scale_fill_gradientn(colours = pal1, trans="log", breaks = my_breaks)
p3 <- p2 + theme(axis.text.x = element_text(angle=45, hjust=1))
print(p3)

#
# 4. store the enrichment results in a table
#
storage_file = 'clusterProfiler_enrichments_Supervenn_RP.tsv'
write.table(ck@compareClusterResult, storage_file, quote=FALSE, sep='\t')
