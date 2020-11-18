# 1. load all necessary packages
library(tidyverse) # package for data wrangling and data visualization
library(DESeq2) # package for analysis of differential expresion
library("org.Hs.eg.db") # package to annotate genes
library(ChIPseeker) # package to annotate ChIP-seq peaks
library(TxDb.Hsapiens.UCSC.hg38.knownGene) # another annotation package to assign genes to genomic regions
library(clusterProfiler) # package to assign funcitonal enrichmets to sets of genes
library(chipenrich)
library(ReactomePA)

###for individual peak Files
library(Gviz)
#library(GenomicInteractions)
#library(rtracklayer)
#library(DOSE)

# 2. define path to and read data-files for RNA-seq and ChIP-seq data


directory = "./data/rna_seq_counts/24h"
sampleFiles = grep("counts",list.files(directory),value=TRUE)
sampleCondition = substr(sampleFiles, 8,14)
sampleTable = data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          experiment = sampleFiles,
                          condition = sampleCondition)

mef2d = readPeakFile("./data/Galaxy_MACS2_MEF2D_OA_repl_q005_dupAll.bed")



# 3. analysis of differential expression

ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
ddsHTSeq

keep = rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq = ddsHTSeq[keep,]


# 4. calculation of significantly changing genes, candidates are sorted by p-value and a p-value cut-off < 0.1 is chosen

ddsHTSeq = DESeq(ddsHTSeq)
res = results(ddsHTSeq)
res = results(ddsHTSeq, name = "condition_siMEF24_vs_siCTR24")
resOrdered = res[order(res$pvalue),]
summary(res)
sum(res$padj < 0.1, na.rm = TRUE)


# 5. post-processing of results from differential expression analysis


# variance-stabilzed transformation
vsd<-vst(ddsHTSeq)
head(assay(vsd), 3)
VSD<-as.data.frame(assay(vsd))
VSD<-rownames_to_column(VSD, var = "ENSEMBL")

plotMA(res, ylim = c(-2,2), main = "significanlty changing genes")
plotCounts(ddsHTSeq, gene=which.min(res$padj), intgroup="condition")
plotPCA(vsd)


# 6. Transformation of data to perform subsequent analysis, like integration with ChIP-seq data, functional annotations and more plotting

MEF2D_si24 = as.data.frame(res)
MEF2D_si24 = rownames_to_column(MEF2D_si24, var = "ENSEMBL")
counts = as.data.frame(counts(ddsHTSeq, normalized = T))
counts = rownames_to_column(counts, var = "ENSEMBL")
Norm_MEF2D_si24 = left_join(MEF2D_si24, counts)
symbols = mapIds(org.Hs.eg.db, keys = Norm_MEF2D_si24$ENSEMBL, keytype = "ENSEMBL", column = "SYMBOL")
symbols = as.data.frame(symbols)
MEF2D_si24_Norm_Symbol = bind_cols(symbols, Norm_MEF2D_si24)


# 7. Use ChIP-seq data for the identification of direct target genes:

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

mef2d_Anno<-annotatePeak(mef2d, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoBar(mef2d_Anno)

tagMatrix <- getTagMatrix(mef2d, windows=promoter)
tagHeatmap(tagMatrix , xlim=c(-3000,3000), color ="red")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")



# 8. identification of genes with binding-sites in their vicinity

mef2d_gene <- seq2gene(mef2d, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(mef2d_gene)
dotplot(pathway2, title="OCI-AML-3 MEF2D")


MEF2D_anno<-as.data.frame(mef2d_Anno)
MEF2D_peaks<-MEF2D_anno[,c(1:7,13,21:23)]
colnames(MEF2D_peaks)<-c("chr","start", "end", "peak_width", "strand", "peakID", "peak_score", "annotation", "distanceToTSS", "ENSEMBL", "SYMBOL")
MEF2D_targets<-MEF2D_peaks%>%left_join(Norm_MEF2D_si24)
#gene<-bitr(MEF2D_targets$ENSEMBL, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
p<-enrichGO(gene = MEF2D_targets$ENSEMBL,OrgDb = 'org.Hs.eg.db',keyType = 'ENSEMBL',ont='MF', pvalueCutoff = 0.05,pAdjustMethod = 'BH')
dotplot(p, title="OCI-AML-3 MEF2D")
