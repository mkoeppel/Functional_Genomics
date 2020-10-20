library(tidyverse) # package for data wrangling and data visualization
library(DESeq2) # package for analysis of differential expresion
library("org.Hs.eg.db") # package to annotate genes 
library(ChIPseeker) # package to annotate ChIP-seq peaks
library(TxDb.Hsapiens.UCSC.hg38.knownGene) # another annotation package to assign genes to genomic regions
library(clusterProfiler) # package to assign funcitonal enrichmets to sets of genes

directory = "./data/rna_seq_counts/48h"
sampleFiles = grep("counts",list.files(directory),value=TRUE) 
sampleCondition = substr(sampleFiles, 8,14) 
sampleTable = data.frame(sampleName = sampleFiles,
                         fileName = sampleFiles,
                         experiment = sampleFiles,
                         condition = sampleCondition)

mef2d = readPeakFile("./data/Galaxy_MACS2_MEF2D_OA_repl_q005_dupAll.bed")

ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                      directory = directory,
                                      design= ~ condition)
ddsHTSeq

keep = rowSums(counts(ddsHTSeq)) >= 10 
ddsHTSeq = ddsHTSeq[keep,]