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
library(GenomicInteractions)
library(rtracklayer)
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

# instead of the summarizing DESeq-function individual steps can also be taken to obtain intermediate results. Below, both possibilites are given:
# 4a: intermediate steps
ddsHTSeq_size_factor = estimateSizeFactors(ddsHTSeq)
sizeFactors(ddsHTSeq_size_factor)
ddsHTSeq_dispersion = estimateDispersions(ddsHTSeq_size_factor)
ddsHTSeq_GLM = nbinomWaldTest(ddsHTSeq_dispersion)
stepwise_results = results(ddsHTSeq_GLM, name="condition_siMEF24_vs_siCTR24")

# 4b: run in wrapper-function
ddsHTSeq = DESeq(ddsHTSeq)
res = results(ddsHTSeq, name="condition_siMEF24_vs_siCTR24")
resOrdered = res[order(res$pvalue),]
summary(res)
sum(res$padj<0.1, na.rm=TRUE)


# 5. post-processing of results from differential expression analysis

# variance-stabilzed transformation
vsd = vst(ddsHTSeq)
head(assay(vsd), 3)
VSD = as.data.frame(assay(vsd))
VSD = rownames_to_column(VSD, var="ENSEMBL")

plotMA(res,
       ylim=c(-2,2),
       main="significantly changing genes")
plotCounts(ddsHTSeq,
           gene = which.min(res$padj),
           intgroup="condition")
plotPCA(vsd)


# 6. Transformation of data to perform subsequent analysis, like integration with ChIP-seq data, functional annotations and more plotting

MEF2D_si24 = as.data.frame(res)
MEF2D_si24 = rownames_to_column(MEF2D_si24, var="ENSEMBL")
counts = as.data.frame(counts(ddsHTSeq, normalized=T))
counts = rownames_to_column(counts, var="ENSEMBL")
Norm_MEF2D_si24 = left_join(MEF2D_si24, counts)
symbols = mapIds(org.Hs.eg.db,
                 keys=Norm_MEF2D_si24$ENSEMBL,
                 keytype="ENSEMBL",
                 column="SYMBOL")
symbols = as.data.frame(symbols)
MEF2D_si24_Norm_Symbol = bind_cols(symbols, Norm_MEF2D_si24)


# 7. Use ChIP-seq data for the identification of direct target genes:

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
promoter = getPromoters(TxDb=txdb,
                        upstream=3000, downstream=3000)

mef2d_Anno = annotatePeak(mef2d,
                          tssRegion=c(-3000, 3000),
                          TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoBar(mef2d_Anno)

tagMatrix = getTagMatrix(mef2d, windows=promoter)
tagHeatmap(tagMatrix ,
           xlim=c(-3000,3000), color="red")
plotAvgProf(tagMatrix,
            xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab="Read Count Frequency") +
            ggtitle("Factor binding signal at start site of genes")


# 8. identification of genes with binding-sites in their vicinity

mef2d_gene = seq2gene(mef2d,
                      tssRegion=c(-1000, 1000),
                      flankDistance=3000, TxDb=txdb)
pathway2 = enrichPathway(mef2d_gene)
dotplot(pathway2, title="pathway enrichment of MEF2D target genes")


MEF2D_anno = as.data.frame(mef2d_Anno)
MEF2D_peaks = MEF2D_anno[,c(1:7,13,21:23)]
colnames(MEF2D_peaks) = c("chr","start", "end",
                          "peak_width", "strand",
                          "peakID", "peak_score",
                          "annotation", "distanceToTSS",
                          "ENSEMBL", "SYMBOL")
MEF2D_targets = MEF2D_peaks%>%left_join(Norm_MEF2D_si24)

#gene<-bitr(MEF2D_targets$ENSEMBL, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
p = enrichGO(gene=MEF2D_targets$ENSEMBL,
             OrgDb="org.Hs.eg.db",
             keyType="ENSEMBL", ont="MF",
             pvalueCutoff=0.05,
             pAdjustMethod="BH")

dotplot(p, title="OCI-AML-3 MEF2D")


# 9. visualization of differentially regulated genes

subset_targets = subset(MEF2D_targets, MEF2D_targets$log2FoldChange>0.5 &
                          MEF2D_targets$padj<0.1 |
                          MEF2D_targets$log2FoldChange< -0.5 &
                          MEF2D_targets$padj<0.1)

subset_targets = subset_targets[!duplicated(subset_targets["ENSEMBL"]),]



library(gplots)

heatAll = subset_targets[,c(18:23)]
colnames(heatAll) = c("siCTR_1", "siCTR_2",
                      "siCTR_3","siMEF2D_1",
                      "siMEF2D_2", "siMEF2D_3")
rownames(heatAll) = subset_targets[,11]

tab = as.matrix(heatAll, as.is=T, header=T, row.names=1)

RColorBrewer::brewer.pal(n=6,"Set1")
cols   = colorpanel(15,"#2c7bb6","grey96","#d7191c")


heatmap.2(tab,col=cols, margins=c(10,5),
          scale="row", dendrogram="col",
          Colv=T, Rowv=T,
          key=TRUE, main="expression of MEF2D target genes",
          density.info="none", trace="none",
          cexRow=0.5, cexCol=0.8)
dev.off()


# 10. rawdata for example target genes (as raw-data is not deposited so far, this currently does not work)

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

# add gene model, ideotrack etc
gtrack = GenomeAxisTrack()
grtrack = GeneRegionTrack(txdb, name="Gene Model",
                          transcriptAnnotation="symbol",
                          collapseTranscripts="longest")
ideoTrack = IdeogramTrack(genome="hg38")

# add data
MEF2D = DataTrack("/Volumes/macAir_backup/MKO_G_CUTnRUN/CUTnRUN_1/bam/G6_OA_4_MEF2D.chr.bam",
                  name="MEF2D",  col=("#d7191c"),
                  fill.mountain=c("#d7191c", "#d7191c"), ylim=c(0,20))

IgG = DataTrack("/Volumes/macAir_backup/MKO_G_CUTnRUN/CUTnRUN_1//bam/G6_OA_1_IgG.chr.bam",
                name="IgG",  col=("#2c7bb6"),
                fill.mountain=c("#2c7bb6", "#2c7bb6"), ylim=c(0,20))

OA_si_ctr = DataTrack("/Volumes/macAir_backup/MKO_G/MKO_G_AML_ngs/MKO_G_013_Lexogen/bam/G13_OA_siCTR24_1.lexo.s.bam",
                      name="ctr-siRNA",  col=("#d7191c"),
                      fill.mountain=c("#d7191c", "#d7191c"), ylim=c(0,40))

OA_si_MEF = DataTrack("/Volumes/macAir_backup/MKO_G/MKO_G_AML_ngs/MKO_G_013_Lexogen/bam/G13_OA_siMEF24_1.lexo.s.bam",
                      name="MEF2D-siRNA",  col=("#2c7bb6"),
                      fill.mountain=c("#2c7bb6", "#2c7bb6"), ylim=c(0,40))


mef2d = readPeakFile("./data/Galaxy_MACS2_MEF2D_OA_repl_q005_dupAll.bed")


# upstream of HDAC9
plotTracks(list(ideoTrack, gtrack, grtrack, MEF2D, IgG,  OA_si_ctr,  OA_si_MEF),
           chromosome ="chr7",  from =18490019, to= 18683544,
           type="polygon", main ="HDAC9, a directly regulated target gene")
