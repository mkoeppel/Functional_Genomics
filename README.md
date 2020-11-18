# Functional_Genomics
Finding real and direct target genes for a transcription factor is a challenging task in molecular biology. \
This R-script script is an example of an analytical pipeline from high-throughput sequencing (NGS) in functional genomics. 
It was created when I was teaching bioinformatics at the University of Applied Sciences in Bremen (HSB) in 2020.

The data used was generated in my lab at the DSMZ-Leibniz-Institut and, so far, is unpublished. Therefore, the raw-data has not yet been deposited to any nucleic acid archive. Provided is a form of pre-processed data, which is small enough to be distributed with the repo. All preprocessing steps are documented below and are standard first level analysis of NGS-data, that is quality-control and mapping/ counting to a reference genome. 

This repo contains the necessary steps to process NGS data derived from a typical functional genomics experiment:

A typical research question for such an experiment would be:
Which genes are directly regulated by a particular transcription factor?
  - This requires knowledge about those genes that change their expression, when the abundance of this particular transcription factor is changed (typically via siRNA-knockdown or overexpression; here it is a knockdown).
  - As genes can also change their expression due to indirect/ secondary effects, binding sites of the transcription factor within the genomic sequences also need to be identified and assigned to target genes.
In combination these two experiements allow the prediction of direct, real target genes with high probability.

As a factor, the Myocyte Enhancer Factor 2D (MEF2D) was examined. It is a transcription factor that is originally involved in the differentiation of muscle cells. However different mutations, usually involving some kind of translocation, cause it to act as an oncogene involved in the onset of different forms of leukemia. 
The cellular system used to examine the behaviour of this factor was the acute myeloid leukemia cell line OCI-AML3, as it has high mRNA levels and also clearly detectable protein-levels of this factor. 
