Transcriptomic deconvolution in cancer and other heterogeneous tissues remains challenging. Available methods lack the ability to estimate both component-specific proportions and expression profiles for individual samples. We present DeMixT, a new tool to deconvolve high dimensional data from mixtures of more than two components. DeMixT implements an iterated conditional mode algorithm and a novel gene-set-based component merging approach to improve accuracy. In a series of experimental validation studies and application to TCGA data, DeMixT showed high accuracy. Improved deconvolution is an important step towards linking tumor transcriptomic data with clinical outcomes. An R package, scripts and data are available.

# DeMixT
A deconvolution framework for mixed transcriptomes from heterogeneous tumor samples with two or three components using expression data from RNAseq or microarray platforms

DeMixT is a frequentist-based method and fast in yielding accurate estimates of cell proportions and compart-ment-specific expression profiles for two-component and three-component deconvolution problem. Our method promises to provide deeper insight into cancer biomarkers and assist in the development of novel prognostic markers and therapeutic strategies. 

The function DeMixT is designed to finish the whole pipeline of deconvolution for two or three components. DeMixT_DE function is designed to estimate the proportions of all mixed samples for each mixing component.DeMixT_GS function is designed to estimate the proportions of all mixed samples for each mixing component with profile likelihood based gene selection. DeMixT_S2 function is designed to estimate the component-specific deconvolved expressions of individual mixed samples for a given set of genes.

# Installation
DeMixT source files are compatible with windows, linux and mac os.

This version is for users who have OpenMP on the computer. To install this package, start R and enter:

devtools::install_github("wwylab/DeMixT")

For more information, please visit:
http://bioinformatics.mdanderson.org/main/DeMixT

# How to install OpenMP
We provide a brief instruction for installing OpenMP, which is needed to enable the parallel computing for DeMixT. Please check the file "How_to_install_OpenMP.pdf".

# Reference profile data from the GTEx study 
The Genotype-Tissue Expression (GTEx) project provides a comprehensive public resource to study tissue-specific gene expression and regulation [1]. RNA sequencing data from 42 normal prostate samples, 67 normal thyroid samples and 20 normal
lung samples without significant pathology in the corresponding tissue types are downloaded. For more information on how DeMixT uses unmatched reference profile data from the GTEx study, please refer to the vignette and our recent bioRxiv preprint [2].

Raw count matrices of selected samples for prostate, lung and thyroid normal tissues can be downloaded from: 
<a href=“https://github.com/wwylab/DeMixTallmaterials/tree/master/Data/Reference_profile%20data_from_the%20GTEx_study%20”>Download here</a>



# Mixed cell line data
This data set is used in our validation experiment for DeMixT. To generate this dataset in RNA-seq, we performed a mixing experiment, in which we mixed mRNAs from three cell lines: lung adenocarcinoma in humans (H1092), cancer-associated fibroblasts (CAFs) and tumor infiltrating lymphocytes (TIL), at different proportions to generate 32 samples, including 9 samples that correspond to three repeats of a pure cell line sample for three cell lines. The RNA amount of each tissue in the mixture samples was calculated on the basis of real RNA concentrations tested in the biologist’s lab.

Knitr documentation for the DeMixT paper (Wang et al.) can be downloaded from the website:
http://bioinformatics.mdanderson.org/Software/DeMixT/online_methods.html.

# Reference
[1] GTEx Consortium. "The Genotype-Tissue Expression (GTEx) pilot analysis: Multitissue gene regulation in humans." Science 348.6235 (2015): 648-660.

[2] Cao, Shaolong, et al. "Differing total mRNA expression shapes the molecular and clinical phenotype of cancer." bioRxiv (2020).
