Transcriptomic deconvolution in cancer and other heterogeneous tissues remains challenging. Available methods lack the ability to estimate both component-specific proportions and expression profiles for individual samples. We present DeMixT, a new tool to deconvolve high dimensional data from mixtures of more than two components. DeMixT implements an iterated conditional mode algorithm and a novel gene-set-based component merging approach to improve accuracy. In a series of experimental validation studies and application to TCGA data, DeMixT showed high accuracy. Improved deconvolution is an important step towards linking tumor transcriptomic data with clinical outcomes. An R package, scripts and data are available.

# DeMixT
A deconvolution framework for mixed transcriptomes from heterogeneous tumor samples with two or three components using expression data from RNAseq or microarray platforms

DeMixT is a frequentist-based method and fast in yielding accurate estimates of cell proportions and compart-ment-specific expression profiles for two-component and three-component deconvolution problem. Our method promises to provide deeper insight into cancer biomarkers and assist in the development of novel prognostic markers and therapeutic strategies. 

The function DeMixT is designed to finish the whole pipeline of deconvolution for two or three components. DeMixT.S1 function is designed to estimate the proportions of all mixed samples for each mixing component. DeMixT.S2 function is designed to estimate the component-specific deconvolved expressions of individual mixed samples for a given set of genes.

# Installation
To install this package, start R and enter:

devtools:::install_github("wwylab/DeMixT/DeMixT_0.1")

DeMixT enables the feature of parallel computing using OpenMP. If OpenMP has not been installed, you can still install and use DeMixT but without parallel computing feature. Please download and install OpenMP if you want to enable the feature. To learn more about how to install OpenMP, please check the file "How to install OpenMP.docx".

# Mixed cell line data
This data set is used in our validation experiment for DeMixT. To generate this dataset in RNA-seq, we performed a mixing experiment, in which we mixed mRNAs from three cell lines: lung adenocarcinoma in humans (H1092), cancer-associated fibroblasts (CAFs) and tumor infiltrating lymphocytes (TIL), at different proportions to generate 32 samples, including 9 samples that correspond to three repeats of a pure cell line sample for three cell lines. The RNA amount of each tissue in the mixture samples was calculated on the basis of real RNA concentrations tested in the biologist’s lab

# How to install OpenMP
We provide a brief instruction for installing OpenMP, which is needed to enable the parallel computing for DeMixT.
