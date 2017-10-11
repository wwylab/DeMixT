# DeMixT
Deconvolution models for mixed transcriptomes from heterogeneous tumor samples with two or three components using expression data from RNAseq or microarray platforms.

DeMixT is a frequentist-based method and fast in yielding accurate estimates of cell proportions and compart-ment-specific expression profiles for two-component and three-component deconvolution problem. Our method promises to provide deeper insight into cancer biomarkers and assist in the development of novel prognostic markers and therapeutic strategies. 

The function DeMixT is designed to finish the whole pipeline of deconvolution for two or three components. DeMixT.S1 function is designed to estimate the proportions of all mixed samples for each mixture component with or without subsetting gene set. DeMixT.S2 function is designed to estimate the deconvolved expressions of individual mixed tumor samples for unknown component given a subset of genes.

# Installation
To install this package, start R and enter:

devtools:::install_github("wwylab/DeMixT/DeMixT_0.1")

DeMixT enables the feature of parallel computing using OpenMP. If OpenMP has not been installed, you can still install and use DeMixT but without parallel computing feature. Please download and install OpenMP if you want to enable the feature. To learn more about how to install OpenMP, please check the file "How to install OpenMP.docx".
