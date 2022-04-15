# DeMixT
A deconvolution framework for mixed transcriptomes from heterogeneous tumor samples with two or three components using expression data from RNAseq or microarray platforms

DeMixT is a frequentist-based method and fast in yielding accurate estimates of cell proportions and compart-ment-specific expression profiles for two-component and three-component deconvolution problem. Our method promises to provide deeper insight into cancer biomarkers and assist in the development of novel prognostic markers and therapeutic strategies. 

The function DeMixT is designed to finish the whole pipeline of deconvolution for two or three components. ``DeMixT_S1`` function is designed to estimate the proportions of all mixed samples for each mixing component. ``DeMixT_GS`` function is designed to estimate the proportions of all mixed samples for each mixing component with profile likelihood based gene selection. ``DeMixT_S2`` function is designed to estimate the component-specific deconvolved expressions of individual mixed samples for a given set of genes.

# Installation
DeMixT source files are compatible with windows, linux and mac os.

This version is for users who have OpenMP or MPI on the computer. To install this package, start R and enter:
```
# install devtools if necessary
install.packages('devtools')

devtools::install_github("wwylab/DeMixT")

# load package
library(DeMixT)
```

# Use ``DeMixT``
A tutorial is available at [https://wwylab.github.io/DeMixT/](https://wwylab.github.io/DeMixT/).

# Cite ``DeMixT``
Wang, Z. et al. Transcriptome Deconvolution of Heterogeneous Tumor Samples with Immune Infiltration. iScience 9, 451–460 (2018).