# DeMixT (v 1.10.1)
DeMixT is a frequentist-based method and fast in yielding accurate estimates of cell proportions and compartment-specific expression profiles for two- and three-component deconvolution from heterogeneous tumor samples. 

# Installation
The DeMixT source files are compatible with Windows, Linux and MacOS.

This version is for users who have OpenMP on the computer. 
We recommend the user to install ``DeMixT`` (v 1.10.0) from Bioconductor: 
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DeMixT")
```

If needed, the user can install the latest ``DeMixT``  (v 1.10.1) from GitHub:

```
if (!require("devtools", quietly = TRUE))
    install.packages('devtools')

devtools::install_github("wwylab/DeMixT")
```

Check if ``DeMixT`` is installed successfully:
```
# load package
library(DeMixT)
```

# Use ``DeMixT``
A tutorial is available at [https://wwylab.github.io/DeMixT/](https://wwylab.github.io/DeMixT/).

# Cite ``DeMixT``

[1] Ahn, J. et al. DeMix: Deconvolution for mixed cancer transcriptomes using raw measured data. Bioinformatics 29, 1865–1871 (2013).

[2] Wang, Z. et al. Transcriptome Deconvolution of Heterogeneous Tumor Samples with Immune Infiltration. iScience 9, 451–460 (2018).
