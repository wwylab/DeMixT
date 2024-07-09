# DeMixT (v 1.20.1)
DeMixT is a frequentist-based method and fast in yielding accurate estimates of cell proportions and compartment-specific expression profiles for two- and three-component deconvolution from heterogeneous tumor samples. 

# Updates
(06/23/2024) Script DeMixT_preprocessing.R now inside the R package for input data preprocessing, to be consistent with the bioconductor version 1.20.
The DeMixT tutorial at https://wwylab.github.io/DeMixT/ is updated.

(06/18/2022) Fixed a bug in ``R/DeMixT_GS.R:130``: changed ``return(class(try(solve(m), silent = T)) == "matrix")`` to ``return(class(try(solve(m), silent = T))[1] == "matrix")`` to make it work properly under both ``R 3.x`` and ``R 4.x``; since under ``R 4.x``, ``class(matrix object)`` returns ``"matrix", "array"``, instead only ``"matrix"``  under ``R 3.x``. Please be aware of this bug for those who installed ``DeMixT``(v1.10.0) from Bioconductor. 


# Installation
The DeMixT package is compatible with Windows, Linux and MacOS. Specifically, for Linux and MacOS, the user can install the latest ``DeMixT``  (v 1.20.1) from GitHub:

```
if (!require("devtools", quietly = TRUE))
    install.packages('devtools')

devtools::install_github("wwylab/DeMixT")
```

For Windows, we recommend the user to install DeMixT (v 1.20.0) from ``Bioconductor``:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DeMixT")
```

Please note, Linux and MacOS users can also install DeMixT from ``Bioconductor``. 

Check if ``DeMixT`` is installed successfully:
```
# load package
library(DeMixT)
```

**We highly recommend the user to use DeMixT on Linux or Windows machines.**
# Issues of DeMixT on MacOS

1. DeMixT relies on OpenMP for parallel computing. Starting from R 4.00, R no longer supports OpenMP on MacOS, meaning the user can only run DeMixT with one core on MacOS.
2. We noticed there may be installation/running errors for DeMixT on MacOS machines with M1 chips. We are tring to fix it.

# Use ``DeMixT``
A tutorial is available at [https://wwylab.github.io/DeMixT/](https://wwylab.github.io/DeMixT/).

# Cite ``DeMixT``

[1] Ahn, J. et al. DeMix: Deconvolution for mixed cancer transcriptomes using raw measured data. Bioinformatics 29, 1865–1871 (2013).

[2] Wang, Z. et al. Transcriptome Deconvolution of Heterogeneous Tumor Samples with Immune Infiltration. iScience 9, 451–460 (2018).

[3] Cao, S. et al. Estimation of tumor cell total mRNA expression in 15 cancer types predicts disease progression. Nature Biotechnology Published online June 13 2022. doi 10.1038/s41587-022-01342-x.
