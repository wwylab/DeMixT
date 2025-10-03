# DeMixT v2.0

This repository contains code for DeMixT 2.0, which integrates and extends both *DeMixT* and *DeMixNB*.

## Overview

- ***DeMixT*** is a frequentist-based method that provides fast and accurate estimates of cell proportions and compartment-specific expression profiles for two- and three-component deconvolution from heterogeneous tumor samples.

- ***DeMixNB***, now incorporated into DeMixT 2.0, is a variant of DeMixT based on the negative binomial distribution. It enables accurate estimation of transcript proportions in two-component deconvolution, particularly effective for sparse RNA-seq data.


## What’s New in DeMixT 2.0

With this upgrade, DeMixT 2.0 unifies the strengths of DeMixT and DeMixNB, making it a more comprehensive and flexible tool for transcriptomic deconvolution across diverse data types and experimental scenarios.


## Installation
The DeMixT package is compatible with Windows, Linux and MacOS. Specifically, for Linux and MacOS, the user can install the latest DeMixT(v 2.0.0) from GitHub:

```
if (!require("devtools", quietly = TRUE))
    install.packages('devtools')

devtools::install_github("wwylab/DeMixT")
```

For Windows, we recommend the user to install from ``Bioconductor`` (to be released, current version is 1.20.0):
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

#### Issues of running on MacOS

- Both DeMixT and DeMixNB rely on OpenMP for parallel computing. Starting from R 4.00, R no longer supports OpenMP on MacOS, meaning the user can only run DeMixT/DeMixNB with one core on MacOS.

<u>We highly recommend the user to use DeMixT and DeMixNB on Linux or Windows machines.</u>

## Usage
- A DeMixT tutorial is available at [DeMixT](https://wwylab.github.io/DeMixT/).
- A DeMixNB tutorial is available at [DeMixNB](https://odin.mdacc.tmc.edu/~wwang7/DeMixNB_tutorial.html).


## Updates
(10/01/2025) Added implementation of the DeMixNB algorithm.

(06/23/2024) Script DeMixT_preprocessing.R now inside the R package for input data preprocessing, to be consistent with the bioconductor version 1.20.
The DeMixT tutorial at https://wwylab.github.io/DeMixT/ is updated.

(06/18/2022) Fixed a bug in ``R/DeMixT_GS.R:130``: changed ``return(class(try(solve(m), silent = T)) == "matrix")`` to ``return(class(try(solve(m), silent = T))[1] == "matrix")`` to make it work properly under both ``R 3.x`` and ``R 4.x``; since under ``R 4.x``, ``class(matrix object)`` returns ``"matrix", "array"``, instead only ``"matrix"``  under ``R 3.x``. Please be aware of this bug for those who installed ``DeMixT``(v1.10.0) from Bioconductor. 




## Citation

[1] Ahn, J. et al. DeMix: Deconvolution for mixed cancer transcriptomes using raw measured data. Bioinformatics 29, 1865–1871 (2013).

[2] Wang, Z. et al. Transcriptome Deconvolution of Heterogeneous Tumor Samples with Immune Infiltration. iScience 9, 451–460 (2018).

[3] Cao, S. et al. Estimation of tumor cell total mRNA expression in 15 cancer types predicts disease progression. Nature Biotechnology Published online June 13 2022. doi 10.1038/s41587-022-01342-x.
