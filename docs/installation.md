---
layout: page
title: Installation
description: ~
---

`DeMixT` is implemented as an R package, which is compatible with Windows, Linux and MacOS. It can be installed from GitHub by:

#### 1. Install OpenMP
`DeMixT` requires OpenMP to enable the parallel computing. We provide a brief instruction for installing OpenMP. Please check the file [https://github.com/wwylab/DeMixT/HowtoinstallOpenMP.docx](https://github.com/wwylab/DeMixT/HowtoinstallOpenMP.docx).
#### 2. Install `devtools` if necessary
```r
install.packages('devtools')
```
#### 3. Install DeMixT

Two ways to install DeMixT. 
- From Github:
  
```
devtools::install_github("wwylab/DeMixT")
```

- From Bioconductor:
  
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DeMixT")
```

#### 4. Load package
```r
library(DeMixT)
```
