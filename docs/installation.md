---
layout: page
title: Installation
disable_anchors: true
description: ~
---

The `DeMixT` source files are compatible with Windows, Linux, and MacOS. This version is for users who have OpenMP on the computer. 

We recommend the user to install ``DeMixT`` (v 1.20.0) from Bioconductor: 
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DeMixT")
```

If needed, the user can install it from GitHub:

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
