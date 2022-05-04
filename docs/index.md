---
layout: full
homepage: true
title: "DeMixT"
disable_anchors: true
description: Cell type-specific deconvolution of heterogeneous tumor samples with two or three components using expression data from RNAseq or microarray platforms
---

## DeMixT Overview

<!-- <div style="display:table; clear:both"> 
<div style="float: left;width: 50%;"> -->

Transcriptomic deconvolution in cancer and other heterogeneous tissues remains challenging. Available methods lack the ability to estimate both component-specific proportions and expression profiles for individual samples. We present DeMixT, a new tool to deconvolve high dimensional data from mixtures of two or three cellular components (i.e. within heterogenous tissues such as cancers). DeMixT implements an iterated conditional mode algorithm and a gene-set-based component merging approach to improve accuracy. In a series of experimental validation studies and application across large datasets of cancer studies, DeMixT showed high accuracy in inference of cell-type specific proportions[1-2]. Improved deconvolution is an important step towards linking tumor transcriptomic data with phenotypes and clinical outcomes.

An example of how to use DeMixT is available [here](tutorial.html). 

<center>
<img src="./etc/demixt.jpg" alt="demixt" width="60%" />
</center>
<!-- </div> -->
<!-- <div style="float: right; width: 50%"> <img src="./etc/demixt.jpg" alt="demixt" /> </div> -->
<!-- </div> -->

## Reference
[1] Ahn, J. et al. DeMix: Deconvolution for mixed cancer transcriptomes using raw measured data. Bioinformatics 29, 1865–1871 (2013).

[2] Wang, Z. et al. Transcriptome Deconvolution of Heterogeneous Tumor Samples with Immune Infiltration. iScience 9, 451–460 (2018).