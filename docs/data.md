---
layout: page
title: Data
disable_anchors: true
description: ~
---

## Dataset used in this tutorial
### TCGA PRAD
Raw expression counts of 458 samples (406 tumor and 52 normal) were downloaded from TCGA data portal (v14.0) ([https://portal.gdc.cancer.gov/](https://portal.gdc.cancer.gov/)). After filtering out the tumor samples whose expression profiles are highly similar to those of normal samples, and the normal samples whose expression profiles are highly similar to those of tumor samples, 342 samples (295 tumor and 47 normal) are remaining and used in the tutorial. The row counts of these samples contained in the file ([PRAD.RData](./etc/PRAD.RData)). 

## Dataset used in DeMixT paper[1]
### Reference profile data from the GTEx study
The Genotype-Tissue Expression (GTEx) project provides a comprehensive public resource to study tissue-specific gene expression and regulation [2]. RNA sequencing data from 42 normal prostate samples, 67 normal thyroid samples, and 20 normal lung samples without significant pathology in the corresponding tissue types were downloaded. For more information on how DeMixT uses unmatched reference profile data from the GTEx study, please refer to the vignette ([https://www.bioconductor.org/packages/release/bioc/vignettes/DeMixT/inst/doc/demixt.html](https://www.bioconductor.org/packages/release/bioc/vignettes/DeMixT/inst/doc/demixt.html)) and our recent paper [3].

Raw count matrices of selected samples for prostate, lung, and thyroid normal tissues can be downloaded from: [https://github.com/wwylab/DeMixTallmaterials/tree/master/Data/Reference_profile_data_from_the_GTEx_study](https://github.com/wwylab/DeMixTallmaterials/tree/master/Data/Reference_profile_data_from_the_GTEx_study).

### Mixed cell line data
This dataset is used in our validation experiment for DeMixT [1]. We performed a mixing experiment, in which we mixed mRNAs from three cell lines: lung adenocarcinoma in humans (H1092), cancer-associated fibroblasts (CAFs), and tumor infiltrating lymphocytes (TIL), at different proportions to generate 32 samples, including 9 samples that correspond to three repeats of a pure cell line sample for the three cell lines. The RNA amount of each tissue in the mixture samples was calculated on the basis of real RNA concentrations.

Knitr documentation for the DeMixT paper [1] can be downloaded from the website: [http://bioinformatics.mdanderson.org/Software/DeMixT/online_methods.html](http://bioinformatics.mdanderson.org/Software/DeMixT/online_methods.html).

## Acknowledgements
GTEx Project - The Genotype-Tissue Expression (GTEx) Project was supported by the Common Fund of the Office of the Director of the National Institutes of Health, and by NCI, NHGRI, NHLBI, NIDA, NIMH, and NINDS. The data used for the analyses described in this manuscript were obtained from the GTEx Portal ([https://gtexportal.org/home/](https://gtexportal.org/home/)) with dbGaP accession number phs000424.vN.pN.

## Reference

[1] Wang, Z. et al. Transcriptome Deconvolution of Heterogeneous Tumor Samples with Immune Infiltration. iScience 9, 451–460 (2018).

[2] GTEx Consortium. The Genotype-Tissue Expression (GTEx) pilot analysis: Multitissue gene regulation in humans. Science 348, 648–660 (2015).

[3] Cao, S. et al. Estimation of tumor cell total mRNA expression in 6,590 cancers predicts disease progression. Nature Biotechnology (2022) (in press).



