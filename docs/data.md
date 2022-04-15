---
layout: page
title: Data
description: ~
---

## TCGA PRAD
Raw expression counts of 342 samples (295 tumor and 47 normal) were downloaded from TCGA data portal ([https://portal.gdc.cancer.gov/](https://portal.gdc.cancer.gov/)). This data is used in the tutorial. It is available at ([PRAD.RData](./etc/PRAD.RData)). 


## Reference profile data from the GTEx study
The Genotype-Tissue Expression (GTEx) project provides a comprehensive public resource to study tissue-specific gene expression and regulation [1]. RNA sequencing data from 42 normal prostate samples, 67 normal thyroid samples and 20 normal lung samples without significant pathology in the corresponding tissue types are downloaded. For more information on how DeMixT uses unmatched reference profile data from the GTEx study, please refer to the vignette and our recent bioRxiv preprint [2].

Raw count matrices of selected samples for prostate, lung and thyroid normal tissues can be downloaded from: [https://github.com/wwylab/DeMixTallmaterials/tree/master/Data/Reference_profile_data_from_the_GTEx_study](https://github.com/wwylab/DeMixTallmaterials/tree/master/Data/Reference_profile_data_from_the_GTEx_study).

## Mixed cell line data
This data set is used in our validation experiment for DeMixT. To generate this dataset in RNA-seq, we performed a mixing experiment, in which we mixed mRNAs from three cell lines: lung adenocarcinoma in humans (H1092), cancer-associated fibroblasts (CAFs) and tumor infiltrating lymphocytes (TIL), at different proportions to generate 32 samples, including 9 samples that correspond to three repeats of a pure cell line sample for three cell lines. The RNA amount of each tissue in the mixture samples was calculated on the basis of real RNA concentrations tested in the biologist’s lab.

Knitr documentation for the DeMixT paper (Wang et al.) can be downloaded from the website: [http://bioinformatics.mdanderson.org/Software/DeMixT/online_methods.html](http://bioinformatics.mdanderson.org/Software/DeMixT/online_methods.html).

## Reference
[1] GTEx Consortium. "The Genotype-Tissue Expression (GTEx) pilot analysis: Multitissue gene regulation in humans." Science 348.6235 (2015): 648-660.

[2] Cao, Shaolong, et al. "Differing total mRNA expression shapes the molecular and clinical phenotype of cancer." bioRxiv (2020). https://www.biorxiv.org/content/10.1101/2020.09.30.306795v1.

## Acknowledgements
GTEx Project - The Genotype-Tissue Expression (GTEx) Project was supported by the Common Fund of the Office of the Director of the National Institutes of Health, and by NCI, NHGRI, NHLBI, NIDA, NIMH, and NINDS. The data used for the analyses described in this manuscript were obtained from the GTEx Portal ([https://gtexportal.org/home/](https://gtexportal.org/home/)) with dbGaP accession number phs000424.vN.pN.