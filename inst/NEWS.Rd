\name{NEWS}
\title{News for Package \pkg{DeMixT}}

\section{Version 1.2.5}{
\itemize{
\item
Added pi01 and pi02 as input values for users to initialize the
proportion estimation.
\item
Added nspikein as an input value in the DeMixT, DeMixT_S1 and DeMixT_GS 
functions to specify how many spike-in normal reference samples need to be 
generated; Setting nspikein at null as a default value, the number of 
spike-in normal reference samples equal the min(200, 0.3 x My), 
where My is the number of mixed samples; By setting nspikein equals 
0, no spike-in normal reference will be generated; If the input value 
of data.N2 is not null, nspikein will be forced to be 0.
\item
Added DeMixT_GS function, new proposed gene selection method which 
applies profile likelihood, for proportion estimation.
\item
Added simulate_2comp function for users to simulate test data for
2-component de-convolution.
\item
Added simulate_3comp function for users to simulate test data for
3-component de-convolution.
\item
Added row names and column names for all output values.
\item
Added gene.selection.method as an input value for DeMixT function. 
The default is 'GS'.
\item
Added ngene.Profile.selected as an input value for DeMixT function.
The default is NA.}
}

\section{Version 0.99.0}{
\itemize{
\item
New package DeMixT, for cell type-specific deconvolution of heterogeneous 
tumor samples with two or three components using expression data from 
RNAseq or microarray platforms
}}
