# LPEadj
A correction for the LPE gene expression test which removes the bias from the standard error adjustment

## Qualification  
LPE and LPEadj are both rather dated methods for testing of gene expression. I currently recommend 
the R package limma for testing of gene expression for both microarrays and RNASeq data when applicable.

## Installation  
library(devtools)  
install_github("carlmurie/LPEadj")  

## Introduction  
  
LPEadj is a correction for the LPE statistical test [1] which consists of two additions to
the LPE method. The LPE package documentation is still correct with the exception of the
two additions listed below and should be consulted for more information on the LPE method.
The correction is in two parts. See [3] for more information on the correction.  

1. The LPEadj method discontinues the LPE practice of setting all variances below the
maximum variance in the ordered distribution of variances to the maximum variance.
In certain cases this practice can set many variances to the maximum and lower the
performance of this algorithm. If the assumption that there are only a few low variances
to be adjusted is correct then it may be safe to use this procedure. This option is
controlled by the doMax parameter (default is FALSE).  

2. The LPEadj method replaces the π/2 adjustment of the variance with an empirically
estimated adjustment based on the sample size of the data. The empirical adjustment
values have been estimated for replicate sizes up to 10. This option is controlled by
the doAdj parameter (default is TRUE).

## Citation  
Murie C, Nadon R, A correction for estimating error when using the Local Pooled Error
statistical test Bioinformatics. 2008 Aug 1;24(15):1735-6. doi: 10.1093/bioinformatics/btn211
