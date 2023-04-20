# benchdamic
BENCHmarking of Differential Abundance detection methods for MICrobial data. 

This package implements a series of analysis to help users to analyze microbiome data. The theoretical considerations supporting the package are described in:  
["Calgaro, M., Romualdi, C., Waldron, L., Risso, D., and Vitulo, N. Assessment of statistical methods from single cell, bulk RNA-seq, and metagenomics applied to microbiome data. Genome Biol 21, 191 (2020)"](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02104-1) 

The package structure and its features are described in:
["Calgaro, M., Romualdi, C., Risso, D., and Vitulo, N. benchdamic: benchmarking of differential abundance methods on microbiome data. Bioinformatics 39 (2023)"](https://doi.org/10.1093/bioinformatics/btac778) 

Not only does the package structure allow the users to test a variety of commonly used methods for differential abundance analysis, but it also enables them to set benchmarks including custom methods on their datasets. Performances of each method are evaluated with respect to i) suitability of distributional assumptions, ii) ability to control false discoveries, iii) concordance of the findings, and iv) enrichment of differentially abundant microbial species in specific conditions. Each step of the assessment is flexible when it comes to the choice of differential abundance methods, their parameters, and input data types. Various graphic outputs lead the users to an informed decision when evaluating the most suitable method to use for their data.

## Installation

If you want to install the development version of `benchdamic` from GitHub, you can use:

```{r}
library(devtools)
install_github("mcalgaro93/benchdamic")
```

To install the released version of `benchdamic` from Bioconductor, you can use:

```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("benchdamic")
```

Please note that due to the large number of methods used by benchdamic, not all of their dependencies may already be present in some environments (e.g. RcppGSL R-package). Please refer to documentation of missing dependencies on how to install them. If problems persist, a possible solution is using the [Bioconductor docker image](https://www.bioconductor.org/help/docker/).
