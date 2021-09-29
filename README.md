# benchdamic
BENCHmarking of Differential Abundance detection methods for MICrobial data. 

This package implements a series of analysis to help users to analyze microbiome data.The steps of the package are described in details in the paper:
[Calgaro, M., Romualdi, C., Waldron, L. et al. Assessment of statistical methods from single cell, bulk RNA-seq, and metagenomics applied to microbiome data. Genome Biol 21, 191 (2020)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02104-1)

Not only does the package structure allow the users to test a variety of commonly used methods for differential abundance analysis, but it also enables them to set benchmarks including custom methods on their datasets. Performances of each method are evaluated with respect to i) suitability of distributional assumptions, ii) ability to control false discoveries, iii) concordance of the findings, and iv) enrichment of differentially abundant microbial species in specific conditions. Each step of the assessment is flexible when it comes to the choice of differential abundance methods, their parameters, and input data types. Various graphic outputs lead the users to an informed decision when evaluating the most suitable method to use for their data.

## Installation

If you want to install the development version of `benchdamic` from GitHub, you can do so with the following.

```{r}
library(devtools)
install_github("mcalgaro93/benchdamic")
```
