#' @title norm_DESeq2
#'
#' @importFrom DESeq2 estimateSizeFactors sizeFactors
#' @importFrom phyloseq sample_sums taxa_are_rows otu_table
#' @importFrom stats quantile median
#' @export
#' @description
#' Calculate normalization factors from a phyloseq object to scale the raw
#' library sizes. Inherited from DESeq2
#' \code{\link[DESeq2]{estimateSizeFactors}} function.
#'
#' @param object phyloseq object containing the counts to be normalized.
#' @param method Method for estimation: either \code{"ratio"},
#' \code{"poscounts"}, or \code{"iterate"}. \code{"ratio"} uses the standard
#' median ratio method introduced in DESeq. The size factor is the median ratio
#' of the sample over a "pseudosample": for each gene, the geometric mean of all
#' samples. \code{"poscounts"} and \code{"iterate"} offer alternative
#' estimators, which can be used even when all features contain a sample with a
#' zero (a problem for the default method, as the geometric mean becomes zero,
#' and the ratio undefined). The \code{"poscounts"} estimator deals with
#' a feature with some zeros, by calculating a modified geometric mean by taking
#' the n-th root of the product of the non-zero counts. This evolved out of use
#' cases with Paul McMurdie's phyloseq package for metagenomic samples. The
#' \code{"iterate"} estimator iterates between estimating the dispersion with a
#' design of ~1, and finding a size factor vector by numerically optimizing the
#' likelihood of the ~1 model.
#' @param ... other parameters for DESeq2
#' \code{\link[DESeq2]{estimateSizeFactors}} function.
#'
#' @return A new column containing the chosen DESeq2-based normalization factors
#' is added to the phyloseq \code{sample_data} slot.
#'
#' @seealso \code{\link[DESeq2]{estimateSizeFactors}} for details.
#' \code{\link{setNormalizations}} and \code{\link{runNormalizations}} to fastly
#' set and run normalizations.
#'
#' @examples
#' set.seed(1)
#' # Create a very simple phyloseq object
#' counts <- matrix(rnbinom(n = 60, size = 3, prob = 0.5), nrow = 10, ncol = 6)
#' metadata <- data.frame("Sample" = c("S1", "S2", "S3", "S4", "S5", "S6"),
#'                        "group" = as.factor(c("A", "A", "A", "B", "B", "B")))
#' ps <- phyloseq::phyloseq(phyloseq::otu_table(counts, taxa_are_rows = TRUE),
#'                          phyloseq::sample_data(metadata))
#'
#' # Calculate the normalization factors
#' ps_NF <- norm_DESeq2(object = ps, method = "poscounts")
#' # The phyloseq object now contains the normalization factors:
#' normFacts <- phyloseq::sample_data(ps_NF)[, "NF.poscounts"]
#' head(normFacts)
#'
#' # VERY IMPORTANT: to convert normalization factors to scaling factors divide
#' # them by the library sizes and renormalize.
#' scaleFacts = normFacts / phyloseq::sample_sums(ps_stool_16S)
#' # Renormalize: multiply to 1
#' scaleFacts = scaleFacts/exp(mean(log(scaleFacts)))

norm_DESeq2 <- function(object, method = c("ratio", "poscounts", "iterate"),
                        ...){
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    ## Calculate size factors
    obj <- phyloseq::phyloseq_to_deseq2(object, design = ~ 1)
    if(missing(method))
        stop("Please supply a normalization method between 'ratio', 'poscounts'
            or 'iterate'.")
    normFacts <- DESeq2::sizeFactors(DESeq2::estimateSizeFactors(obj,
                                                                 type = method, ...))
    phyloseq::sample_data(object)[,paste("NF", method, sep = ".")] <- normFacts
    return(object)
}# END - function: norm_DESeq2
