#' @title norm_edgeR
#'
#' @importFrom edgeR calcNormFactors
#' @importFrom phyloseq otu_table sample_data
#' @importFrom stats quantile
#' @importFrom SummarizedExperiment colData
#' @export
#' @description
#' Calculate normalization factors from a phyloseq or TreeSummarizedExperiment
#' object. Inherited from edgeR \code{\link{calcNormFactors}} function.
#'
#' @inheritParams get_counts_metadata
#' @param method normalization method to be used. Choose between \code{TMM},
#' \code{TMMwsp}, \code{RLE}, \code{upperquartile}, \code{posupperquartile} or
#' \code{none}.
#' @param verbose an optional logical value. If \code{TRUE}, information about
#' the steps of the algorithm is printed. Default \code{verbose = TRUE}.
#' @inheritParams edgeR::calcNormFactors
#'
#' @return A new column containing the chosen edgeR-based normalization factors 
#' is added to the \code{sample_data} slot of the phyloseq object or the 
#' \code{colData} slot of the TreeSummarizedExperiment object.
#' @seealso \code{\link[edgeR]{calcNormFactors}} for details.
#'
#' @seealso \code{\link{setNormalizations}} and \code{\link{runNormalizations}}
#' to fastly set and run normalizations.
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
#' ps_NF <- norm_edgeR(object = ps, method = "TMM")
#'
#' # The phyloseq object now contains the normalization factors:
#' normFacts <- phyloseq::sample_data(ps_NF)[, "NF.TMM"]
#' head(normFacts)
#' 
#' # VERY IMPORTANT: edgeR uses normalization factors to normalize library sizes
#' # not counts. They are used internally by a regression model. To make edgeR 
#' # normalization factors available for other methods, such as DESeq2 or other 
#' # DA methods based on scaling or size factors, we need to transform them into
#' # size factors. This is possible by multiplying the factors for the library 
#' # sizes and renormalizing. 
#' normLibSize = normFacts * colSums(phyloseq::otu_table(ps_stool_16S))
#' # Renormalize: multiply to 1
#' sizeFacts = normLibSize/exp(colMeans(log(normLibSize)))

norm_edgeR <- function(object, assay_name = "counts", method = c("TMM", 
    "TMMwsp", "RLE", "upperquartile", "posupperquartile", "none"), 
    refColumn = NULL, logratioTrim = .3, sumTrim = 0.05, doWeighting = TRUE, 
    Acutoff = -1e10, p = 0.75, verbose = TRUE, ...){
    counts_and_metadata <- get_counts_metadata(object, assay_name = assay_name)
    counts <- counts_and_metadata[[1]]
    is_phyloseq <- counts_and_metadata[[3]]
    if (is.na(method))
        stop("Please, supply a valid method between 'TMM', 'TMMwsp', 'RLE',
            'upperquartile', 'posupperquartile' or 'none'")
    else if (method == "posupperquartile"){
        scaledCounts <- t(counts) / colSums(counts)
        tmpNF <- apply(scaledCounts, MARGIN = 1L, FUN = function(x)
            quantile(x[x != 0], probs = .75))
        normFacts <- tmpNF/exp(mean(log(tmpNF)))
    } else {
        normFacts <- edgeR::calcNormFactors(counts, method = method, refColumn =
            refColumn, logratioTrim = logratioTrim, sumTrim = sumTrim,
            doWeighting = doWeighting, Acutoff = Acutoff, p = p, ...) }
    # END - ifelse: posupperquartile = upperquartile only of non-zero counts
    if (all(is.na(normFacts)))
        stop("Failed to compute normalization factors!")
    NF.col <- paste("NF", method, sep = ".")
    if(is_phyloseq){
        phyloseq::sample_data(object)[, NF.col] <- normFacts
    } else {
        SummarizedExperiment::colData(object)[, NF.col] <- normFacts
    }
    if(verbose)
        message(NF.col, " column has been added.")
    return(object)
}# END - function: norm_edgeR
