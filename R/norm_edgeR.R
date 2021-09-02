#' @title norm_edgeR
#'
#' @importFrom edgeR calcNormFactors
#' @importFrom phyloseq sample_sums taxa_are_rows otu_table
#' @importFrom stats quantile
#' @export
#' @description
#' Calculate scaling factors from a phyloseq object to scale the raw library
#' sizes. Inherited from edgeR \code{\link{calcNormFactors}} function.
#'
#' @param object a phyloseq object containing the counts to be normalized.
#' @param method normalization method to be used. Choose between \code{TMM},
#' \code{TMMwsp}, \code{RLE}, \code{upperquartile}, \code{posupperquartile} or
#' \code{none}.
#' @inheritParams edgeR::calcNormFactors
#'
#' @return A new column containing the chosen edgeR-based scaling factors is
#' added to the phyloseq \code{sample_data} slot. The effective library sizes
#' to use in downstream analysis must be multiplied by the normalization factors.
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
#' # Calculate the scaling factors
#' ps_NF <- norm_edgeR(object = ps, method = "TMM")
#'
#' # The phyloseq object now contains the scaling factors:
#' scaleFacts <- phyloseq::sample_data(ps_NF)[, "NF.TMM"]
#' head(scaleFacts)
#'
#' # VERY IMPORTANT: to convert scaling factors to normalization factors
#' # multiply them by the library sizes and renormalize.
#' normFacts = scaleFacts * phyloseq::sample_sums(ps_stool_16S)
#' # Renormalize: multiply to 1
#' normFacts = normFacts/exp(colMeans(log(normFacts)))

norm_edgeR <- function(object, method = c("TMM", "TMMwsp", "RLE",
                                          "upperquartile", "posupperquartile", "none"), refColumn = NULL,
                       logratioTrim = .3, sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e10,
                       p = 0.75, ...){
    counts <- as(phyloseq::otu_table(object), "matrix")
    if (!phyloseq::taxa_are_rows(object))
        counts <- t(counts)
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
    phyloseq::sample_data(object)[,paste("NF", method, sep = ".")] <- normFacts
    return(object)
}# END - function: norm_edgeR
