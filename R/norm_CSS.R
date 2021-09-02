#' @title norm_CSS
#'
#' @importFrom metagenomeSeq newMRexperiment
#' @importFrom phyloseq taxa_are_rows otu_table
#' @importFrom stats median
#' @export
#' @description
#' Calculate scaling factors from a phyloseq object to scale the raw library
#' sizes. Inherited from metagenomeSeq `calcNormFactors` function which performs
#' the Cumulative Sum Scaling normalization.
#'
#' @param object phyloseq object containing the counts to be normalized.
#' @param method  normalization scaling parameter (default
#' \code{method = "default"}). If \code{"median"}, the median of the
#' normalization factors is used as scaling (Paulson et al. 2013).
#'
#' @return A new column containing the CSS scaling factors is added to the
#' phyloseq \code{sample_data} slot.
#'
#' @seealso \code{\link[metagenomeSeq]{calcNormFactors}} for details.
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
#' # Calculate the scaling factors
#' ps_NF <- norm_CSS(object = ps, method = "median")
#' # The phyloseq object now contains the scaling factors:
#' scaleFacts <- phyloseq::sample_data(ps_NF)[, "NF.CSSmedian"]
#' head(scaleFacts)
#'
#' # VERY IMPORTANT: to convert scaling factors to normalization factors
#' # multiply them by the library sizes and renormalize.
#' normFacts = scaleFacts * phyloseq::sample_sums(ps_stool_16S)
#' # Renormalize: multiply to 1
#' normFacts = normFacts/exp(colMeans(log(normFacts)))

norm_CSS <- function(object, method = "default")
{
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    counts <- as(phyloseq::otu_table(object), "matrix")
    obj <- metagenomeSeq::newMRexperiment(counts = counts)
    normFacts <- metagenomeSeq::calcNormFactors(obj = obj)
    normFacts <- drop(as(normFacts, "matrix"))
    # Default: log2(normFacts/1000 + 1)
    # Original metagenomeSeq paper: log2(normFacts/median(libsize) +1)
    if(method == "default")
        normFacts <- log2(normFacts/1000 + 1)
    else if (method == "median")
        normFacts <- log2(normFacts/stats::median(normFacts) + 1)
    else stop("Please choose a scaling method between 'default' or 'median'.")
    # Remember to useCSSoffset = FALSE in fitZig function
    phyloseq::sample_data(object)[,paste0("NF.CSS", method)] <-
        normFacts
    return(object)
}# END - function: norm_CSS
