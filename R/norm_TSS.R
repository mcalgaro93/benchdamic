#' @title norm_TSS
#'
#' @importFrom phyloseq taxa_are_rows sample_sums
#' @export
#' @description
#' Calculate the raw library sizes from a phyloseq object. If used to divide
#' counts, known as Total Sum Scaling normalization (TSS).
#'
#' @param object phyloseq object containing the counts to be normalized.
#' @param method normalization method to be used.
#'
#' @return A new column containing the TSS scaling factors is added to the
#' phyloseq \code{sample_data} slot.
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
#' ps_NF <- norm_TSS(object = ps)
#' # The phyloseq object now contains the scaling factors:
#' scaleFacts <- phyloseq::sample_data(ps_NF)[, "NF.TSS"]
#' head(scaleFacts)
#'
#' # VERY IMPORTANT: to convert scaling factors to normalization factors
#' # multiply them by the library sizes and renormalize.
#' normFacts = scaleFacts * phyloseq::sample_sums(ps_stool_16S)
#' # Renormalize: multiply to 1
#' normFacts = normFacts/exp(colMeans(log(normFacts)))

norm_TSS <- function(object, method = "TSS")
{
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    normFacts <- 1/phyloseq::sample_sums(object)
    phyloseq::sample_data(object)[,"NF.TSS"] <- normFacts
    return(object)
}# END - function: norm_TSS
