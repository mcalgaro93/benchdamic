#' @title norm_CSS
#'
#' @importFrom metagenomeSeq newMRexperiment
#' @importFrom phyloseq otu_table sample_data
#' @importFrom stats median
#' @importFrom SummarizedExperiment colData
#' @export
#' @description
#' Calculate normalization factors from a phyloseq or TreeSummarizedExperiment
#' object. Inherited from metagenomeSeq 
#' \code{\link[metagenomeSeq]{calcNormFactors}} function which performs the 
#' Cumulative Sum Scaling normalization.
#'
#' @inheritParams get_counts_metadata
#' @param method normalization method to be used (only CSS).
#' @param verbose an optional logical value. If \code{TRUE}, information about
#' the steps of the algorithm is printed. Default \code{verbose = TRUE}.
#'
#' @return A new column containing the CSS normalization factors is added to 
#' the \code{sample_data} slot of the phyloseq object or the \code{colData} 
#' slot of the TreeSummarizedExperiment object.
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
#' # Calculate the normalization factors
#' ps_NF <- norm_CSS(object = ps, method = "CSS")
#' # The phyloseq object now contains the normalization factors:
#' CSSFacts <- phyloseq::sample_data(ps_NF)[, "NF.CSS"]
#' head(CSSFacts)
#' 
#' # VERY IMPORTANT: metagenomeSeq uses scaling factors to normalize counts 
#' # (even though they are called normalization factors). These factors are used
#' # internally by a regression model. To make CSS size factors available for 
#' # edgeR, we need to transform them into normalization factors. This is 
#' # possible by dividing the factors for the library sizes and renormalizing. 
#' normCSSFacts = CSSFacts / colSums(phyloseq::otu_table(ps_stool_16S))
#' # Renormalize: multiply to 1
#' normFacts = normCSSFacts/exp(colMeans(log(normCSSFacts)))

norm_CSS <- function(object, assay_name = "counts", method = "CSS", 
    verbose = TRUE)
{
    if(method != "CSS")
        stop("The only accepted method is 'CSS'.")
    counts_and_metadata <- get_counts_metadata(object, assay_name = assay_name)
    counts <- counts_and_metadata[[1]]
    is_phyloseq <- counts_and_metadata[[3]]
    obj <- metagenomeSeq::newMRexperiment(counts = counts)
    if(verbose){
        normFacts <- metagenomeSeq::calcNormFactors(obj = obj)
    } else {
        normFacts <- suppressMessages(metagenomeSeq::calcNormFactors(obj = obj))
    }
    normFacts <- drop(as(normFacts, "matrix"))
    NF.col <- paste("NF", method, sep = ".")
    if(is_phyloseq){
        phyloseq::sample_data(object)[, NF.col] <- normFacts
    } else {
        SummarizedExperiment::colData(object)[, NF.col] <- normFacts
    }
    if(verbose)
        message(NF.col, " column has been added.")
    return(object)
}# END - function: norm_CSS
