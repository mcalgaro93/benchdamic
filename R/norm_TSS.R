#' @title norm_TSS
#'
#' @importFrom phyloseq otu_table sample_data
#' @importFrom SummarizedExperiment colData
#' @export
#' @description
#' Calculate the Total Sum Scaling (TSS) factors for a phyloseq or a 
#' TreeSummarizedExperiment object, i.e. the library sizes. If the counts are 
#' divided by the scaling factors, a relative abundance is obtained.
#'
#' @inheritParams get_counts_metadata
#' @param method normalization method to be used.
#' @param verbose an optional logical value. If \code{TRUE}, information about
#' the steps of the algorithm is printed. Default \code{verbose = TRUE}.
#'
#' @return A new column containing the TSS scaling factors is added to the 
#' \code{sample_data} slot of the phyloseq object or the \code{colData} slot of 
#' the TreeSummarizedExperiment object.
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

norm_TSS <- function(object, assay_name = "counts", method = "TSS", 
    verbose = TRUE) {
    counts_and_metadata <- get_counts_metadata(object, assay_name = assay_name)
    counts <- counts_and_metadata[[1]]
    is_phyloseq <- counts_and_metadata[[3]]
    normFacts <- colSums(counts)
    NF.col <- paste("NF", method, sep = ".")
    if(is_phyloseq){
        phyloseq::sample_data(object)[,NF.col] <- normFacts
    } else {
        SummarizedExperiment::colData(object)[, NF.col] <- normFacts
    }
    if(verbose)
        message(NF.col, " column has been added.")
    return(object)
}# END - function: norm_TSS
