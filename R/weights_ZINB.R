#' @title weights_ZINB
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data
#' @importFrom stats model.matrix
#' @importFrom zinbwave zinbFit computeObservationalWeights
#' @importFrom BiocParallel SerialParam
#' @export
#' @description
#' Computes the observational weights of the counts under a zero-inflated
#' negative binomial (ZINB) model. For each count, the ZINB distribution is
#' parametrized by three parameters: the mean value and the dispersion of the
#' negative binomial distribution, and the probability of the zero component.
#'
#' @param object phyloseq object containing the counts and the sample data.
#' @param design character name of the metadata columns, formula, or design
#' matrix with rows corresponding to samples and columns to coefficients to be
#' estimated (the user needs to explicitly include the intercept in the design).
#' @inheritParams zinbwave::zinbFit
#'
#' @return A matrix of weights.
#'
#' @seealso \code{\link[zinbwave]{zinbFit}} for zero-inflated negative binomial
#' parameters' estimation and
#' \code{\link[zinbwave]{computeObservationalWeights}} for weights extraction.
#'
#' @examples
#' set.seed(1)
#' # Create a very simple phyloseq object
#' counts <- matrix(rnbinom(n = 60, size = 3, prob = 0.5), nrow = 10, ncol = 6)
#' metadata <- data.frame("Sample" = c("S1", "S2", "S3", "S4", "S5", "S6"))
#' ps <- phyloseq::phyloseq(phyloseq::otu_table(counts, taxa_are_rows = TRUE),
#'     phyloseq::sample_data(metadata))
#' # Calculate the ZINB weights
#' zinbweights <- weights_ZINB(object = ps, K = 0, design = "~ 1")

weights_ZINB <- function(object, design, K = 0, commondispersion = TRUE,
    zeroinflation = TRUE, verbose = FALSE, ...){
    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)
    if(is.character(design)){
        if(grepl("~", design))
            design <- as.formula(design)
        else design <- as.formula(paste0("~", design))}
    if(is(design, "formula"))
        design <- stats::model.matrix(object = design,
            data = data.frame(metadata))
    zinbmodel <- zinbwave::zinbFit(Y = counts, X = design, K = K,
        commondispersion = commondispersion, zeroinflation = TRUE, verbose =
        verbose, BPPARAM = BiocParallel::SerialParam(), ... = ...)
    w <- zinbwave::computeObservationalWeights(model = zinbmodel, x = counts)
    return(w)
}# END - function: weights_ZINB
