#' @title DA_dearseq
#'
#' @importFrom dearseq dear_seq
#' @importFrom SummarizedExperiment assays
#' @export
#' @description
#' Fast run for dearseq differential abundance detection method.
#'
#' @inheritParams DA_edgeR
#' @param covariates a character vector containing the colnames of the 
#' covariates to include in the model.
#' @param variables2test a character vector containing the colnames of the 
#' variable of interest.
#' @param test a character string indicating which method to use to approximate 
#' the variance component score test, either 'permutation' or 'asymptotic'
#' (default \code{test = "permutation"}).
#' @inheritParams dearseq::dear_seq
#'
#' @return A list object containing the matrix of p-values `pValMat`,
#' a matrix of summary statistics for each tag `statInfo` which are still the 
#' p-values as this method does not produce other values, and a suggested 
#' `name` of the final object considering the parameters passed to the 
#' function.
#'
#' @seealso \code{\link[dearseq]{dear_seq}} for analysis of differential 
#' expression/abundance through a variance component test.
#'
#' @examples
#' set.seed(1)
#' # Create a very simple phyloseq object
#' counts <- matrix(rnbinom(n = 60, size = 3, prob = 0.5), nrow = 10, ncol = 6)
#' metadata <- data.frame("Sample" = c("S1", "S2", "S3", "S4", "S5", "S6"),
#'     "group" = as.factor(c("A", "A", "A", "B", "B", "B")))
#' ps <- phyloseq::phyloseq(phyloseq::otu_table(counts, taxa_are_rows = TRUE),
#'      phyloseq::sample_data(metadata))
#' # Differential abundance
#' DA_dearseq(object = ps, pseudo_count = FALSE, covariates = NULL, 
#'     variables2test = "group", sample_group = NULL, test = "asymptotic",
#'     preprocessed = FALSE, verbose = TRUE)

DA_dearseq <- function(object, assay_name = "counts", pseudo_count = FALSE, 
    covariates = NULL, variables2test = NULL, sample_group = NULL,
    test = c("permutation", "asymptotic"), preprocessed = FALSE, 
    n_perm = 1000, verbose = TRUE){
    counts_and_metadata <- get_counts_metadata(object, assay_name = assay_name)
    counts <- counts_and_metadata[[1]]
    metadata <- counts_and_metadata[[2]]
    is_phyloseq <- counts_and_metadata[[3]]
    # Name building
    name <- "dearseq"
    method <- "DA_dearseq"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        message("Adding a pseudo count... \n")
        counts <- counts + 1
        name <- paste(name, ".pseudo", sep = "")
    }
    # Check the assay
    if (!is_phyloseq){
        if(verbose)
            message("Using the ", assay_name, " assay.")
        name <- paste(name, ".", assay_name, sep = "")
    } 
    # Check the variable2test
    if(is.null(variables2test)){
        stop(method, "\n", 
            "variables2test: please choose a variable to test.")
    }
    # Check the sample_group
    if(!is.null(sample_group)){
        sample_group <- metadata[, sample_group]
        name <- paste(name, ".", "repeated", sep = "")
    }
    # Check the test
    if(length(test) != 1){
        stop(method, "\n", 
           "test: please choose one test for this istance of", 
           " differential abundance analysis.")
    } 
    name <- paste(name, ".", test, sep = "")
    if(test == "permutation"){
        name <- paste(name, ".", n_perm, sep = "")
    }
    # Check the preprocessed
    if(!preprocessed){
        if(assay_name != "counts"){
            warning(method, "\n", 
                "The 'assay_name' = ", assay_name, " but 'preprocessed' = ", 
                " FALSE. Be sure the data are not preprocessed.")
        }
    } else {
        name <- paste(name, ".", "preprocessed", sep = "")
    }
    input <- SummarizedExperiment::SummarizedExperiment(
        assays = list("counts" = counts),
        colData = metadata)
    if(verbose){
        
        res <- dearseq::dear_seq(object = input, covariates = covariates,
            variables2test = variables2test, 
            sample_group = sample_group, 
            which_test = test, parallel_comp = FALSE, n_perm = n_perm, 
            progressbar = TRUE, preprocessed = preprocessed)
    } else {
        res <- suppressWarnings(suppressMessages(
            dearseq::dear_seq(object = input, covariates = covariates,
            variables2test = variables2test, 
            sample_group = sample_group,
            which_test = test, parallel_comp = FALSE, n_perm = n_perm, 
            progressbar = FALSE, preprocessed = preprocessed)))
    }
    statInfo <- pValMat <- as.data.frame(res["pvals"])
    colnames(statInfo) <- colnames(pValMat) <- c("rawP", "adjP")
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
}# END - function: DA_dearseq

#' @title set_dearseq
#'
#' @export
#' @description
#' Set the parameters for dearseq differential abundance detection method.
#'
#' @inheritParams DA_dearseq
#' @param expand logical, if TRUE create all combinations of input parameters
#' (default \code{expand = TRUE}).
#'
#' @return A named list containing the set of parameters for \code{DA_dearseq}
#' method.
#'
#' @seealso \code{\link{DA_dearseq}}
#'
#' @examples
#' # Set some basic combinations of parameters for dearseq 
#' base_dearseq <- set_dearseq(pseudo_count = FALSE, variables2test = "group",
#'     test = c("permutation", "asymptotic"), expand = TRUE)

set_dearseq <- function(assay_name = "counts", pseudo_count = FALSE, 
    covariates = NULL, variables2test = NULL, sample_group = NULL,
    test = c("permutation", "asymptotic"), preprocessed = FALSE, 
    n_perm = 1000, expand = TRUE) {
    method <- "DA_dearseq"
    if (is.null(assay_name)) {
        stop(method, "\n", "'assay_name' is required (default = 'counts').")
    }
    if (!is.logical(pseudo_count)) {
        stop(method, "\n", "'pseudo_count' must be logical.")
    }
    if (!is.logical(preprocessed)) {
        stop(method, "\n", "'preprocessed' must be logical.")
    }
    if (is.null(variables2test)) {
        stop(method, "\n", "'variables2test' must be specified.")
    }
    if (sum(!is.element(test, c("permutation", "asymptotic"))) > 0) {
        stop(method, "\n", 
             "One or more 'test' are not available for dearseq",
             " Please choose between 'permutation' or 'asymptotic'.")
    }
    if (expand) {
        parameters <- expand.grid(method = method, assay_name = assay_name, 
            pseudo_count = pseudo_count, test = test, 
            preprocessed = preprocessed, n_perm = n_perm,
            stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, assay_name = assay_name, 
            pseudo_count = pseudo_count, test = test, 
            preprocessed = preprocessed, n_perm = n_perm)
    }
    # data.frame to list
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        x <- append(x = x, values = list("covariates" = covariates,
            "variables2test" = variables2test, "sample_group" = sample_group), 
            after = 3)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}
