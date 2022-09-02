#' @title DA_NOISeq
#'
#' @importFrom NOISeq noiseqbio
#' @importFrom SummarizedExperiment assays
#' @export
#' @description
#' Fast run for NOISeqBIO differential abundance detection method. It computes 
#' differential expression between two experimental conditions.
#'
#' @inheritParams DA_edgeR
#' @inheritParams NOISeq::noiseqbio
#' @param contrast character vector with exactly, three elements: a string 
#' indicating the name of factor whose levels are the conditions to be 
#' compared, the name of the level of interest, and the name of the other 
#' level. 
#'
#' @return A list object containing the matrix of p-values `pValMat`,
#' a matrix of summary statistics for each tag `statInfo`, and a suggested 
#' `name` of the final object considering the parameters passed to the 
#' function. NOISeq does not produce p-values but an estimated probability of 
#' differential expression for each feature. Note that these probabilities are 
#' not equivalent to p-values. The higher the probability, the more likely that 
#' the difference in abundance is due to the change in the experimental 
#' condition and not to chance... Hence, `pValMat` matrix is filled with 
#' \code{1 - prob} values which can be interpreted as 1 - FDR. Where FDR can 
#' be considered as an adjusted p-value (see NOISeq vignette).
#'
#' @seealso \code{\link[NOISeq]{noiseqbio}} for analysis of differential 
#' expression/abundance between two experimental conditions from read count 
#' data.
#'
#' @examples
#' set.seed(1)
#' # Create a very simple phyloseq object
#' counts <- matrix(rnbinom(n = 60, size = 3, prob = 0.5), nrow = 10, ncol = 6)
#' metadata <- data.frame("Sample" = c("S1", "S2", "S3", "S4", "S5", "S6"),
#'                        "group" = as.factor(c("A", "A", "A", "B", "B", "B")))
#' ps <- phyloseq::phyloseq(phyloseq::otu_table(counts, taxa_are_rows = TRUE),
#'                          phyloseq::sample_data(metadata))
#' # Differential abundance
#' DA_NOISeq(object = ps, pseudo_count = FALSE, contrast = c("group", "B", "A"),
#'     norm = "tmm", verbose = FALSE)

DA_NOISeq <- function(object, assay_name = "counts", pseudo_count = FALSE, 
    contrast = NULL, norm = c("rpkm", "uqua", "tmm", "n"), verbose = TRUE){
    counts_and_metadata <- get_counts_metadata(object, assay_name = assay_name)
    counts <- counts_and_metadata[[1]]
    metadata <- counts_and_metadata[[2]]
    is_phyloseq <- counts_and_metadata[[3]]
    # Name building
    name <- "NOISeq"
    method <- "DA_NOISeq"
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
    if(!is.character(contrast) | length(contrast) != 3){
        stop(method, "\n", 
            "contrast: please supply a character vector with exactly", 
            " three elements: a string indicating the name of factor whose",
            " levels are the conditions to be compared,",  
            " the name of the level of interest, and the", 
            " name of the other level.")
    }
    if(is.element(contrast[1], colnames(metadata))){
        if(!is.factor(metadata[, contrast[1]])){
            if(verbose){
                message("Converting variable ", contrast[1], " to factor.")
            }
            metadata[, contrast[1]] <- as.factor(metadata[, contrast[1]])
        }
        if(!is.element(contrast[2], levels(metadata[, contrast[1]])) | 
           !is.element(contrast[3], levels(metadata[, contrast[1]]))){
            stop(method, "\n", 
                 "contrast: ", contrast[2], " and/or ", contrast[3], 
                 " are not levels of ", contrast[1], " variable.")
        }
        if(verbose){
            message("Setting ", contrast[3], " the reference level for ", 
                    contrast[1], " variable.")
        }
        metadata[, contrast[1]] <- stats::relevel(metadata[, contrast[1]], 
            ref = contrast[3])
    }
    if (length(norm) > 1) {
        stop(method, "\n",
            "norm: please choose one normalization for this istance of", 
            " differential abundance analysis.")
    }
    this_norms <- c("rpkm", "uqua", "tmm", "n")
    if(!is.element(norm, this_norms)) {
        stop(method, "\n",
            norm, " normalization is not an available NOISeqBIO normalization.",
            " Choose between 'rpkm', 'uqua', 'tmm', or 'n'. To choose a custom",
            " normalization set norm to 'n' and provide already normalized",
            " data")
    }
    name <- paste(name, ".", norm, sep = "")
    input <- NOISeq::readData(data = counts, factors = metadata)
    if(verbose){
        res <- NOISeq::noiseqbio(input = input, k = 0.5, norm = norm,
            nclust = min(c(nrow(counts) - 1, 15)), factor = contrast[1], 
            conditions = contrast[2:3], filter = 0)
    } else { # We need to use invisible with capture.output to suppress cat()
        invisible(capture.output(
            res <- NOISeq::noiseqbio(input = input, k = 0.5, norm = norm,
                nclust = min(c(nrow(counts) - 1, 15)), factor = contrast[1], 
                conditions = contrast[2:3], filter = 0)))
    }
    statInfo <- as.data.frame(res@results[[1]])
    pValMat <- 1 - statInfo[,c("prob", "prob")]
    colnames(pValMat) <- c("rawP", "adjP")
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
}# END - function: DA_NOISeq

#' @title set_NOISeq
#'
#' @export
#' @description
#' Set the parameters for NOISeq differential abundance detection method.
#'
#' @inheritParams DA_NOISeq
#' @param expand logical, if TRUE create all combinations of input parameters
#' (default \code{expand = TRUE}).
#'
#' @return A named list containing the set of parameters for \code{DA_NOISeq}
#' method.
#'
#' @seealso \code{\link{DA_NOISeq}}
#'
#' @examples
#' # Set a basic combination of parameters for NOISeq with 'tmm' normalization
#' base_NOISeq <- set_NOISeq(pseudo_count = FALSE, norm = "tmm",
#'     contrast = c("group", "B", "A"), expand = FALSE)
#' # try many normalizations
#' many_NOISeq <- set_NOISeq(pseudo_count = FALSE, 
#'     norm = c("tmm", "uqua", "rpkm", "n"), contrast = c("group", "B", "A"))
set_NOISeq <- function(assay_name = "counts", pseudo_count = FALSE, 
    contrast = NULL, norm = c("rpkm", "uqua", "tmm", "n"), expand = TRUE) {
    method <- "DA_NOISeq"
    if (is.null(assay_name)) {
        stop(method, "\n", "'assay_name' is required (default = 'counts').")
    }
    if (!is.logical(pseudo_count)) {
        stop(method, "\n", "'pseudo_count' must be logical.")
    }
    if (is.null(contrast)) {
        stop(method, "\n", "'contrast' must be specified.")
    }
    if (!is.character(contrast) & length(contrast) != 3){
        stop(method, "\n", 
             "contrast: please supply a character vector with exactly", 
             " three elements: a string indicating the name of factor whose",
             " levels are the conditions to be compared,",  
             " the name of the level of interest, and the", 
             " name of the other level.")
    }
    if (sum(!is.element(norm, c("rpkm", "uqua", "tmm", "n"))) > 0) {
        stop(method, "\n", 
            "One or more elements into 'norm' are not available for NOISeqBIO.",
            " Please choose between 'rpkm', 'uqua', 'tmm', or 'n'. To choose a",
            " custom normalization set norm to 'n' and provide already", 
            " normalized data")
    }
    if (expand) {
        parameters <- expand.grid(method = method, assay_name = assay_name, 
            pseudo_count = pseudo_count, norm = norm, stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, assay_name = assay_name, 
            pseudo_count = pseudo_count, norm = norm)
    }
    # data.frame to list
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        x <- append(x = x, values = list("contrast" = contrast), after = 3)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}
