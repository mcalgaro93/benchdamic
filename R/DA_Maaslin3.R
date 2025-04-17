#' @title DA_maaslin3
#'
#' @importFrom maaslin3 maaslin3
#' @importFrom SummarizedExperiment assays
#' @importFrom phyloseq otu_table sample_data phyloseq taxa_are_rows
#' @export
#' @description
#' Fast run for maaslin3 differential abundance detection method.
#'
#' @inheritParams DA_edgeR
#' @param contrast character vector with exactly, three elements: a string 
#' indicating the name of factor whose levels are the conditions to be 
#' compared, the name of the level of interest, and the name of the other 
#' level. 
#' @inheritParams maaslin3::maaslin3
#'
#' @return A list object containing the matrix of p-values `pValMat`,
#' a matrix of summary statistics for each tag `statInfo`, and a suggested 
#' `name` of the final object considering the parameters passed to the 
#' function.
#'
#' @seealso \code{\link[maaslin3]{maaslin3}}.
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
#' DA_Maaslin2(object = ps, normalization = "CLR", transform = "NONE",
#'     analysis_method = "LM", correction = "BH", random_effects = NULL,
#'     fixed_effects = "group", contrast = c("group", "B", "A"),
#'     verbose = FALSE)

DA_maaslin3 <- function(object, assay_name = "counts", 
    normalization = c("TSS", "CLR", "NONE"), 
    transform = c("LOG", "NONE"), 
    median_comparison_abundance = c(TRUE, FALSE),
    subtract_median = c(TRUE, FALSE),
    correction = "BH", random_effects = NULL, fixed_effects = NULL,
    small_random_effects = FALSE,
    contrast = NULL, reference = NULL, verbose = TRUE){
    counts_and_metadata <- get_counts_metadata(object, assay_name = assay_name)
    counts <- counts_and_metadata[[1]]
    metadata <- counts_and_metadata[[2]]
    is_phyloseq <- counts_and_metadata[[3]]
    # Name building
    name <- "maaslin3"
    method <- "DA_maaslin3"
    # Check the assay
    if (!is_phyloseq){
        if(verbose)
            message("Using the ", assay_name, " assay.")
        name <- paste(name, ".", assay_name, sep = "")
    } 
    # Check normalization
    if(sum(!is.element(normalization, 
        c("TSS", "CLR", "NONE"))) > 0){
        stop(method, "\n", 
             "normalization: please choose one normalization between 'TSS',",
             " 'CLR', or 'NONE'.")
    }
    # Check transform
    if(sum(!is.element(transform, c("LOG", "NONE"))) > 0){
        stop(method, "\n", 
             "transform: please choose one transfomation between 'LOG'",
             " or 'NONE'.")
    }
    # Remove senseless combinations
    if(!(median_comparison_abundance == subtract_median))
        stop(method, "\n", 
             "median_comparison_abundance must equal subtract_median")
    if(normalization == 'CLR' & transform != 'NONE')
        stop(method, "\n", 
             "if normalization is CLR, transform must be NONE")
    name <- paste(name, ".", normalization, "norm.", transform, "trans.", 
        "medCompare", median_comparison_abundance, sep = "")
    if(!is.character(contrast) | length(contrast) != 3)
        stop(method, "\n", 
             "contrast: please supply a character vector with exactly", 
             " three elements: the name of a variable used in",  
             " 'design', the name of the level of interest, and the", 
             " name of the reference level.")
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
    if(verbose){
        res <- maaslin3(input_data = t(counts), input_metadata = metadata, 
            output = tempdir(), normalization = normalization, 
            transform = transform, standardize = TRUE, 
            median_comparison_abundance = median_comparison_abundance,
            subtract_median = subtract_median,
            max_significance = 0, small_random_effects = small_random_effects,
            random_effects = random_effects, fixed_effects = fixed_effects, 
            correction = correction, plot_summary_plot = FALSE, 
            plot_associations = FALSE)
    } else {
        utils::capture.output(file = tempfile(),
        res <- maaslin3(input_data = t(counts), input_metadata = metadata, 
            output = tempdir(), normalization = normalization, 
            transform = transform, standardize = TRUE, 
            median_comparison_abundance = median_comparison_abundance,
            subtract_median = subtract_median,
            max_significance = 0, small_random_effects = small_random_effects,
            random_effects = random_effects, fixed_effects = fixed_effects, 
            correction = correction, plot_summary_plot = FALSE, 
            plot_associations = FALSE))
    }
    results <- as.data.frame(res[['fit_data_abundance']][["results"]])
    statInfo <- results[results[, "metadata"] == contrast[1] &
        results[, "value"] == contrast[2], ]
    ord <- match(rownames(counts), statInfo[, "feature"])
    statInfo <- statInfo[ord, ]
    pValMat <- statInfo[, c("pval_joint", "qval_joint")] 
    colnames(pValMat) <- c("rawP", "adjP")
    rownames(statInfo) <- statInfo[, "feature"] <- rownames(pValMat) <- 
        rownames(counts)
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
}# END - function: DA_maaslin3

#' @title set_maaslin3
#'
#' @export
#' @description
#' Set the parameters for maaslin3 differential abundance detection method.
#'
#' @inheritParams DA_maaslin3
#' @param expand logical, if TRUE create all combinations of input parameters
#' (default \code{expand = TRUE}).
#'
#' @return A named list containing the set of parameters for \code{DA_maaslin3}
#' method.
#'
#' @seealso \code{\link{DA_maaslin3}}
#'
#' @examples
#' # Set some basic combinations of parameters for maaslin3
#' base_maaslin3 <- set_maaslin3(normalization = "TSS", transform = "LOG",
#'     median_comparison_abundance = TRUE, subtract_median = TRUE,
#'     fixed_effects = "group",
#'     contrast = c("group", "B", "A"))
#' many_maaslin3 <- set_maaslin3(normalization = c("TSS", "CLR", "NONE"), 
#'     transform = c("LOG", "NONE"),
#'     median_comparison_abundance = c(TRUE, FALSE),
#'     subtract_median = c(TRUE, FALSE),
#'     fixed_effects = "group",
#'     contrast = c("group", "B", "A"))
set_maaslin3 <- function(assay_name = "counts",
                         normalization = c("TSS", "CLR", "NONE"), 
                         transform = c("LOG", "NONE"), 
                         median_comparison_abundance = c(TRUE, FALSE),
                         subtract_median = c(TRUE, FALSE),
                         correction = "BH", random_effects = NULL, fixed_effects = NULL,
                         contrast = NULL, reference = NULL, expand = TRUE) {
    
    method <- "DA_maaslin3"
    
    # Check required parameters
    if (is.null(assay_name)) {
        stop(method, "\n", "'assay_name' is required (default = 'counts').")
    }
    if (is.null(fixed_effects)) {
        stop(method, "\n", "'fixed_effects' are missing.")
    }
    if (is.null(contrast)) {
        stop(method, "\n", "'contrast' must be specified.")
    }
    if (!is.character(contrast) || length(contrast) != 3) {
        stop(method, "\n", 
             "contrast: please supply a character vector with exactly",
             " three elements: a string indicating the name of a factor,",
             " the level of interest, and the reference level.")
    }
    if(sum(!is.element(normalization, c("TSS", "CLR", "NONE"))) > 0) {
        stop(method, "\n", 
             "normalization: please choose one normalization between",
             " 'TSS', 'CLR', or 'NONE'.")
    }
    if(sum(!is.element(transform, c("LOG", "NONE"))) > 0) {
        stop(method, "\n", 
             "transform: please choose one transformation between",
             " 'LOG' or 'NONE'.")
    }
    
    # Create a grid of parameter combinations.
    if (expand) {
        parameters <- expand.grid(method = method, assay_name = assay_name,
            normalization = normalization, transform = transform,
            median_comparison_abundance = median_comparison_abundance,
            subtract_median = subtract_median,
            correction = correction,
            stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, assay_name = assay_name,
            normalization = normalization, transform = transform,
            median_comparison_abundance = median_comparison_abundance,
            subtract_median = subtract_median,
            correction = correction,
            stringsAsFactors = FALSE)
    }
    
    # Remove senseless combinations:
    wrong_index <- c(which(parameters[, "median_comparison_abundance"] != 
                        parameters[, "subtract_median"]),
                    which(parameters[, "normalization"] == "CLR" & 
                        parameters[, "transform"] != "NONE"))
    if(length(wrong_index) > 0){
        message("Removing incompatible sets.")
        parameters <- parameters[-wrong_index, ]
    }
    
    # Convert data.frame to list of parameter sets
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        # Append additional parameters not included in the expansion grid.
        x <- append(x = x, values = list("random_effects" = random_effects, 
            "fixed_effects" = fixed_effects, "contrast" = contrast, 
            "reference" = reference), after = 6)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}
