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
#' @details
#' Some maaslin3 parameters are not available for customization in this 
#' implementation. For this reason they assume default values or are internally
#' assigned. The latter case is represented by:
#' \itemize{
#' \item \code{warn_prevalence} which is internally set to \code{FALSE};
#' \item \code{subtract_median} which is internally set to the same 
#' \code{median_comparison_abundance} value;
#' \item \code{zero_threshold} which is automatically set to -1 when 
#' \code{transform = "PLOG"};
#' \item \code{evaluate_only} is automatically set to \code{"abundance"} when
#' \code{transform = "PLOG"}.
#' }
#' Please refer to maaslin3's guide to choose proper parameter combinations.
#' 
#' @return A list object containing the matrix of abundance models related 
#' p-values \code{pValMat}, prevalence models related p-values 
#' \code{pValMat_prev}, joint abundance and prevalence models p-values 
#' \code{pValMat_joint}, a matrix of summary statistics for each tag based on 
#' abundance models \code{statInfo} or prevalence models \code{statInfo_prev}, 
#' and a suggested \code{name} of the final object considering the parameters 
#' passed to the function. 
#' In case of \code{transform = "PLOG"} only \code{pValMat}, \code{statInfo}, 
#' and \code{name} are returned.
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
#' DA_maaslin3(object = ps, formula = "~ group", normalization = "CLR", 
#'     transform = "NONE", correction = "BH", contrast = c("group", "B", "A"), 
#'     verbose = FALSE)

DA_maaslin3 <- function(object, assay_name = "counts", 
    formula = NULL, contrast = NULL,
    normalization = c("TSS", "CLR", "NONE"), 
    transform = c("LOG", "PLOG", "NONE"), 
    median_comparison_abundance = TRUE,
    correction = "BH", 
    verbose = TRUE){
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
    if(length(normalization) > 1)
        stop(method, "\n", 
             "normalization: please choose one normalization for this istance",
             " of differential abundance analysis.")
    if(sum(!is.element(normalization, c("TSS", "CLR", "NONE"))) > 0){
        stop(method, "\n", 
             "normalization: please choose one normalization between 'TSS',",
             " 'CLR', or 'NONE'.")
    }
    # Check transform
    if(length(transform) > 1)
        stop(method, "\n", 
             "transform: please choose one transform for this istance",
             " of differential abundance analysis.")
    if(sum(!is.element(transform, c("LOG", "PLOG", "NONE"))) > 0){
        stop(method, "\n", 
             "transform: please choose one transfomation between 'LOG',",
             " 'PLOG', or 'NONE'.")
    }
    # PLOG: automatically sets zero_threshold and abundance models only
    zero_threshold <- 0
    evaluate_only <- NULL
    if(transform == "PLOG"){
        zero_threshold <- -1
        evaluate_only <- "abundance"
    }
    # Check compatibility between normalization and transform
    if(normalization == 'CLR' & transform != 'NONE')
        stop(method, "\n", 
             "if normalization is CLR, transform must be NONE.")
    name <- paste(name, ".", normalization, "norm.", transform, "trans", 
        ifelse(median_comparison_abundance, ".medCompare", ""), sep = "")
    if(!is.character(contrast) | length(contrast) != 3)
        stop(method, "\n", 
             "contrast: please supply a character vector with exactly", 
             " three elements: the name of a variable used in",  
             " 'fixed_effects', the name of the level of interest, and the", 
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
            subtract_median = median_comparison_abundance,
            warn_prevalence = FALSE, zero_threshold = zero_threshold,
            evaluate_only = evaluate_only, max_significance = 0,
            formula = formula, correction = correction, 
            plot_summary_plot = FALSE, plot_associations = FALSE,
            verbosity = "FINEST")
    } else {
        utils::capture.output(file = tempfile(),
        res <- maaslin3(input_data = t(counts), input_metadata = metadata, 
            output = tempdir(), normalization = normalization, 
            transform = transform, standardize = TRUE, 
            median_comparison_abundance = median_comparison_abundance,
            subtract_median = median_comparison_abundance,
            warn_prevalence = FALSE, evaluate_only = evaluate_only,
            zero_threshold = zero_threshold,
            max_significance = 0,
            formula = formula, correction = correction, 
            plot_summary_plot = FALSE, plot_associations = FALSE,
            verbosity = "WARN"))
    }
    # Results for abundance
    results <- as.data.frame(res[['fit_data_abundance']][["results"]])
    statInfo <- results[results[, "metadata"] == contrast[1] &
        results[, "value"] == contrast[2], ]
    ord <- match(rownames(counts), statInfo[, "feature"])
    statInfo <- statInfo[ord, ]
    pValMat <- statInfo[, c("pval_individual", "qval_individual")] 
    colnames(pValMat) <- c("rawP", "adjP")
    # When transform = "PLOG" only abundance models are fit
    if(transform != "PLOG"){
        # Results for prevalence
        results_prev <- as.data.frame(res[['fit_data_prevalence']][["results"]])
        statInfo_prev <- results_prev[
            results_prev[, "metadata"] == contrast[1] &
            results_prev[, "value"] == contrast[2], ]
        ord_prev <- match(rownames(counts), statInfo_prev[, "feature"])
        statInfo_prev <- statInfo_prev[ord_prev, ]
        pValMat_prev <- statInfo_prev[, c("pval_individual", "qval_individual")] 
        colnames(pValMat_prev) <- c("rawP", "adjP")
        # Abundance and Prevalence joint results
        pValMat_joint <- statInfo_prev[, c("pval_joint", "qval_joint")]
        colnames(pValMat_joint) <- c("rawP", "adjP")
        rownames(statInfo) <- rownames(statInfo_prev) <- 
            statInfo[, "feature"] <- statInfo_prev[, "feature"] <- 
            rownames(pValMat) <- rownames(pValMat_prev) <- 
            rownames(pValMat_joint) <- rownames(counts)
        return(list("pValMat" = pValMat, "statInfo" = statInfo, 
            "pValMat_prev" = pValMat_prev, "statInfo_prev" = statInfo_prev,
            "pValMat_joint" = pValMat_joint, "name" = name))
    } else {
        rownames(statInfo) <- statInfo[, "feature"] <- rownames(pValMat) <- 
            rownames(counts)
        return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
    }
    
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
#' @inherit DA_maaslin3
#'
#' @seealso \code{\link{DA_maaslin3}}
#'
#' @examples
#' # Set some basic combinations of parameters for maaslin3
#' base_maaslin3 <- set_maaslin3(normalization = "TSS", transform = "LOG",
#'     median_comparison_abundance = TRUE, formula = ~ group,
#'     contrast = c("group", "B", "A"))
#' many_maaslin3 <- set_maaslin3(normalization = c("TSS", "CLR", "NONE"), 
#'     transform = c("LOG", "NONE"),
#'     median_comparison_abundance = c(TRUE, FALSE),
#'     formula = ~ group, contrast = c("group", "B", "A"))
set_maaslin3 <- function(assay_name = "counts",
    normalization = c("TSS", "CLR", "NONE"), 
    transform = c("LOG", "PLOG", "NONE"), 
    median_comparison_abundance = c(TRUE, FALSE),
    correction = "BH", formula = NULL, contrast = NULL,
    expand = TRUE) {
    
    method <- "DA_maaslin3"
    
    # Check required parameters
    if (is.null(assay_name)) {
        stop(method, "\n", "'assay_name' is required (default = 'counts').")
    }
    if (is.null(formula)) {
        stop(method, "\n", "'formula' is missing.")
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
             "normalization: please choose normalizations between",
             " 'TSS', 'CLR', or 'NONE'.")
    }
    if(sum(!is.element(transform, c("LOG", "PLOG", "NONE"))) > 0) {
        stop(method, "\n", 
             "transform: please choose transformations between",
             " 'LOG', 'PLOG', or 'NONE'.")
    }
    
    # Create a grid of parameter combinations.
    if (expand) {
        parameters <- expand.grid(method = method, assay_name = assay_name,
            normalization = normalization, transform = transform,
            median_comparison_abundance = median_comparison_abundance,
            correction = correction,
            stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, assay_name = assay_name,
            normalization = normalization, transform = transform,
            median_comparison_abundance = median_comparison_abundance,
            correction = correction, stringsAsFactors = FALSE)
    }
    
    # Remove senseless combinations:
    wrong_index <- which(parameters[, "normalization"] == "CLR" & 
                         parameters[, "transform"] != "NONE")
    if(length(wrong_index) > 0){
        message("Removing incompatible sets.")
        parameters <- parameters[-wrong_index, ]
    }
    
    # Convert data.frame to list of parameter sets
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        # Append additional parameters not included in the expansion grid.
        x <- append(x = x, values = list("formula" = formula, 
            "contrast" = contrast), after = 6)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}
