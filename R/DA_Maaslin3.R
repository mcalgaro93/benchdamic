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
#' @param stat_type Whether to return statistics based on abundance 
#' ("abundance") or prevalence ("prevalence") models.
#' @param pvalue_type Whether to return p-values based on abundance 
#' ("abundance") models, prevalence ("prevalence") models, or joint 
#' ("joint") p-values. Choose "abundance" or "joint" when \code{stat_type} is
#' set to "abundance", choose "prevalence" when \code{stat_type} is set to
#' "prevalence".
#' 
#' @details
#' Some maaslin3 parameters are not available for customization in this 
#' implementation. For this reason they assume default values or are internally
#' assigned. The latter case is represented by:
#' \itemize{
#' \item \code{warn_prevalence} which is internally set to \code{TRUE};
#' \item \code{subtract_median} which is internally set to the same 
#' \code{median_comparison_abundance} value;
#' \item \code{zero_threshold} which is automatically set to -1 when 
#' \code{transform = "PLOG"};
#' \item \code{evaluate_only} is automatically set to \code{"abundance"} when
#' \code{transform = "PLOG"}.
#' }
#' 
#' MaAsLin 3 produces both abundance and prevalence associations with 
#' individual p and adjusted p-values (specific to abundance or prevalence) as 
#' well as joint p and adjusted p-values for testing whether a metadatum is 
#' associated with either the abundance or prevalence. To avoid issues with 
#' having twice as many associations as other tools (from both abundance and 
#' prevalence), \code{stat_type} can be set to report the desired abundance or 
#' prevalence associations. When the abundance and prevalence associations are 
#' expected to go in the same direction, \code{pvalue_type = "joint"} allows to
#' return p-values and adjusted p-values taken from the joint p-values and 
#' adjusted p-values. 
#' Please refer to maaslin3's guide to choose proper parameter combinations.
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
#' DA_maaslin3(object = ps, formula = "~ group", normalization = "CLR", 
#'     transform = "NONE", correction = "BH", contrast = c("group", "B", "A"), 
#'     verbose = FALSE, stat_type = "abundance", pvalue_type = "joint")

DA_maaslin3 <- function(object, assay_name = "counts", 
    formula = NULL, contrast = NULL,
    normalization = c("TSS", "CLR", "NONE"), 
    transform = c("LOG", "PLOG", "NONE"), 
    median_comparison_abundance = TRUE,
    stat_type = c("abundance", "prevalence"),
    pvalue_type = c("abundance", "prevalence", "joint"),
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
    # Check stat_type
    if(length(stat_type) > 1)
        stop(method, "\n", 
             "stat_type: please choose one stat_type for this istance",
             " of differential abundance analysis.")
    if(sum(!is.element(stat_type, c("abundance", "prevalence"))) > 0){
        stop(method, "\n", 
             "stat_type: please choose one stat_type between 'abundance',",
             " or 'prevalence'.")
    }
    # Check pvalue_type
    if(length(pvalue_type) > 1)
        stop(method, "\n", 
             "pvalue_type: please choose one pvalue_type for this istance",
             " of differential abundance analysis.")
    if(sum(!is.element(pvalue_type, 
        c("abundance", "prevalence", "joint"))) > 0){
        stop(method, "\n", 
             "pvalue_type: please choose one pvalue_type between 'abundance',",
             " 'prevalence', or 'joint'.")
    }
    if(stat_type == "abundance" & pvalue_type == "prevalence"){
        stop(method, "\n", 
             "pvalue_type: please choose one pvalue_type between 'abundance',",
             " or 'joint'.")
    }
    if(stat_type == "prevalence" & pvalue_type != "prevalence"){
        stop(method, "\n", 
             "pvalue_type: please choose 'prevalence' pvalue_type.")
    }
    # PLOG: automatically sets zero_threshold and abundance models only
    zero_threshold <- 0
    evaluate_only <- NULL
    if(transform == "PLOG"){
        zero_threshold <- -1
        evaluate_only <- "abundance"
        # Check compatibility between PLOG and stat_type, pvalue_type
        if(is.element("prevalence", c(stat_type, pvalue_type))){
            stop(method, "\n", 
                 "if transform is 'PLOG', stat_type and pvalue_type must be",
                 " 'abundance'.")
        }
    }
    # Check compatibility between normalization and transform
    if(normalization == 'CLR' & transform != 'NONE')
        stop(method, "\n", 
             "if normalization is CLR, transform must be NONE.")
    name <- paste(name, ".", normalization, "norm.", transform, "trans", 
        ifelse(median_comparison_abundance, ".med", ""), sep = "")
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
            warn_prevalence = TRUE, zero_threshold = zero_threshold,
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
            warn_prevalence = TRUE, evaluate_only = evaluate_only,
            zero_threshold = zero_threshold,
            max_significance = 0,
            formula = formula, correction = correction, 
            plot_summary_plot = FALSE, plot_associations = FALSE,
            verbosity = "WARN"))
    }
    # Results for abundance
    results_ab <- as.data.frame(res[['fit_data_abundance']][["results"]])
    statInfo_ab <- results_ab[results_ab[, "metadata"] == contrast[1] &
        results_ab[, "value"] == contrast[2], ]
    ord <- match(rownames(counts), statInfo_ab[, "feature"])
    statInfo_ab <- statInfo_ab[ord, ]
    pValMat_ab <- statInfo_ab[, c("pval_individual", "qval_individual")] 
    colnames(pValMat_ab) <- c("rawP", "adjP")
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
        rownames(statInfo_ab) <- rownames(statInfo_prev) <- 
            statInfo_ab[, "feature"] <- statInfo_prev[, "feature"] <- 
            rownames(pValMat_ab) <- rownames(pValMat_prev) <- 
            rownames(pValMat_joint) <- rownames(counts)
        # Build name
        # pvalue_type name
        if(pvalue_type == "abundance"){
            pValMat <- pValMat_ab
            name_p <- "ab"
        } else if(pvalue_type == "prevalence"){
            pValMat <- pValMat_prev
            name_p <- "prev"
        } else {
            pValMat <- pValMat_joint
            name_p <- "joint"
        }
        # stat_type name
        if(stat_type == "abundance"){
            statInfo <- statInfo_ab
            name_stat <- "ab"
        } else {
            statInfo <- statInfo_prev
            name_stat <- "prev"
        }
        # Check if stat_type and pvalue_type are the same
        if(name_p != name_stat){
            name <- paste(name, ".", name_stat, "S.", name_p, "P", sep = "")
        } else {
            name <- paste(name, ".", name_stat, "SP", sep = "")
        }
        return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
    } else {
        rownames(statInfo_ab) <- statInfo_ab[, "feature"] <- 
            rownames(pValMat_ab) <- rownames(counts)
        name <- paste(name, ".", "abSP", sep = "")
        return(list("pValMat" = pValMat_ab, "statInfo" = statInfo_ab, 
            "name" = name))
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
#'     median_comparison_abundance = TRUE, stat_type = "abundance",
#'     pvalue_type = "abundance", formula = ~ group,
#'     contrast = c("group", "B", "A"))
#' many_maaslin3 <- set_maaslin3(normalization = c("TSS", "CLR"), 
#'     transform = c("LOG", "NONE"),
#'     median_comparison_abundance = c(TRUE, FALSE),
#'     stat_type = "abundance", pvalue_type = c("abundance", "joint"),
#'     formula = ~ group, contrast = c("group", "B", "A"))
set_maaslin3 <- function(assay_name = "counts",
    normalization = c("TSS", "CLR", "NONE"), 
    transform = c("LOG", "PLOG", "NONE"), 
    median_comparison_abundance = c(TRUE, FALSE),
    stat_type = c("abundance", "prevalence"),
    pvalue_type = c("abundance", "prevalence", "joint"),
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
    if(sum(!is.element(stat_type, c("abundance", "prevalence"))) > 0) {
        stop(method, "\n", 
             "stat_type: please choose stat_type between",
             " 'abundance' or 'prevalence'.")
    }
    if(sum(!is.element(pvalue_type, c("abundance", "prevalence", "joint"))) > 0) {
        stop(method, "\n", 
             "pvalue_type: please choose pvalue_type between",
             " 'abundance', 'prevalence', or 'joint'.")
    }
    # Create a grid of parameter combinations.
    if (expand) {
        parameters <- expand.grid(method = method, assay_name = assay_name,
            normalization = normalization, transform = transform,
            median_comparison_abundance = median_comparison_abundance,
            stat_type = stat_type, pvalue_type = pvalue_type,
            correction = correction, stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, assay_name = assay_name,
            normalization = normalization, transform = transform,
            median_comparison_abundance = median_comparison_abundance,
            stat_type = stat_type, pvalue_type = pvalue_type,
            correction = correction, stringsAsFactors = FALSE)
    }
    
    # Remove senseless combinations:
    wrong_index <- c(which(parameters[, "normalization"] == "CLR" & 
                         parameters[, "transform"] != "NONE"),
                     which(parameters[, "stat_type"] == "abundance" &
                           parameters[, "pvalue_type"] == "prevalence"),
                     which(parameters[, "stat_type"] == "prevalence" &
                           parameters[, "pvalue_type"] != "prevalence"))
    if(length(wrong_index) > 0){
        message("Removing incompatible sets.")
        parameters <- parameters[-wrong_index, ]
    }
    
    # Convert data.frame to list of parameter sets
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        # Append additional parameters not included in the expansion grid.
        x <- append(x = x, values = list("formula" = formula, 
            "contrast" = contrast), after = 8)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}
