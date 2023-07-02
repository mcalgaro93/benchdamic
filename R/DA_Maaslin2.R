#' @title DA_Maaslin2
#'
#' @importFrom Maaslin2 Maaslin2
#' @importFrom SummarizedExperiment assays
#' @importFrom phyloseq otu_table sample_data phyloseq taxa_are_rows
#' @export
#' @description
#' Fast run for Maaslin2 differential abundance detection method.
#'
#' @inheritParams DA_edgeR
#' @param contrast character vector with exactly, three elements: a string 
#' indicating the name of factor whose levels are the conditions to be 
#' compared, the name of the level of interest, and the name of the other 
#' level. 
#' @inheritParams Maaslin2::Maaslin2
#'
#' @return A list object containing the matrix of p-values `pValMat`,
#' a matrix of summary statistics for each tag `statInfo`, and a suggested 
#' `name` of the final object considering the parameters passed to the 
#' function.
#'
#' @seealso \code{\link[Maaslin2]{Maaslin2}}.
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

DA_Maaslin2 <- function(object, assay_name = "counts", 
    normalization = c("TSS", "CLR", "CSS", "NONE", "TMM"), 
    transform = c("LOG", "LOGIT", "AST", "NONE"), 
    analysis_method = c("LM", "CPLM", "ZICP", "NEGBIN", "ZINB"), 
    correction = "BH", random_effects = NULL, fixed_effects = NULL,
    contrast = NULL, reference = NULL, verbose = TRUE){
    counts_and_metadata <- get_counts_metadata(object, assay_name = assay_name)
    counts <- counts_and_metadata[[1]]
    metadata <- counts_and_metadata[[2]]
    is_phyloseq <- counts_and_metadata[[3]]
    # Name building
    name <- "Maaslin2"
    method <- "DA_Maaslin2"
    # Check the assay
    if (!is_phyloseq){
        if(verbose)
            message("Using the ", assay_name, " assay.")
        name <- paste(name, ".", assay_name, sep = "")
    } 
    # Check normalization
    if(sum(!is.element(normalization, 
        c("TSS", "CLR", "CSS", "NONE", "TMM"))) > 0){
        stop(method, "\n", 
             "normalization: please choose one normalization between 'TSS',",
             " 'CLR', 'CSS', 'NONE', or 'TMM'.")
    }
    # Check transform
    if(sum(!is.element(transform, c("LOG", "LOGIT", "AST", "NONE"))) > 0){
        stop(method, "\n", 
             "transform: please choose one transfomation between 'LOG',",
             " 'LOGIT', 'AST', or 'NONE'.")
    }
    # Check analysis method
    if(sum(!is.element(analysis_method, 
        c("LM", "CPLM", "ZICP", "NEGBIN", "ZINB"))) > 0){
        stop(method, "\n", 
             "analysis_method: please choose one analysis method between 'LM',",
             " 'CPLM', 'ZICP', 'NEGBIN', or 'ZINB'.")
    }
    # Remove senseless combinations
    # analysis_method == "LM", no CSS or TMM with non-counts data
    
    if(analysis_method == "LM" & is.element(normalization,  c("CSS", "TMM")))
        stop(method, "\n", 
             "LM should be used with non-counts values. Your normalization",
             " method should be one between 'CLR', 'TSS', or 'NONE'.")
    if(analysis_method == "LM" & normalization == "CLR" & transform != "NONE")
        stop(method, "\n", 
             "LM can be used with 'CLR' values. Please do not add further",
             " transformations, use transform = 'NONE'.")
    if(is.element(analysis_method, c("CPLM", "ZICP")) & 
        (normalization == "CLR" | transform != "NONE"))
        stop(method, "\n", 
            "CPLM and ZICP can be used on positive values. 'CLR' normalization",
            " and some transformations could generate negative values. Use",
            " transform = 'NONE'.")
    if(is.element(analysis_method, c("NEGBIN", "ZINB")) & 
       (is.element(normalization, c("CLR", "TSS")) | transform != "NONE"))
        stop(method, "\n", 
             "NEGBIN and ZINB models should be used on count data. Do not use",
             " 'CLR' and 'TSS' normalizations or any transformation.")
    name <- paste(name, ".", normalization, "norm.", transform, "trans.", 
        analysis_method, sep = "")
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
        res <- Maaslin2(input_data = t(counts), input_metadata = metadata, 
            output = tempdir(), min_abundance = 0, min_prevalence = 0, 
            min_variance = 0, normalization = normalization, 
            transform = transform, standardize = TRUE, 
            analysis_method = analysis_method, max_significance = 0, 
            random_effects = random_effects, fixed_effects = fixed_effects, 
            correction = correction, plot_heatmap = FALSE, plot_scatter = FALSE)
    } else {
        utils::capture.output(file = tempfile(),
            res <- Maaslin2(input_data = t(counts), input_metadata = metadata, 
                output = tempdir(), min_abundance = 0, min_prevalence = 0, 
                min_variance = 0, normalization = normalization, 
                transform = transform, standardize = TRUE, 
                analysis_method = analysis_method, max_significance = 0, 
                random_effects = random_effects, fixed_effects = fixed_effects, 
                correction = correction, plot_heatmap = FALSE, 
                plot_scatter = FALSE))
    }
    results <- as.data.frame(res[["results"]])
    statInfo <- results[results[, "metadata"] == contrast[1] &
        results[, "value"] == contrast[2], ]
    ord <- match(rownames(counts), statInfo[, "feature"])
    statInfo <- statInfo[ord, ]
    pValMat <- statInfo[, c("pval", "qval")] 
    colnames(pValMat) <- c("rawP", "adjP")
    rownames(statInfo) <- rownames(pValMat) <- rownames(counts)
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
}# END - function: DA_Maaslin2

#' @title set_Maaslin2
#'
#' @export
#' @description
#' Set the parameters for Maaslin2 differential abundance detection method.
#'
#' @inheritParams DA_Maaslin2
#' @param expand logical, if TRUE create all combinations of input parameters
#' (default \code{expand = TRUE}).
#'
#' @return A named list containing the set of parameters for \code{DA_Maaslin2}
#' method.
#'
#' @seealso \code{\link{DA_Maaslin2}}
#'
#' @examples
#' # Set some basic combinations of parameters for Maaslin2
#' base_Maaslin2 <- set_Maaslin2(normalization = "TSS", transform = "LOG",
#'     analysis_method = "LM", fixed_effects = "group",
#'     contrast = c("group", "B", "A"))
#' many_Maaslin2 <- set_Maaslin2(normalization = c("TSS", "CLR", "CSS", "TMM", 
#'     "NONE"), transform = c("LOG", "NONE"),
#'     analysis_method = c("LM", "NEGBIN"), fixed_effects = "group",
#'     contrast = c("group", "B", "A"))
set_Maaslin2 <- function(assay_name = "counts",
    normalization = c("TSS", "CLR", "CSS", "NONE", "TMM"), 
    transform = c("LOG", "LOGIT", "AST", "NONE"), 
    analysis_method = c("LM", "CPLM", "ZICP", "NEGBIN", "ZINB"), 
    correction = "BH", random_effects = NULL, fixed_effects = NULL,
    contrast = NULL, reference = NULL, expand = TRUE) {
    method <- "DA_Maaslin2"
    if (is.null(assay_name)) {
        stop(method, "\n", "'assay_name' is required (default = 'counts').")
    }
    if(is.null(fixed_effects)){
        stop(method, "\n", "'fixed_effects' are missing.")
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
    if(sum(!is.element(normalization, 
        c("TSS", "CLR", "CSS", "NONE", "TMM"))) > 0){
        stop(method, "\n", 
             "normalization: please choose one normalization between",
             " 'TSS', 'CLR', 'CSS', 'NONE', or 'TMM'.")
    }
    if(sum(!is.element(transform, 
                       c("LOG", "LOGIT", "AST", "NONE"))) > 0){
        stop(method, "\n", 
             "transform: please choose one transform method between",
             " 'LOG', 'LOGIT', 'AST', or 'NONE'.")
    }
    if(sum(!is.element(analysis_method, 
                       c("LM", "CPLM", "ZICP", "NEGBIN", "ZINB"))) > 0){
        stop(method, "\n", 
             "analysis_method: please choose one analysis method between",
             " 'LM', 'CPLM', 'ZICP', 'NEGBIN', or 'ZINB'.")
    }
    if (expand) {
        parameters <- expand.grid(method = method, assay_name = assay_name,
            normalization = normalization, transform = transform, 
            analysis_method = analysis_method, correction = correction,
            stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, assay_name = assay_name,
             normalization = normalization, transform = transform, 
             analysis_method = analysis_method, correction = correction)
    }
    # Remove senseless combinations
    # analysis_method == "LM", no CSS or TMM with non-counts data
    wrong_LM_index <- which(parameters[, "analysis_method"] == "LM" & 
        is.element(parameters[, "normalization"],  c("CSS", "TMM")))
    wrong_LM_CLR <- which(parameters[, "analysis_method"] == "LM" & 
        parameters[, "normalization"] == "CLR" &
        parameters[, "transform"] != "NONE")
    # analysis_method == "CPLM" or "ZICP", only for positive values
    wrong_CPLM_index <- which(
        is.element(parameters[, "analysis_method"], c("CPLM", "ZICP"))  & 
        (parameters[, "normalization"] == "CLR" | 
        parameters[, "transform"] != "NONE"))
    # analysis_method == "NEGBIN" or "ZINB", only count data
    wrong_NB_index <- which(
        is.element(parameters[, "analysis_method"], c("NEGBIN", "ZINB")) & 
        (is.element(parameters[, "normalization"],  c("CLR", "TSS")) |
        parameters[, "transform"] != "NONE"))
    wrong_index <- c(wrong_LM_index, wrong_LM_CLR, wrong_CPLM_index, 
        wrong_NB_index)
    if(length(wrong_index) > 0){
        message("Removing incompatible sets of normalization, transformation,",
            " and analysis method.")
        parameters <- parameters[-wrong_index, ]
    }
    # data.frame to list
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        x <- append(x = x, values = list("random_effects" = random_effects, 
            "fixed_effects" = fixed_effects, "contrast" = contrast, 
            "reference" = reference), after = 6)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}
