#' @title DA_linda
#'
#' @importFrom MicrobiomeStat linda
#' @importFrom SummarizedExperiment assays
#' @importFrom phyloseq otu_table sample_data phyloseq taxa_are_rows
#' @export
#' @description
#' Fast run for linda differential abundance detection method.
#'
#' @inheritParams DA_edgeR
#' @param contrast character vector with exactly, three elements: a string 
#' indicating the name of factor whose levels are the conditions to be 
#' compared, the name of the level of interest, and the name of the other 
#' level. 
#' @inheritParams MicrobiomeStat::linda
#'
#' @return A list object containing the matrix of p-values `pValMat`,
#' a matrix of summary statistics for each tag `statInfo`, and a suggested 
#' `name` of the final object considering the parameters passed to the 
#' function.
#'
#' @seealso \code{\link[MicrobiomeStat]{linda}}.
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
#' DA_linda(object = ps, formula = "~ group", contrast = c("group", "B", "A"), 
#'     is.winsor = TRUE, zero.handling = "pseudo-count", verbose = FALSE)

DA_linda <- function(object, assay_name = "counts", formula = NULL, 
    contrast = NULL, is.winsor = TRUE, outlier.pct = 0.03,
    zero.handling = c("pseudo-count", "imputation"), pseudo.cnt = 0.5,
    alpha = 0.05, p.adj.method = "BH", verbose = TRUE){
    counts_and_metadata <- get_counts_metadata(object, assay_name = assay_name)
    counts <- counts_and_metadata[[1]]
    metadata <- counts_and_metadata[[2]]
    is_phyloseq <- counts_and_metadata[[3]]
    # Name building
    name <- "linda"
    method <- "DA_linda"
    # Check the assay
    if (!is_phyloseq){
        if(verbose)
            message("Using the ", assay_name, " assay.")
        name <- paste(name, ".", assay_name, sep = "")
    } 
    # Check winsor
    if(is.winsor)
        name <- paste(name, ".", "win", outlier.pct, sep = "")
    # Check zero.handling method
    # Adaptive not included, there is a bug in original code and is always 
    # pseudo-count
    if(length(zero.handling) != 1){
        stop(method, "\n", 
             "zero.handling: please choose only one zero.handling method", 
             " between 'pseudo-count' or 'imputation'.")
    }
    if(!is.element(zero.handling, c("pseudo-count", "imputation")))
        stop(method, "\n", 
             "zero.handling: please choose one test between 'pseudo-count'",
             " or 'imputation'.")
    if(zero.handling == "pseudo-count"){
        name <- paste(name, ".", "pc", pseudo.cnt, sep = "")
    } else {
        name <- paste(name, ".", zero.handling, sep = "")
    }
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
    res <- MicrobiomeStat::linda(feature.dat = counts, 
        meta.dat = metadata, formula = formula, alpha = alpha, 
        is.winsor = is.winsor, outlier.pct = outlier.pct, 
        zero.handling = zero.handling, pseudo.cnt = pseudo.cnt,
        p.adj.method = p.adj.method, prev.filter = 0, 
        mean.abund.filter = 0, verbose = verbose)
    statInfo <- as.data.frame(res[["output"]]
        [[paste(contrast[c(1,2)], collapse = "")]])
    pValMat <- statInfo[, c("pvalue", "padj")] 
    colnames(pValMat) <- c("rawP", "adjP")
    rownames(statInfo) <- rownames(pValMat) <- rownames(counts)
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
}# END - function: DA_linda

#' @title set_linda
#'
#' @export
#' @description
#' Set the parameters for linda differential abundance detection method.
#'
#' @inheritParams DA_linda
#' @param expand logical, if TRUE create all combinations of input parameters
#' (default \code{expand = TRUE}).
#'
#' @return A named list containing the set of parameters for \code{DA_linda}
#' method.
#'
#' @seealso \code{\link{DA_linda}}
#'
#' @examples
#' # Set some basic combinations of parameters for ANCOM with bias correction
#' base_linda <- set_linda(formula = "~ group", contrast = c("group", "B", "A"), 
#'     zero.handling = "pseudo-count", expand = TRUE)
#' many_linda <- set_linda(formula = "~ group", contrast = c("group", "B", "A"), 
#'     is.winsor = c(TRUE, FALSE), 
#'     zero.handling = c("pseudo-count", "imputation"), expand = TRUE)
set_linda <- function(assay_name = "counts", formula = NULL, 
    contrast = NULL, is.winsor = TRUE, outlier.pct = 0.03,
    zero.handling = c("pseudo-count", "imputation"), pseudo.cnt = 0.5,
    alpha = 0.05, p.adj.method = "BH", expand = TRUE) {
    method <- "DA_linda"
    if (is.null(assay_name)) {
        stop(method, "\n", "'assay_name' is required (default = 'counts').")
    }
    if (!is.logical(is.winsor)) {
        stop(method, "\n", 
            "'is.winsor' must be logical.")
    }
    if(!is.null(formula)){
        if (!is.character(formula)){
            stop(method, "\n", "'formula' should be a character formula.")
        }
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
    if(sum(!is.element(zero.handling, c("pseudo-count", "imputation"))) > 0){
        stop(method, "\n", 
             "zero.handling: please choose one zero handling method between",
             " 'pseudo-count' or 'imputation'.")
    }
    if (expand) {
        parameters <- expand.grid(method = method, assay_name = assay_name,
            is.winsor = is.winsor, outlier.pct = outlier.pct, 
            zero.handling = zero.handling, pseudo.cnt = pseudo.cnt, 
            alpha = alpha, p.adj.method = p.adj.method, 
            stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, assay_name = assay_name,
            is.winsor = is.winsor, outlier.pct = outlier.pct, 
            zero.handling = zero.handling, pseudo.cnt = pseudo.cnt, 
            alpha = alpha, p.adj.method = p.adj.method)
    }
    # Remove senseless combinations
    # is.winsor = FALSE with many outlier.pct
    not_winsor <- which(!parameters[, "is.winsor"])
    unique_not_winsor <- duplicated(parameters[not_winsor, -4])
    if(sum(unique_not_winsor) > 0){
        message("Removing duplicated instances without winsorisation.")
        parameters <- parameters[-not_winsor[unique_not_winsor], ]
    }
    # zero.handling = "imputation" with many pseudo.cnt
    imputation <- which(parameters[, "zero.handling"] == "imputation")
    unique_imputation <- duplicated(parameters[imputation, -6])
    if(sum(unique_imputation) > 0){
        message("Removing duplicated instances with zero imputation.")
        parameters <- parameters[-imputation[unique_imputation], ]
    }
    # data.frame to list
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        x <- append(x = x, values = list("formula" = formula, 
            "contrast" = contrast), after = 2)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}
