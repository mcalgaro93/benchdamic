#' @title DA_basic
#'
#' @importFrom SummarizedExperiment assays
#' @export
#' @description
#' Fast run for basic differential abundance detection methods such as wilcox 
#' and t tests.
#'
#' @inheritParams DA_edgeR
#' @param test name of the test to perform. Choose between "t" or "wilcox".
#' @param paired boolean. Choose whether the test is paired or not (default 
#' \code{paired = FALSE}). If \code{paired = TRUE} be sure to provide the 
#' object properly ordered (by the grouping variable).
#' @param contrast character vector with exactly, three elements: a string 
#' indicating the name of factor whose levels are the conditions to be 
#' compared, the name of the level of interest, and the name of the other 
#' level. 
#'
#' @return A list object containing the matrix of p-values `pValMat`,
#' a matrix of summary statistics for each tag `statInfo`, and a suggested 
#' `name` of the final object considering the parameters passed to the 
#' function.
#'
#' @seealso \code{\link{DA_Seurat}} for a similar implementation of basic 
#' tests.
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
#' DA_basic(object = ps, pseudo_count = FALSE, contrast = c("group", "B", "A"),
#'     test = "t", verbose = FALSE)

DA_basic <- function(object, assay_name = "counts", pseudo_count = FALSE, 
    contrast = NULL, test = c("t", "wilcox"), paired = FALSE, verbose = TRUE){
    counts_and_metadata <- get_counts_metadata(object, assay_name = assay_name)
    counts <- counts_and_metadata[[1]]
    metadata <- counts_and_metadata[[2]]
    is_phyloseq <- counts_and_metadata[[3]]
    # Name building
    name <- "basic"
    method <- "DA_basic"
    # Check test
    if (length(test) > 1) {
        stop(method, "\n",
             "test: please choose one test for this istance of", 
             " differential abundance analysis.")
    }
    this_tests <- c("t", "wilcox")
    if(!is.element(test, this_tests)) {
        stop(method, "\n",
             test, " test is not available. Choose between 't' or 'wilcox'.")
    }
    name <- paste(name, ".", test, sep = "")
    # check paired
    if(paired){
        name <- paste(name, ".", "paired", sep = "")
        if(verbose){
            message(method, "\n",
                 "'paired' = TRUE. Be sure to have ordered samples in both",
                 " groups.")
        }
    }
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
    group1 <- which(metadata[, contrast[1]] == contrast[2])
    group2 <- which(metadata[, contrast[1]] == contrast[3])
    if(paired & length(group1) != length(group2)){
        stop(method, "\n", 
             "paired = TRUE but ", contrast[2], " and ", contrast[3], 
             " have ", length(group1), " and ", length(group2), " samples",
             " respectively.")
    }
    statInfo <- plyr::ldply(apply(X = counts, MARGIN = 1, 
        FUN = function(feature){
            if(!paired){
                # Compute the average values
                avg_group1 <- mean(feature[group1]) 
                avg_group2 <- mean(feature[group2])
                # Compute the logFC
                logFC <- log1p(avg_group1) - log1p(avg_group2)
            } else {
                logFC <- mean(log1p(feature[group1]) - log1p(feature[group2]))
            }
        
        if(test == "t"){
            results <- stats::t.test(x = feature[group1], y = feature[group2], 
                paired = paired)
            out <- data.frame(
                "statistic" = results[["statistic"]], 
                "pvalue" = results[["p.value"]],
                "lowCI" = results[["conf.int"]][1],
                "uppCI" = results[["conf.int"]][2])
        } 
        if(test == "wilcox"){
            results <- stats::wilcox.test(x = feature[group1], 
                y = feature[group2], paired = paired, exact = FALSE)
            out <- data.frame(
                "statistic" = results[["statistic"]], 
                "pvalue" = results[["p.value"]])
        }
        out <- data.frame(out, "logFC" = logFC)
        return(out)
    }))
    colnames(statInfo)[1] <- "taxon"
    # Creating pValMat
    padj <- stats::p.adjust(p = statInfo$pvalue, method = "BH")
    pValMat <- data.frame("rawP" = statInfo[, "pvalue"], "adjP" = padj)
    rownames(statInfo) <- rownames(pValMat) <- statInfo[, 1]
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
}# END - function: DA_basic

#' @title set_basic
#'
#' @export
#' @description
#' Set the parameters for basic differential abundance detection methods such 
#' as t and wilcox.
#'
#' @inheritParams DA_basic
#' @param expand logical, if TRUE create all combinations of input parameters
#' (default \code{expand = TRUE}).
#'
#' @return A named list containing the set of parameters for \code{DA_basic}
#' method.
#'
#' @seealso \code{\link{DA_basic}}
#'
#' @examples
#' # Set some basic methods
#' basic_methods <- set_basic(pseudo_count = FALSE, test = c("t", "wilcox"),
#'     contrast = c("group", "B", "A"), expand = TRUE)

set_basic <- function(assay_name = "counts", pseudo_count = FALSE, 
    contrast = NULL, test = c("t", "wilcox"), paired = FALSE, expand = TRUE) {
    method <- "DA_basic"
    if (is.null(assay_name)) {
        stop(method, "\n", "'assay_name' is required (default = 'counts').")
    }
    if (!is.logical(pseudo_count) | !is.logical(paired)) {
        stop(method, "\n", "'pseudo_count' and 'paired' must be logical.")
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
    if (sum(!is.element(test, c("t", "wilcox"))) > 0) {
        stop(method, "\n", 
            "One or more elements into 'test' are not available.",
            " Please choose between 't' or 'wilcox'.")
    }
    if (expand) {
        parameters <- expand.grid(method = method, assay_name = assay_name, 
            pseudo_count = pseudo_count, test = test, paired = paired,
            stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, assay_name = assay_name, 
            pseudo_count = pseudo_count, test = test, paired = paired)
    }
    # data.frame to list
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        x <- append(x = x, values = list("contrast" = contrast), after = 3)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}
