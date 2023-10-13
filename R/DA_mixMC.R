#' @title DA_mixMC
#'
#' @importFrom mixOmics tune.splsda splsda perf plotLoadings
#' @importFrom plyr join
#' @importFrom SummarizedExperiment assays
#' @importFrom phyloseq otu_table sample_data phyloseq taxa_are_rows
#' @importFrom grDevices pdf dev.off
#' @export
#' @description
#' Fast run for mixMC sPLS-DA method for biomarker identification. It performs
#' a CLR transformation on the `counts + pseudo_counts` values. Then the 
#' sPLS-DA is tuned through a leave-one-out cross validation procedure.
#'
#' @inheritParams DA_edgeR
#' @param pseudo_count a positive numeric value for the pseudo-count to be 
#' added. Default is 1.
#' @param contrast character vector with exactly, three elements: a string 
#' indicating the name of factor whose levels are the conditions to be 
#' compared, the name of the level of interest, and the name of the other 
#' level. 
#' @param ID_variable a character string indicating the name of the variable 
#' name corresponding to the repeated measures units (e.g., the subject ID).
#'
#' @return A list object containing the matrix of p-values `pValMat`,
#' a matrix of summary statistics for each tag `statInfo`, and a suggested 
#' `name` of the final object considering the parameters passed to the 
#' function. mixMC does not produce p-values. The frequency and the importance
#' values are produced instead. The frequency indicates the stability of 
#' the features across the folds of the cross validation. The importance 
#' indicates the magnitude of the discrimination for the features and their 
#' direction. Hence, `pValMat` matrix is filled with \code{1 - frequency} values 
#' which are not p-values. To find discriminant features a threshold on this 
#' statistic can be used (liberal < 1, < 0.5, < 0.1 conservative).
#'
#' @seealso \code{\link[mixOmics]{splsda}}, \code{\link[mixOmics]{perf}}, 
#' \code{\link[mixOmics]{tune.splsda}}.
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
#' DA_mixMC(object = ps, pseudo_count = 1, contrast = c("group", "B", "A"),
#'     verbose = FALSE)

DA_mixMC <- function(object, pseudo_count = 1, assay_name = "counts", 
    contrast = NULL, ID_variable = NULL, verbose = TRUE){
    counts_and_metadata <- get_counts_metadata(object, assay_name = assay_name)
    counts <- counts_and_metadata[[1]]
    metadata <- counts_and_metadata[[2]]
    is_phyloseq <- counts_and_metadata[[3]]
    # Name building
    name <- "mixMC"
    method <- "DA_mixMC"
    # Check the assay
    if (!is_phyloseq){
        if(verbose){
            message("Using the ", assay_name, " assay.")
        }  
        name <- paste(name, ".", assay_name, sep = "")
    } 
    # add a pseudo-count
    if(verbose){
        message("Adding a pseudo count... \n")
    }
    X <- t(counts) + pseudo_count
    name <- paste(name, ".pc", pseudo_count, sep = "")
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
    Y <- metadata[, contrast[1]]
    grid.keepX = c(seq(2,9,1), seq(10, ncol(X), 10))
    if(is.null(ID_variable)){
        if(verbose)
            message("Tuning: sPLS-DA, 1 component, leave-one-out cross", 
                " validation.")
        splsda_tune = mixOmics::tune.splsda(X = X, Y = Y, logratio = 'CLR', 
            ncomp = 1, test.keepX = grid.keepX, validation = "loo", 
            progressBar = FALSE)
        optimal.keepX = splsda_tune[["choice.keepX"]]
        if(verbose)
            message("Computing: sPLS-DA, 1 component.")
        cor_splsda = mixOmics::splsda(X = X, Y = Y, logratio= "CLR", ncomp = 1,
            keepX = optimal.keepX)
    } else{
        if(verbose){
            message("Using ", ID_variable, " to identify repeated measures.",
                " Make sure you have complete cases and your data are ordered",
                " by ", contrast[1], " and ", ID_variable, ".")
        }
        multilevel <- as.factor(metadata[, ID_variable])
        if(verbose)
            message("Tuning: sPLS-DA, 1 component, leave-one-out cross", 
                    " validation, repeated measures.")
        splsda_tune = mixOmics::tune.splsda(X = X, Y = Y, logratio = 'CLR', 
            ncomp = 1, test.keepX = grid.keepX, validation = "loo", 
            multilevel = multilevel, progressBar = FALSE)
        optimal.keepX = splsda_tune[["choice.keepX"]]
        if(verbose)
            message("Computing: sPLS-DA, 1 component, repeated measures.")
        cor_splsda = mixOmics::splsda(X = X, Y = Y, logratio= "CLR", ncomp = 1,
            keepX = optimal.keepX, multilevel = multilevel)
    }
    perf.splsda = mixOmics::perf(cor_splsda, validation = "loo", 
        progressBar = FALSE)
    # Extract feature stability
    if(verbose)
        message("Extracting feature stability.")
    stability <- data.frame(perf.splsda[["features"]][["stable"]])
    colnames(stability) = c("feature", "frequency")
    # Extract direction
    grDevices::pdf(NULL)
    loadings_plot <- mixOmics::plotLoadings(cor_splsda, contrib = "max", 
        plot = TRUE)
    invisible(grDevices::dev.off())
    importance <- loadings_plot[["X"]]["importance"]
    importance[, "feature"] <- rownames(importance)
    # Create the statInfo
    statInfo <- data.frame("feature" = rownames(counts))
    statInfo <- plyr::join(x = statInfo, y = stability, by = "feature", 
        type = "left")
    statInfo <- plyr::join(x = statInfo, y = importance, by = "feature", 
        type = "left")
    # Fill with zero all NA values
    statInfo[is.na(statInfo)] <- 0
    pValMat <- data.frame("rawP" = 1 - statInfo[, "frequency"], 
        "adjP" = 1 - statInfo[, "frequency"])
    rownames(statInfo) <- rownames(pValMat) <- statInfo[, "feature"]
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
}# END - function: DA_mixMC

#' @title set_mixMC
#'
#' @export
#' @description
#' Set the parameters for mixMC sPLS-DA.
#'
#' @inheritParams DA_mixMC
#' @param expand logical, if TRUE create all combinations of input parameters
#' (default \code{expand = TRUE}).
#'
#' @return A named list contaning the set of parameters for \code{DA_mixMC}
#' method.
#'
#' @seealso \code{\link{DA_mixMC}}
#'
#' @examples
#' # Set some basic combinations of parameters for mixMC
#' base_mixMC <- set_mixMC(pseudo_count = 1, contrast = c("group", "B", "A"))
#' many_mixMC <- set_mixMC(pseudo_count = c(0.1, 0.5, 1), 
#'     contrast = c("group", "B", "A"))
set_mixMC <- function(assay_name = "counts", pseudo_count = 1,
    contrast = NULL, ID_variable = NULL, expand = TRUE) {
    method <- "DA_mixMC"
    if (is.null(assay_name)) {
        stop(method, "\n", "'assay_name' is required (default = 'counts').")
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
    if (expand) {
        parameters <- expand.grid(method = method, assay_name = assay_name,
            pseudo_count = pseudo_count, stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, assay_name = assay_name,
            pseudo_count = pseudo_count)
    }
    # data.frame to list
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        x <- append(x = x, values = list("contrast" = contrast, 
            "ID_variable" = ID_variable), after = 3)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}
