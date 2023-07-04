#' @title DA_ZicoSeq
#'
#' @importFrom GUniFrac ZicoSeq
#' @importFrom SummarizedExperiment assays
#' @importFrom phyloseq otu_table sample_data phyloseq taxa_are_rows
#' @export
#' @description
#' Fast run for ZicoSeq differential abundance detection method.
#'
#' @inheritParams DA_edgeR
#' @param contrast character vector with exactly, three elements: a string 
#' indicating the name of factor whose levels are the conditions to be 
#' compared, the name of the level of interest, and the name of the other 
#' level. 
#' @inheritParams GUniFrac::ZicoSeq
#'
#' @return A list object containing the matrix of p-values `pValMat`,
#' a matrix of summary statistics for each tag `statInfo`, and a suggested 
#' `name` of the final object considering the parameters passed to the 
#' function.
#'
#' @seealso \code{\link[GUniFrac]{ZicoSeq}}.
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
#' DA_ZicoSeq(object = ps, normalization = "CLR", transform = "NONE",
#'     analysis_method = "LM", correction = "BH", random_effects = NULL,
#'     fixed_effects = "group", contrast = c("group", "B", "A"),
#'     verbose = FALSE)

DA_ZicoSeq <- function(object, assay_name = "counts", 
    contrast = NULL, strata = NULL, adj.name = NULL,
    feature.dat.type = c("count", "proportion", "other"), 
    is.winsor = TRUE, outlier.pct = 0.03, 
    winsor.end = c("top", "bottom", "both"),
    is.post.sample = TRUE, post.sample.no = 25, 
    perm.no = 99, link.func = list(function(x) sign(x) * (abs(x))^0.5),
    ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
    verbose = TRUE){
    counts_and_metadata <- get_counts_metadata(object, assay_name = assay_name)
    counts <- counts_and_metadata[[1]]
    metadata <- counts_and_metadata[[2]]
    is_phyloseq <- counts_and_metadata[[3]]
    # Name building
    name <- "ZicoSeq"
    method <- "DA_ZicoSeq"
    # Check the assay
    if (!is_phyloseq){
        if(verbose)
            message("Using the ", assay_name, " assay.")
        name <- paste(name, ".", assay_name, sep = "")
    } 
    # Check feature.dat.type
    if(length(feature.dat.type) > 1)
        stop(method, "\n", 
             "feature.dat.type: please choose one data type for this",
             " istance of differential abundance analysis.")
    if(sum(!is.element(feature.dat.type, 
        c("count", "proportion", "other"))) > 0){
        stop(method, "\n", 
             "feature.dat.type: please choose one data type between 'count',",
             " 'proportion', or 'other'.")
    }
    # Check winsor.end
    if(length(winsor.end) > 1)
        stop(method, "\n", 
             "winsor.end: please choose one winsor end for this",
             " istance of differential abundance analysis.")
    if(sum(!is.element(winsor.end, c("top", "bottom", "both"))) > 0){
        stop(winsor.end, "\n", 
             "winsor.end: please choose one winsor end between 'top',",
             " 'bottom', or 'both'.")
    }
    # Warn for unsuggested combinations
    # feature.dat.type == "count", winsor.end != "top"
    if(feature.dat.type == "count" & winsor.end != "top")
        warning(method, "\n", 
             "With count data, winsorisation is suggested at the top end")
    # feature.dat.type != "count", winsor.end != "both"
    if(feature.dat.type != "count" & winsor.end != "both")
        warning(method, "\n", 
                "With non-count data, winsorisation may be useful at both",
                " ends.")
    if(is.winsor){
        name <- paste(name, ".winsor", outlier.pct, winsor.end, sep = "")
    }
    # Check posterior sampling
    if(is.post.sample == FALSE){
        
    } else if (is.post.sample & feature.dat.type != "count"){
        warning(method, "\n", 
                "With non-count data, posterior sampling is not performed.")
    } else {
        if(ncol(counts) >= 40)
            name <- paste(name, ".post", post.sample.no, sep = "")
    }
    # Update name with reference taxa
    name <- paste(name, ".ref", ref.pct, ".excl", excl.pct, sep = "")
    # Check contrast
    if(!is.character(contrast) | length(contrast) != 3)
        stop(method, "\n", 
             "contrast: please supply a character vector with exactly", 
             " three elements: a string indicating the name of factor whose",
             " levels are the conditions to be compared,",  
             " the name of the level of interest, and the", 
             " name of the other level.")
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
    # Check strata
    if(!is.null(strata)){
        if(!is.element(strata, colnames(metadata))){
            stop(method, "\n", 
                 "strata: ", strata, " is not a variable present in metadata.")
        } else strata <- metadata[, strata]
    }
    if(verbose){
        res <- ZicoSeq(meta.dat = metadata, feature.dat = counts, 
            grp.name = contrast[1], adj.name = adj.name, 
            feature.dat.type = feature.dat.type, 
            prev.filter = 0, mean.abund.filter = 0, max.abund.filter = 0, 
            min.prop = 0, is.winsor = is.winsor, outlier.pct = outlier.pct, 
            winsor.end = winsor.end, is.post.sample = is.post.sample, 
            post.sample.no = post.sample.no, link.func = link.func, 
            stats.combine.func = max, perm.no = perm.no, strata = strata, 
            ref.pct = ref.pct, stage.no = stage.no, excl.pct = excl.pct, 
            is.fwer = TRUE, verbose = verbose, return.feature.dat = FALSE)
    } else {
        utils::capture.output(file = tempfile(),
            suppressMessages(
            res <- GUniFrac::ZicoSeq(meta.dat = metadata, feature.dat = counts, 
               grp.name = contrast[1], adj.name = adj.name, 
               feature.dat.type = feature.dat.type, 
               prev.filter = 0, mean.abund.filter = 0, max.abund.filter = 0, 
               min.prop = 0, is.winsor = is.winsor, outlier.pct = outlier.pct, 
               winsor.end = winsor.end, is.post.sample = is.post.sample, 
               post.sample.no = post.sample.no, link.func = link.func, 
               stats.combine.func = max, perm.no = perm.no, strata = strata,
               ref.pct = ref.pct, stage.no = stage.no, excl.pct = excl.pct, 
               is.fwer = TRUE, verbose = verbose, return.feature.dat = FALSE)))
    }
    statInfo <- as.data.frame(t(res[["coef.list"]][[1]]))
    statInfo[, "effect"] <- statInfo[, paste0(contrast[1], contrast[2])]
    statInfo[, c("R2", "F0", "RSS", "p.raw", "p.adj.fdr", "p.adj.fwer")] <- 
        cbind(res[["R2"]], res[["F0"]], res[["RSS"]], res[["p.raw"]], 
            res[["p.adj.fdr"]], res[["p.adj.fwer"]])
    pValMat <- statInfo[, c("p.raw", "p.adj.fdr")] 
    colnames(pValMat) <- c("rawP", "adjP")
    rownames(statInfo) <- rownames(pValMat) <- rownames(counts)
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
}# END - function: DA_ZicoSeq

#' @title set_ZicoSeq
#'
#' @export
#' @description
#' Set the parameters for ZicoSeq differential abundance detection method.
#'
#' @inheritParams DA_ZicoSeq
#' @param expand logical, if TRUE create all combinations of input parameters
#' (default \code{expand = TRUE}).
#'
#' @return A named list containing the set of parameters for \code{DA_ZicoSeq}
#' method.
#'
#' @seealso \code{\link{DA_ZicoSeq}}
#'
#' @examples
#' # Set some basic combinations of parameters for ZicoSeq
#' base_ZicoSeq <- set_ZicoSeq(contrast = c("group", "B", "A"), 
#'     feature.dat.type = "count", winsor.end = "top")
#' many_ZicoSeq <- set_ZicoSeq(contrast = c("group", "B", "A"), 
#'     feature.dat.type = "count", outlier.pct = c(0.03, 0.05),
#'     winsor.end = "top", is.post.sample = c(TRUE, FALSE))
set_ZicoSeq <- function(assay_name = "counts", contrast = NULL, strata = NULL, 
    adj.name = NULL, feature.dat.type = c("count", "proportion", "other"), 
    is.winsor = TRUE, outlier.pct = 0.03, 
    winsor.end = c("top", "bottom", "both"), is.post.sample = TRUE, 
    post.sample.no = 25, perm.no = 99,
    link.func = list(function(x) sign(x) * (abs(x))^0.5), ref.pct = 0.5,
    stage.no = 6, excl.pct = 0.2, expand = TRUE) {
    method <- "DA_ZicoSeq"
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
    # Check feature.dat.type
    if(length(feature.dat.type) > 1)
        stop(method, "\n", 
             "feature.dat.type: you can set only one data type.")
    if(sum(!is.element(feature.dat.type, 
                       c("count", "proportion", "other"))) > 0){
        stop(method, "\n", 
             "feature.dat.type: please choose one data type between 'count',",
             " 'proportion', or 'other'.")
    }
    # Check winsor.end
    if(sum(!is.element(winsor.end, c("top", "bottom", "both"))) > 0){
        stop(winsor.end, "\n", 
             "winsor.end: please choose one winsor end between 'top',",
             " 'bottom', or 'both'.")
    }
    if (expand) {
        parameters <- expand.grid(method = method, assay_name = assay_name,
            feature.dat.type = feature.dat.type, is.winsor = is.winsor, 
            outlier.pct = outlier.pct, winsor.end = winsor.end, 
            is.post.sample = is.post.sample, post.sample.no = post.sample.no, 
            perm.no = perm.no, ref.pct = ref.pct, stage.no = stage.no, 
            excl.pct = excl.pct, stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, assay_name = assay_name,
             feature.dat.type = feature.dat.type, is.winsor = is.winsor, 
             outlier.pct = outlier.pct, winsor.end = winsor.end, 
             is.post.sample = is.post.sample, post.sample.no = post.sample.no, 
             perm.no = perm.no, ref.pct = ref.pct, stage.no = stage.no, 
             excl.pct = excl.pct)
    }
    # Remove senseless combinations
    # is.winsor = FALSE with many outlier.pct or different winsor.end
    not_winsor <- which(!parameters[, "is.winsor"])
    unique_not_winsor <- duplicated(parameters[not_winsor, -c(5, 6)])
    if(sum(unique_not_winsor) > 0){
        message("Removing duplicated instances without winsorisation.")
        parameters <- parameters[-not_winsor[unique_not_winsor], ]
    }
    # is.post.sample = FALSE with many post.sample.no
    not_post.sample <- which(!parameters[, "is.post.sample"])
    unique_not_post.sample <- duplicated(parameters[not_post.sample, -8])
    if(sum(unique_not_post.sample) > 0){
        message("Removing duplicated instances without posterior sampling.")
        parameters <- parameters[-not_post.sample[unique_not_post.sample], ]
    }
    # is.post.sample = TRUE with wrong data type
    wrong_post.sample <- which(parameters[, "is.post.sample"] & 
        parameters[, "feature.dat.type"] != "count")
    if(length(wrong_post.sample) > 0){
        message("Removing posterior sampling istances for non-count data.")
        parameters <- parameters[-wrong_post.sample, ]
    }
    # data.frame to list
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        x <- append(x = x, values = list("contrast" = contrast, 
            "strata" = strata, "adj.name" = adj.name), after = 2)
        x <- append(x = x, values = list("link.func" = link.func), after = 11)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}
