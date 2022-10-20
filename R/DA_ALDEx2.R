#' @title DA_ALDEx2
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data
#' @importFrom ALDEx2 aldex
#' @importFrom stats p.adjust
#' @export
#' @description
#' Fast run for the ALDEx2's differential abundance detection method.
#' Support for Welch's t, Wilcoxon, Kruskal-Wallace, Kruskal-Wallace 
#' glm ANOVA-like, and glm tests.
#'
#' @inheritParams DA_DESeq2
#' @param design a character with the name of a variable to group samples and 
#' compare them or a formula to compute a model.matrix (when 
#' \code{test = "glm"}).
#' @param mc.samples an integer. The number of Monte Carlo samples to use when
#' estimating the underlying distributions. Since we are estimating central
#' tendencies, 128 is usually sufficient.
#' @inheritParams ALDEx2::aldex.clr
#' @param test a character string. Indicates which tests to perform. "t" runs
#' Welch's t test while "wilcox" runs Wilcoxon test. "kw" runs 
#' Kruskal-Wallace test while "kw_glm" runs glm ANOVA-like test. "glm" runs a 
#' generalized linear model.
#' @param contrast character vector with exactly three elements: the name of a
#' variable used in "design", the name of the level of interest, and the 
#' name of the reference level. If "kw" or "kw_glm" as test, contrast vector is
#' not used.
#' @inheritParams ALDEx2::aldex
#'
#' @return A list object containing the matrix of p-values `pValMat`, the matrix
#' of summary statistics for each tag `statInfo`, and a suggested `name` of the
#' final object considering the parameters passed to the function.
#'
#' @seealso \code{\link[ALDEx2]{aldex}} for the Dirichlet-Multinomial model
#' estimation. Several and more complex tests are present in the ALDEx2
#' framework.
#'
#' @examples
#' set.seed(1)
#' # Create a very simple phyloseq object
#' counts <- matrix(rnbinom(n = 300, size = 3, prob = 0.5), nrow = 50, ncol = 6)
#' metadata <- data.frame("Sample" = c("S1", "S2", "S3", "S4", "S5", "S6"),
#'                        "group" = as.factor(c("A", "A", "A", "B", "B", "B")))
#' ps <- phyloseq::phyloseq(phyloseq::otu_table(counts, taxa_are_rows = TRUE),
#'                          phyloseq::sample_data(metadata))
#' # Differential abundance with t test and denom defined by the user
#' DA_t <- DA_ALDEx2(ps, design = "group", test = "t", denom = c(1,2,3), 
#'     paired.test = FALSE, contrast = c("group", "B", "A"))
#' # Differential abundance with wilcox test and denom = "iqlr"
#' DA_w <- DA_ALDEx2(ps, design = "group", test = "wilcox", denom = "iqlr", 
#'     paired.test = FALSE, contrast = c("group", "B", "A"))
#' # Differential abundance with kw test and denom = "zero"
#' # mc.samples = 2 to speed up (128 suggested)
#' DA_kw <- DA_ALDEx2(ps, design = "group", test = "kw", denom = "zero", 
#'     mc.samples = 2)
#' # Differential abundance with kw_glm test and denom = "median"
#' DA_kw_glm <- DA_ALDEx2(ps, design = "group", test = "kw", denom = "median", 
#'     mc.samples = 2)
#' # Differential abundance with glm test and denom = "all"
#' DA_glm <- DA_ALDEx2(ps, design = ~ group, test = "glm", denom = "all", 
#'     mc.samples = 2, contrast = c("group", "B", "A"))

DA_ALDEx2 <- function(object, assay_name = "counts", pseudo_count = FALSE, 
    design = NULL, mc.samples = 128, test = c("t", "wilcox", "kw", "kw_glm", 
    "glm"), paired.test = FALSE, denom = "all", contrast = NULL, 
    verbose = TRUE){
    counts_and_metadata <- get_counts_metadata(object, assay_name = assay_name)
    counts <- counts_and_metadata[[1]]
    metadata <- counts_and_metadata[[2]]
    is_phyloseq <- counts_and_metadata[[3]]
    # Name building
    name <- "ALDEx2"
    method <- "DA_ALDEx2"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        if(verbose)
            message("Adding a pseudo count... \n")
        counts <- counts + 1
        name <- paste(name,".pseudo",sep = "")
    }
    # Check the assay_name
    if (!is_phyloseq){
        if(verbose)
            message("Using the ", assay_name, " assay.")
        name <- paste(name, ".", assay_name, sep = "")
    } 
    # Check test
    if(length(test) != 1){
        stop(method, "\n", 
             "test: please choose only one test between 't', 'wilcox',", 
             " 'kw', 'kw_glm', or 'glm' for this instance of differential",  
             " abundance analysis.")
    }
    if(!is.element(test, c("t", "wilcox", "kw", "kw_glm", "glm")))
        stop(method, "\n", 
             "test: please choose one test between 't', 'wilcox', 'kw',",
             " 'kw_glm', or 'glm'.")
    # Check contrast
    if(!is.element(test, c("kw", "kw_glm"))){
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
    }
    # Check design
    conds <- NULL
    if (is.character(design)) { # When it is a character
        if (grepl("~", design)) { # It could be a formula
            if(test != "glm")
                stop(method, "\n", 
                     "design: formula type is accepted only with test = 'glm'.",
                    " Please supply the name of a grouping factor instead.")
            design <- as.formula(design)
        } else { # Or it is just the name of a variable in metadata
            if(test == "glm")
                stop(method, "\n", 
                     "test: 'glm' test requires a formula, not a variable", 
                    " name in 'design'.")
            if(!is.element(design, colnames(metadata))){
                stop(method, "\n", 
                     "'design' = ", design, " but there is not a", 
                     " metadata column with this name.")
            } # No need to check for contrast if test is kw or kw_glm
            if(!is.element(test, c("kw", "kw_glm"))){
                if(contrast[1] != design){
                    stop(method, "\n", 
                         "'design' = ", design, " but 'contrast' is based on", 
                        " the ", contrast[1], " variable. They should match.")
                }
            }
            conds <- unlist(metadata[, design])
        }
    } 
    if (is(design, "formula")) {
        if(test != "glm")
            stop(method, "\n", 
                 "design: formula type is accepted only with test = 'glm'.",
                 " Please supply the name of a grouping factor instead.")
        conds <- stats::model.matrix(object = design, data = data.frame(
            metadata))
        if(!is.element(paste0(contrast[1], contrast[2]), colnames(conds))){
            stop(method, "\n", 
                 "contrast: 'contrast' is based on the ", contrast[1], 
                " variable but it is not included into the design formula or", 
                " the ", contrast[2], " level is not the reference.")
        }
    }
    if(is.null(design) | is.null(conds)) 
        stop(method, "\n", 
            "design: Please supply a character with the name of a variable to", 
            " group samples and compare them or a formula to compute a", 
            " model.matrix (if test = 'glm').")
    # Check denom
    denomName <- NULL
    if(is.numeric(denom)){ # they are user defined indexes
        denomName <- "user_denom" # we can't use all of them for naming!
    } else { 
        if(length(denom) > 1)
            stop(method, "\n", 
                 "denom: Please choose only one denom for this instance of",
                " differential abundance analysis.")
        denomName <- denom
        if(!is.element(denom, c("all", "iqlr", "zero", "lvha", "median"))){
            stop(method, "\n", 
                 "denom: Please choose one denom between 'all', 'iqlr',",  
                 " 'zero' 'lvha', or 'median'. Otherwise supply a vector", 
                 " of row indices to use as denominator.")
        } else {
            if(test == "glm" & denom != "all")
                stop(method, "\n", 
                     "denom: when test = 'glm' denom must be 'all'.")
        }
    }
    name <- paste(name, ".", denomName, sep = "")
    name <- paste(name, ".", test, sep = "")
    # Check paired.test and add it to name
    test_to_do <- NULL
    if(is.element(test, c("t", "wilcox"))){
        if(paired.test){
            name <- paste(name, ".", "paired", sep = "")
        } else name <- paste(name, ".", "unpaired", sep = "") 
        test_to_do <- "t"
    } else if (is.element(test, c("kw", "kw_glm"))){
        test_to_do <- "kw"
    } else test_to_do <- "glm"
    if(verbose & !is.element(test, c("kw", "kw_glm"))){
        message("Extracting results for ", contrast[1], " variable, ",
                contrast[2], " (reference level = ", contrast[3], ")")
    }
    # Perform the analysis
    if(verbose){
        statInfo <- ALDEx2::aldex(reads = counts, conditions = conds,
            mc.samples = mc.samples, test = test_to_do, 
            paired.test = paired.test, denom = denom, verbose = verbose)
    } else {
        statInfo <- suppressMessages(ALDEx2::aldex(reads = counts, 
            conditions = conds, mc.samples = mc.samples, test = test_to_do, 
            paired.test = paired.test, denom = denom, verbose = verbose))
    }
    if(test == "t"){
        pValMat <- data.frame(statInfo[, c("we.ep", "we.eBH")])
    } else if(test == "wilcox"){
        pValMat <- data.frame(statInfo[, c("wi.ep", "wi.eBH")])
    } else if(test == "kw"){
        pValMat <- data.frame(statInfo[, c("kw.ep", "kw.eBH")])
    } else if(test == "kw_glm"){
        pValMat <- data.frame(statInfo[, c("glm.ep", "glm.eBH")])
    } else if(test == "glm"){
        p_val_col <- paste0(contrast[1], contrast[2], ".pval")
        if(verbose){
            message("Extracting p-values using the '", p_val_col, "' column.")
        }
        pValMat <- data.frame(statInfo[, c(p_val_col, p_val_col)])
        pValMat[, 2] <- stats::p.adjust(pValMat[, 2], method = "BH")
    }
    colnames(pValMat) <- c("rawP", "adjP")
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
}# END - function: DA_ALDEx2

#' @title set_ALDEx2
#'
#' @export
#' @description
#' Set the parameters for ALDEx2 differential abundance detection method.
#'
#' @inheritParams DA_ALDEx2
#' @param expand logical, if TRUE create all combinations of input parameters
#' (default \code{expand = TRUE})
#'
#' @return A named list containing the set of parameters for
#' \code{DA_ALDEx2} method.
#'
#' @seealso \code{\link{DA_ALDEx2}}
#'
#' @examples
#' # Set some basic combinations of parameters for ALDEx2
#' base_ALDEx2 <- set_ALDEx2(design = "group", 
#'     contrast = c("group", "grp2", "grp1"))
#' # Set a specific set of normalization for ALDEx2 (even of other
#' # packages!)
#' setNorm_ALDEx2 <- set_ALDEx2(design = "group", 
#'     contrast = c("group", "grp2", "grp1"))
#' # Set many possible combinations of parameters for ALDEx2
#' all_ALDEx2 <- set_ALDEx2(design = "group", denom = c("iqlr", "zero"),
#'     test = c("t", "wilcox"), contrast = c("group", "grp2", "grp1"))

set_ALDEx2 <- function(assay_name = "counts", pseudo_count = FALSE, 
    design = NULL, mc.samples = 128, test = "t", paired.test = FALSE, 
    denom = "all", contrast = NULL, expand = TRUE) {
    method <- "DA_ALDEx2"
    if (is.null(assay_name)) {
        stop(method, "\n", 
             "'assay_name' is required (default = 'counts'). ")
    }
    if (!is.logical(pseudo_count)) {
        stop(method, "\n", "'pseudo_count' must be logical.")
    }
    if(sum(!is.element(test, c("t", "wilcox", "kw", "kw_glm", "glm"))) > 0){
        stop(method, "\n", 
             "test: please choose one test between 't', 'wilcox', 'kw',",
             " 'kw_glm', or 'glm'.")
    }
    if (!is.logical(paired.test)) {
        stop(method, "\n", "'paired.test' must be logical.")
    }
    custom_denom <- FALSE
    if(sum(!is.element(denom, c("all", "iqlr", "zero", "lvha", "median"))) > 0){
        if(is.numeric(denom))
            custom_denom <- TRUE
        else{
            stop(method, "\n", 
                "denom: please choose denom between 'all', 'iqlr', 'zero',",
                " 'lvha', 'median' or a numeric vector of custom features to",
                " use as denom.")
        }
    }
    if (is.null(design) | is.null(contrast)) {
        stop(method, "\n", "'design' and 'contrast' must be specified.")
    }
    if (!is.character(design) & !is(design, "formula")){
        stop(method, "\n", "'design' should be a character or a formula.")
    }
    if (!is.character(contrast) & length(contrast) != 3){
        stop(method, "\n", 
            "contrast: Please supply a character vector with exactly", 
            " three elements: the name of a variable used in",  
            " 'design', the name of the level of interest, and the", 
            " name of the reference level.")
    }
    if (expand) {
        if(!custom_denom){
            parameters <- expand.grid(method = method, assay_name = assay_name,
                pseudo_count = pseudo_count, mc.samples = mc.samples, 
                test = test, paired.test = paired.test, 
                denom = denom, stringsAsFactors = FALSE)
        } else {
            parameters <- expand.grid(method = method, assay_name = assay_name,
                pseudo_count = pseudo_count, mc.samples = mc.samples, 
                test = test, paired.test = paired.test,
                stringsAsFactors = FALSE)
        }
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        if(!custom_denom){
            parameters <- data.frame(method = method, assay_name = assay_name,
                pseudo_count = pseudo_count, mc.samples = mc.samples, 
                test = test, paired.test = paired.test, 
                denom = denom)
        } else {
            parameters <- data.frame(method = method, assay_name = assay_name,
                pseudo_count = pseudo_count, mc.samples = mc.samples, 
                test = test, paired.test = paired.test)
        }
    }
    # Remove senseless combinations
    fake_paired <- which(parameters[, "paired.test"] == TRUE & 
        !is.element(parameters[, "test"], c("t", "wilcox")))
    if(length(fake_paired) > 0){
        message("Removing paired.test with test not equal to 't' or 'wilcox'.")
        parameters <- parameters[-fake_paired, ]
    }
    glm_not_denom <- which(parameters[, "test"] == "glm" & 
        parameters[, "test"] != "all")
    if(length(glm_not_denom) > 0){
        message("Removing combinations with denom not equal to 'all' when test",
            " is equal to 'glm'. Those options are not supported.")
        parameters <- parameters[- glm_not_denom, ]
    }
    # data.frame to list
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        x <- append(x = x, values = list("design" = design), after = 3)
        x <- append(x = x, values = list("contrast" = contrast), after = 8)
        if(custom_denom){
            x <- append(x = x, values = list("denom" = denom), after = 8)
        }
        return(x)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}
