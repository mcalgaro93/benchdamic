#' @title DA_metagenomeSeq
#'
#' @importFrom phyloseq otu_table sample_data taxa_names
#' @importFrom phyloseq phyloseq_to_metagenomeSeq
#' @importFrom metagenomeSeq fitZig fitFeatureModel MRfulltable
#' @importFrom stats model.matrix median
#' @export
#' @description
#' Fast run for the metagenomeSeq's differential abundance detection method.
#'
#' @inheritParams DA_edgeR
#' @param design the model for the count distribution. Can be the variable name,
#' or a character similar to "~ 1 + group", or a formula.
#' @param coef coefficient of interest to grab log fold-changes.
#' @param norm name of the normalization method to use in the differential 
#' abundance analysis. Choose the native metagenomeSeq normalization method 
#' \code{CSS}. Alternatively (only for advanced users), if \code{norm} is equal 
#' to "TMM", "TMMwsp", "RLE", "upperquartile", "posupperquartile", or "none" 
#' from \code{\link{norm_edgeR}}, "ratio", "poscounts", or "iterate" from 
#' \code{\link{norm_DESeq2}}, or "TSS" from \code{\link{norm_TSS}}, the 
#' factors are automatically transformed into scaling factors. If custom 
#' factors are supplied, make sure they are compatible with metagenomeSeq 
#' normalization factors.
#' @param model character equal to "fitFeatureModel" for differential abundance 
#' analysis using a zero-inflated log-normal model, "fitZig" for a complex 
#' mathematical optimization routine to estimate probabilities that a zero for 
#' a particular feature in a sample is a technical zero or not. The latter model
#' relies heavily on the limma package (default 
#' \code{model = "fitFeatureModel"}).
#' @inheritParams metagenomeSeq::fitZig
#' @inheritParams metagenomeSeq::MRcoefs
#'
#' @return A list object containing the matrix of p-values `pValMat`, the matrix
#' of summary statistics for each tag `statInfo`, and a suggested `name` of the
#' final object considering the parameters passed to the function.
#'
#' @seealso \code{\link[metagenomeSeq]{fitZig}} for the Zero-Inflated Gaussian
#' regression model estimation and \code{\link[metagenomeSeq]{MRfulltable}}
#' for results extraction.
#'
#' @examples
#' set.seed(1)
#' # Create a very simple phyloseq object
#' counts <- matrix(rnbinom(n = 60, size = 3, prob = 0.5), nrow = 10, ncol = 6)
#' metadata <- data.frame("Sample" = c("S1", "S2", "S3", "S4", "S5", "S6"),
#'                        "group" = as.factor(c("A", "A", "A", "B", "B", "B")))
#' ps <- phyloseq::phyloseq(phyloseq::otu_table(counts, taxa_are_rows = TRUE),
#'                          phyloseq::sample_data(metadata))
#' # Calculate the CSS normalization factors
#' ps_NF <- norm_CSS(object = ps, method = "CSS")
#' # The phyloseq object now contains the normalization factors:
#' normFacts <- phyloseq::sample_data(ps_NF)[, "NF.CSS"]
#' head(normFacts)
#' # Differential abundance
#' DA_metagenomeSeq(object = ps_NF, pseudo_count = FALSE, design = ~ group,
#'     coef = 2, norm = "CSS")

DA_metagenomeSeq <- function(object, assay_name = "counts", 
    pseudo_count = FALSE, design = NULL, coef = 2, norm = "CSS", 
    model = "fitFeatureModel", verbose = TRUE){
    counts_and_metadata <- get_counts_metadata(object, assay_name = assay_name)
    counts <- counts_and_metadata[[1]]
    taxa_names <- rownames(counts)
    libSizes <- colSums(counts)
    metadata <- counts_and_metadata[[2]]
    is_phyloseq <- counts_and_metadata[[3]]
    # Name building
    name <- "metagenomeSeq"
    method <- "DA_metagenomeSeq"
    if(!is_phyloseq){
        object <- phyloseq::phyloseq(
            otu_table = phyloseq::otu_table(counts, taxa_are_rows = TRUE), 
            sample_data = phyloseq::sample_data(object = data.frame(metadata)))
        if(verbose)
            message("Using the ", assay_name, " assay.")
        name <- paste(name, ".", assay_name, sep = "")
    } 
    obj <- phyloseq::phyloseq_to_metagenomeSeq(object)
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        if(verbose)
            message("Adding a pseudo count...")
        counts <- counts + 1
        name <- paste(name,".pseudo",sep = "")}
    if(length(norm) > 1)
        stop(method, "\n", 
            "norm: please choose one normalization for this istance of",
            " differential abundance analysis.")
    # Check if the column with the normalization factors is present
    NF.col <- paste("NF", norm, sep = ".")
    if (!any(colnames(metadata) == NF.col)) {
        stop(method, "\n", 
            "norm: can't find the ", NF.col," column in your object.",
            " Make sure to add the normalization factors column in your",
            " object first.")
    }
    NFs <- unlist(metadata[, NF.col])
    # Check for implemented normalizations
    this_norms <- "CSS"
    scaling_norms <- c("ratio", "iterate", "poscounts")
    noscaling_norms <- c("TMM", "TMMwsp", "RLE", "upperquartile", 
    "posupperquartile", "none", "TSS")
    if(!is.element(norm, this_norms)) {
        if(verbose){
            warning(method, "\n", norm, 
                " normalization is not a native metagenomeSeq",
                " normalization. Make sure you know what you are doing,", 
                " otherwise choose 'CSS'.")
        } 
        if(is.element(norm, noscaling_norms)){
            if(verbose)
                message("Automatically converting NF.", norm, 
                        " factors to normalization factors", 
                        " compatible with metagenomeSeq.")
            if(norm == "TSS"){
                NFs <- 1
                if(verbose)
                    message("'TSS' normalization: automatically setting NFs", 
                        " to ones.")
            }
            if(!is.element(norm, c("CSS", "TSS"))){ 
                NFs <- (NFs*libSizes)/exp(mean(log(NFs*libSizes)))
            }
        } else {
            if(verbose){
                warning(method, "\n", 
                        "Make sure that the provided NF.", norm, " factors", 
                        " are compatible with metagenomeSeq normalization", 
                        " factors.")
            } 
        }
    }
    name <- paste(name, ".", norm, sep = "")
    models <- c("fitZig", "fitFeatureModel")
    if(!is.element(model, models)) {
        stop(method, "\n", 
            "model: ", model, " is not supported. Please, choose between",
            " 'fitZig' or 'fitFeatureModel'.")
    }
    if(verbose)
        message("Differential abundance on ", norm, " normalized data using ",
            model, " method.")
    metadata[, NF.col] <- log2(NFs/stats::median(NFs) + 1) # used in fitZig
    metagenomeSeq::normFactors(object = obj) <- NFs # used in fitFeatureModel
    if(is.character(design)){
        if(grepl("~", design))
            design <- as.formula(design)
        else design <- as.formula(paste0("~", design))
    } else if(is(design, "formula")){
        if(model == "fitZig"){
            design <- as.formula(paste0(paste0(design,collapse = " "), " + ",
                NF.col))
        } 
    } else {
        stop(method, "\n", "'design' should be a character or a formula.")
    }
    design <- stats::model.matrix(object = design, 
        data = data.frame(metadata))
    control <- metagenomeSeq::zigControl(maxit = 1000, verbose = verbose)
    if(model == "fitZig"){
        fit <- try(metagenomeSeq::fitZig(obj = obj, mod = design,
            verbose = FALSE, useCSSoffset = FALSE, control = control))
    } else {
        if(verbose){
            fit <- try(metagenomeSeq::fitFeatureModel(obj = obj, mod = design))
        } else {
            fit <- suppressWarnings(try(metagenomeSeq::fitFeatureModel(
                obj = obj, mod = design)))
        }
    }
    name <- paste(name, ".", model, sep = "")
    if(is(fit, "try-error")){
        res <- matrix(NA, ncol = 2, nrow = nrow(counts))
        stop(method, "\n", 
            "Something went wrong during model estimation. Make sure you are",
            " using the correct 'model' (only two groups are allowed with",
            " 'fitFeatureModel').")
    } else {
        statInfo <- metagenomeSeq::MRcoefs(obj = fit, by = coef, number = nrow(
            counts))
        statInfo <- statInfo[taxa_names, ]
        if(verbose){
            if(model == "fitZig"){
                message("Extracting results for ", colnames(statInfo[coef]),
                    " coefficient")
            } else {
                coef_names <- colnames(design)
                coef_index <- which(coef_names == coef)
                message("Extracting results for ", coef_names[coef_index], 
                    " coefficient")
            }
        }
    }
    pValMat <- statInfo[, c("pvalues", "adjPvalues")]
    colnames(pValMat) <- c("rawP", "adjP")
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
}# END - function: DA_metagenomeSeq

#' @title set_metagenomeSeq
#'
#' @export
#' @description
#' Set the parameters for metagenomeSeq differential abundance detection method.
#'
#' @inheritParams DA_metagenomeSeq
#' @param expand logical, if TRUE create all combinations of input parameters
#' (default \code{expand = TRUE})
#'
#' @return A named list containing the set of parameters for
#' \code{DA_metagenomeSeq} method.
#'
#' @seealso \code{\link{DA_metagenomeSeq}}
#'
#' @examples
#' # Set a basic combination of parameters for metagenomeSeq
#' base_mgs <- set_metagenomeSeq(design = ~ group, coef = 2)
#' # Set a specific model for metagenomeSeq
#' setModel_mgs <- set_metagenomeSeq(design = ~ group, coef = 2, 
#'     model = "fitZig")
#' # Set many possible combinations of parameters for metagenomeSeq
#' all_mgs <- set_metagenomeSeq(pseudo_count = c(TRUE, FALSE), design = ~ group,
#'     coef = 2, model = c("fitFeatureModel", "fitZig"), norm = "CSS")

set_metagenomeSeq <- function(assay_name = "counts", pseudo_count = FALSE, 
    design = NULL, coef = 2, norm = "CSS", model = "fitFeatureModel", 
    expand = TRUE) {
    method <- "DA_metagenomeSeq"
    if (is.null(assay_name)) {
        stop(method, "\n", "'assay_name' is required (default = 'counts').")
    }
    if (!is.logical(pseudo_count)) {
        stop(method, "\n", "'pseudo_count' must be logical.")
    }
    if (is.null(design) | is.null(coef) | is.null(model)) {
        stop(method, "\n", "'design', 'coef', and 'model' are required.")
    }
    if (!is.character(design) & !is(design, "formula")){
        stop(method, "\n", "'design' should be a character or a formula.")
    } else design <- as.formula(design)
    if (sum(!is.element(model, c("fitFeatureModel", "fitZig"))) > 0) {
        stop(method, "\n", 
            "model: One or more elements into 'model' are not compatible with",
            " metagenomeSeq.")
    }
    if (sum(!is.element(norm, "CSS")) > 0) {
        warning(method, "\n", 
            "One or more elements into 'norm' are not native to",
            " metagenomeSeq.")
    }
    if (expand) {
        parameters <- expand.grid(method = method, assay_name = assay_name,
            pseudo_count = pseudo_count, norm = norm, model = model, 
            stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, assay_name = assay_name,
            pseudo_count = pseudo_count, norm = norm, model = model)
    }
    # data.frame to list
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        x <- append(x = x, values = list("design" = deparse(design),
            "coef" = coef), after = 3)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}

