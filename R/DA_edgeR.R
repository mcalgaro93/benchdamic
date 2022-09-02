#' @title DA_edgeR
#'
#' @importFrom edgeR DGEList estimateDisp estimateGLMRobustDisp glmQLFit
#' glmQLFTest getDispersion
#' @importFrom stats model.matrix p.adjust
#' @export
#' @description
#' Fast run for edgeR differential abundance detection method.
#'
#' @inheritParams get_counts_metadata
#' @param pseudo_count add 1 to all counts if TRUE (default
#' \code{pseudo_count = FALSE}).
#' @param group_name character giving the name of the column containing
#' information about experimental group/condition for each sample/library.
#' @param design character or formula to specify the model matrix.
#' @param robust logical, should the estimation of \code{prior.df} be
#' robustified against outliers?
#' @param coef integer or character index vector indicating which coefficients
#' of the linear model are to be tested equal to zero.
#' @param norm name of the normalization method to use in the differential 
#' abundance analysis. Choose between the native edgeR normalization methods, 
#' such as \code{TMM}, \code{TMMwsp}, \code{RLE}, \code{upperquartile}, 
#' \code{posupperquartile}, or \code{none}. Alternatively (only for advanced 
#' users), if \code{norm} is equal to "ratio", "poscounts", or "iterate" from 
#' \code{\link{norm_DESeq2}}, "CSS" from \code{\link{norm_CSS}}, or "TSS" from 
#' \code{\link{norm_TSS}}, the scaling factors are automatically transformed 
#' into normalization factors. If custom factors are supplied, make sure they 
#' are compatible with edgeR normalization factors.
#' @param weights an optional numeric matrix giving observational weights.
#' @param verbose an optional logical value. If \code{TRUE}, information about
#' the steps of the algorithm is printed. Default \code{verbose = TRUE}.
#'
#' @return A list object containing the matrix of p-values \code{pValMat},
#' the dispersion estimates \code{dispEsts}, the matrix of summary statistics
#' for each tag \code{statInfo}, and a suggested \code{name} of the final object
#' considering the parameters passed to the function.
#'
#' @seealso \code{\link[edgeR]{DGEList}} for the edgeR DEG object creation,
#' \code{\link[edgeR]{estimateDisp}} and
#' \code{\link[edgeR]{estimateGLMRobustDisp}} for dispersion estimation, and
#' \code{\link[edgeR]{glmQLFit}} and \code{\link[edgeR]{glmQLFTest}} for the
#' quasi-likelihood negative binomial model fit.
#'
#' @examples
#' set.seed(1)
#' # Create a very simple phyloseq object
#' counts <- matrix(rnbinom(n = 60, size = 3, prob = 0.5), nrow = 10, ncol = 6)
#' metadata <- data.frame("Sample" = c("S1", "S2", "S3", "S4", "S5", "S6"),
#'                        "group" = as.factor(c("A", "A", "A", "B", "B", "B")))
#' ps <- phyloseq::phyloseq(phyloseq::otu_table(counts, taxa_are_rows = TRUE),
#'                          phyloseq::sample_data(metadata))
#'
#' # Calculate the TMM normalization factors
#' ps_NF <- norm_edgeR(object = ps, method = "TMM")
#' # The phyloseq object now contains the normalization factors:
#' normFacts <- phyloseq::sample_data(ps_NF)[, "NF.TMM"]
#' head(normFacts)
#'
#' # Differential abundance
#' DA_edgeR(object = ps_NF, pseudo_count = FALSE, group_name = "group",
#'          design = ~ group, coef = 2, robust = FALSE, norm = "TMM")

DA_edgeR <- function(object, assay_name = "counts", pseudo_count = FALSE, 
    group_name = NULL, design = NULL, robust = FALSE, coef = 2, norm = c("TMM", 
    "TMMwsp", "RLE", "upperquartile", "posupperquartile", "none"), weights, 
    verbose = TRUE) {
    counts_and_metadata <- get_counts_metadata(object, assay_name = assay_name)
    counts <- counts_and_metadata[[1]]
    libSizes <- colSums(counts)
    metadata <- counts_and_metadata[[2]]
    is_phyloseq <- counts_and_metadata[[3]]
    # Name building
    name <- "edgeR"
    method <- "DA_edgeR"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count) {
        if(verbose)
            message("Adding a pseudo count... \n")
        counts <- counts + 1
        name <- paste(name, ".pseudo", sep = "")
    }
    # Check the assay_name
    if (!is_phyloseq){
        if(verbose)
            message("Using the ", assay_name, " assay.")
        name <- paste(name, ".", assay_name, sep = "")
    } 
    if (length(norm) > 1) {
        stop(method, "\n",
            "norm: please choose one normalization for this istance of", 
            " differential abundance analysis.")
    }
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
    this_norms <- c("TMM", "TMMwsp", "RLE", "upperquartile", 
        "posupperquartile", "none")
    other_norms <- c("poscounts", "ratio", "iterate", "CSS", "TSS")
    if(!is.element(norm, this_norms)) {
        if(verbose){
            warning(method, "\n",
                norm, " normalization is not a native edgeR",
                " normalization. Make sure you know what you are doing,", 
                " otherwise choose between 'TMM', 'TMMwsp', 'RLE',", 
                " 'upperquartile', 'posupperquartile', or 'none'.")
        } 
        if(is.element(norm, other_norms)){
            if(verbose)
                message("Automatically converting NF.", norm, 
                    " scaling/size factors to normalization factors.")
            if(norm == "TSS"){
                NFs <- 1
                if(verbose)
                    message("'TSS' normalization: automatically setting NFs",
                        " to ones.")
            } else {
                NFs <- (NFs/libSizes)/exp(mean(log(NFs/libSizes)))
            }
        } else {
            if(prod(NFs) != 1){
                stop(method, "\n",
                    "norm: the custom normalization factors do not multiply",
                    " to 1. Make sure that the provided NF.", norm, " factors",
                    " are compatible with edgeR normalization factors.")
            } else {
                if(verbose){
                    warning(method, "\n",
                        "Make sure that the provided NF.", norm, " factors", 
                        " are compatible with edgeR normalization factors.")
                } 
            }
        }
    }
    name <- paste(name, ".", norm, sep = "")
    if(verbose)
        message("Differential abundance on ", norm, " normalized data")
    group <- unlist(metadata[, group_name])
    dge <- edgeR::DGEList(counts = counts, norm.factors = NFs, group = group)
    if (missing(weights)) {
        if(verbose)
            message("Estimating Differential Abundance without weighting")
    } else {
        if(is.null(weights)){
            if(verbose)
                message("Estimating Differential Abundance without weighting")
        } else {
            if(verbose)
                message("Estimating Differential Abundance with weights")
            dge$weights <- weights
            name <- paste(name, ".weighted", sep = "")
        }
    }
    if (is.character(design)) {
        if (grepl("~", design)) {
            design <- as.formula(design)
        } else {
            design <- as.formula(paste0("~", design))
        }
    }
    if (is(design, "formula")) {
        design <- stats::model.matrix(object = design, data = data.frame(
            metadata
        ))
    }
    if (!robust) {
        dge <- edgeR::estimateDisp(y = dge, design = design)
    } else {
        if(verbose)
            message("Estimating robust dispersions")
        dge <- edgeR::estimateGLMRobustDisp(y = dge, design = design)
        name <- paste(name, ".robust", sep = "")
    }
    if(verbose)
        message("Extracting results")
    dispEsts <- edgeR::getDispersion(dge)
    glmFit <- edgeR::glmQLFit(
        y = dge, dispersion = dispEsts,
        robust = robust, design = design
    )
    coefNames <- colnames(coef(glmFit))
    # coef check
    if(is.character(coef)){
        coef <- which(coefNames == coef)
        if(length(coef) == 0){
            stop(method, "\n", 
                 "Wrong 'coef' argument. Please choose numbers from 1 to ",
                 length(coefNames), " or names among: ",
                 paste(coefNames, collapse = ", "))
        }
    } else if(is.numeric(coef)){
        if(sum(coef > length(coefNames)) > 0){
            stop(method, "\n", 
                 "Wrong 'coef' argument. Please choose numbers from 1 to ",
                 length(coefNames), " or names among: ",
                 paste(coefNames, collapse = ", "))
        }
    }
    if(verbose){
        message("Extracting results for ", 
            paste(coefNames[coef], collapse = ","))
    }
    glmRes <- edgeR::glmQLFTest(glmFit, coef = coef)
    statInfo <- glmRes[["table"]]
    pval <- statInfo[, "PValue"]
    pValMat <- data.frame("rawP" = pval, "adjP" = stats::p.adjust(pval, "BH"))
    rownames(pValMat) <- rownames(statInfo)
    return(list(
        "pValMat" = pValMat, "statInfo" = statInfo, "dispEsts" =
            dispEsts, "name" = name
    ))
} # END - function: DA_edgeR

#' @title set_edgeR
#'
#' @export
#' @description
#' Set the parameters for edgeR differential abundance detection method.
#'
#' @inheritParams DA_edgeR
#' @param weights_logical logical vector, if true a matrix of observation
#' weights must be supplied (default \code{weights_logical = FALSE}).
#' @param expand logical, if TRUE create all combinations of input parameters
#' (default \code{expand = TRUE}).
#'
#' @return A named list containing the set of parameters for \code{DA_edgeR}
#' method.
#'
#' @seealso \code{\link{DA_edgeR}}
#'
#' @examples
#' # Set some basic combinations of parameters for edgeR
#' base_edgeR <- set_edgeR(group_name = "group", design = ~ group, coef = 2)
#'
#' # Set a specific set of normalization for edgeR
#' setNorm_edgeR <- set_edgeR(group_name = "group", design = ~ group, coef = 2,
#'     norm = c("TMM", "RLE"))
#'
#' # Set many possible combinations of parameters for edgeR
#' all_edgeR <- set_edgeR(pseudo_count = c(TRUE, FALSE), group_name = "group",
#'     design = ~ group, robust = c(TRUE, FALSE), coef = 2,
#'     weights_logical = c(TRUE, FALSE))

set_edgeR <- function(assay_name = "counts", pseudo_count = FALSE, 
    group_name = NULL, design = NULL, robust = FALSE, coef = 2, norm = c("TMM", 
    "TMMwsp", "RLE", "upperquartile", "posupperquartile", "none"), 
    weights_logical = FALSE, expand = TRUE) {
    method <- "DA_edgeR"
    if (is.null(assay_name)) {
        stop(method, "\n",             
            "'assay_name' is required (default = 'counts'). ")
    }
    if (!is.logical(pseudo_count) | !is.logical(robust) |
        !is.logical(weights_logical)) {
        stop(method, "\n",
            "'pseudo_count', 'robust', and 'weights_logical' must be logical.")
    }
    if (is.null(group_name) | is.null(design) | is.null(coef)) {
        stop(method, "\n",
            "'group_name', 'design', and 'coef' are required.")
    }
    if (!is.character(group_name)){
        stop(method, "\n", "'group_name' should be a character.")
    }
    if (!is.character(design) & !is(design, "formula")){
        stop(method, "\n", "'design' should be a character or a formula.")
    } else design <- as.formula(design)
    if (sum(!is.element(norm, c("TMM", "TMMwsp", "RLE", "upperquartile", 
        "posupperquartile", "none"))) > 0) {
        warning(method, "\n",
            "One or more elements into 'norm' are not native to edgeR.")
    }
    if (expand) {
        parameters <- expand.grid(method = method, assay_name = assay_name, 
            pseudo_count = pseudo_count, robust = robust, norm = norm, 
            weights = weights_logical, stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, assay_name = assay_name, 
            pseudo_count = pseudo_count, robust = robust, norm = norm, 
            weights = weights_logical)
    }
    # data.frame to list
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        x <- append(x = x, values = list("group_name" = group_name,
                                         "design" = deparse(design),
                                         "coef" = coef), after = 3)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}
