#' @title DA_limma
#'
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom stats model.matrix weights
#' @export
#' @description
#' Fast run for limma voom differential abundance detection method.
#'
#' @inheritParams DA_edgeR
#' @param design character name of the metadata columns, formula, or design
#' matrix with rows corresponding to samples and columns to coefficients to be
#' estimated.
#' @inheritParams limma::topTable
#'
#' @return A list object containing the matrix of p-values `pValMat`, the matrix
#' of summary statistics for each tag `statInfo`, and a suggested `name` of the
#' final object considering the parameters passed to the function.
#'
#' @seealso \code{\link[limma]{voom}} for the mean-variance relationship
#' estimation, \code{\link[limma]{lmFit}} for the linear model framework.
#'
#' @examples
#' set.seed(1)
#' # Create a very simple phyloseq object
#' counts <- matrix(rnbinom(n = 60, size = 3, prob = 0.5), nrow = 10, ncol = 6)
#' metadata <- data.frame("Sample" = c("S1", "S2", "S3", "S4", "S5", "S6"),
#'                        "group" = as.factor(c("A", "A", "A", "B", "B", "B")))
#' ps <- phyloseq::phyloseq(phyloseq::otu_table(counts, taxa_are_rows = TRUE),
#'                          phyloseq::sample_data(metadata))
#' # Calculate the TMM normalization factors
#' ps_NF <- norm_edgeR(object = ps, method = "TMM")
#' # The phyloseq object now contains the normalization factors:
#' normFacts <- phyloseq::sample_data(ps_NF)[, "NF.TMM"]
#' head(normFacts)
#' # Differential abundance
#' DA_limma(object = ps_NF, pseudo_count = FALSE, design = ~ group, coef = 2,
#'     norm = "TMM")

DA_limma <- function(object, assay_name = "counts", pseudo_count = FALSE, 
    design = NULL, coef = 2, norm = c("TMM", "TMMwsp", "RLE", "upperquartile", 
    "posupperquartile", "none"), weights, verbose = TRUE){
    counts_and_metadata <- get_counts_metadata(object, assay_name = assay_name)
    counts <- counts_and_metadata[[1]]
    libSizes <- colSums(counts)
    metadata <- counts_and_metadata[[2]]
    is_phyloseq <- counts_and_metadata[[3]]
    # Name building
    name <- "limma"
    method <- "DA_limma"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        if(verbose)
            message("Adding a pseudo count... \n")
        counts <- counts + 1
        name <- paste(name,".pseudo",sep = "")}
    # Check the assay
    if (!is_phyloseq){
        if(verbose)
            message("Using the ", assay_name, " assay.")
        name <- paste(name, ".", assay_name, sep = "")
    } 
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
    this_norms <- c("TMM", "TMMwsp", "RLE", "upperquartile", 
        "posupperquartile", "none")
    other_norms <- c("poscounts", "ratio", "iterate", "CSS", "TSS")
    if(!is.element(norm, this_norms)) {
        if(verbose){
            warning(method, "\n", norm, " normalization is not a native edgeR",
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
        message("Differential abundance on ", norm," normalized data")
    if(is.character(design)){
        if(grepl("~", design))
            design <- as.formula(design)
        else design <- as.formula(paste0("~", design))}
    if(is(design, "formula"))
        design <- stats::model.matrix(object = design, data = data.frame(
            metadata))
    v <- limma::voom(counts = counts, design = design, 
        lib.size = libSizes * NFs, plot = FALSE)
    if(missing(weights)){
        if(verbose)
            message("Estimating Differential Abundance without weighting")
        fit <- limma::lmFit(object = v, design = design)
    } else {
        if(is.null(weights)){
            if(verbose)
                message("Estimating Differential Abundance without weighting")
            fit <- limma::lmFit(object = v, design = design)
        } else {
            if(verbose)
                message("Estimating Differential Abundance with weights")
            name <- paste(name,".weighted",sep = "")
            fit <- limma::lmFit(object = v, design = design, weights =
                stats::weights(v) * weights)
        }
    }
    fit <- limma::eBayes(fit)
    if(verbose)
        message("Extracting results for ", colnames(coef(fit[, coef])),
            " coefficient")
    statInfo <- limma::topTable(fit, coef = coef, n = nrow(counts), sort.by =
        "none")
    pValMat <- statInfo[, c("P.Value", "adj.P.Val")]
    colnames(pValMat) <- c("rawP", "adjP")
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
}# END - function: DA_limma

#' @title set_limma
#'
#' @export
#' @description
#' Set the parameters for limma differential abundance detection method.
#'
#' @inheritParams DA_limma
#' @param weights_logical logical vector, if TRUE a matrix of observational
#' weights will be used for differential abundance analysis (default
#' \code{weights_logical = FALSE}).
#' @param expand logical, if TRUE create all combinations of input parameters
#' (default \code{expand = TRUE}).
#'
#' @return A named list containing the set of parameters for \code{DA_limma}
#' method.
#'
#' @seealso \code{\link{DA_limma}}
#'
#' @examples
#' # Set some basic combinations of parameters for limma
#' base_limma <- set_limma(design = ~ group, coef = 2)
#' # Set a specific set of normalization for limma (even of other packages!)
#' setNorm_limma <- set_limma(design = ~ group, coef = 2,
#'     norm = c("TMM", "upperquartile"))
#' # Set many possible combinations of parameters for limma
#' all_limma <- set_limma(pseudo_count = c(TRUE, FALSE), design = ~ group,
#'     coef = 2, weights_logical = c(TRUE, FALSE))

set_limma <- function(assay_name = "counts", pseudo_count = FALSE, 
    design = NULL, coef = 2, norm = c("TMM", "TMMwsp", "RLE", "upperquartile", 
    "posupperquartile", "none"), weights_logical = FALSE, expand = TRUE) {
    method <- "DA_limma"
    if (is.null(assay_name)) {
        stop(method, "\n", "'assay_name' is required (default = 'counts').")
    }
    if (!is.logical(pseudo_count) | !is.logical(weights_logical)) {
        stop(method, "\n", 
            "'pseudo_count' and 'weights_logical' must be logical.")
    }
    if (is.null(design) | is.null(coef)) {
        stop(method, "\n", "'design', and 'coef' are required.")
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
            pseudo_count = pseudo_count, norm = norm, weights = weights_logical,
            stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, assay_name = assay_name,
            pseudo_count = pseudo_count, norm = norm, weights = weights_logical)
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
