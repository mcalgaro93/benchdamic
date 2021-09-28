#' @title DA_limma
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data
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
#' # Calculate the TMM scaling factors
#' ps_NF <- norm_edgeR(object = ps, method = "TMM")
#' # The phyloseq object now contains the scaling factors:
#' scaleFacts <- phyloseq::sample_data(ps_NF)[, "NF.TMM"]
#' head(scaleFacts)
#' # Differential abundance
#' DA_limma(object = ps_NF, pseudo_count = FALSE, design = ~ group, coef = 2,
#'     norm = "TMM")

DA_limma <- function(object, pseudo_count = FALSE, design = NULL, coef = 2,
    norm = c("TMM", "TMMwsp", "RLE", "upperquartile", "posupperquartile",
        "none", "ratio", "poscounts", "iterate", "TSS", "CSSmedian",
        "CSSdefault"), weights, verbose = TRUE){
    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)
    # Name building
    name <- "limma"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        if(verbose)
            message("Adding a pseudo count... \n")
        counts <- counts + 1
        name <- paste(name,".pseudo",sep = "")}
    if(length(norm) > 1)
        stop("Please choose one normalization for this istance of differential
        abundance analysis.")
    if(norm == "TSS")
        NFs = 1
    else {
        # Check if the column with the normalization factors is present
        NF.col <- paste("NF", norm, sep = ".")
        if(!any(colnames(metadata) == NF.col))
            stop("Can't find the ", NF.col," column in your object.",
                " Make sure to add the normalization factors column in your",
                " object first.")
        NFs = unlist(metadata[,NF.col])
        # DESeq2 NFs are supplied -> make them scaling factors!
        if(is.element(norm, c("ratio", "poscounts", "iterate")))
            NFs <- NFs/colSums(counts)
        NFs = NFs/exp(mean(log(NFs)))} # Make NFs multiply to 1
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
    v <- limma::voom(counts = counts, design = design, lib.size = colSums(
        counts) * NFs, plot = FALSE)
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
#'     norm = c("TMM", "poscounts"))
#' # Set many possible combinations of parameters for limma
#' all_limma <- set_limma(pseudo_count = c(TRUE, FALSE), design = ~ group,
#'     coef = 2, weights_logical = c(TRUE,FALSE))

set_limma <- function(pseudo_count = FALSE, design = NULL, coef = 2,
    norm = c("TMM", "TMMwsp", "RLE", "upperquartile", "posupperquartile",
        "none"), weights_logical = FALSE, expand = TRUE) {
    method <- "DA_limma"
    if (!is.logical(pseudo_count) | !is.logical(weights_logical)) {
        stop("'pseudo_count' and 'weights_logical' must be logical.")
    }
    if (is.null(design) | is.null(coef)) {
        stop("'design', and 'coef' are required.")
    }
    if (!is.character(design) & !is(design, "formula")){
        stop("'design' should be a character or a formula.")
    } else design <- as.formula(design)
    if (sum(!is.element(norm, c(
        "TMM", "TMMwsp", "RLE", "upperquartile",
        "posupperquartile", "none", "TSS"
    ))) > 0) {
        warning("One or more elements into 'norm' are not native to edgeR.")
    }
    if (expand) {
        parameters <- expand.grid(method = method, pseudo_count = pseudo_count,
            norm = norm, weights = weights_logical, stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, pseudo_count = pseudo_count,
            norm = norm, weights = weights_logical)
    }
    # data.frame to list
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        x <- append(x = x, values = list("design" = deparse(design),
                                         "coef" = coef), after = 2)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}
