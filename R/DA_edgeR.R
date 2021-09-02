#' @title DA_edgeR
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data
#' @importFrom edgeR DGEList estimateDisp estimateGLMRobustDisp glmQLFit
#' glmQLFTest getDispersion
#' @importFrom stats model.matrix p.adjust
#' @export
#' @description
#' Fast run for edgeR differential abundance detection method.
#'
#' @param object phyloseq object.
#' @param pseudo_count add 1 to all counts if TRUE (default
#' \code{pseudo_count = FALSE}).
#' @param group_name character giving the name of the column containing
#' information about experimental group/condition for each sample/library.
#' @param design character or formula to specify the model matrix.
#' @param robust logical, should the estimation of \code{prior.df} be
#' robustified against outliers?
#' @param coef integer or character index vector indicating which coefficients
#' of the linear model are to be tested equal to zero.
#' @param norm name of the normalization method used to compute the
#' scaling factors to use in the differential abundance analysis. If \code{norm}
#' is equal to "ratio", "poscounts", or "iterate" the normalization factors are
#' automatically transformed into scaling factors.
#' @param weights an optional numeric matrix giving observational weights.
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
#' # Calculate the TMM scaling factors
#' ps_NF <- norm_edgeR(object = ps, method = "TMM")
#' # The phyloseq object now contains the scaling factors:
#' scaleFacts <- phyloseq::sample_data(ps_NF)[, "NF.TMM"]
#' head(scaleFacts)
#'
#' # Differential abundance
#' DA_edgeR(object = ps_NF, pseudo_count = FALSE, group_name = "group",
#'          design = ~ group, coef = 2, robust = FALSE, norm = "TMM")

DA_edgeR <- function(object, pseudo_count = FALSE, group_name = NULL,
    design = NULL, robust = FALSE, coef = 2, norm = c( "TMM", "TMMwsp", "RLE",
    "upperquartile", "posupperquartile", "none", "ratio", "poscounts",
    "iterate", "TSS", "CSSmedian", "CSSdefault"), weights) {
    # Check the orientation
    if (!phyloseq::taxa_are_rows(object)) {
        object <- t(object)
    }
    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)
    # Name building
    name <- "edgeR"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count) {
        message("Adding a pseudo count... \n")
        counts <- counts + 1
        name <- paste(name, ".pseudo", sep = "")
    }
    if (length(norm) > 1) {
        stop("Please choose one normalization for this istance of differential
            abundance analysis.")
    }
    if (norm == "TSS") {
        NFs <- 1
    } else {
        # Check if the column with the normalization factors is present
        NF.col <- paste("NF", norm, sep = ".")
        if (!any(colnames(metadata) == NF.col)) {
            stop(paste0("Can't find the ", NF.col, " column in your object. Be
            sure to add the normalization factors column in your object first."))
        }
        NFs <- unlist(metadata[, NF.col])
        # DESeq2 NFs are size factors. To obtain normalized counts
        # we need to divide the raw counts by the NFs
        # DESeq2 NFs are supplied -> make them scaling factors!
        if (is.element(norm, c("ratio", "poscounts", "iterate"))) {
            NFs <- NFs / colSums(counts)
        }
        NFs <- NFs / exp(mean(log(NFs)))
    } # Make NFs multiply to 1
    name <- paste(name, ".", norm, sep = "")
    message(paste0("Differential abundance on ", norm, " normalized data"))
    group <- unlist(metadata[, group_name])
    dge <- edgeR::DGEList(counts = counts, norm.factors = NFs, group = group)
    if (missing(weights)) {
        message("Estimating Differential Abundance without weighting")
    } else {
        if(is.null(weights)){
            message("Estimating Differential Abundance without weighting")
        } else {
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
        message(paste0("Estimating robust dispersions"))
        dge <- edgeR::estimateGLMRobustDisp(y = dge, design = design)
        name <- paste(name, ".robust", sep = "")
    }
    message(paste0("Extracting results"))
    dispEsts <- edgeR::getDispersion(dge)
    glmFit <- edgeR::glmQLFit(
        y = dge, dispersion = dispEsts,
        robust = robust, design = design
    )
    glmRes <- edgeR::glmQLFTest(glmFit, coef = coef)
    message(paste0(
        "Extracting results for ", colnames(coef(glmRes))[coef],
        " coefficient"
    ))
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
#' # Set a specific set of normalization for edgeR (even of other packages!)
#' setNorm_edgeR <- set_edgeR(group_name = "group", design = ~ group, coef = 2,
#'     norm = c("TMM", "poscounts"))
#'
#' # Set many possible combinations of parameters for edgeR
#' all_edgeR <- set_edgeR(pseudo_count = c(TRUE, FALSE), group_name = "group",
#'     design = ~ group, robust = c(TRUE, FALSE), coef = 2,
#'     weights_logical = c(TRUE,FALSE))

set_edgeR <- function(pseudo_count = FALSE, group_name = NULL,
    design = NULL, robust = FALSE, coef = 2, norm = c("TMM", "TMMwsp", "RLE",
    "upperquartile", "posupperquartile", "none"), weights_logical = FALSE,
    expand = TRUE) {
    method <- "DA_edgeR"
    if (!is.logical(pseudo_count) | !is.logical(robust) |
        !is.logical(weights_logical)) {
        stop("'pseudo_count', 'robust', and 'weights_logical' must be logical.")
    }
    if (is.null(group_name) | is.null(design) | is.null(coef)) {
        stop("'group_name', 'design', and 'coef' are required.")
    }
    if (!is.character(group_name)){
        stop("'group_name' should be a character.")
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
            robust = robust, norm = norm, weights = weights_logical,
            stringsAsFactors = FALSE
        )
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, pseudo_count = pseudo_count,
            robust = robust, norm = norm, weights = weights_logical)
    }
    # data.frame to list
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        x <- append(x = x, values = list("group_name" = group_name,
                                         "design" = deparse(design),
                                         "coef" = coef), after = 2)
    })
    names(out) <- paste0(method, ".", 1:length(out))
    return(out)
}
