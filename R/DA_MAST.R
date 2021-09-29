#' @title DA_MAST
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data
#' @importFrom MAST zlm summary FromMatrix
#' @importFrom SummarizedExperiment colData assay
#' @export
#' @description
#' Fast run for MAST differential abundance detection method.
#'
#' @inheritParams DA_corncob
#' @param rescale Rescale count data, per million if 'default', or per median
#' library size if 'median' ('median' is suggested for metagenomics data).
#' @param design The model for the count distribution. Can be the variable name,
#' or a character similar to "~ 1 + group", or a formula, or a `model.matrix`
#' object.
#'
#' @return A list object containing the matrix of p-values `pValMat`, the matrix
#' of summary statistics for each tag `statInfo`, and a suggested `name` of the
#' final object considering the parameters passed to the function.
#'
#' @seealso \code{\link[MAST]{zlm}} for the Truncated Gaussian Hurdle model
#' estimation.
#'
#' @examples
#' set.seed(1)
#' # Create a very simple phyloseq object
#' counts <- matrix(rnbinom(n = 60, size = 3, prob = 0.5), nrow = 10, ncol = 6)
#' metadata <- data.frame("Sample" = c("S1", "S2", "S3", "S4", "S5", "S6"),
#'                        "group" = as.factor(c("A", "A", "A", "B", "B", "B")))
#' ps <- phyloseq::phyloseq(phyloseq::otu_table(counts, taxa_are_rows = TRUE),
#'                          phyloseq::sample_data(metadata))
#' # No use of scaling factors
#' ps_NF <- norm_edgeR(object = ps, method = "none")
#' # The phyloseq object now contains the scaling factors:
#' scaleFacts <- phyloseq::sample_data(ps_NF)[, "NF.none"]
#' head(scaleFacts)
#' # Differential abundance
#' DA_MAST(object = ps_NF, pseudo_count = FALSE, rescale = "median",
#'     design = ~ group, norm = "none", coefficient = "groupB")

DA_MAST <- function(object, pseudo_count = FALSE,
    rescale = c("median", "default"), design, coefficient = NULL,
    norm = c("TMM", "TMMwsp", "RLE",
        "upperquartile", "posupperquartile", "none", "ratio", "poscounts",
        "iterate", "TSS", "CSSmedian", "CSSdefault"),
    verbose = TRUE){
    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- data.frame(phyloseq::sample_data(object))
    # Name building
    name <- "MAST"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        if(verbose)
            message("Adding a pseudo count...")
        counts <- counts + 1
        name <- paste(name,".pseudo",sep = "")}
    if(length(norm) > 1)
        stop("Please choose one normalization for this istance of differential",
            " abundance analysis.")
    NF.col <- paste("NF", norm, sep = ".")
    # Check if the column with the normalization factors is present
    if(!any(colnames(metadata) == NF.col))
        stop("Can't find the ", NF.col," column in your object. Make sure to",
            " add the normalization factors column in your object first.")
    name <- paste(name, ".", norm, sep = "")
    NFs = unlist(metadata[, NF.col])
    # Check if the NFs are scaling factors. If so, make them norm. factors
    if(is.element(norm, c("TMM", "TMMwsp", "RLE", "upperquartile",
                          "posupperquartile", "CSSmedian", "CSSdefault", "TSS")))
        NFs <- NFs * colSums(counts)
    NFs <- NFs/exp(mean(log(NFs)))
    norm_counts <- round(counts %*% diag(1/NFs), digits = 0)
    colnames(norm_counts) <- colnames(counts)
    if(verbose)
        message("Differential abundance on ", norm, " normalized data")
    if(length(rescale) > 1 | !is.element(rescale,c("default","median"))){
        stop("Please choose between 'default' or 'median' for the rescale",
            " parameter. 'median' is suggested for metagenomics data.")
    } else if(rescale == "median"){
        if(verbose)
            message("per ",rescale, "-lib.size rescaled data")
        tpm <- norm_counts * median(colSums(norm_counts)) / colSums(norm_counts)
    } else {
        if(verbose)
            message(rescale," (per million) rescaled data")
        tpm <- norm_counts * 10^6 / colSums(norm_counts)
    }
    name <- paste(name, ".", rescale, sep = "")
    tpm <- log2(tpm + 1)
    if(verbose){
        sca <- MAST::FromMatrix(tpm, cData = metadata)
    } else {
        sca <- suppressMessages(MAST::FromMatrix(tpm, cData = metadata))
    }
    # here, we keep all OTUs so that we can fairly compare MAST and the other
    # methods. So, no adaptive thresholding or filtering by gene expression.
    SummarizedExperiment::assays(sca) <- list(tpm =
        SummarizedExperiment::assay(sca))
    ngeneson <- apply(norm_counts, 2, function(x) mean(x>0))
    metadata[, "cngeneson"] <- ngeneson - mean(ngeneson)
    SummarizedExperiment::colData(sca)[, "cngeneson"] <- metadata[, "cngeneson"]
    if(is.character(design)){
        if(grepl("~", design))
            design <- as.formula(design)
        else design <- as.formula(paste0("~", design))}
    if(is(design, "formula"))
        design <- as.formula(paste0(paste0(design, collapse = " "),
            " + cngeneson"))
    ## differential expression
    if(verbose){
        fit <- MAST::zlm(formula = design, sca = sca, method = "bayesglm",
            ebayes = TRUE)
    } else {
        fit <- suppressMessages(
            MAST::zlm(formula = design, sca = sca, method = "bayesglm",
                ebayes = TRUE))
    }
    if(!is.element(coefficient, colnames(coef(fit, "C"))))
        stop("Please supply the coefficient of interest as a single word",
            " formed by the variable name and the non reference level. (e.g.:",
            " 'ConditionDisease' if the reference level for the variable",
            " 'Condition' is 'control')")
    if(verbose){
        summaryDt = data.frame(MAST::summary(fit, doLRT = coefficient)[[
            "datatable"]])
    } else {
        summaryDt = suppressMessages(
            data.frame(MAST::summary(fit, doLRT = coefficient)[[
                "datatable"]]))
    }
    contrast <- component <- NULL
    fcHurdle <- merge(x = summaryDt[summaryDt[,"contrast"] == coefficient &
        summaryDt[,"component"] == 'H', c("primerid", "Pr..Chisq.")], y =
        summaryDt[summaryDt[,"contrast"] == coefficient & summaryDt[,
        "component"] == "logFC", c("primerid", "coef", "ci.hi", "ci.lo")],
        by = "primerid")
    statInfo = data.frame(logFC = fcHurdle[, "coef"], logFC.lo = fcHurdle[,
        "ci.lo"], logFC.hi = fcHurdle[, "ci.hi"], rawP = fcHurdle[,
        "Pr..Chisq."], adjP = stats::p.adjust(fcHurdle[, "Pr..Chisq."], 'BH'))
    rownames(statInfo) <- fcHurdle[, "primerid"]
    pValMat <- statInfo[, c("rawP", "adjP")]
    statInfo <- statInfo[rownames(counts),]
    pValMat <- pValMat[rownames(counts),]
    list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name)
}# END - function: MAST

#' @title set_MAST
#'
#' @export
#' @description
#' Set the parameters for MAST differential abundance detection method.
#'
#' @inheritParams DA_MAST
#' @param expand logical, if TRUE create all combinations of input parameters
#' (default \code{expand = TRUE})
#'
#' @return A named list containing the set of parameters for \code{DA_MAST}
#' method.
#'
#' @seealso \code{\link{DA_MAST}}
#'
#' @examples
#' # Set some basic combinations of parameters for MAST
#' base_MAST <- set_MAST(design = ~ group, coefficient = "groupB")
#' # Set a specific set of normalization for MAST (even of other packages!)
#' setNorm_MAST <- set_MAST(design = ~ group, coefficient = "groupB",
#'     norm = c("TSS", "poscounts", "TMM"))
#' # Set many possible combinations of parameters for MAST
#' all_MAST <- set_MAST(pseudo_count = c(TRUE, FALSE), rescale = c("median",
#'     "default"), design = ~ group, coefficient = "groupB", norm = c("TSS",
#'     "poscounts"))

set_MAST <- function(pseudo_count = FALSE, rescale = c("median", "default"),
    design = NULL, coefficient = NULL, norm = "TSS", expand = TRUE) {
    method <- "DA_MAST"
    if (!is.logical(pseudo_count)) {
        stop("'pseudo_count' must be logical.")
    }
    if (is.null(design) | is.null(coefficient)) {
        stop("'design' and 'coefficient' are required.")
    }
    if (sum(!is.element(norm, c("TSS", "none"))) > 0) {
        warning("One or more elements into 'norm' are not native to MAST")
    }
    if(sum(!is.element(rescale, c("median","default"))) > 0){
        stop("Please choose rescale between 'median' and/or 'default'.")
    }
    if (expand) {
        parameters <- expand.grid(method = method, pseudo_count = pseudo_count,
                                  rescale = rescale, norm = norm,
                                  stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, pseudo_count = pseudo_count,
                                 rescale = rescale, norm = norm)
    }
    # data.frame to list
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        x <- append(x = x, values = list("design" = design,
                                         "coefficient" = coefficient), after = 2)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}
