#' @title DA_MAST
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data
#' @importFrom MAST zlm summary FromMatrix
#' @importFrom SummarizedExperiment colData assay
#' @importFrom stats4 coef
#' @export
#' @description
#' Fast run for MAST differential abundance detection method.
#'
#' @inheritParams DA_DESeq2
#' @param rescale Rescale count data, per million if 'default', or per median
#' library size if 'median' ('median' is suggested for metagenomics data).
#' @param design The model for the count distribution. Can be the variable name,
#' or a character similar to "~ 1 + group", or a formula, or a `model.matrix`
#' object.
#' @param coefficient The coefficient of interest as a single word formed by 
#' the variable name and the non reference level. (e.g.: 'ConditionDisease' if 
#' the reference level for the variable 'Condition' is 'control').
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
#'                          
#' # Differential abundance
#' DA_MAST(object = ps, pseudo_count = FALSE, rescale = "median",
#'     design = ~ group, coefficient = "groupB")

DA_MAST <- function(object, assay_name = "counts", pseudo_count = FALSE,
    rescale = c("median", "default"), design, coefficient = NULL, 
    verbose = TRUE){
    counts_and_metadata <- get_counts_metadata(object, assay_name = assay_name)
    counts <- counts_and_metadata[[1]]
    metadata <- counts_and_metadata[[2]]
    is_phyloseq <- counts_and_metadata[[3]]
    # Name building
    name <- "MAST"
    method <- "DA_MAST"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        if(verbose)
            message("Adding a pseudo count...")
        counts <- counts + 1
        name <- paste(name,".pseudo",sep = "")}
    # Check the assay
    if (!is_phyloseq){
        if(verbose)
            message("Using the ", assay_name, " assay.")
        name <- paste(name, ".", assay_name, sep = "")
    } 
    if(length(rescale) > 1 | !is.element(rescale,c("default","median"))){
        stop(method, "\n", 
            "rescale: please choose between 'default' or 'median' for the", 
            " rescale parameter. 'median' is suggested for metagenomics data.")
    } else if(rescale == "median"){
        if(verbose)
            message("per ",rescale, "-lib.size rescaled data")
        tpm <- counts * median(colSums(counts)) / colSums(counts)
    } else {
        if(verbose)
            message(rescale," (per million) rescaled data")
        tpm <- counts * 10^6 / colSums(counts)
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
    ngeneson <- apply(counts, 2, function(x) mean(x>0))
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
    if(!is.element(coefficient, colnames(stats4::coef(fit, "C"))))
        stop(method, "\n", 
            "coefficient: please supply the coefficient of interest as a",
            " single word formed by the variable name and the non reference", 
            " level. (e.g.: 'ConditionDisease' if the reference level for the",
            " variable 'Condition' is 'control')")
    if(verbose){
        summaryDt <- data.frame(MAST::summary(fit, doLRT = coefficient)[[
            "datatable"]])
    } else {
        summaryDt <- suppressMessages(
            data.frame(MAST::summary(fit, doLRT = coefficient)[[
                "datatable"]]))
    }
    contrast <- component <- NULL
    fcHurdle <- merge(x = summaryDt[summaryDt[,"contrast"] == coefficient &
        summaryDt[,"component"] == 'H', c("primerid", "Pr..Chisq.")], y =
        summaryDt[summaryDt[,"contrast"] == coefficient & summaryDt[,
        "component"] == "logFC", c("primerid", "coef", "ci.hi", "ci.lo")],
        by = "primerid")
    statInfo <- data.frame(logFC = fcHurdle[, "coef"], logFC.lo = fcHurdle[,
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
#' # Set many possible combinations of parameters for MAST
#' all_MAST <- set_MAST(pseudo_count = c(TRUE, FALSE), rescale = c("median",
#'     "default"), design = ~ group, coefficient = "groupB")

set_MAST <- function(assay_name = "counts", pseudo_count = FALSE, 
    rescale = c("median", "default"), design = NULL, coefficient = NULL, 
    expand = TRUE) {
    method <- "DA_MAST"
    if (is.null(assay_name)) {
        stop(method, "\n", "'assay_name' is required (default = 'counts').")
    }
    if (!is.logical(pseudo_count)) {
        stop(method, "\n", "'pseudo_count' must be logical.")
    }
    if (is.null(design) | is.null(coefficient)) {
        stop(method, "\n", "'design' and 'coefficient' are required.")
    }
    if(sum(!is.element(rescale, c("median","default"))) > 0){
        stop(method, "\n", 
            "rescale: please choose rescale between 'median' and/or",
            " 'default'.")
    }
    if (expand) {
        parameters <- expand.grid(method = method, assay_name = assay_name,
            pseudo_count = pseudo_count, rescale = rescale, 
            stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, assay_name = assay_name,
            pseudo_count = pseudo_count, rescale = rescale)
    }
    # data.frame to list
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        x <- append(x = x, values = list("design" = design,
            "coefficient" = coefficient), after = 3)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}
