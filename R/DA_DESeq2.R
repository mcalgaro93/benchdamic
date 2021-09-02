#' @title DA_DESeq2
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data phyloseq_to_deseq2
#' @importFrom DESeq2 sizeFactors DESeq dispersions results
#' @importFrom SummarizedExperiment assays
#' @export
#' @description
#' Fast run for DESeq2 differential abundance detection method.
#'
#' @inheritParams phyloseq::phyloseq_to_deseq2
#' @inheritParams DA_edgeR
#' @param contrast character vector with exactly three elements: the name of a
#' factor in the design formula, the name of the numerator level for the fold
#' change, and the name of the denominator level for the fold change.
#' @param alpha the significance cutoff used for optimizing the independent
#' filtering (by default 0.05). If the adjusted p-value cutoff (FDR) will be a
#' value other than 0.05, alpha should be set to that value.
#' @param norm name of the normalization method used to compute the
#' normalization factors to use in the differential abundance analysis. If
#' \code{norm} is equal to "TMM", "TMMwsp", "RLE", "upperquartile",
#' "posupperquartile", "CSSmedian", "CSSdefault", "TSS" the scaling factors are
#' automatically transformed into normalization factors.
#'
#' @return A list object containing the matrix of p-values `pValMat`,
#' the dispersion estimates `dispEsts`, the matrix of summary statistics for
#' each tag `statInfo`, and a suggested `name` of the final object considering
#' the parameters passed to the function.
#'
#' @seealso \code{\link[phyloseq]{phyloseq_to_deseq2}} for phyloseq to DESeq2
#' object conversion, \code{\link[DESeq2]{DESeq}} and
#' \code{\link[DESeq2]{results}} for the differential abundance method.
#'
#' @examples
#' set.seed(1)
#' # Create a very simple phyloseq object
#' counts <- matrix(rnbinom(n = 60, size = 3, prob = 0.5), nrow = 10, ncol = 6)
#' metadata <- data.frame("Sample" = c("S1", "S2", "S3", "S4", "S5", "S6"),
#'                        "group" = as.factor(c("A", "A", "A", "B", "B", "B")))
#' ps <- phyloseq::phyloseq(phyloseq::otu_table(counts, taxa_are_rows = TRUE),
#'                          phyloseq::sample_data(metadata))
#' # Calculate the poscounts normalization factors
#' ps_NF <- norm_DESeq2(object = ps, method = "poscounts")
#' # The phyloseq object now contains the normalization factors:
#' scaleFacts <- phyloseq::sample_data(ps_NF)[, "NF.poscounts"]
#' head(scaleFacts)
#' # Differential abundance
#' DA_DESeq2(object = ps_NF, pseudo_count = FALSE, design = ~ group, contrast =
#'               c("group", "B", "A"), norm = "poscounts")

DA_DESeq2 <- function(object, pseudo_count = FALSE, design = NULL, contrast =
    NULL, alpha = 0.05, norm = c("TMM", "TMMwsp", "RLE", "upperquartile",
    "posupperquartile", "none", "ratio", "poscounts", "iterate", "TSS",
    "CSSmedian", "CSSdefault"), weights){
    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)
    # Name building
    name <- "DESeq2"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        message("Adding a pseudo count... \n")
        counts <- counts + 1
        phyloseq::otu_table(object) <- counts
        name <- paste(name,".pseudo",sep = "")}
    if (is.character(design)){
        design <- as.formula(design)
    }
    dds <- phyloseq::phyloseq_to_deseq2(object, design = design)
    if(length(norm) > 1)
        stop("Please choose one normalization for this istance of differential
        abundance analysis.")
    # Check if the column with the normalization factors is present
    NF.col <- paste("NF", norm, sep = ".")
    if(!any(colnames(metadata) == NF.col))
        stop(paste0("Can't find the ", NF.col," column in your object. Make sure
            to add the normalization factors column in your object first."))
    NFs = unlist(metadata[,NF.col])
    # edgeR, TSS, and CSS NFs supplied -> make them normalization factors!
    if(is.element(norm, c("TMM", "TMMwsp", "RLE", "upperquartile",
                          "posupperquartile", "CSSmedian", "CSSdefault", "TSS"))){
        NFs <- NFs*colSums(counts)}
    DESeq2::sizeFactors(dds) = NFs/exp(mean(log(NFs))) # Make NFs multiply to 1
    name <- paste(name, ".", norm, sep = "")
    message(paste0("Differential abundance on ", norm," normalized data"))
    if(missing(weights))
        message("Estimating Differential Abundance without weighting")
    else {
        if(is.null(weights)){
            message("Estimating Differential Abundance without weighting")
        } else {
            message("Estimating Differential Abundance with weights")
            weights[which(weights < 1e-6)] <- 1e-06
            SummarizedExperiment::assays(dds, withDimnames = FALSE)$weights <-
                weights
            name <- paste(name,".weighted",sep = "")
        }
    }
    ### Run DESeq
    ddsRes <- DESeq2::DESeq(object = dds, test = "LRT", reduced = ~ 1,
                            parallel = FALSE)
    dispEsts <- DESeq2::dispersions(ddsRes)
    if(missing(contrast) | (!is.character(contrast) & length(contrast) != 3))
        stop(paste0("Please supply a character vector with exactly three
            elements: the name of a factor in the design formula, the name of
            the numerator level for the fold change, and the name of the
            denominator level for the fold change."))
    else message(paste0("Extracting results for ", contrast[1]," variable, ",
                        contrast[2], " / ", contrast[3]))
    res <- DESeq2::results(ddsRes, alpha = alpha, contrast = contrast)
    statInfo <- as(res, "data.frame")
    pValMat <- statInfo[,c("pvalue", "padj")]
    colnames(pValMat) <- c("rawP", "adjP")
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "dispEsts" =
                    dispEsts, "name" = name))
}# END - function: DA_DESeq2

#' @title set_DESeq2
#'
#' @export
#' @description
#' Set the parameters for DESeq2 differential abundance detection method.
#'
#' @inheritParams DA_DESeq2
#' @param weights_logical logical vector, if TRUE a matrix of observational
#' weights will be used for differential abundance analysis (default
#' \code{weights_logical = FALSE}).
#' @param expand logical, if TRUE create all combinations of input parameters
#' (default \code{expand = TRUE}).
#'
#' @return A named list containing the set of parameters for \code{DA_DESeq2}
#' method.
#'
#' @seealso \code{\link{DA_DESeq2}}
#'
#' @examples
#' # Set some basic combinations of parameters for DESeq2
#' base_DESeq2 <- set_DESeq2(design = ~ group, contrast = c("group", "B", "A"))
#' # Set a specific set of normalization for DESeq2 (even of other packages!)
#' setNorm_DESeq2 <- set_DESeq2(design = ~ group, contrast =
#'     c("group", "B", "A"), norm = c("TMM", "poscounts"))
#' # Set many possible combinations of parameters for edgeR
#' all_DESeq2 <- set_DESeq2(pseudo_count = c(TRUE, FALSE), design = ~ group,
#'     contrast = c("group", "B", "A"), weights_logical = c(TRUE,FALSE))
set_DESeq2 <- function(pseudo_count = FALSE, design = NULL,
    contrast = NULL, alpha = 0.05, norm = c("ratio", "poscounts", "iterate"),
    weights_logical = FALSE, expand = TRUE) {
    method <- "DA_DESeq2"
    if (!is.logical(pseudo_count) | !is.logical(weights_logical)) {
        stop("'pseudo_count' and 'weights_logical' must be logical.")
    }
    if (is.null(design) | is.null(contrast)) {
        stop("'design' and 'contrast' must be specified.")
    }
    if (!is.character(design) & !is(design, "formula")){
        stop("'design' should be a character or a formula.")
    } else design <- as.formula(design)
    if (!is.character(contrast) & length(contrast) != 3){
        stop(paste0("Please supply a character vector with exactly three
            elements: the name of a factor in the design formula, the name of
            the numerator level for the fold change, and the name of the
            denominator level for the fold change."))
    }
    if (sum(!is.element(norm, c(
        "ratio", "poscounts", "iterate", "none", "TSS"
    ))) > 0) {
        warning("One or more elements into 'norm' are not native to DESeq2")
    }
    if (expand) {
        parameters <- expand.grid(method = method, pseudo_count = pseudo_count,
            alpha = alpha, norm = norm, weights = weights_logical,
            stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, pseudo_count = pseudo_count,
            alpha = alpha, norm = norm, weights = weights_logical)
    }
    # data.frame to list
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        x <- append(x = x, values = list("design" = deparse(design),
                                         "contrast" = contrast), after = 2)
    })
    names(out) <- paste0(method, ".", 1:length(out))
    return(out)
}
