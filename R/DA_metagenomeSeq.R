#' @title DA_metagenomeSeq
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data taxa_names
#' @importFrom phyloseq phyloseq_to_metagenomeSeq
#' @importFrom metagenomeSeq fitZig MRfulltable
#' @importFrom stats model.matrix
#' @export
#' @description
#' Fast run for the metagenomeSeq's differential abundance detection method.
#'
#' @inheritParams DA_edgeR
#' @param design The model for the count distribution. Can be the variable name,
#' or a character similar to "~ 1 + group", or a formula, or a `model.matrix`
#' object.
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
#' # Calculate the CSSdefault scaling factors
#' ps_NF <- norm_CSS(object = ps, method = "default")
#' # The phyloseq object now contains the scaling factors:
#' scaleFacts <- phyloseq::sample_data(ps_NF)[, "NF.CSSdefault"]
#' head(scaleFacts)
#' # Differential abundance
#' DA_metagenomeSeq(object = ps_NF, pseudo_count = FALSE, design = ~ group,
#'     coef = 2, norm = "CSSdefault")

DA_metagenomeSeq <- function(object, pseudo_count = FALSE, design = NULL,
    coef = 2, norm = c("TMM", "TMMwsp", "RLE", "upperquartile",
        "posupperquartile", "none", "ratio", "poscounts", "iterate", "TSS",
        "CSSmedian", "CSSdefault")){
    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)
    obj <- phyloseq::phyloseq_to_metagenomeSeq(object)
    # Name building
    name <- "metagenomeSeq"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        message("Adding a pseudo count... \n")
        counts <- counts + 1
        name <- paste(name,".pseudo",sep = "")}
    if(length(norm) > 1)
        stop("Please choose one normalization for this istance of differential
        abundance analysis.")
    NF.col <- paste("NF", norm, sep = ".")
    if(norm == "TSS"){
        phyloseq::sample_data(obj)[NF.col] <- NFs <- 1L
    } else {
        # Check if the column with the normalization factors is present
        if(!any(colnames(metadata) == NF.col))
            stop(paste0("Can't find the ", NF.col," column in your object. Make
            sure to add the normalization factors column in your object first."
            ))
        NFs = unlist(metadata[,NF.col])
        # DESeq2 NFs are supplied -> make them scaling factors!
        if(is.element(norm, c("ratio", "poscounts", "iterate")))
            NFs <- NFs/colSums(counts)
        NFs = NFs/exp(mean(log(NFs)))
    } # Make NFs multiply to 1
    name <- paste(name, ".", norm, sep = "")
    message(paste0("Differential abundance on ", norm," normalized data"))
    metagenomeSeq::normFactors(object = obj) <- NFs
    if(is.character(design)){
        if(grepl("~", design))
            design <- as.formula(design)
        else design <- as.formula(paste0("~", design))}
    if(is(design, "formula")){
        design <- as.formula(paste0(paste0(design,collapse = " "), " + ",
            NF.col))
        design <- stats::model.matrix(object = design,
            data = data.frame(metadata))
    }
    suppressWarnings(fit <- try(metagenomeSeq::fitZig(obj = obj, mod = design,
        verbose = FALSE, useCSSoffset = FALSE,
        control = metagenomeSeq::zigControl(maxit = 1000)), silent = TRUE))
    if(is(fit, "try-error")){
        res = matrix(NA, ncol = 2, nrow = nrow(counts))
        stop("Something went wrong during fitZig estimation.")
    } else {
        statInfo <- metagenomeSeq::MRcoefs(obj = fit, by = coef, number = nrow(
            counts))
        statInfo <- statInfo[phyloseq::taxa_names(object),]
        message(paste0("Extracting results for ", colnames(statInfo[coef]),
                       " coefficient"))}
    pValMat <- statInfo[, c("pvalues", "adjPvalues")]
    colnames(pValMat) = c("rawP", "adjP")
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
#' # Set some basic combinations of parameters for metagenomeSeq
#' base_mgs <- set_metagenomeSeq(design = ~ group, coef = 2)
#' # Set a specific set of normalization for metagenomeSeq (even of other
#' # packages!)
#' setNorm_mgs <- set_metagenomeSeq(design = ~ group, coef = 2,
#'     norm = c("CSSmedian", "TMM"))
#' # Set many possible combinations of parameters for metagenomeSeq
#' all_mgs <- set_metagenomeSeq(pseudo_count = c(TRUE, FALSE), design = ~ group,
#'     coef = 2, norm = c("CSSmedian", "CSSdefault", "TMM", "TSS"))

set_metagenomeSeq <- function(pseudo_count = FALSE, design = NULL, coef = 2,
                      norm = c("CSSmedian", "CSSdefault"), expand = TRUE) {
    method <- "DA_metagenomeSeq"
    if (!is.logical(pseudo_count)) {
        stop("'pseudo_count' must be logical.")
    }
    if (is.null(design) | is.null(coef)) {
        stop("'design', and 'coef' are required.")
    }
    if (!is.character(design) & !is(design, "formula")){
        stop("'design' should be a character or a formula.")
    } else design <- as.formula(design)
    if (sum(!is.element(norm, c("CSSmedian", "CSSdefault"))) > 0) {
        warning(paste("One or more elements into 'norm' are not native to
                       metagenomeSeq."))
    }
    if (expand) {
        parameters <- expand.grid(method = method, pseudo_count = pseudo_count,
            norm = norm, stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, pseudo_count = pseudo_count,
            norm = norm)
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

