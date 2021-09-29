#' @title DA_ALDEx2
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data
#' @importFrom ALDEx2 aldex
#' @export
#' @description
#' Fast run for the ALDEx2's differential abundance detection method.
#'
#' @inheritParams DA_DESeq2
#' @param mc.samples An integer. The number of Monte Carlo samples to use when
#' estimating the underlying distributions. Since we are estimating central
#' tendencies, 128 is usually sufficient.
#' @param denom A character string. Indicates which features to retain as the
#' denominator for the Geometric Mean calculation. Using "iqlr" accounts for
#' data with systematic variation and centers the features on the set features
#' that have variance that is between the lower and upper quartile of variance.
#' Using "zero" is a more extreme case where there are many non-zero features in
#' one condition but many zeros in another. In this case the geometric mean of
#' each group is calculated using the set of per-group non-zero features.
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
#' # No use of scaling factors
#' ps_NF <- norm_edgeR(object = ps, method = "none")
#' # The phyloseq object now contains the scaling factors:
#' scaleFacts <- phyloseq::sample_data(ps_NF)[, "NF.none"]
#' head(scaleFacts)
#' # Differential abundance
#' DA_ALDEx2(ps_NF, conditions = "group", test = "t", denom = "iqlr",
#'     norm = "none")

DA_ALDEx2 <- function(object, pseudo_count = FALSE, conditions = NULL,
    mc.samples = 128, test = c("t","wilcox"), denom = "iqlr", norm = c("TMM",
    "TMMwsp", "RLE", "upperquartile", "posupperquartile", "none", "ratio",
    "poscounts", "iterate", "TSS", "CSSmedian", "CSSdefault"), verbose = TRUE){
    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)
    # Name building
    name <- "ALDEx2"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        if(verbose)
            message("Adding a pseudo count... \n")
        counts <- counts + 1
        name <- paste(name,".pseudo",sep = "")}
    if(length(norm) > 1)
        stop("Please choose one normalization for this istance of differential",
            " abundance analysis.")
    NF.col <- paste("NF", norm, sep = ".")
    # Check if the column with the normalization factors is present
    if(!any(colnames(metadata) == NF.col)){
        stop("Can't find the ", NF.col," column in your object.",
             " Make sure to add the normalization factors column in your",
             " object first.")}
    name <- paste(name, ".", norm, sep = "")
    NFs = unlist(metadata[, NF.col])
    # Check if the NFs are scaling factors. If so, make them norm. factors
    if(is.element(norm, c("TMM", "TMMwsp", "RLE", "upperquartile",
        "posupperquartile", "CSSmedian", "CSSdefault", "TSS")))
        NFs <- NFs * colSums(counts)
    NFs <- NFs/exp(mean(log(NFs)))
    norm_counts <- round(counts %*% diag(1/NFs), digits = 0)
    colnames(norm_counts) <- colnames(counts)
    if(is.null(conditions))
        stop("Please supply the name of the variable of interest or the",
            " entire character vector.")
    else if(length(conditions) == 1)
        conditions = unlist(metadata[, conditions])
    name <- paste(name, ".", denom, sep = "")
    if(!is.element(test, c("t","wilcox")) | length(test) != 1)
        stop("Please choose between p-values produced by Welch t-test (t) or",
             " by the Wilcoxon test (wilcox).")
    name <- paste(name, ".", test, sep = "")
    if(verbose){
        statInfo <- ALDEx2::aldex(reads = norm_counts, conditions = conditions,
            mc.samples = mc.samples, test = test, effect = TRUE,
            include.sample.summary = FALSE, denom = denom, verbose = verbose)
    } else {
        statInfo <- suppressMessages(ALDEx2::aldex(reads = norm_counts,
            conditions = conditions, mc.samples = mc.samples, test = test,
            effect = TRUE, include.sample.summary = FALSE, denom = denom,
            verbose = verbose))
    }
    if(test == "t")
        pValMat <- data.frame(statInfo[, c("we.ep", "we.eBH")])
    else pValMat <- data.frame(statInfo[, c("wi.ep", "wi.eBH")])
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
#' base_ALDEx2 <- set_ALDEx2(conditions = "group")
#' # Set a specific set of normalization for ALDEx2 (even of other
#' # packages!)
#' setNorm_ALDEx2 <- set_ALDEx2(conditions = "group", norm = c("TSS", "TMM"))
#' # Set many possible combinations of parameters for ALDEx2
#' all_ALDEx2 <- set_ALDEx2(conditions = "group", denom = c("iqlr", "zero"),
#'     test = c("t", "wilcox"))

set_ALDEx2 <- function(pseudo_count = FALSE, conditions = NULL,
    mc.samples = 128, test = "t", denom = "iqlr", norm = "TSS", expand = TRUE) {
    method <- "DA_ALDEx2"
    if (!is.logical(pseudo_count)) {
        stop("'pseudo_count' must be logical.")
    }
    if (is.null(conditions)) {
        stop("'conditions' is required.")
    }
    if (sum(!is.element(norm, c("TSS", "none"))) > 0) {
        warning("One or more elements into 'norm' are not native to ALDEx2.")
    }
    if(sum(!is.element(test, c("t","wilcox"))) > 0){
        stop("Please choose between p-values produced by Welch t-test 't'",
            "or by the Wilcoxon test 'wilcox'.")
    }
    if(sum(!is.element(denom, c("iqlr","zero"))) > 0){
        stop("Please choose denom between 'iqlr' and 'zero'.")
    }
    if (expand) {
        parameters <- expand.grid(method = method, pseudo_count = pseudo_count,
            mc.samples = mc.samples, test = test, denom = denom, norm = norm,
            stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, pseudo_count = pseudo_count,
            mc.samples = mc.samples, test = test, denom = denom, norm = norm)
    }
    # data.frame to list
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        x <- append(x = x, values = list("conditions" = conditions), after = 2)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}
