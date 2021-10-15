#' @title DA_corncob
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data
#' @importFrom corncob differentialTest
#' @importFrom stats4 coef
#' @export
#' @description
#' Fast run for corncob differential abundance detection method.
#'
#' @inheritParams DA_DESeq2
#' @param formula an object of class \code{formula} without the response: a
#' symbolic description of the model to be fitted to the abundance.
#' @param phi.formula an object of class \code{formula} without the response: a
#' symbolic description of the model to be fitted to the dispersion.
#' @param coefficient The coefficient of interest as a single word formed by the
#' variable name and the non reference level. (e.g.: 'ConditionDisease' if the
#' reference level for the variable 'Condition' is 'control').
#' @param test Character. Hypothesis testing procedure to use. One of
#' \code{"Wald"} or \code{"LRT"} (likelihood ratio test).
#' @param boot Boolean. Defaults to \code{FALSE}. Indicator of whether or not to
#' use parametric bootstrap algorithm. (See \code{\link[corncob]{pbWald}} and
#' \code{\link[corncob]{pbLRT}}).
#' @inheritParams corncob::differentialTest
#'
#' @return A list object containing the matrix of p-values `pValMat`, the matrix
#' of summary statistics for each tag `statInfo`, and a suggested `name` of the
#' final object considering the parameters passed to the function.
#'
#' @seealso \code{\link[corncob]{bbdml}} and
#' \code{\link[corncob]{differentialTest}} for differential abundance and
#' differential variance evaluation.
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
#' DA_corncob(object = ps_NF, formula = ~ group, phi.formula = ~ group,
#'     formula_null = ~ 1, phi.formula_null = ~ group, coefficient = "groupB",
#'     norm = "none", test = "Wald")

DA_corncob <- function(object, pseudo_count = FALSE, formula, phi.formula,
    formula_null, phi.formula_null, test, boot = FALSE, coefficient = NULL,
    norm = c("TMM", "TMMwsp", "RLE", "upperquartile", "posupperquartile",
    "none", "ratio", "poscounts", "iterate", "TSS", "CSSmedian", "CSSdefault"),
    verbose = TRUE){
    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)
    # Name building
    name <- "corncob"
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
    # Check if the NFs are scaling factors. If so, make them norm.factors
    if(is.element(norm, c("TMM", "TMMwsp", "RLE", "upperquartile",
        "posupperquartile", "CSSmedian", "CSSdefault", "TSS")))
        NFs <- NFs * colSums(counts)
    NFs <- NFs/exp(mean(log(NFs)))
    norm_counts <- round(counts %*% diag(1/NFs), digits = 0)
    colnames(norm_counts) <- colnames(counts)
    if(verbose)
        message("Differential abundance on ", norm," normalized data")
    if(missing(test))
        stop("Please supply the test to perform, 'Wald' or 'LRT'.")
    else name <- paste(name, ".", test, sep = "")
    if(boot)
        name <- paste(name, ".", "boot", sep = "")
    ## differential expression
    requireNamespace("phyloseq")
    fit <- corncob::differentialTest(formula = formula, phi.formula =
        phi.formula, formula_null = formula_null, phi.formula_null =
        phi.formula_null, data = norm_counts, sample_data = metadata, test =
        test, boot = boot)
    # differentialTest's output has a complex structure:
    # extraction of summary table for each estimated model
    # mu.(Intercept), mu.condition, and phi.(Intercept), phi.condition
    # only the mu.condition is of interest.
    if(is.null(coefficient) | !is.element(coefficient,fit[["restrictions_DA"]]))
        stop("Please supply the coefficient of interest as a single word",
            " formed by the variable name and the non reference level. (e.g.:",
            " 'ConditionDisease' if the reference level for the variable",
            " 'Condition' is 'control')")
    statInfo <- plyr::ldply(.data = fit[["all_models"]], .fun = function(model){
        stats4::coef(model)[paste0("mu.",coefficient),]})
    if(length(fit[["restrictions_DA"]])>0)
        if(verbose)
            message("Differential abundance across",
                paste0(fit[["restrictions_DA"]], collapse = " and "))
    if(length(fit[["restrictions_DV"]])>0)
        if(verbose)
            message("Differential variability across",
                paste0(fit[["restrictions_DV"]], collapse = " and "))
    pValMat <- data.frame("rawP" = fit[["p"]], "adjP" = fit[["p_fdr"]])
    rownames(statInfo) <- rownames(pValMat) <- names(fit[["p"]])
    list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name)
}# END - function: corncob

#' @title set_corncob
#'
#' @export
#' @description
#' Set the parameters for corncob differential abundance detection method.
#'
#' @inheritParams DA_corncob
#' @param expand logical, if TRUE create all combinations of input parameters
#' (default \code{expand = TRUE})
#'
#' @return A named list containing the set of parameters for \code{DA_corncob}
#' method.
#'
#' @seealso \code{\link{DA_corncob}}
#'
#' @examples
#' # Set some basic combinations of parameters for corncob
#' base_corncob <- set_corncob(formula = ~ group, phi.formula = ~ group,
#'     formula_null = ~ 1, phi.formula_null = ~ group, coefficient = "groupB")
#' # Set a specific set of normalization for corncob (even of other packages!)
#' setNorm_corncob <- set_corncob(formula = ~ group, phi.formula = ~ group,
#'     formula_null = ~ 1, phi.formula_null = ~ group, coefficient = "groupB",
#'     norm = c("TMM", "TSS", "poscounts"))
#' # Set many possible combinations of parameters for corncob
#' all_corncob <- set_corncob(pseudo_count = c(TRUE, FALSE), formula = ~ group,
#'     phi.formula = ~ group, formula_null = ~ 1, phi.formula_null = ~ group,
#'     coefficient = "groupB", boot = c(TRUE, FALSE))

set_corncob <- function(pseudo_count = FALSE, formula = NULL, phi.formula = NULL,
    formula_null = NULL, phi.formula_null = NULL, test = c("Wald", "LRT"),
    boot = FALSE, coefficient = NULL, norm = "TSS", expand = TRUE) {
    method <- "DA_corncob"
    if (!is.logical(pseudo_count) | !is.logical(boot)) {
        stop("'pseudo_count' and 'boot' must be logical.")
    }
    if (is.null(coefficient)) {
        stop("'coefficient' is required.")
    }
    if (is.null(formula) | is.null(formula_null)) {
        stop("Please specify 'formula' and 'formula_null'.")
    }
    if (is.null(phi.formula) | is.null(phi.formula_null)) {
        stop("Please specify 'phi.formula' and 'phi.formula_null'.")
    }
    if (sum(!is.element(norm, c("TSS", "none"))) > 0) {
        warning("One or more elements into 'norm' are not native to corncob.")
    }
    if(sum(!is.element(test, c("Wald","LRT"))) > 0){
        stop("Please choose the test between 'Wald' and 'LRT'.")
    }
    if (expand) {
        parameters <- expand.grid(method = method, pseudo_count = pseudo_count,
                                  test = test, boot = boot, norm = norm,
                                  stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, pseudo_count = pseudo_count,
                                 test = test, boot = boot, norm = norm)
    }
    # data.frame to list
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        x <- append(x = x, values = list("formula" = formula,
                                         "formula_null" = formula_null,
                                         "phi.formula" = phi.formula,
                                         "phi.formula_null" = phi.formula_null,
                                         "coefficient" = coefficient), after = 2)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}
