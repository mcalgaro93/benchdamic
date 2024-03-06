#' @title DA_corncob
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data
#' @importFrom corncob differentialTest
#' @importFrom stats4 coef
#' @importFrom plyr is.formula
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
#'                          
#' # Differential abundance
#' DA_corncob(object = ps, formula = ~ group, phi.formula = ~ group,
#'     formula_null = ~ 1, phi.formula_null = ~ group, coefficient = "groupB",
#'     test = "Wald")

DA_corncob <- function(object, assay_name = "counts", pseudo_count = FALSE, 
    formula, phi.formula, formula_null, phi.formula_null, test, boot = FALSE, 
    coefficient = NULL, verbose = TRUE){
    counts_and_metadata <- get_counts_metadata(object, assay_name = assay_name)
    counts <- counts_and_metadata[[1]]
    metadata <- counts_and_metadata[[2]]
    is_phyloseq <- counts_and_metadata[[3]]
    # Name building
    name <- "corncob"
    method <- "DA_corncob"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        if(verbose)
            message("Adding a pseudo count...")
        counts <- counts + 1
        name <- paste(name,".pseudo",sep = "")
    }
    # Check the assay
    if (!is_phyloseq){
        if(verbose)
            message("Using the ", assay_name, " assay.")
        name <- paste(name, ".", assay_name, sep = "")
    } 
    if (!plyr::is.formula(formula) | !plyr::is.formula(formula_null)) {
        stop(method, "\n", 
            "Please specify 'formula' and 'formula_null' as formula objects.")
    }
    if (!plyr::is.formula(phi.formula) | !plyr::is.formula(phi.formula_null)) {
        stop(method, "\n", 
            "Please specify 'phi.formula' and 'phi.formula_null' as formula", 
            " objects.")
    }
    if(missing(test))
        stop(method, "\n", 
            "test: please supply the test to perform, 'Wald' or 'LRT'.")
    else name <- paste(name, ".", test, sep = "")
    if(boot)
        name <- paste(name, ".", "boot", sep = "")
    ## differential expression
    requireNamespace("phyloseq")
    fit <- corncob::differentialTest(formula = formula, phi.formula =
        phi.formula, formula_null = formula_null, phi.formula_null =
        phi.formula_null, data = counts, sample_data = metadata, test =
        test, boot = boot)
    # differentialTest's output has a complex structure:
    # extraction of summary table for each estimated model
    # mu.(Intercept), mu.condition, and phi.(Intercept), phi.condition
    # only the mu.condition is of interest.
    if(is.null(coefficient) | !is.element(coefficient,fit[["restrictions_DA"]]))
        stop(method, "\n", 
            "coefficient: please supply the coefficient of interest as a", 
            " single word formed by the variable name and the non reference",
            " level. (e.g.: 'ConditionDisease' if the reference level for the",
            " variable 'Condition' is 'control')")
    statInfo <- plyr::ldply(.data = fit[["all_models"]], .fun = function(model){
        stats4::coef(model)[paste0("mu.",coefficient),]})
    if(length(fit[["restrictions_DA"]])>0)
        if(verbose)
            message("Differential abundance across ",
                paste0(fit[["restrictions_DA"]], collapse = " and "))
    if(length(fit[["restrictions_DV"]])>0)
        if(verbose)
            message("Differential variability across ",
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
#' # Set many possible combinations of parameters for corncob
#' all_corncob <- set_corncob(pseudo_count = c(TRUE, FALSE), formula = ~ group,
#'     phi.formula = ~ group, formula_null = ~ 1, phi.formula_null = ~ group,
#'     coefficient = "groupB", boot = c(TRUE, FALSE))

set_corncob <- function(assay_name = "counts", pseudo_count = FALSE, 
    formula = NULL, phi.formula = NULL, formula_null = NULL, 
    phi.formula_null = NULL, test = c("Wald", "LRT"), boot = FALSE, 
    coefficient = NULL, expand = TRUE) {
    method <- "DA_corncob"
    if (is.null(assay_name)) {
        stop(method, "\n", "'assay_name' is required (default = 'counts').")
    }
    if (!is.logical(pseudo_count) | !is.logical(boot)) {
        stop(method, "\n", "'pseudo_count' and 'boot' must be logical.")
    }
    if (is.null(coefficient)) {
        stop(method, "\n", "'coefficient' is required.")
    }
    if (!plyr::is.formula(formula) | !plyr::is.formula(formula_null)) {
        stop(method, "\n", 
            "Please specify 'formula' and 'formula_null' as formula objects.")
    }
    if (!plyr::is.formula(phi.formula) | !plyr::is.formula(phi.formula_null)) {
        stop(method, "\n", 
            "Please specify 'phi.formula' and 'phi.formula_null' as formula", 
            " objects.")
    }
    if(sum(!is.element(test, c("Wald","LRT"))) > 0){
        stop(method, "\n", 
            "test: please choose the test between 'Wald' and 'LRT'.")
    }
    if (expand) {
        parameters <- expand.grid(method = method, assay_name = assay_name,
            pseudo_count = pseudo_count, test = test, boot = boot, 
            stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, assay_name = assay_name,
            pseudo_count = pseudo_count, test = test, boot = boot)
    }
    # data.frame to list
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        x <- append(x = x, values = list("formula" = formula,
            "formula_null" = formula_null,
            "phi.formula" = phi.formula,
            "phi.formula_null" = phi.formula_null,
            "coefficient" = coefficient), after = 3)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}
