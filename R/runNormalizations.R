#' @title checkNormalization
#'
#' @export
#' @description
#' Check if the normalization function's name and the method's name to compute
#' normalization/scaling factors are correctly matched.
#'
#' @param fun a character with the name of normalization function (e.g.
#' "norm_edgeR", "norm_DESeq2", "norm_CSS"...).
#' @param method a character with the normalization method (e.g.
#' "TMM", "upperquartile"... if the \code{fun} is "norm_edgeR").
#' @param ... other arguments if needed (e.g. for \code{\link{norm_edgeR}}
#' normalizations).
#'
#' @return a list object containing the normalization method and its
#' parameters.
#'
#' @seealso \code{\link{setNormalizations}}, \code{\link{norm_edgeR}},
#' \code{\link{norm_DESeq2}}, \code{\link{norm_CSS}}, \code{\link{norm_TSS}}
#'
#' @examples
#' # Check if TMM normalization belong to "norm_edgeR"
#' check_TMM_normalization <- checkNormalization(fun = "norm_edgeR",
#'     method = "TMM")
checkNormalization <- function(fun, method, ...){
    normalization_list <- list(
        norm_edgeR = c("TMM", "TMMwsp", "RLE", "upperquartile",
            "posupperquartile", "none"),
        norm_DESeq2 = c("ratio", "poscounts", "iterate"),
        norm_CSS = c("median", "default"),
        norm_TSS = "TSS"
    )
    if(length(fun) > 1 | length(method) > 1){
        stop("Plase supply only one normalization.")
    }
    if(is.element(el = method, set = normalization_list[[fun]])){
        return(list(fun = fun, method = method, ...))
    } else {
        stop(method, " normalization doesn't belong to ", fun, " function.")
    }
}

#' @title setNormalizations
#'
#' @export
#' @description
#' Set the methods and parameters to compute normalization/scaling factors.
#'
#' @inheritParams checkNormalization
#'
#' @return a list object containing the normalization methods and their
#' parameters.
#'
#' @seealso \code{\link{runNormalizations}}, \code{\link{norm_edgeR}},
#' \code{\link{norm_DESeq2}}, \code{\link{norm_CSS}}, \code{\link{norm_TSS}}
#'
#' @examples
#' # Set a TMM normalization
#' my_TMM_normalization <- setNormalizations(fun = "norm_edgeR", method = "TMM")
#'
#' # Set some simple normalizations
#' my_normalizations <- setNormalizations()
#'
#' # Add a custom normalization
#' my_normalizations <- c(my_normalizations,
#'     myNormMethod1 = list("myNormMethod", "parameter1", "parameter2"))

setNormalizations <- function(fun = c("norm_edgeR", "norm_DESeq2", "norm_CSS",
    "norm_edgeR"), method = c("TMM", "poscounts", "median", "none")){
    if(length(fun) != length(method)){
        stop("Numbers of methods and functions are different.")
    } else {
        mapply(checkNormalization, fun, method, SIMPLIFY = FALSE)
    }
}

#' @title runNormalizations
#'
#' @export
#' @description
#' Add normalization/scaling factors to a phyloseq object
#'
#' @param normalization_list a list object containing the normalization methods
#' and their parameters.
#' @param object a phyloseq object.
#' @param verbose an optional logical value. If \code{TRUE}, information about
#' the steps of the algorithm is printed. Default \code{verbose = TRUE}.
#'
#' @return A phyloseq object containing the normalization/scaling factors.
#'
#' @seealso \code{\link{setNormalizations}}
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
#' # Set some simple normalizations
#' my_normalizations <- setNormalizations()
#'
#' # Add them to the phyloseq object
#' ps <- runNormalizations(normalization_list = my_normalizations, object = ps)
runNormalizations <- function(normalization_list, object, verbose = TRUE) {
    tryCatch(
        expr = {
            for(x in normalization_list){
                fun <- as.character(x[["fun"]])
                if(verbose)
                    cat("      + Running now:", fun, "\n")
                params <- unlist(lapply(x[-1], paste, collapse = "."))
                param_names <- paste(names(x[-1]))
                if(verbose)
                    cat("        Parameters:", paste(param_names, "=", params,
                        sep = "", collapse = ", "), "\n")
                args_list <- append(x = x[-1], values = list("object" = object,
                    "verbose" = verbose), after = 0)
                object <- do.call(what = fun, args = args_list)
            }
            return(object)
        },
        error = function(e) {
            message(conditionMessage(e))
        }
    )
}
