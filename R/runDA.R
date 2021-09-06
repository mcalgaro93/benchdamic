#' @title runDA
#'
#' @export
#' @description
#' Run the differential abundance detection methods.
#' @param method_list a list object containing the methods and their parameters.
#' @param object a phyloseq object.
#' @param weights an optional numeric matrix giving observational weights.
#'
#' @return A named list containing the results for each method.
#'
#' @examples
#' set.seed(1)
#' # Create a very simple phyloseq object
#' counts <- matrix(rnbinom(n = 60, size = 3, prob = 0.5), nrow = 10, ncol = 6)
#' metadata <- data.frame("Sample" = c("S1", "S2", "S3", "S4", "S5", "S6"),
#'     "group" = as.factor(c("A", "A", "A", "B", "B", "B")))
#' ps <- phyloseq::phyloseq(phyloseq::otu_table(counts, taxa_are_rows = TRUE),
#'     phyloseq::sample_data(metadata))
#'
#' # Set some simple normalizations
#' my_norm <- setNormalizations()
#'
#' # Add them to the phyloseq object
#' ps <- runNormalizations(normalization_list = my_norm, object = ps)
#'
#' # Set some limma instances
#' my_methods <- set_limma(design = ~ group, coef = 2,
#'     norm = c("TMM", "poscounts", "CSSmedian"))
#'
#' # Run the methods
#' results <- runDA(method_list = my_methods, object = ps)
runDA <- function(method_list, object, weights) {
    tryCatch(
        expr = {
            out <- lapply(X = method_list, FUN = function(x) {
                method <- as.character(x[["method"]])
                cat("Running now:", method, "\n")
                params <- unlist(lapply(x[-1], paste, collapse = "."))
                param_names <- paste(names(x[-1]))
                cat("Parameters:", paste(param_names, "=", params, sep = "", collapse = ", "), "\n")
                if(is.element(el = "weights", set = names(x)))
                    if(x[["weights"]]){
                        x["weights"] <- NULL
                        x <- append(x = x, values = list("weights" = weights))
                    } else x["weights"] <- NULL
                args_list <- append(x = x[-1], values = list("object" = object), after = 0)
                do.call(what = method, args = args_list)
            })
            names(out) <- unlist(lapply(out, FUN = function(x) x[["name"]]))
            return(out)
        },
        error = function(e) {
            message(conditionMessage(e))
        }
    )
}
