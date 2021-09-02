#' @title createSplits
#'
#' @export
#' @importFrom phyloseq sample_data
#' @description
#' Given the phyloseq object from which the random splits should be created,
#' this function produces a list of 2 \code{data.frame} objects: \code{Subset1}
#' and \code{Subset2} with as many rows as the number of splits and as many
#' columns as the half of the number of samples.
#' @param object a \code{phyloseq} object.
#' @param varName name of a factor variable with 2 levels.
#' @param paired name of the unique subject identifier variable. If specified,
#' paired samples will remain in the same split. (default = NULL).
#' @param balanced If \code{TRUE} a balanced design will be created for the
#' splits. (Ignored if paired is supplied).
#' @param N number of splits to generate.
#'
#' @return A list of 2 \code{data.frame} objects: \code{Subset1} and
#' \code{Subset2} containing \code{N} rows and half of the total number of
#' samples columns. Each cell contains a unique sample identifier.
#' @examples
#' data(ps_plaque_16S)
#' set.seed(123)
#'
#' # Balanced design for repeated measures
#' splits_df <- createSplits(
#'     object = ps_plaque_16S, varName =
#'         "HMP_BODY_SUBSITE", paired = "RSID", balanced = TRUE, N = 100
#' )
#'
#' # Balanced design for independent samples
#' splits_df <- createSplits(
#'     object = ps_plaque_16S, varName =
#'         "HMP_BODY_SUBSITE", balanced = TRUE, N = 100
#' )
#'
#' # Unbalanced design
#' splits_df <- createSplits(
#'     object = ps_plaque_16S, varName =
#'         "HMP_BODY_SUBSITE", balanced = FALSE, N = 100
#' )

createSplits <- function(object, varName = NULL, paired = NULL, balanced = TRUE,
                         N = 1000) {
    metadata <- data.frame(phyloseq::sample_data(object))
    # Take variable names and levels
    if (is.null(varName)) {
        stop("Please supply the name of the variable to perform the splitting")
    }
    else if (!is.element(varName, colnames(metadata))) {
        stop("Variable not found")
    }
    else {
        variable <- metadata[, varName]
        if (!is.factor(variable)) {
            warning(paste(
                "The variable", varName, "is not a factor.",
                "Coercing to factor."
            ))
            variable <- as.factor(variable)
        }
        var_levels <- levels(variable)
        if (length(var_levels) != 2) {
              stop(paste("The variable", varName, "has not 2 levels."))
          } else {
            grp1_name <- var_levels[1]
            grp2_name <- var_levels[2]
        }
    }
    if (!is.null(paired)) {
        # 1 ID is enough to identify the two paired samples
        ID <- unique(metadata[, paired])
        num_ID <- length(ID)
        half <- num_ID / 2
        Subset1 <- matrix(NA, nrow = N, ncol = 2 * half)
        Subset2 <- matrix(NA, nrow = N, ncol = 2 * half)
        for (i in seq_len(N)) {
            chosenID <- sample(x = ID, size = half, replace = FALSE)
            Subset1[i, ] <- rownames(metadata[metadata[, paired] %in% chosenID,
                ])
            Subset2[i, ] <- rownames(metadata[metadata[, paired] %in% setdiff(
                ID,
                chosenID
            ), ])
        }
    }
    else {
        # find indexes for the 2 levels
        ID1 <- rownames(metadata[metadata[, varName] == grp1_name, ])
        ID2 <- rownames(metadata[metadata[, varName] == grp2_name, ])
        num_ID1 <- length(ID1)
        num_ID2 <- length(ID2)
        if (balanced) {
            min_length <- min(num_ID1, num_ID2)
            ID1 <- ID1[seq_len(min_length)]
            ID2 <- ID2[seq_len(min_length)]
            num_ID1 <- num_ID2 <- min_length
            half <- rep(min_length %/% 2, 2)
        } else {
            half <- c(num_ID1, num_ID2) %/% 2
        }
        Subset1 <- matrix(NA, nrow = N, ncol = sum(half))
        Subset2 <- matrix(NA, nrow = N, ncol = sum(c(num_ID1, num_ID2)) - sum(
            half
        ))
        for (i in seq_len(N)) {
            chosenID1 <- sample(x = ID1, size = half[1], replace = FALSE)
            chosenID2 <- sample(x = ID2, size = half[2], replace = FALSE)
            Subset1[i, ] <- c(chosenID1, chosenID2)
            Subset2[i, ] <- c(setdiff(ID1, chosenID1), setdiff(ID2, chosenID2))
        }
    }
    rownames(Subset1) <- rownames(Subset2) <- paste0("Comparison", seq_len(N))
    return(list("Subset1" = Subset1, "Subset2" = Subset2))
} # END - function: createSplits

#' @title runSplits
#'
#' @export
#' @importFrom phyloseq prune_samples sample_names filter_taxa
#' @importFrom phyloseq filter_taxa prune_samples
#' @description
#' Run the differential abundance detection methods on split datasets.
#'
#' @param split_list A list of 2 \code{data.frame} objects: \code{Subset1} and
#' \code{Subset2} produced by the \code{\link{createSplits}} function.
#' @inheritParams runDA
#' @param normalization_list a list object containing the normalization method
#' names and their parameters produced by \code{\link{setNormalizations}}.
#'
#' @return A named list containing the results for each method.
#'
#' @examples
#' data(ps_plaque_16S)
#'
#' # Balanced design for independent samples
#' my_splits <- createSplits(
#'     object = ps_plaque_16S, varName =
#'     "HMP_BODY_SUBSITE", balanced = TRUE, N = 10 # N = 100 suggested
#' )
#'
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ HMP_BODY_SUBSITE, coef = 2,
#'     norm = c("TMM", "CSSmedian"))
#'
#' # Set the normalization methods according to the DA methods
#' my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
#'     method = c("TMM", "median"))
#'
#' # Run methods on split datasets
#' results <- runSplits(split_list = my_splits, method_list = my_limma,
#'     normalization_list = my_norm, object = ps_plaque_16S)

runSplits <- function(split_list, method_list, normalization_list, object){
    out <- lapply(X = split_list, FUN = function(subset) {
        # apply -> Comparison1, Comparison2, ..., ComparisonN
        apply(X = subset, MARGIN = 1, FUN = function(splits) {
            # Splitting
            cat("Splitting the samples...\n")
            ps <- phyloseq::prune_samples(phyloseq::sample_names(object) %in%
                splits, object)
            # Keep only present taxa
            cat("Removing not present taxa...\n")
            ps <- phyloseq::filter_taxa(ps, function(x) sum(x > 0) > 0, 1)
            # Adding scaling and normalization factors
            cat("Computing normalizations...\n")
            ps <- runNormalizations(normalization_list = normalization_list,
                object = ps)
            # Compute weights if necessary
            weights_info <- unlist(lapply(X = method_list, FUN = function(x){
                x[["weights"]]
            }))
            if(sum(weights_info) > 0){
                cat("Computing ZINB weights...\n")
                weights <- weights_ZINB(ps, design = ~ 1)
            } else weights <- NULL
            ### DA analysis ###
            cat("Differential abundance:\n")
            runDA(method_list = method_list, object = ps, weights = weights)
        })
    })
    return(out)
}

#' @title createConcordance
#'
#' @export
#' @importFrom ffpe CATplot
#' @importFrom plyr ddply
#' @description
#' Compute the between and within method concordances comparing the lists of
#' extracted statistics from the outputs of the differential abundance detection
#' methods.
#'
#' @inheritParams extractStatistics
#'
#' @return A long format \code{data.frame} object with several columns:
#' \itemize{
#'     \item{\code{comparison}}{ which indicates the comparison number;}
#'     \item{\code{n_features}}{ which indicates the total number of taxa in
#'     the comparison dataset;}
#'     \item{\code{method1}}{ which contains the first method name;}
#'     \item{\code{method2}}{ which contains the first method name;}
#'     \item{\code{rank}}{;}
#'     \item{\code{concordance}}{ which is defined as the cardinality of the
#'     intersection of the top rank elements of each list, divided by rank, i.e.
#'     , \eqn{(L_{1:rank} \bigcap M_{1:rank})/(rank)}, where L and M represent
#'     the lists of the extracted statistics of method1 and method2
#'     respectively (averaged values between subset1 and subset2).}}
#'
#' @seealso \code{\link{extractStatistics}} and \code{\link{areaCAT}}.
#'
#' @examples
#' data(ps_plaque_16S)
#'
#' # Balanced design for independent samples
#' my_splits <- createSplits(
#'     object = ps_plaque_16S, varName =
#'     "HMP_BODY_SUBSITE", balanced = TRUE, N = 10 # N = 100 suggested
#' )
#'
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ HMP_BODY_SUBSITE, coef = 2,
#'     norm = c("TMM", "CSSmedian"))
#'
#' # Set the normalization methods according to the DA methods
#' my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
#'     method = c("TMM", "median"))
#'
#' # Run methods on split datasets
#' results <- runSplits(split_list = my_splits, method_list = my_limma,
#'     normalization_list = my_norm, object = ps_plaque_16S)
#'
#' # Concordance for p-values
#' concordance_pvalues <- createConcordance(
#'     object = results, slot = "pValMat", colName = "rawP", type = "pvalue"
#' )
#'
#' # Concordance for log fold changes
#' concordance_logfc <- createConcordance(
#'     object = results, slot = "statInfo", colName = "logFC", type = "logfc"
#' )
#'
#' # Concordance for log fold changes in the first method and p-values in the
#' # other
#' concordance_logfc_pvalues <- createConcordance(
#'     object = results, slot = c("statInfo", "pValMat"),
#'     colName = c("logFC", "rawP"), type = c("logfc", "pvalue")
#' )

createConcordance <- function(object, slot = "pValMat", colName = "rawP",
                              type = "pvalue") {
    data <- lapply(X = object, FUN = function(subset) {
        lapply(X = subset, FUN = function(comparison) {
            c(extractStatistics(
                object = comparison, slot = slot,
                colName = colName, type = type
            ),
            n_features =
                nrow(comparison[[1]][["pValMat"]])
            )
        })
    })
    subset1 <- data[["Subset1"]]
    subset2 <- data[["Subset2"]]
    conc_df <- NULL
    for (i in seq_len(length(subset1))) {
        comparison1 <- subset1[[i]]
        comparison2 <- subset2[[i]]
        n_features1 <- comparison1[["n_features"]]
        n_features2 <- comparison2[["n_features"]]
        n_features <- min(n_features1, n_features2)
        n_methods <- length(comparison1) - 1
        method_names <- names(comparison1)[seq_len(n_methods)]
        for (j in seq_len(n_methods)) {
            for (k in seq_len(n_methods)) {
                if (j != k) { # Compute the Between Methods concordance
                    # For subset1
                    vec1 <- comparison1[[method_names[j]]]
                    vec2 <- comparison1[[method_names[k]]]
                    conc1 <- ffpe::CATplot(
                        vec1 = vec1, vec2 = vec2, make.plot =
                            FALSE
                    )
                    conc1 <- data.frame(conc1, "n_features" = n_features1)
                    # For subset2
                    vec1 <- comparison2[[method_names[j]]]
                    vec2 <- comparison2[[method_names[k]]]
                    conc2 <- ffpe::CATplot(
                        vec1 = vec1, vec2 = vec2, make.plot =
                            FALSE
                    )
                    conc2 <- data.frame(conc2, "n_features" = n_features2)
                    # Together
                    conc <- rbind(conc1, conc2)
                } else { # Compute the Within Method concordance
                    vec1 <- comparison1[[method_names[j]]]
                    vec2 <- comparison2[[method_names[k]]]
                    conc <- ffpe::CATplot(
                        vec1 = vec1, vec2 = vec2, make.plot =
                            FALSE
                    )
                    conc <- data.frame(conc, "n_features" = n_features)
                }
                conc[, "method1"] <- method_names[j]
                conc[, "method2"] <- method_names[k]
                conc[, "comparison"] <- paste0("Comparison", i)
                conc_df <- rbind(conc_df, conc)
            }
        }
    }
    concordance <- plyr::ddply(.data = conc_df, .variables = ~ comparison +
        n_features + method1 + method2 + rank, .fun = function(x) {
        data.frame(
            "concordance" = mean(x[, "concordance"])
        )
    })
    return(concordance)
}

#' @title areaCAT
#'
#' @export
#' @importFrom plyr ddply
#' @importFrom graphics plot abline
#' @description
#' Compute the area between the bisector and the concordance curve.
#'
#' @param concordance A long format \code{data.frame} produced by
#' \link{createConcordance} function.
#' @param plotIt Plot the concordance (default \code{plotIt = FALSE}).
#'
#' @return A long format \code{data.frame} object with several columns:
#' \itemize{
#'     \item{\code{comparison}}{ which indicates the comparison number;}
#'     \item{\code{n_features }}{ which indicates the total number of taxa in
#'     the comparison dataset;}
#'     \item{\code{method1}}{ which contains the first method name;}
#'     \item{\code{method2}}{ which contains the first method name;}
#'     \item{\code{rank}}{;}
#'     \item{\code{concordance}}{ which is defined as the cardinality of the
#'     intersection of the top rank elements of each list, divided by rank, i.e.
#'     , \eqn{(L_{1:rank} \bigcap M_{1:rank})/(rank)}, where L and M represent
#'     the lists of the extracted statistics of method1 and method2
#'     respectively;}
#'     \item{\code{heightOver}}{ which is the distance between the bisector and
#'     the concordance value;}
#'     \item{\code{areaOver}}{ which is the cumulative sum of the
#'     \code{heightOver} value.}}
#'
#' @seealso \code{\link{createConcordance}} and \code{\link{plotConcordance}}
#'
#' @examples
#' data(ps_plaque_16S)
#'
#' # Balanced design for independent samples
#' my_splits <- createSplits(
#'     object = ps_plaque_16S, varName =
#'     "HMP_BODY_SUBSITE", balanced = TRUE, N = 10 # N = 100 suggested
#' )
#'
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ HMP_BODY_SUBSITE, coef = 2,
#'     norm = c("TMM", "CSSmedian"))
#'
#' # Set the normalization methods according to the DA methods
#' my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
#'     method = c("TMM", "median"))
#'
#' # Run methods on split datasets
#' results <- runSplits(split_list = my_splits, method_list = my_limma,
#'     normalization_list = my_norm, object = ps_plaque_16S)
#'
#' # Concordance for p-values
#' concordance_pvalues <- createConcordance(
#'     object = results, slot = "pValMat", colName = "rawP", type = "pvalue"
#' )
#'
#' # Add area over the concordance curve
#' concordance_area <- areaCAT(concordance = concordance_pvalues)

areaCAT <- function(concordance, plotIt = FALSE) {
    plyr::ddply(.data = concordance, .variables = ~ comparison + n_features +
        method1 + method2, .fun = function(conc) {
        MaxArea <- unique(conc[, "n_features"])
        estimated <- conc[, "concordance"]
        # y = x values -> bisector
        theoretical <- seq_along(estimated) / unique(conc[, "n_features"])
        if (plotIt) {
            graphics::plot(x = theoretical, y = estimated, type = "l", xlim = c(
                0, 1
            ), ylim = c(0, 1), main = paste0("CAT for ", unique(conc[
                ,
                "method1"
            ]), " & ", unique(conc[, "method2"])))
            graphics::abline(0, 1)
        }
        # Consider the y = x line
        difference <- estimated - theoretical
        HeightOver <- (estimated - theoretical) / MaxArea # rescaling
        HeightOver[difference <= 0] <- 0 # removing negative values
        return(data.frame(
            "rank" = conc[, "rank"],
            "concordance" = conc[, "concordance"], "heightOver" = HeightOver,
            "areaOver" = cumsum(HeightOver)
        ))
    })
}
