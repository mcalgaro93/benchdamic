#' @title prepareObserved
#'
#' @importFrom stats median
#' @importFrom phyloseq otu_table taxa_are_rows
#' @export
#' @description
#' Continuity corrected logarithms of the average counts and fraction of zeroes
#' by feature.
#'
#' @inheritParams fitNB
#' @param scale If specified it refers to the character vector used in
#' \code{\link{fitHURDLE}} function. Either \code{median} or \code{default} to
#' choose between the median library size or one million as scaling factors for
#' raw counts.
#'
#' @return A data frame containing the continuity corrected logarithm for the
#' raw count mean values for each taxon of the matrix of counts in the \code{Y}
#' column and the observed zero rate in the \code{Y0} column. If \code{scale} is
#' specified the continuity corrected logarithm for the mean CPM
#' (\code{scale = "default"}) or the mean counts per median library size
#' (\code{scale = "median"}) is computed instead.
#'
#' @seealso \code{\link{meanDifferences}}
#'
#' @examples
#' # Generate some random counts
#' counts <- matrix(rnbinom(n = 60, size = 3, prob = 0.5), nrow = 10, ncol = 6)
#'
#' observed1 <- prepareObserved(counts)
#' # For the comparison with HURDLE model
#' observed2 <- prepareObserved(counts, scale = "median")
prepareObserved <- function(object, assay_name = "counts", scale = NULL) {
    if(is(object, "phyloseq") | is(object, "TreeSummarizedExperiment")){
        counts_and_metadata <- get_counts_metadata(object,
            assay_name = assay_name)
        counts <- counts_and_metadata[[1]]
    } else if (is.matrix(object)) {
        NULL
    } else {
        stop("Please supply a phyloseq object, a TreeSummarizedExperiment,", 
             " or a matrix of counts.")
    }
    if (!is.null(scale)) {
        if (scale == "median") {
            counts <- counts * stats::median(colSums(counts)) / colSums(counts)
        } else if (scale == "default") {
            counts <- counts * 1e6 / colSums(counts)
        } else {
            stop("When specified, 'scale' must be 'median' or 'default'")
        }
    }
    Y <- log1p(rowMeans(counts))
    Y0 <- rowMeans(counts == 0)
    return(data.frame("Y" = Y, "Y0" = Y0))
}

#' @title meanDifferences
#'
#' @export
#' @description
#' Compute the differences between the estimated and the observed continuity
#' corrected logarithms of the average count values (MD), and between the
#' estimated average probability to observe a zero and the the observed zero
#' rate (ZPD).
#'
#' @param estimated a two column data.frame, output of \code{\link{fitNB}},
#' \code{\link{fitZINB}}, \code{\link{fitDM}}, \code{\link{fitZIG}}, or
#' \code{\link{fitHURDLE}} functions. More in general, a data frame
#' containing the continuity corrected logarithm for the average of the fitted
#' values for each row of a matrix of counts in the \code{Y} column, and the
#' estimated probability to observe a zero in the \code{Y0} column.
#' @param observed a two column data.frame, output of
#' \code{\link{prepareObserved}} function. More in general, a data frame
#' containing the continuity corrected logarithm for the average of the observed
#' values for each row of a matrix of counts in the \code{Y} column, and the
#' estimated proportion of zeroes in the \code{Y0} column.
#'
#' @return a \code{data.frame} containing the differences between the estimated
#' and the observed continuity corrected logarithms of the average count values
#' in the \code{MD} column, and between the estimated average probability to
#' observe a zero and the the observed zero rate in the \code{ZPD} column.
#'
#' @seealso \code{\link{prepareObserved}}.
#'
#' @examples
#' # Randomly generate the observed and estimated data.frames
#' observed <- data.frame(Y = rpois(10, 5), Y0 = runif(10, 0, 1))
#' estimated <- data.frame(Y = rpois(10, 5), Y0 = runif(10, 0, 1))
#'
#' # Compute the mean differences between estimated and observed data.frames
#' meanDifferences(estimated, observed)
meanDifferences <- function(estimated, observed) {
    if (sum(colnames(estimated) != colnames(observed)) > 0) {
        stop("Estimated and Observed data.frames have different colnames")
    }
    if (sum(colnames(estimated) != c("Y", "Y0")) > 0 | sum(colnames(observed) !=
        c("Y", "Y0")) > 0) {
        stop("Please rename the colnames of the data.frames as Y and Y0")
    }
    if (nrow(estimated) != nrow(observed)) {
        stop("Estimated and Observed data.frames have different number of rows")
    }
    MD <- estimated[, "Y"] - observed[, "Y"]
    ZPD <- estimated[, "Y0"] - observed[, "Y0"]
    return(data.frame("MD" = MD, "ZPD" = ZPD))
}

#' @title RMSE
#'
#' @export
#' @description
#' Computes the Root Mean Square Error (RMSE) from a vector of differences.
#'
#' @param differences a vector of differences.
#'
#' @return RMSE value
#'
#' @seealso \code{\link{prepareObserved}} and \code{\link{meanDifferences}}.
#'
#' @examples
#' # Generate the data.frame of Mean Differences and Zero Probability Difference
#' MD_df <- data.frame(MD = rpois(10, 5), ZPD = runif(10, -1, 1))
#'
#' # Calculate RMSE for MD and ZPD values
#' RMSE(MD_df[, "MD"])
#' RMSE(MD_df[, "ZPD"])
RMSE <- function(differences) {
    sqrt(mean(differences^2, na.rm = TRUE))
}

#' @title fitModels
#'
#' @export
#' @description
#' A wrapper function that fits the specified models for each taxon of the count
#' data and computes the mean difference (MD) and zero probability difference
#' (ZPD) between estimated and observed values.
#'
#' @inheritParams fitNB
#' @param models character vector which assumes the values \code{NB},
#' \code{ZINB}, \code{DM}, \code{ZIG}, and \code{HURDLE}.
#' @param scale_HURDLE character vector, either \code{median} or \code{default}
#' to choose between the median of the library size or one million to scale raw
#' counts for the truncated gaussian hurdle model.
#'
#' @return list of \code{data.frame} objects for each model. The first two
#' columns contain the properly transformed observed values for mean and zero
#' proportion, while the third and the fourth columns contain the estimated
#' values for the mean and the zero rate respectively.
#'
#' @seealso \code{\link{fitNB}}, \code{\link{fitZINB}}, \code{\link{fitDM}},
#' \code{\link{fitZIG}}, and \code{\link{fitHURDLE}} for the model estimations,
#' \code{\link{prepareObserved}} for raw counts preparation, and
#' \code{\link{meanDifferences}} for the Mean Difference (MD) and Zero
#' Probability Difference (ZPD) computations.
#'
#' @examples
#' # Generate some random counts
#' counts <- matrix(rnbinom(n = 60, size = 3, prob = 0.5), nrow = 10, ncol = 6)
#' # Estimate the counts assuming several distributions
#' GOF <- fitModels(
#'     object = counts, models = c(
#'         "NB", "ZINB",
#'         "DM", "ZIG", "HURDLE"
#'     ), scale_HURDLE = c("median", "default")
#' )
#'
#' head(GOF)
fitModels <- function(object, assay_name = "counts", models = c("NB", "ZINB", 
    "DM", "ZIG", "HURDLE"), scale_HURDLE = c("default", "median"), 
    verbose = TRUE) {
    fittedModels <- list()
    observed <- prepareObserved(object, assay_name = assay_name)
    if ("NB" %in% models) {
        fitted <- fitNB(object = object, assay_name = assay_name, 
            verbose = verbose)
        MD <- meanDifferences(estimated = fitted, observed = observed)
        fittedModels[["NB"]] <- data.frame(observed, MD)
    }
    if ("ZINB" %in% models) {
        fitted <- fitZINB(object = object, assay_name = assay_name, 
            verbose = verbose)
        MD <- meanDifferences(estimated = fitted, observed = observed)
        fittedModels[["ZINB"]] <- data.frame(observed, MD)
    }
    if ("DM" %in% models) {
        fitted <- fitDM(object = object, assay_name = assay_name, 
            verbose = verbose)
        MD <- meanDifferences(estimated = fitted, observed = observed)
        fittedModels[["DM"]] <- data.frame(observed, MD)
    }
    if ("ZIG" %in% models) {
        fitted <- fitZIG(object = object, assay_name = assay_name, 
            verbose = verbose)
        MD <- meanDifferences(estimated = fitted, observed = observed)
        fittedModels[["ZIG"]] <- data.frame(observed, MD)
    }
    if ("HURDLE" %in% models) {
        fitted <- fitHURDLE(object = object, assay_name = assay_name, 
            scale = scale_HURDLE[1], verbose = verbose)
        observed <- prepareObserved(object = object, assay_name = assay_name, 
            scale = scale_HURDLE[1])
        MD <- meanDifferences(estimated = fitted, observed = observed)
        name <- paste0("HURDLE_", scale_HURDLE[1])
        fittedModels[[name]] <- data.frame(observed, MD)
        if (length(scale_HURDLE) == 2) {
            fitted <- fitHURDLE(object = object, assay_name = assay_name, 
                scale = scale_HURDLE[2], verbose = verbose)
            observed <- prepareObserved(object = object, 
                assay_name = assay_name, scale = scale_HURDLE[2])
            MD <- meanDifferences(estimated = fitted, observed = observed)
            name <- paste0("HURDLE_", scale_HURDLE[2])
            fittedModels[[name]] <- data.frame(observed, MD)
        }
    }
    return(fittedModels)
}
