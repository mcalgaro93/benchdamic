#' @title plotMD
#'
#' @importFrom plyr ddply
#' @import ggplot2
#' @export
#' @description
#' A function to plot mean difference (MD) and zero probability difference (ZPD)
#' values between estimated and observed values.
#'
#' @param data a list, output of the \code{\link{fitModels}} function or a
#' `data.frame` object with \code{Model}, \code{Y}, \code{Y0}, \code{MD}, and
#' \code{ZPD} columns containing the model name, the observed values for the
#' mean and the zero proportion and the differences between observed and
#' estimated values.
#' @param difference character vector, either \code{MD} or \code{ZPD} to plot
#' the differences between estimated and observed mean counts or the differences
#' between estimated zero probability and observed zero proportion.
#' @param split Display each model mean differences in different facets
#' (default \code{split = TRUE}). If \code{FALSE}, points are not displayed for
#' more clear representation.
#'
#' @return a \code{ggplot} object.
#'
#' @seealso \code{\link{fitModels}} and \code{\link{RMSE}} for the model
#' estimations and the RMSE computations respectively. \code{\link{plotRMSE}}
#' for the graphical evaluation of the RMSE values.
#'
#' @examples
#' # Generate some random counts
#' counts = matrix(rnbinom(n = 60, size = 3, prob = 0.5), nrow = 10, ncol = 6)
#'
#' # Estimate the counts assuming several distributions
#' GOF <- fitModels(
#'     counts = counts, models = c(
#'         "NB", "ZINB",
#'         "DM", "ZIG", "HURDLE"
#'     ), scale_ZIG = c("median", "default"), scale_HURDLE =
#'         c("median", "default")
#' )
#'
#' # Plot the results
#' plotMD(data = GOF, difference = "MD", split = TRUE)
#' plotMD(data = GOF, difference = "ZPD", split = TRUE)
plotMD <- function(data, difference = NULL, split = TRUE) {
    Y <- MD <- Model <- Y0 <- ZPD <- NULL
    if (is.list(data)) {
        data <- plyr::ldply(data, .id = "Model")
    }
    if (difference == "MD") {
        if ("Model" %in% colnames(data) & "Y" %in% colnames(data) & "MD" %in%
            colnames(data)) {
            RMSE_MD <- plyr::ddply(
                .data = data, .variables = ~Model, .fun =
                    function(m) cbind("RMSE" = RMSE(m[, "MD"]))
            )
            gobj <- ggplot(data = data, aes(x = Y, y = MD, color = Model)) +
                ggtitle(
                    label = "Mean Differences plot", subtitle =
                        paste0(
                            "Observed = log(mean(counts*)+1)", "\n",
                            "Estimated = log(mean(fitted*)+1)"
                        )
                )
            if (split) {
                gobj <- gobj + geom_text(
                    data = RMSE_MD, color = "black",
                    aes(x = mean(data[, "Y"]), y = max(data[, "MD"],
                        na.rm =
                            TRUE
                    ), label = paste0("RMSE:", round(RMSE, 2)))
                )
            }
        } else {
            stop("data should contain 'Model', 'Y', and 'MD' columns for
            model name, observed values and mean difference values respectively.
            ")
        }
    } else if (difference == "ZPD") {
        if ("Model" %in% colnames(data) & "Y0" %in% colnames(data) & "ZPD" %in%
            colnames(data)) {
            RMSE_ZPD <- plyr::ddply(
                .data = data, .variables = ~Model, .fun =
                    function(m) cbind("RMSE" = RMSE(m[, "ZPD"]))
            )
            gobj <- ggplot(data = data, aes(x = Y0, y = ZPD, color = Model)) +
                ggtitle(
                    label = "Zero Probability Differences plot", subtitle =
                        paste0(
                            "Observed = mean(counts=0)", "\n",
                            "Estimated = mean(P(Y=0))"
                        )
                )
            if (split) {
                gobj <- gobj + geom_text(
                    data = RMSE_ZPD, color = "black",
                    aes(x = mean(data[, "Y0"]), y = max(data[, "ZPD"],
                        na.rm = TRUE
                    ), label = paste0("RMSE:", round(RMSE, 4)))
                )
            }
        } else {
            stop("data should contain 'Model', 'Y0', and 'ZPD' columns for
            model name, zero rate observed values and zero probability
            difference values respectively.")
        }
    } else {
        stop("Difference must be 'MD' or 'ZPD'")
    }
    gobj <- gobj +
        geom_hline(aes(yintercept = 0), lty = 2, color = "black") +
        theme(legend.position = "bottom") + xlab("Observed") +
        ylab("Estimated-Observed")
    if (length(unique(data[, "Model"])) > 1) {
        if (split) {
            gobj <- gobj + facet_grid(~Model, labeller = labeller(
                .cols =
                    label_both
            )) + geom_point(pch = 21) + geom_smooth(
                color =
                    "black"
            )
        } else {
            gobj <- gobj + geom_smooth()
        }
    }
    return(gobj)
}

#' @title plotRMSE
#'
#' @importFrom plyr ddply
#' @import ggplot2
#' @export
#' @description
#' A function to plot RMSE values computed for mean difference (MD) and zero
#' probability difference (ZPD) values between estimated and observed values.
#'
#' @inheritParams plotMD
#'
#' @return a \code{ggplot} object.
#'
#' @seealso \code{\link{fitModels}} and \code{\link{RMSE}} for the model
#' estimations and the RMSE computations respectively. \code{\link{plotMD}} for
#' the graphical evaluation.
#'
#' @examples
#' # Generate some random counts
#' counts = matrix(rnbinom(n = 60, size = 3, prob = 0.5), nrow = 10, ncol = 6)
#'
#' # Estimate the counts assuming several distributions
#' GOF <- fitModels(
#'     counts = counts, models = c(
#'         "NB", "ZINB",
#'         "DM", "ZIG", "HURDLE"
#'     ), scale_ZIG = c("median", "default"), scale_HURDLE =
#'         c("median", "default")
#' )
#'
#' # Plot the RMSE results
#' plotRMSE(data = GOF, difference = "MD")
#' plotRMSE(data = GOF, difference = "ZPD")
plotRMSE <- function(data, difference = NULL) {
    Model <- NULL
    if (is.list(data)) {
          data <- plyr::ldply(data, .id = "Model")
      }
    if (difference == "MD") {
        if ("Model" %in% colnames(data) & "Y" %in% colnames(data) & "MD" %in%
            colnames(data)) {
            RMSE <- plyr::ddply(
                .data = data, .variables = ~Model, .fun =
                    function(m) cbind("RMSE" = RMSE(m[, "MD"]))
            )
            gobj <- ggplot2::ggplot(data = RMSE, aes(
                x = Model, y = RMSE, fill =
                    Model
            )) +
                geom_col() +
                geom_label(aes(
                    x = Model, y = RMSE,
                    label = round(RMSE, 2)
                ), fill = "white") +
                ggtitle(
                    label =
                        "RMSE", subtitle = "Mean differences"
                )
        } else {
            stop("data should contains 'Model', 'Y', and 'MD' columns for
            model name, observed values and mean difference values respectively.
            ")
        }
    } else if (difference == "ZPD") {
        if ("Model" %in% colnames(data) & "Y0" %in% colnames(data) & "ZPD" %in%
            colnames(data)) {
            RMSE <- plyr::ddply(
                .data = data, .variables = ~Model, .fun =
                    function(m) cbind("RMSE" = RMSE(m[, "ZPD"]))
            )
            gobj <- ggplot2::ggplot(data = RMSE, aes(
                x = Model, y = RMSE, fill =
                    Model
            )) +
                geom_col() +
                geom_label(aes(
                    x = Model, y = RMSE,
                    label = round(RMSE, 4)
                ), fill = "white") +
                ggtitle(
                    label =
                        "RMSE", subtitle = "Zero probability difference"
                )
        } else {
            stop("df should contains 'Model', 'Y0', and 'ZPD' columns for
            model name, zero rate observed values and zero probability
            difference values respectively.")
        }
    } else {
        stop("Difference must be 'MD' or 'ZPD'")
    }
    gobj <- gobj + theme(legend.position = "bottom", axis.text.x = element_text(
        angle = 90, hjust = 1, vjust = 0.5
    )) + scale_x_discrete(
        limits =
            RMSE[order(RMSE[, "RMSE"]), "Model"]
    )
    return(gobj)
}
