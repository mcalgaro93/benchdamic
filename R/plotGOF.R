#' @title plotMD
#'
#' @importFrom plyr ddply
#' @import ggplot2
#' @export
#' @description
#' A function to plot mean difference (MD) and zero probability difference (ZPD)
#' values between estimated and observed values.
#'
#' @param data a list, output of the \code{\link{fitModels}} function. Each
#' element of the list is a `data.frame` object with \code{Model}, \code{Y}, 
#' \code{Y0}, \code{MD}, and \code{ZPD} columns containing the model name, the 
#' observed values for the mean and the zero proportion and the differences 
#' between observed and estimated values.
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
#'     object = counts, models = c(
#'         "NB", "ZINB",
#'         "DM", "ZIG", "HURDLE"
#'     ), scale_HURDLE = c("median", "default")
#' )
#'
#' # Plot the results
#' plotMD(data = GOF, difference = "MD", split = TRUE)
#' plotMD(data = GOF, difference = "ZPD", split = TRUE)
plotMD <- function(data, difference = NULL, split = TRUE) {
    Y <- MD <- Model <- Y0 <- ZPD <- NULL
    if (!is.element(difference, c("MD", "ZPD"))) {
        stop("difference: the parameter must be 'MD' or 'ZPD'.")
    }
    if (is(data, "list")) {
        RMSE <- plotRMSE(data, difference = difference, plotIt = FALSE)
        data <- plyr::ldply(data, .id = "Model")
    } else {
        stop("'data' must be a list produced by the fitModels() method.")
    }
    gobj <- ggplot(data = data, aes_string(
        x = ifelse(difference == "MD", "Y", "Y0"),
        y = difference, color = "Model")) +
        ggtitle(label = ifelse(difference == "MD",
            "Mean Differences plot", "Zero Probability Differences plot"),
            subtitle = ifelse(difference == "MD",
                paste0(
                    "Observed = log(mean(counts*)+1)", "\n",
                    "Estimated = log(mean(fitted*)+1)"
                ),
                paste0(
                    "Observed = mean(counts=0)", "\n",
                    "Estimated = mean(P(Y=0))"
                )))
    if (split) {
        gobj <- gobj + geom_text(data = RMSE, color = "black",
            aes(x = mean(data[, ifelse(difference == "MD", "Y", "Y0")]),
                y = max(data[, difference], na.rm = TRUE ),
                label = paste0("RMSE:", round(RMSE,
                ifelse(difference == "MD", 2, 4)))))
    }
    gobj <- gobj +
        geom_hline(aes(yintercept = 0), lty = 2, color = "black") +
        theme(legend.position = "bottom") + xlab("Observed") +
        ylab("Estimated-Observed")
    if (length(unique(data[, "Model"])) >= 1) {
        if (split) {
            gobj <- gobj + facet_grid(~ Model, labeller = labeller(
                .cols = label_both)) +
                geom_point(pch = 21) +
                geom_smooth(method = "loess", formula = "y ~ x", 
                    color = "black")
        } else {
            gobj <- gobj + geom_smooth(method = "loess", formula = "y ~ x")
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
#' @param plotIt logical. Should plotting be done? (default
#' \code{plotIt = TRUE})
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
#'     object = counts, models = c(
#'         "NB", "ZINB",
#'         "DM", "ZIG", "HURDLE"
#'     ), scale_HURDLE = c("median", "default")
#' )
#'
#' # Plot the RMSE results
#' plotRMSE(data = GOF, difference = "MD")
#' plotRMSE(data = GOF, difference = "ZPD")
plotRMSE <- function(data, difference = NULL, plotIt = TRUE) {
    Model <- NULL
    if (!is.element(difference, c("MD", "ZPD"))) {
        stop("difference: the parameter must be 'MD' or 'ZPD'.")
    }
    if (is(data, "list")) {
        data <- plyr::ldply(data, .id = "Model")
    } else {
        stop("'data' must be a list produced by the fitModels() method.")
    }
    RMSE <- plyr::ddply(.data = data, .variables = ~ Model,
        .fun = function(m) cbind("RMSE" = RMSE(m[, difference])))
    gobj <- ggplot2::ggplot(data = RMSE, aes(x = Model, y = RMSE,
        fill = Model)) +
        geom_col() +
        geom_label(aes(x = Model, y = RMSE, label = round(RMSE, 2)),
            fill = "white") +
        ggtitle(label = "RMSE", subtitle = ifelse(difference == "MD",
            "Mean Difference", "Zero Probability Difference")) +
        theme(legend.position = "bottom", axis.text.x = element_text(
            angle = 90, hjust = 1, vjust = 0.5)) +
        scale_x_discrete(limits = RMSE[order(RMSE[, "RMSE"]), "Model"])
    if(plotIt){
        return(gobj)
    } else return(RMSE)
}
