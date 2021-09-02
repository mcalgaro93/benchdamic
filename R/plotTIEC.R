#' @title plotFPR
#'
#' @export
#' @importFrom plyr ddply
#' @import ggplot2
#' @description
#' Draw the boxplots of the proportions of p-values lower than 0.01, 0.05, and
#' 0.1 thresholds for each method.
#'
#' @param df_FPR a \code{data.frame} produced by the \code{\link{createTIEC}}
#' function, containing the FPR values.
#' @param cols named vector of colors.
#'
#' @return A ggplot object.
#'
#' @examples
#' # Load some data
#' data(ps_stool_16S)
#'
#' # Generate the patterns for 10 mock comparison for an experiment
#' # (N = 1000 is suggested)
#' mocks <- createMocks(nsamples = phyloseq::nsamples(ps_stool_16S), N = 10)
#' head(mocks)
#'
#' # Add some normalization/scaling factors to the phyloseq object
#' my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
#'     method = c("TMM", "median"))
#' ps_stool_16S <- runNormalizations(normalization_list = my_norm,
#'     object = ps_stool_16S)
#'
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ group, coef = 2,
#'     norm = c("TMM", "CSSmedian"))
#'
#' # Run methods on mock datasets
#' results <- runMocks(mocks = mocks, method_list = my_limma,
#'     object = ps_stool_16S)
#'
#' # Prepare results for Type I Error Control
#' TIEC_summary <- createTIEC(results)
#'
#' # Plot the results
#' plotFPR(df_FPR = TIEC_summary$df_FPR)
#' plotQQ(df_QQ = TIEC_summary$df_QQ, zoom = c(0, 0.1))
#' plotKS(df_KS = TIEC_summary$df_KS)
plotFPR <- function(df_FPR, cols = NULL) {
    # Create the vector of colors
    if (is.null(cols)) {
          cols <- createColors(variable = df_FPR$Method)
      }
    if (is.null(names(cols))) {
          stop("'cols' vector is not a named vector of colors.")
      }
    # Methods' ordering for FPR levels
    ord001 <- order(plyr::ddply(
        .data = df_FPR, .variables = ~Method, .fun =
            function(x) mean(x$FPR_obs001)
    )$V1)
    ord005 <- order(plyr::ddply(
        .data = df_FPR, .variables = ~Method, .fun =
            function(x) mean(x$FPR_obs005)
    )$V1)
    ord01 <- order(plyr::ddply(
        .data = df_FPR, .variables = ~Method, .fun =
            function(x) mean(x$FPR_obs01)
    )$V1)
    df_FPR001 <- df_FPR005 <- df_FPR01 <- df_FPR
    df_FPR001$Method <- factor(df_FPR$Method, levels = levels(df_FPR$Method)[
        ord001
    ], ordered = TRUE)
    df_FPR005$Method <- factor(df_FPR$Method, levels = levels(df_FPR$Method)[
        ord005
    ], ordered = TRUE)
    df_FPR01$Method <- factor(df_FPR$Method, levels = levels(df_FPR$Method)[
        ord01
    ], ordered = TRUE)
    ### Plot ###
    # to avoid notes during the check
    Method <- FPR_obs001 <- FPR_obs005 <- FPR_obs01 <- NULL
    g <- ggplot(data = df_FPR, aes(color = Method)) +
        geom_boxplot(data = df_FPR001, aes(
            x = "0.01", y =
                FPR_obs001
        )) +
        geom_boxplot(data = df_FPR005, aes(
            x = "0.05", y =
                FPR_obs005
        )) +
        geom_boxplot(data = df_FPR01, aes(
            x = "0.1", y =
                FPR_obs01
        )) +
        geom_segment(aes(
            x = 1 - 0.5, xend = 1 + 0.5, y =
                0.01, yend = 0.01
        ), color = "red", lty = 2) +
        geom_segment(aes(
            x = 2 - 0.5, xend = 2 + 0.5, y =
                0.05, yend = 0.05
        ), color = "red", lty = 2) +
        geom_segment(aes(
            x = 3 - 0.5, xend = 3 + 0.5, y = 0.1,
            yend = 0.1
        ), color = "red", lty = 2) +
        xlab(expression(Nominal ~ alpha)) +
        ylab(expression(Observed ~ alpha)) +
        ggtitle(
            label = "False Positive Rate", subtitle =
                "For each method"
        ) +
        scale_color_manual(values = cols)
    return(g)
} # END - function: plotFPR

#' @title plotQQ
#'
#' @export
#' @importFrom plyr ddply
#' @import ggplot2
#' @description
#' Draw the average QQ-plots across the mock comparisons.
#'
#' @param df_QQ Coordinates to draw the QQ-plot to compare the mean observed
#' p-value distribution across comparisons, with the theoretical uniform
#' distribution.
#' @param cols named vector of colors.
#' @param zoom 2-dimesional vector containing the starting and the
#' final coordinates (default: \code{c(0, 0.1)})
#'
#' @return A ggplot object.
#'
#' @examples
#' # Load some data
#' data(ps_stool_16S)
#'
#' # Generate the patterns for 10 mock comparison for an experiment
#' # (N = 1000 is suggested)
#' mocks <- createMocks(nsamples = phyloseq::nsamples(ps_stool_16S), N = 10)
#' head(mocks)
#'
#' # Add some normalization/scaling factors to the phyloseq object
#' my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
#'     method = c("TMM", "median"))
#' ps_stool_16S <- runNormalizations(normalization_list = my_norm,
#'     object = ps_stool_16S)
#'
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ group, coef = 2,
#'     norm = c("TMM", "CSSmedian"))
#'
#' # Run methods on mock datasets
#' results <- runMocks(mocks = mocks, method_list = my_limma,
#'     object = ps_stool_16S)
#'
#' # Prepare results for Type I Error Control
#' TIEC_summary <- createTIEC(results)
#'
#' # Plot the results
#' plotFPR(df_FPR = TIEC_summary$df_FPR)
#' plotQQ(df_QQ = TIEC_summary$df_QQ, zoom = c(0, 0.1))
#' plotKS(df_KS = TIEC_summary$df_KS)
plotQQ <- function(df_QQ, cols = NULL, zoom = c(0, 0.1)) {
    # Create the vector of colors
    if (is.null(cols)) {
          cols <- createColors(variable = df_QQ$Method)
      }
    if (is.null(names(cols))) {
          stop("'cols' vector is not a named vector of colors.")
      }
    # Create segment for the 0 to 0.1 zoom of the plot
    df_segments <- data.frame(x1 = c(0.01, 0, 0.05, 0, 0.1, 0), x2 = c(
        0.01,
        0.01, 0.05, 0.05, 0.1, 0.1
    ), y1 = c(0, 0.01, 0, 0.05, 0, 0.1), y2 = c(
        0.01, 0.01, 0.05, 0.05, 0.1, 0.1
    ))
    # to avoid notes during the check
    Method <- pval_theoretical <- pval_observed <- x1 <- x2 <- y1 <- y2 <- NULL
    g <- ggplot(data = df_QQ, aes(
        x = pval_theoretical,
        y = pval_observed, color = Method
    )) +
        geom_line() +
        geom_abline() +
        coord_cartesian(xlim = zoom, ylim = zoom) +
        xlab("Theoretical p-value") +
        ylab("Observed p-value") +
        ggtitle(label = "Average QQ-plot", subtitle = paste(
            "From",
            zoom[1], "to", zoom[2]
        )) +
        geom_segment(aes(
            x = x1, xend = x2, y = y1,
            yend = y2
        ), color = "red", lty = 2, data = df_segments) +
        scale_color_manual(values = cols)
    return(g)
}

#' @title plotKS
#'
#' @export
#' @importFrom plyr ddply
#' @import ggplot2
#' @description
#' Draw the boxplots of the Kolmogorov-Smirnov test statistics for the p-value
#' distributions across the mock comparisons.
#'
#' @param df_KS a \code{data.frame} produced by the \code{\link{createTIEC}}
#' function containing the KS statistics and their p-values.
#' @param cols named vector of colors.
#'
#' @return A ggplot object.
#'
#' @examples
#' # Load some data
#' data(ps_stool_16S)
#'
#' # Generate the patterns for 10 mock comparison for an experiment
#' # (N = 1000 is suggested)
#' mocks <- createMocks(nsamples = phyloseq::nsamples(ps_stool_16S), N = 10)
#' head(mocks)
#'
#' # Add some normalization/scaling factors to the phyloseq object
#' my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
#'     method = c("TMM", "median"))
#' ps_stool_16S <- runNormalizations(normalization_list = my_norm,
#'     object = ps_stool_16S)
#'
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ group, coef = 2,
#'     norm = c("TMM", "CSSmedian"))
#'
#' # Run methods on mock datasets
#' results <- runMocks(mocks = mocks, method_list = my_limma,
#'     object = ps_stool_16S)
#'
#' # Prepare results for Type I Error Control
#' TIEC_summary <- createTIEC(results)
#'
#' # Plot the results
#' plotFPR(df_FPR = TIEC_summary$df_FPR)
#' plotQQ(df_QQ = TIEC_summary$df_QQ, zoom = c(0, 0.1))
#' plotKS(df_KS = TIEC_summary$df_KS)
plotKS <- function(df_KS, cols = NULL) {
    # Create the vector of colors
    if (is.null(cols)) {
          cols <- createColors(variable = df_KS$Method)
      }
    if (is.null(names(cols))) {
          stop("'cols' vector is not a named vector of colors.")
      }
    ord <- order(plyr::ddply(
        .data = df_KS, .variables = ~Method, .fun =
            function(x) mean(x$KS)
    )$V1)
    df_KS$Method <- factor(df_KS$Method,
        levels = levels(df_KS$Method)[ord],
        ordered = TRUE
    )
    # to avoid notes during the check
    Method <- KS <- NULL
    g <- ggplot(df_KS, aes(
        x = Method, y = KS, color =
            Method
    )) +
        geom_boxplot() +
        coord_flip() +
        ylab("KS statistics") +
        xlab("Methods") +
        ggtitle(
            label = "K-S statistics", subtitle =
                "Ordered methods"
        ) +
        scale_x_discrete(
            limits = rev(levels(df_KS$Method)), position =
                "top"
        ) +
        theme(legend.position = "none", panel.spacing = unit(
            0,
            "lines"
        ), plot.margin = unit(c(0, 0, 0, 1), "cm")) +
        scale_color_manual(values = cols)
    return(g)
}
