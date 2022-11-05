#' @title plotFPR
#'
#' @export
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @importFrom tidytext reorder_within scale_x_reordered
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
#'     method = c("TMM", "CSS"))
#' ps_stool_16S <- runNormalizations(normalization_list = my_norm,
#'     object = ps_stool_16S)
#'
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ group, coef = 2,
#'     norm = c("TMM", "CSS"))
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
#' plotFDR(df_FDR = TIEC_summary$df_FDR)
#' plotQQ(df_QQ = TIEC_summary$df_QQ, zoom = c(0, 0.1))
#' plotKS(df_KS = TIEC_summary$df_KS)
#' plotLogP(df_QQ = TIEC_summary$df_QQ)
plotFPR <- function(df_FPR, cols = NULL) {
    # Create the vector of colors
    if (is.null(cols)) {
        cols <- createColors(variable = df_FPR$Method)
    }
    if (is.null(names(cols))) {
        stop("'cols' vector is not a named vector of colors.")
    }
    colnames(df_FPR) <- c("Comparison", "Method", "FPR 0.01", "FPR 0.05", 
        "FPR 0.1")
    df_FPR_melted <- reshape2::melt(df_FPR)
    ### Plot ###
    # to avoid notes during the check
    Method <- value <- variable <- y <- NULL
    data_segm <- data.frame("y" = c(0.01, 0.05, 0.1), 
        "variable" = c("FPR 0.01", "FPR 0.05", "FPR 0.1"))
    g <- ggplot(data = df_FPR_melted, 
        aes(color = Method, 
            x = tidytext::reorder_within(Method, value, variable, 
            fun = stats::median), y = value)) +
        geom_boxplot() +
        facet_grid(~ variable, scales = "free_x") +
        tidytext::scale_x_reordered() +
        geom_hline(data = data_segm, 
            aes(yintercept = y), 
            color = "red", lty = 2) +
        xlab("Method") +
        ylab(expression(Observed ~ alpha)) +
        ggtitle(label = "False Positive Rate", subtitle = "by FPR level") +
        scale_color_manual(values = cols) + 
        theme(legend.position = "none", axis.text.x = element_text(
            angle = 90, hjust = 1, vjust = 0.5))
    return(g)
} # END - function: plotFPR

#' @title plotFDR
#'
#' @export
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @importFrom tidytext reorder_within scale_x_reordered
#' @import ggplot2
#' @description
#' Draw the nominal false discovery rates for the 0.01, 0.05, and 0.1 levels.
#'
#' @param df_FDR a \code{data.frame} produced by the \code{\link{createTIEC}}
#' function, containing the FDR values.
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
#'     method = c("TMM", "CSS"))
#' ps_stool_16S <- runNormalizations(normalization_list = my_norm,
#'     object = ps_stool_16S)
#'
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ group, coef = 2,
#'     norm = c("TMM", "CSS"))
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
#' plotFDR(df_FDR = TIEC_summary$df_FDR)
#' plotQQ(df_QQ = TIEC_summary$df_QQ, zoom = c(0, 0.1))
#' plotKS(df_KS = TIEC_summary$df_KS)
#' plotLogP(df_QQ = TIEC_summary$df_QQ)
plotFDR <- function(df_FDR, cols = NULL) {
    # Create the vector of colors
    if (is.null(cols)) {
        cols <- createColors(variable = df_FDR$Method)
    }
    if (is.null(names(cols))) {
        stop("'cols' vector is not a named vector of colors.")
    }
    colnames(df_FDR) <- c("Method", "FDR 0.01", "FDR 0.05", "FDR 0.1")
    df_FDR_melted <- reshape2::melt(df_FDR)
    ### Plot ###
    # to avoid notes during the check
    Method <- value <- variable <- y <- NULL
    data_segm <- data.frame("y" = c(0.01, 0.05, 0.1), 
        "variable" = c("FDR 0.01", "FDR 0.05", "FDR 0.1"))
    g <- ggplot(data = df_FDR_melted, aes(color = Method, 
            x = tidytext::reorder_within(Method, value, variable, 
                fun = stats::median), 
            y = value)) +
        geom_point() +
        facet_grid(~ variable, scales = "free_x") +
        tidytext::scale_x_reordered() +
        geom_hline(data = data_segm, 
                   aes(yintercept = y), 
                   color = "red", lty = 2) +
        xlab("Method") +
        ylab("Average FDR") +
        ggtitle(label = "False Discovery Rate", 
            subtitle = "by significance thresholds") +
        scale_color_manual(values = cols) + 
        theme(legend.position = "none", axis.text.x = element_text(
            angle = 90, hjust = 1, vjust = 0.5))
    return(g)
} # END - function: plotFDR

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
#' @param split boolean value. If \code{TRUE}, the qq-plots are 
#' reported separately for each method (default \code{split = FALSE}). Setting 
#' it to \code{TRUE} is hardly suggested when the number of methods is high or 
#' when their colors are similar. 
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
#'     method = c("TMM", "CSS"))
#' ps_stool_16S <- runNormalizations(normalization_list = my_norm,
#'     object = ps_stool_16S)
#'
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ group, coef = 2,
#'     norm = c("TMM", "CSS"))
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
#' plotFDR(df_FDR = TIEC_summary$df_FDR)
#' plotQQ(df_QQ = TIEC_summary$df_QQ, zoom = c(0, 0.1))
#' plotKS(df_KS = TIEC_summary$df_KS)
#' plotLogP(df_QQ = TIEC_summary$df_QQ)
plotQQ <- function(df_QQ, cols = NULL, zoom = c(0, 0.1), split = FALSE) {
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
    g <- ggplot(data = df_QQ, aes(x = pval_theoretical, y = pval_observed, 
        color = Method)) +
        geom_smooth(se = FALSE, na.rm = TRUE, span = 0.3, method = "loess", 
            formula = y ~ x) +
        geom_abline() +
        coord_equal(xlim = zoom, ylim = zoom) +
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
        if(split){
            g <- g + 
                facet_wrap(~ Method) +
                theme(legend.position = "none")
        }
    return(g)
}

#' @title plotLogP
#'
#' @export
#' @importFrom plyr ddply
#' @import ggplot2
#' @importFrom ggridges geom_density_ridges2 position_raincloud
#' @importFrom stats runif median p.adjust
#' @description
#' Draw the p-values or the average p-values distribution across the mock 
#' comparisons in logarithmic scale.
#'
#' @param df_pval a \code{data.frame} produced by the \code{\link{createTIEC}}
#' function, containing the p-values for each taxon, method, and mock 
#' comparison. It is used to draw the negative log10 p-values distribution.
#' If df_pval is supplied, let \code{df_QQ = NULL}.
#' @param df_QQ a \code{data.frame} produced by the \code{\link{createTIEC}}
#' function, containing the average p-values for each quantile and method. It is
#' used to draw the negative log10 average p-values distribution. If df_QQ is 
#' supplied, let \code{df_pval = NULL}.
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
#'     method = c("TMM", "CSS"))
#' ps_stool_16S <- runNormalizations(normalization_list = my_norm,
#'     object = ps_stool_16S)
#'
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ group, coef = 2,
#'     norm = c("TMM", "CSS"))
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
#' plotFDR(df_FDR = TIEC_summary$df_FDR)
#' plotQQ(df_QQ = TIEC_summary$df_QQ, zoom = c(0, 0.1))
#' plotKS(df_KS = TIEC_summary$df_KS)
#' plotLogP(df_QQ = TIEC_summary$df_QQ)
plotLogP <- function(df_pval = NULL, df_QQ = NULL, cols = NULL) {
    if(is.null(df_pval) & is.null(df_QQ)){
        stop("Please supply 'df_pval' for all p-values analysis or ", 
            " 'df_QQ' for average p-values analysis.")
    }
    if(!is.null(df_pval) & !is.null(df_QQ)){
        stop("Please supply one between 'df_pval' for all p-values analysis", 
            " or 'df_QQ' for average p-values analysis.")
    }
    # Add IDEAL distribution
    p <- stats::runif(10000)
    if(!is.null(df_pval)){
        gtitle <- "All p-values log distribution"
        df_pval_ideal <- data.frame(Comparison = "Comparison IDEAL", 
            Method = "IDEAL", variable = "IDEAL", pval = p, 
            padj = stats::p.adjust(p, method = "BH"))
        df_to_plot <- rbind(df_pval, df_pval_ideal)
        # Create the vector of colors
        if (is.null(cols)) {
            cols <- createColors(variable = df_pval$Method)
        }
    }
    if(!is.null(df_QQ)){
        gtitle <- "Average p-values log distribution"
        df_pval_ideal <- data.frame(Method = "IDEAL", 
            pval_theoretical = p, pval_observed = p)
        df_to_plot <- rbind(df_QQ, df_pval_ideal)
        colnames(df_to_plot) <- c("Method", "pval_theoretical", "pval")
        # Create the vector of colors
        if (is.null(cols)) {
            cols <- createColors(variable = df_QQ$Method)
        }
    }
    if (is.null(names(cols))) {
        stop("'cols' vector is not a named vector of colors.")
    }
    # Order method by closeness to IDEAL method
    ord <- order(plyr::ddply(.data = df_to_plot, .variables = ~ Method, .fun =
        function(x){
            sum(stats::quantile(-log10(x$pval), c(0.99, 0.95, 0.9), 
                na.rm = TRUE) + log10(c(0.01, 0.05, 0.1)))
        })$V1, decreasing = TRUE)
    ordered_levels <- levels(df_to_plot$Method)[ord]
    IDEAL_index <- which(ordered_levels == "IDEAL")
    IDEAL_first <- c(ordered_levels[-IDEAL_index], "IDEAL")
    df_to_plot$Method <- factor(df_to_plot$Method,
        levels = IDEAL_first,
        ordered = TRUE)
    
    bold <- "IDEAL"
    bold.labels <- ifelse(is.element(levels(df_to_plot[, "Method"]), bold), 
        yes = "bold", no = "plain")
    color.labels <- ifelse(is.element(levels(df_to_plot[, "Method"]), bold), 
        yes = "red", no = "black")
    
    # to avoid notes during the check
    Method <- pval <-  NULL
    g <- ggplot(data = df_to_plot, 
        aes(x = -log10(pval), y = Method, fill = Method)) +
        # Add histogram
        ggridges::geom_density_ridges2(stat = "binline", binwidth = 0.1, 
            scale = 0.9, alpha = 0.5, draw_baseline = FALSE, na.rm = TRUE, 
            boundary = TRUE) +
        # Add density
        ggridges::geom_density_ridges(aes(vline_color = stat(quantile)), 
            alpha = 0.2, scale = 0.9, na.rm = TRUE, quantile_lines = TRUE, 
            quantiles = c(0.9, 0.95, 0.99), vline_size = 1, from = 0, 
            bandwidth = 0.1,
            position = ggridges::position_raincloud(adjust_vlines = TRUE)) +
        geom_vline(xintercept = -log10(c(0.1, 0.05, 0.01)), lty = 2) +
        ylab("Method") +
        xlab(expression(-log[10](p-value))) +
        scale_x_continuous(sec.axis = sec_axis(trans = ~., name = "p-value",
            breaks = -log10(c(0.1, 0.05, 0.01)), labels = c(0.1, 0.05, 0.01))) +
        ggtitle(label = gtitle, subtitle = 
            paste0("IDEAL method as reference, red bars represents the 90,", 
                " 95, and 99 percentiles")) +
        scale_fill_manual(values = c(cols, "IDEAL" = "#ff0000")) +
        scale_color_manual(aesthetics = "vline_color",
            values = c("#ffbaba", "#ff7b7b", "#ff5252", "#ff0000")) +
        theme(legend.position = "none",
            axis.text.y = suppressWarnings(element_text(face = bold.labels, 
                colour = color.labels, vjust = 0)),
            axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5))
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
#'     method = c("TMM", "CSS"))
#' ps_stool_16S <- runNormalizations(normalization_list = my_norm,
#'     object = ps_stool_16S)
#'
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ group, coef = 2,
#'     norm = c("TMM", "CSS"))
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
#' plotFDR(df_FDR = TIEC_summary$df_FDR)
#' plotQQ(df_QQ = TIEC_summary$df_QQ, zoom = c(0, 0.1))
#' plotKS(df_KS = TIEC_summary$df_KS)
#' plotLogP(df_QQ = TIEC_summary$df_QQ)
plotKS <- function(df_KS, cols = NULL) {
    # Create the vector of colors
    if (is.null(cols)) {
          cols <- createColors(variable = df_KS$Method)
      }
    if (is.null(names(cols))) {
          stop("'cols' vector is not a named vector of colors.")
      }
    ord <- order(plyr::ddply(
        .data = df_KS, .variables = ~ Method, .fun =
            function(x) median(x$KS)
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
                "Methods ordered by median K-S"
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
