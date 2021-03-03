#' @title createTIEC
#'
#' @export
#' @importFrom plyr ldply ddply
#' @importFrom stats ks.test
#' @description
#' Extract the list of p-values from the outputs of the differential abundance
#' detection methods to compute several statistics to study the ability to
#' control the type I error.
#'
#' @param object Output of the differential abundance tests on mock comparisons.
#' Must follow a specific structure with comparison, method, matrix of
#' p-values, and method's name (See vignette for detailed information).
#'
#' @return A \code{list} of \code{data.frame}s:
#' \itemize{
#'  \item{\code{df_pval}}{3 columns per number_of_features x methods x
#'  comparisons rows data.frame. The three columns are called Comparison, pval,
#'  and method;}
#'  \item{\code{df_FPR}}{5 columns per methods x comparisons rows data.frame.
#'  For each set of method and comparison, the proportion of false discoveries,
#'  considering 3 threshold (0.01, 0.05, 0.1) are reported;}
#'  \item{\code{df_QQ}}{contains the coordinates to draw the QQ-plot to compare
#'  the mean observed p-value distribution across comparisons, with the
#'  theoretical uniform distribution;}
#'  \item{\code{df_KS}}{5 columns and methods x comparisons rows data.frame. For
#'  each set of method and comparison, the Kolmogorov-Smirnov test statistics
#'  and p-values are reported in KS and KS_pval columns respectively.}
#' }
#'
#' @seealso [benchdamic::createMocks()]

createTIEC <- function(object){
    # Create a list of data frames with two columns: pval and method name
    # One data frame for each comparison
    df_list_pval <- lapply(X = object, FUN = function(Comparison){
        plyr::ldply(Comparison, function(method){

          # Check objects and column names
          if(!is.element("pValMat", names(method)))
            stop("'pValMat' matrix not found for one of the methods in input.")
          if(!is.element("name", names(method)))
            stop("'name' character not found for one of the methods in input.")
          if(!is.element("rawP", colnames(method$pValMat)))
            stop("'rawP' column not found in 'pValMat' matrix.")

            data.frame(pval = method$pValMat[, "rawP"],
                       Method = factor(method$name))

        },.id = "Method")
    })
    # Melt down to a single data frame with the Comparison column added
    df_pval <- plyr::ldply(df_list_pval, .id = "Comparison")

    ### FPR ###
    # Count the p-values which are lower than a selected threshold
    df_pval_FPR <- plyr::ddply(.data = df_pval,
                               .variables = ~ Comparison + Method,
                               .fun = function(x){
        # index for not NA p-values
        k_pval <- !is.na(x$pval)
        x$FPR_obs001 <- sum(x$pval < 0.01, na.rm = TRUE) / sum(k_pval)
        x$FPR_obs005 <- sum(x$pval < 0.05, na.rm = TRUE) / sum(k_pval)
        x$FPR_obs01 <- sum(x$pval < 0.1, na.rm = TRUE) / sum(k_pval)
        return(x)
    })

    # Compute the mean for each FPR threshold
    df_FPR <- plyr::ddply(.data = df_pval_FPR,
                          .variables = ~ Comparison + Method,
                          .fun = function(x){
        colMeans(x[,c("FPR_obs001", "FPR_obs005", "FPR_obs01")])
    })

    ### QQ and KS ###
    # Create data frame for empirical and theoretical p-values
    # Kolmogorov-Smirnov test
    df_QQ_KS <- suppressWarnings(plyr::ddply(.data = df_pval,
                                       .variables = ~ Comparison + Method,
                                       .fun = function(x){
        # Ordered list of p-values
        x <- x[order(x$pval),]
        # Index of not NA p-values
        k_pval <- !is.na(x$pval)
        # Theoretical uniform distribution
        x$pval_theoretical[k_pval] <- (seq_len(sum(k_pval))) / (sum(k_pval))
        x$pval_theoretical_rounded[k_pval] <- round((seq_len(sum(k_pval))) /
                                                      (sum(k_pval)), digits = 2)
        x$KS <- stats::ks.test(x = x$pval_theoretical[k_pval],
                               y = x$pval[k_pval])$statistic
        x$KS_pval <- stats::ks.test(x = x$pval_theoretical[k_pval],
                                    y = x$pval[k_pval])$p.value
        return(x)
    }))

    ### QQ ###
    # For each theoretical p-value, the mean of observed one is computed
    df_QQ <- plyr::ddply(.data = df_QQ_KS,
                         .variables = ~ Method + pval_theoretical_rounded,
                         .fun = function(x){
        pval_mean <- mean(x$pval)
        return(pval_mean)
    })
    names(df_QQ) <- c("Method", "pval_theoretical", "pval_observed")

    ### KS ###
    # Compute the mean for KS statistics and p-values
    df_KS <- plyr::ddply(.data = df_QQ_KS,
                         .variables = ~ Comparison + Method,
                         .fun = function(x){
        colMeans(x[, c("KS", "KS_pval")])})

    return(list(df_pval = df_pval,
                df_FPR = df_FPR,
                df_QQ = df_QQ,
                df_KS = df_KS))
}# END - function: createTIEC

#' @title createColors
#'
#' @export
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @description
#' Produce a qualitative set of colors.
#'
#' @param variable character vector or factor variable.
#'
#' @return A named vector containing the color codes.

createColors <- function(variable){
    if(is.factor(variable)){
        levels <- levels(variable)
    } else levels <- unique(variable)

    n_levels <- length(levels)

    if(n_levels > 74)
        stop("To many levels to color them differently.")

    pal.info <- RColorBrewer::brewer.pal.info
    qual_col_pals = pal.info[pal.info$category == 'qual',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal,
                               qual_col_pals$maxcolors,
                               rownames(qual_col_pals)))
    cols <- col_vector[seq_len(n_levels)]
    names(cols) <- levels
    return(cols)
}# END - function: createColors

#' @title plotFPR
#'
#' @export
#' @importFrom plyr ddply
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_segment xlab ylab ggtitle
#' @importFrom ggplot2 scale_color_manual
#' @description
#' Draw the boxplots of the proportions of p-values lower than 0.01, 0.05, and
#' 0.1 thresholds for each method.
#'
#' @param df_FPR 5 columns per methods x comparisons rows data.frame produced by
#' the [bechdamic::createTIEC()] function, containing the FPR values.
#' @param cols named vector of colors.
#'
#' @return A ggplot object.

plotFPR <- function(df_FPR,
                    cols = NULL){
    # Create the vector of colors
    if(is.null(cols))
        cols <- createColors(variable = df_FPR$Method)
    if(is.null(names(cols)))
        stop("'cols' vector is not a named vector of colors.")

    # Methods' ordering for FPR levels
    ord001 <- order(plyr::ddply(df_FPR,
                                .variables = ~ Method,
                                function(x) mean(x$FPR_obs001))$V1)
    ord005 <- order(plyr::ddply(df_FPR,
                                .variables = ~ Method,
                                function(x) mean(x$FPR_obs005))$V1)
    ord01 <- order(plyr::ddply(df_FPR,
                               .variables = ~ Method,
                               function(x) mean(x$FPR_obs01))$V1)

    df_FPR001 <- df_FPR005 <- df_FPR01 <- df_FPR

    df_FPR001$Method <- factor(df_FPR$Method,
        levels = levels(df_FPR$Method)[ord001],ordered = TRUE)

    df_FPR005$Method <- factor(df_FPR$Method,
        levels = levels(df_FPR$Method)[ord005],ordered = TRUE)

    df_FPR01$Method <- factor(df_FPR$Method,
        levels = levels(df_FPR$Method)[ord01],ordered = TRUE)

    ### Plot ###
    # to avoid notes during the check
    Method <- FPR_obs001 <- FPR_obs005 <- FPR_obs01 <- NULL

    g <- ggplot2::ggplot(data = df_FPR, ggplot2::aes(color = Method)) +
            ggplot2::geom_boxplot(data = df_FPR001,
                                  ggplot2::aes(x = "0.01",y = FPR_obs001)) +
            ggplot2::geom_boxplot(data = df_FPR005,
                                  ggplot2::aes(x = "0.05",y = FPR_obs005)) +
            ggplot2::geom_boxplot(data = df_FPR01,
                                  ggplot2::aes(x = "0.1",y = FPR_obs01)) +
            ggplot2::geom_segment(ggplot2::aes(x = 1 - 0.5,
                                               xend = 1 + 0.5,
                                               y = 0.01,
                                               yend = 0.01),
                                  color = "red", lty = 2) +
            ggplot2::geom_segment(ggplot2::aes(x = 2 - 0.5,
                                               xend = 2 + 0.5,
                                               y = 0.05,
                                               yend = 0.05),
                                  color = "red", lty = 2) +
            ggplot2::geom_segment(ggplot2::aes(x = 3 - 0.5,
                                               xend = 3 + 0.5,
                                               y = 0.1,
                                               yend = 0.1),
                                  color = "red", lty = 2) +
            ggplot2::xlab(expression(Nominal ~ alpha)) +
            ggplot2::ylab(expression(Observed ~ alpha)) +
            ggplot2::ggtitle(label = "False Positive Rate",
                             subtitle = "For each method") +
            ggplot2::scale_color_manual(values = cols)

    return(g)
} # END - function: plotFPR

#' @title plotQQ
#'
#' @export
#' @importFrom plyr ddply
#' @importFrom ggplot2 ggplot aes geom_line geom_segment xlab ylab ggtitle
#' @importFrom ggplot2 geom_abline coord_cartesian scale_color_manual
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

plotQQ <- function(df_QQ,
                   cols = NULL,
                   zoom = c(0,0.1)){
    # Create the vector of colors
    if(is.null(cols))
        cols <- createColors(variable = df_QQ$Method)
    if(is.null(names(cols)))
        stop("'cols' vector is not a named vector of colors.")

    # Create segment for the 0 to 0.1 zoom of the plot
    df_segments <- data.frame(x1 = c(0.01, 0, 0.05, 0, 0.1, 0),
                              x2 = c(0.01, 0.01, 0.05, 0.05, 0.1, 0.1),
                              y1 = c(0, 0.01, 0, 0.05, 0, 0.1),
                              y2 = c(0.01, 0.01, 0.05, 0.05, 0.1, 0.1))

    # to avoid notes during the check
    Method <- pval_theoretical <- pval_observed <- x1 <- x2 <- y1 <- y2 <- NULL

    g <- ggplot2::ggplot(data = df_QQ, ggplot2::aes(x = pval_theoretical,
                                                    y = pval_observed,
                                                    color = Method)) +
        ggplot2::geom_line() +
        ggplot2::geom_abline() +
        ggplot2::coord_cartesian(xlim = zoom, ylim = zoom) +
        ggplot2::xlab("Theoretical p-value") +
        ggplot2::ylab("Observed p-value") +
        ggplot2::ggtitle(label = "Average QQ-plot",
                         subtitle = paste("From", zoom[1], "to", zoom[2])) +
        ggplot2::geom_segment(ggplot2::aes(x = x1, xend = x2,
                                           y = y1, yend = y2),
                              color = "red",
                              lty = 2,
                              data = df_segments) +
        ggplot2::scale_color_manual(values = cols)
    return(g)
}

#' @title plotKS
#'
#' @export
#' @importFrom plyr ddply
#' @importFrom ggplot2 ggplot aes geom_boxplot coord_flip xlab ylab ggtitle
#' @importFrom ggplot2 scale_x_discrete theme scale_color_manual
#' @description
#' Draw the boxplots of the Kolmogorov-Smirnov test statistics for the p-value
#' distributions across the mock comparisons.
#'
#' @param df_KS 5 columns per methods x comparisons rows data.frame produced by
#' the [benchdamic::createTIEC] function.
#' @param cols named vector of colors.
#'
#' @return A ggplot object.

plotKS <- function(df_KS,
                    cols = NULL){

    # Create the vector of colors
    if(is.null(cols))
        cols <- createColors(variable = df_KS$Method)
    if(is.null(names(cols)))
        stop("'cols' vector is not a named vector of colors.")

    ord <- order(plyr::ddply(.data = df_KS,
                             .variables = ~ Method,
                             .fun = function(x) mean(x$KS))$V1)

    df_KS$Method <- factor(df_KS$Method,
                           levels = levels(df_KS$Method)[ord],
                           ordered = TRUE)

    # to avoid notes during the check
    Method <- KS <- NULL

    g <- ggplot2::ggplot(df_KS,
                         ggplot2::aes(x = Method,y = KS, color = Method)) +
        ggplot2::geom_boxplot() +
        ggplot2::coord_flip() +
        ggplot2::ylab("KS statistics") + ggplot2::xlab("Methods") +
        ggplot2::ggtitle(label = "K-S statistics", subtitle =
                             "Ordered methods") +
        ggplot2::scale_x_discrete(limits = rev(levels(df_KS$Method)),
                                  position = "top") +
        ggplot2::theme(legend.position = "none",
                       panel.spacing = unit(0, "lines"),
                       plot.margin = unit(c(0,0,0,1), "cm")) +
        ggplot2::scale_color_manual(values = cols)

    return(g)
}
