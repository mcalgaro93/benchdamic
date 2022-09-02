#' @title plotConcordanceHeatmap
#'
#' @keywords internal
#' @import ggplot2
#' @description
#' Plots the heatmap of concordances.
#'
#' @param c_df A simplified concordance \code{data.frame} produced in
#' \code{\link{plotConcordance}} function.
#' @param threshold The threshold for rank (x-axis upper limit if all methods
#' have a higher number of computed statistics).
#' @param cols A named vector containing the color hex codes.
#'
#' @return a \code{ggplot2} object
#'
#' @seealso \code{\link{createConcordance}} and \code{\link{plotConcordance}}

plotConcordanceHeatmap <- function(c_df, threshold, cols) {
    concordance <- max_areaOver <- method1 <- method2 <- NULL
    # If we consider the trapezoid with bases B = 1, b = 1-threshold/n_features
    # and height h = threshold, its area will be (B+b)*h/2.
    # Then rescale the area by n_features in order to make it from 0 to 0.5
    n_features <- max(c_df[, "n_features"])
    rescaled_area <- (2 - threshold / n_features) * threshold / 2 / n_features
    ggplot(data = c_df, mapping = aes(
        x = rank, y =
            concordance
    )) +
        facet_grid(method1 ~ method2, scales = "free_x", switch = "y") +
        geom_ribbon(aes(
            ymin = rank / n_features, ymax = concordance,
            fill = max_areaOver
        )) +
        geom_line() +
        geom_line(
            mapping = aes(x = rank, y = rank / n_features),
            lty = 2
        ) +
        theme(
            legend.position = "bottom",
            axis.line.x.bottom = element_blank(),
            axis.line.y.right = element_blank(),
            strip.text = element_blank(),
            strip.background = element_blank(),
            panel.spacing = unit(0.1, "cm"),
            plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
        ) +
        coord_cartesian(
            xlim = c(1, threshold), ylim = c(0, 1),
            expand = FALSE
        ) +
        geom_segment(data = c_df[c_df[, "method1"] == levels(c_df[, "method1"])[
            1
        ], ], mapping = aes(
            x = 1, xend = threshold, y = 1,
            yend = 1, size = 1, color = method2
        ), show.legend = FALSE) +
        geom_segment(data = c_df[c_df[, "method2"] == levels(c_df[, "method1"])[
            1
        ], ], mapping = aes(
            x = 1, xend = 1, y = 0,
            yend = 1, size = 1, color = method1
        ), show.legend = FALSE) +
        geom_rect(
            data = c_df[c_df[, "method1"] == c_df[, "method2"], ],
            mapping = aes(
                xmin = 1, xmax = threshold, ymin = 0, ymax =
                    1
            ), color = "red", fill = NA, show.legend = FALSE
        ) +
        scale_fill_distiller(palette = "RdYlBu", limits = c(
            0,
            rescaled_area
        ), direction = 1) +
        scale_color_manual(values = cols) +
        scale_y_continuous(breaks = c(0.5, 1), position = "right") +
        scale_x_continuous(limits = c(1, threshold)) +
        xlab("Rank") +
        ylab("Average concordance") +
        guides(
            fill = guide_colorbar(
                title =
                    "Average area over the bisector:", direction = "horizontal",
                ticks.colour = "black", frame.colour = "black"
            ),
            colour = guide_legend(NULL)
        )
}

#' @title plotConcordanceDendrogram
#'
#' @keywords internal
#' @importFrom ggdendro dendro_data label
#' @import ggplot2
#' @description
#' Plots the method's dendrogram of concordances.
#'
#' @param hc Hierarchical clustering results produced in
#' \code{\link{plotConcordance}} function.
#' @param direction vertical (default \code{direction = "v"}) or horizontal
#' (\code{direction = "h"}).
#' @param cols A named vector containing the color hex codes.
#'
#' @return a \code{ggplot2} object
#'
#' @seealso \code{\link{createConcordance}} and \code{\link{plotConcordance}}

plotConcordanceDendrogram <- function(hc, direction = "v", cols) {
    x <- y <- xend <- yend <- NULL
    g_dendro <- ggplot() +
        geom_segment(
            data = ggdendro::dendro_data(hc)[["segments"]],
            mapping = aes(x = x, y = y, xend = xend, yend = yend)
        ) +
        theme(
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            panel.background = element_rect(fill = "white"),
            panel.grid = element_blank(),
            legend.position = "none",
            panel.spacing = unit(0, "lines"),
            plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
        )
    if (direction == "v") {
        g_dendro <- g_dendro +
            geom_label(
                data = ggdendro::dendro_data(hc)[["labels"]],
                mapping = aes(
                    x = x, y = y, label = label, hjust = 1,
                    fill = label
                ), nudge_y = 0
            ) +
            coord_flip() +
            scale_y_reverse(expand = c(0, 0, 0, 0)) +
            scale_x_reverse() +
            scale_fill_manual(values = cols)
    } else if (direction == "h") {
        g_dendro <- g_dendro +
            geom_point(
                data = ggdendro::dendro_data(hc)[["labels"]],
                mapping = aes(x = x, y = y, color = label), nudge_y =
                    0.1, size = 5
            ) +
            scale_y_continuous(expand = c(0, 0, 0, 0)) +
            scale_x_continuous(expand = c(0, 0, 0, 0)) +
            scale_color_manual(values = cols)
    } else {
        stop("Please supply a direction value between 'h' for horizontal or
                'v' for vertical.")
    }
    return(g_dendro)
}

#' @title plotConcordance
#'
#' @export
#' @importFrom plyr ddply
#' @importFrom reshape2 dcast
#' @importFrom stats hclust as.dist
#' @importFrom cowplot plot_grid
#' @description
#' Produce a list of graphical outputs summarizing the between and
#' within method concordance.
#'
#' @param concordance A long format \code{data.frame} produced by
#' \code{\link{createConcordance}} function.
#' @inheritParams plotConcordanceHeatmap
#'
#' @return A 2 elements list of \code{ggplot2} class objects:
#' \itemize{
#'     \item{\code{concordanceDendrogram}}{ which contains the
#'     vertically directioned dendrogram for the methods involved in the
#'     concordance analysis;}
#'     \item{\code{concordanceHeatmap}}{ which contains the heatmap of between
#'     and within method concordances.}}
#'
#' @seealso \code{\link{createConcordance}}
#'
#' @examples
#' data(ps_plaque_16S)
#' 
#' # Balanced design
#' my_splits <- createSplits(
#'     object = ps_plaque_16S, varName = "HMP_BODY_SUBSITE", balanced = TRUE,
#'     paired = "RSID", N = 10 # N = 100 suggested
#' )
#' 
#' # Make sure the subject ID variable is a factor
#' phyloseq::sample_data(ps_plaque_16S)[, "RSID"] <- as.factor(
#'     phyloseq::sample_data(ps_plaque_16S)[["RSID"]])
#'
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ RSID + HMP_BODY_SUBSITE, 
#'     coef = "HMP_BODY_SUBSITESupragingival Plaque",
#'     norm = c("TMM", "CSS"))
#'
#' # Set the normalization methods according to the DA methods
#' my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
#'     method = c("TMM", "CSS"))
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
#' # plot concordances from rank 1 to 50.
#' plotConcordance(
#'     concordance = concordance_pvalues,
#'     threshold = 50
#' )
plotConcordance <- function(concordance, threshold = NULL, cols = NULL) {
    area_df <- areaCAT(concordance = concordance, plotIt = FALSE)
    # Maximum rank for each method and comparison
    maxRank <- plyr::ddply(.data = area_df, .variables = ~ comparison +
        method1 + method2, .fun = function(x) {
        c(maxRank = max(x[, "rank"]))
    })
    # The maximum threshold allowed for each pair of methods is the minimum
    # maxRank
    threshold_df <- plyr::ddply(.data = maxRank, .variables = ~ method1 +
        method2, .fun = function(x) {
        c(threshold = min(x[, "maxRank"]))
    })
    # Common threshold
    common_threshold <- min(threshold_df[, "threshold"])
    # cat(common_threshold,"\n")
    if (is.null(threshold)) {
        threshold <- common_threshold
    } else {
        if (threshold > common_threshold) {
            message(threshold, "\n")
            threshold <- common_threshold
            message(threshold, "\n")
        }
    }
    # Keep only ranks that are lower than threshold
    concordance_df <- area_df[area_df[, "rank"] <= threshold, ]
    # Average values for plotting
    c_df <- plyr::ddply(.data = concordance_df, .variables = ~ method1 +
        method2 + rank, .fun = function(x) {
        colMeans(x[, c("n_features", "concordance", "areaOver")])
    })
    # Extract the area at the maximum rank (= threshold)
    max_areaOver <- plyr::ddply(c_df,
        .variables = ~ method1 + method2, .fun =
            function(x) c(areaOver = max(x[, "areaOver"]))
    )
    # Repeat the max_areaOver for all ranks in pair of method1 and method2
    c_df[, "max_areaOver"] <- rep(max_areaOver[, "areaOver"], each = threshold)
    # Clustering of methods
    # Organize area values in a n_methods x n_methods matrix
    dist_df <- reshape2::dcast(
        data = max_areaOver, formula = method1 ~ method2,
        value.var = "areaOver"
    )
    dist_df <- dist_df[, 2:ncol(dist_df)]
    rownames(dist_df) <- colnames(dist_df)
    # Compute distances between methods
    distances <- stats::as.dist(1 - dist_df)
    hc <- stats::hclust(d = distances)
    # Ordering the levels of method1 and method2 variables
    c_df[, "method1"] <- factor(c_df[, "method1"], levels = colnames(dist_df)[hc
    [["order"]]])
    c_df[, "method2"] <- factor(c_df[, "method2"], levels = colnames(dist_df)[hc
    [["order"]]])
    if (is.null(cols)) {
          cols <- createColors(hc[["labels"]])
      }
    g_list <- list(
        "concordanceDendrogram" = plotConcordanceDendrogram(
            hc = hc,
            direction = "v", cols = cols
        ),
        "concordanceHeatmap" = plotConcordanceHeatmap(
            c_df = c_df, threshold =
                threshold, cols = cols
        )
    )
    return(g_list)
}
