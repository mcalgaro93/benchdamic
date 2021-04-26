#' @title getStatistics
#'
#' @export
#' @description
#' Extract the list of p-values or/and log fold changes from the output of a
#' differential abundance detection method.
#'
#' @param method Output of a differential abundance detection method.
#' \code{pValMat}, \code{statInfo} matrices, and method's \code{name} must be
#' present (See vignette for detailed information).
#' @param slot The slot name where to extract values
#' (default \code{slot = "pValMat"}).
#' @param colName The column name of the slot where to extract values
#' (default \code{colName = "rawP"}).
#' @param type The value type of the column selected where to extract values.
#' Two values are possible: \code{"pvalue"} or \code{"logfc"}
#' (default \code{type = "pvalue"}).
#' @param direction \code{statInfo}'s column name containing information about
#' the signs of differential abundance (usually log fold changes)
#' (default \code{direction = NULL}).
#' @param verbose Boolean to display the kind of extracted values
#' (default \code{verbose = FALSE}).
#'
#' @return A vector or a \code{data.frame}. If \code{direction = NULL},
#' the \code{colname} column values, transformed according to \code{type} (not
#' tranformed if \code{type = "pvalue"}, \code{-abs(value)} if
#' \code{type = "logfc"}), of the \code{slot} are reported, otherwise the
#' \code{direction} column of the \code{statInfo} matrix is added to the output.
#'
#' @seealso \code{\link{extractStatistics}}
#'
#' @examples
#' data("ps_plaque_16S")
#'
#' # Add scaling factors
#' ps_plaque_16S <- norm_edgeR(object = ps_plaque_16S, method = "TMM")
#' # DA analysis
#' da.limma <- DA_limma(
#'     object = ps_plaque_16S,
#'     design = ~ 1 + HMP_BODY_SUBSITE,
#'     coef = 2,
#'     norm = "TMM"
#' )
#' # get p-values
#' getStatistics(method = da.limma, slot = "pValMat", colName = "rawP",
#'     type = "pvalue", direction = NULL)
#' # get negative abs(logFC) values
#' getStatistics(method = da.limma, slot = "statInfo", colName = "logFC",
#'     type = "logfc", direction = NULL)
#' # get p-values and logFC
#' getStatistics(method = da.limma, slot = "pValMat", colName = "rawP",
#'     type = "pvalue", direction = "logFC")

getStatistics <- function(method, slot = "pValMat", colName = "rawP",
    type = "pvalue", direction = NULL, verbose = FALSE){
    # Info extraction
    method_name <- method[["name"]]
    if(!is.element(slot, names(method))){
        stop(paste0("'", slot,"' slot not found for ", method_name))
    } else info <- method[[slot]]
    if(!is.element(colName, colnames(info))){
        stop(paste0("'", colName, "' column not found in '", slot,
            "' slot for ", method_name))
    } else info_col <- info[, colName]
    names(info_col) <- rownames(info)
    msg <- paste0("\nMethod: ", method_name)
    # Output vector
    if(type == "pvalue"){
        msg <- paste0(msg, "\n * '", colName, "'")
        out <- info_col#[!is.na(info_col) & info_col < 1]
    } else if(type == "logfc"){
        msg <- paste0(msg, "\n * -|", colName, "|")
        out <- -abs(info_col)#[!is.na(info_col)])
    } else stop("Please choose between type: pvalue or logfc.")
        msg <- paste0(msg, " column, as ", type, " type, of ", slot,
            " matrix.")
    # Check if direction is not null
    if(!is.null(direction)){
        msg <- paste0(msg, "\n * '", direction,
            "' column, as direction, of statInfo matrix.")
        # Extract statInfo
        statInfo <- method[["statInfo"]]
        if(!is.element(direction, colnames(statInfo))){
            stop(paste0(direction, " column not found for ", method_name))
        } else {
            out <- data.frame(out, statInfo[, direction])
            colnames(out) <- c(colName, direction)
        }
    }
    if(verbose){
        message(msg, appendLF = FALSE)
    }
    return(out)
}

#' @title extractStatistics
#'
#' @export
#' @description
#' Extract the list of p-values or/and log fold changes from the outputs of the
#' differential abundance detection methods.
#'
#' @param object Output of differential abundance detection methods.
#' \code{pValMat}, \code{statInfo} matrices, and method's \code{name} must be
#' present (See vignette for detailed information).
#' @param slot A character vector with 1 or number-of-methods-times repeats of
#' the slot names where to extract values for each method
#' (default \code{slot = "pValMat"}).
#' @param colName A character vector with 1 or number-of-methods-times repeats
#' of the column name of the slot where to extract values for each method
#' (default \code{colName = "rawP"}).
#' @param type A character vector with 1 or number-of-methods-times repeats
#' of the value type of the column selected where to extract values for each
#' method. Two values are possible: \code{"pvalue"} or \code{"logfc"}
#' (default \code{type = "pvalue"}).
#' @param direction A character vector with 1 or number-of-methods-times repeats
#' of the \code{statInfo}'s column name containing information about the signs
#' of differential abundance (usually log fold changes) for each method
#' (default \code{direction = NULL}).
#' @param verbose Boolean to display the kind of extracted values
#' (default \code{verbose = FALSE}).
#'
#' @return A vector or a \code{data.frame} for each method. If
#' \code{direction = NULL}, the \code{colname} column values, transformed
#' according to \code{type} (not tranformed if \code{type = "pvalue"},
#' \code{-abs(value)} if \code{type = "logfc"}), of the \code{slot} are reported
#' , otherwise the \code{direction} column of the \code{statInfo} matrix is
#' added to the output.
#'
#' @seealso \code{\link{getStatistics}}
#'
#' @examples
#' data("ps_plaque_16S")
#'
#' # Add scaling factors
#' ps_plaque_16S <- norm_edgeR(object = ps_plaque_16S, method = "TMM")
#' ps_plaque_16S <- norm_CSS(object = ps_plaque_16S, method = "median")
#'
#' # Perform DA analysis
#' Plaque_16S_DA <- list()
#' Plaque_16S_DA <- within(Plaque_16S_DA, {
#'     # DA analysis
#'     da.limma <- DA_limma(
#'         object = ps_plaque_16S,
#'         design = ~ 1 + HMP_BODY_SUBSITE,
#'         coef = 2,
#'         norm = "TMM"
#'     )
#'     da.limma.css <- DA_limma(
#'         object = ps_plaque_16S,
#'         design = ~ 1 + HMP_BODY_SUBSITE,
#'         coef = 2,
#'         norm = "CSSmedian"
#'     )
#' })
#'
#' # Extract statistics for concordance analysis:
#' # Only p-values
#' extracted_pvalues = extractStatistics(object = Plaque_16S_DA, slot =
#'     "pValMat", colName = "rawP", type = "pvalue")
#' # Only transformed log fold changes -abs(logFC)
#' extracted_abslfc = extractStatistics(object = Plaque_16S_DA, slot =
#'     "statInfo", colName = "logFC", type = "logfc")
#' # Only transformed log fold changes for a method and p-values for the other
#' extracted_abslfc_pvalues = extractStatistics(object = Plaque_16S_DA,
#'     slot = c("statInfo", "pValMat"), colName = c("logFC", "rawP"), type =
#'     c("logfc","pvalue"))
#'
#' # Extract statistics for enrichment analysis:
#' # p-values and log fold changes
#' extracted_pvalues_and_lfc = extractStatistics(object = Plaque_16S_DA,
#'     slot = "pValMat", colName = "rawP", type = "pvalue", direction = "logFC")
#' # transformed log fold changes and untouched log fold changes
#' extracted_abslfc_and_lfc = extractStatistics(object = Plaque_16S_DA,
#'     slot = "statInfo", colName = "logFC", type = "logfc", direction =
#'     "logFC")
#' # Only transformed log fold changes for a method and p-values for the other
#' extracted_mix = extractStatistics(object = Plaque_16S_DA,
#'     slot = c("statInfo", "pValMat"), colName = c("logFC", "rawP"), type =
#'     c("logfc","pvalue"), direction = "logFC")

extractStatistics <- function(object, slot = "pValMat", colName = "rawP",
    type = "pvalue", direction = NULL, verbose = FALSE){
    n_methods = length(object)
    # Check the dimension of slot, colName, and type.
    if(length(slot) == 1)
        slot <- rep(slot, n_methods)
    if(length(colName) == 1)
        colName <- rep(colName, n_methods)
    if(length(type) == 1)
        type <- rep(type, n_methods)
    # Error if they have unequal lengths
    if(length(slot) != n_methods | length(colName) != n_methods |
        length(type) != n_methods)
        stop("Unequal lengths for slot, colName, or type arguments.")
    # Check if direction is defined
    if(!is.null(direction)){
        if(length(direction) == 1)
            direction <- rep(direction, n_methods)
        if(length(direction) != n_methods)
            stop("Wrong length for direction argument.")
    }
    # Rename method names
    names(object) <- unlist(lapply(object, function(method) method[["name"]]))
    if(is.null(direction)){
        out <- mapply(getStatistics, method = object, slot = slot,
            colName = colName, type = type, MoreArgs = list(direction = NULL,
            verbose = verbose), SIMPLIFY = FALSE)
    } else {
        out <- mapply(getStatistics, method = object, slot = slot,
            colName = colName, type = type, direction = direction, MoreArgs =
            list(verbose = verbose), SIMPLIFY = FALSE)
    }
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
#' data("ps_plaque_16S")
#'
#' set.seed(123)
#' phyloseq::sample_data(ps_plaque_16S)$HMP_BODY_SUBSITE <- factor(
#'     phyloseq::sample_data(ps_plaque_16S)$HMP_BODY_SUBSITE)
#' # At least N = 100 is suggested
#' split_list <- createSplits(object = ps_plaque_16S, varName =
#'     "HMP_BODY_SUBSITE", paired = "RSID", balanced = TRUE, N = 10)
#' # Run some methods
#' # lapply -> Subset1 and Subset2
#' Plaque_16S_splitsDA <- lapply(X = split_list, FUN = function(subset){
#'     # apply -> Comparison1, Comparison2, ..., Comparison10
#'     apply(X = subset, MARGIN = 1, FUN = function(splits) {
#'         # Splitting
#'         ps <- phyloseq::prune_samples(phyloseq::sample_names(ps_plaque_16S)
#'             %in% splits, ps_plaque_16S)
#'         # Keep only present taxa
#'         ps <- phyloseq::filter_taxa(ps, function(x) sum(x>0)>0, 1)
#'         # Adding scaling and normalization factors
#'         ps <- norm_edgeR(object = ps, method = "TMM")
#'         ps <- norm_CSS(object = ps, method = "median")
#'         # DA analysis
#'         returnList = list()
#'         returnList = within(returnList, {
#'             da.limma <- DA_limma(object = ps, design = ~ 1 +
#'                 HMP_BODY_SUBSITE, coef = 2, norm = "TMM")
#'             da.limma.css <- DA_limma(object = ps, design = ~ 1 +
#'                 HMP_BODY_SUBSITE, coef = 2, norm = "CSSmedian")
#'         })
#'         return(returnList)
#'     })
#' })
#'
#' # Concordance for p-values
#' concordance_pvalues = createConcordance(object = Plaque_16S_splitsDA, slot =
#'     "pValMat", colName = "rawP", type = "pvalue")
#' # Concordance for log fold changes
#' concordance_logfc = createConcordance(object = Plaque_16S_splitsDA, slot =
#'     "statInfo", colName = "logFC", type = "logfc")
#' # Concordance for log fold changes in the first method and p-values in the
#' # other
#' concordance_logfc_pvalues = createConcordance(object = Plaque_16S_splitsDA,
#'     slot = c("statInfo", "pValMat"), colName = c("logFC", "rawP"), type =
#'     c("logfc","pvalue"))

createConcordance <- function(object, slot = "pValMat", colName = "rawP",
    type = "pvalue"){
    data <- lapply(X = object, FUN = function(subset){
        lapply(X = subset, FUN = function(comparison){
            c(extractStatistics(object = comparison, slot = slot,
                colName = colName, type = type), n_features =
                nrow(comparison[[1]][["pValMat"]]))
        })
    })
    subset1 <- data[["Subset1"]]
    subset2 <- data[["Subset2"]]
    conc_df <- NULL
    for(i in seq_len(length(subset1))){
        comparison1 = subset1[[i]]
        comparison2 = subset2[[i]]
        n_features1 <- comparison1[["n_features"]]
        n_features2 <- comparison2[["n_features"]]
        n_features <- min(n_features1, n_features2)
        n_methods = length(comparison1) - 1
        method_names <- names(comparison1)[seq_len(n_methods)]
        for(j in seq_len(n_methods)){
            for(k in seq_len(n_methods)){
                if(j != k){ # Compute the Between Methods concordance
                    # For subset1
                    vec1 <- comparison1[[method_names[j]]]
                    vec2 <- comparison1[[method_names[k]]]
                    conc1 <- ffpe::CATplot(vec1 = vec1, vec2 = vec2, make.plot =
                        FALSE)
                    conc1 <- data.frame(conc1, "n_features" = n_features1)
                    # For subset2
                    vec1 <- comparison2[[method_names[j]]]
                    vec2 <- comparison2[[method_names[k]]]
                    conc2 <- ffpe::CATplot(vec1 = vec1, vec2 = vec2, make.plot =
                        FALSE)
                    conc2 <- data.frame(conc2, "n_features" = n_features2)
                    # Together
                    conc <- rbind(conc1, conc2)
                } else { # Compute the Within Method concordance
                    vec1 <- comparison1[[method_names[j]]]
                    vec2 <- comparison2[[method_names[k]]]
                    conc <- ffpe::CATplot(vec1 = vec1, vec2 = vec2, make.plot =
                        FALSE)
                    conc <- data.frame(conc, "n_features" = n_features)
                }
                conc[, "method1"] <- method_names[j]
                conc[, "method2"] <- method_names[k]
                conc[, "comparison"] <- paste0("Comparison",i)
                conc_df <- rbind(conc_df, conc)
            }
        }
    }
    concordance <- plyr::ddply(.data = conc_df, .variables = ~ comparison +
        n_features + method1 + method2 + rank, .fun = function(x) data.frame(
        "concordance" = mean(x[, "concordance"])))
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
#' data("ps_plaque_16S")
#'
#' set.seed(123)
#' phyloseq::sample_data(ps_plaque_16S)$HMP_BODY_SUBSITE <- factor(
#'     phyloseq::sample_data(ps_plaque_16S)$HMP_BODY_SUBSITE)
#' split_list <- createSplits(object = ps_plaque_16S,
#'                            varName = "HMP_BODY_SUBSITE",
#'                            paired = "RSID",
#'                            balanced = TRUE,
#'                            N = 10) # At least 100 is suggested
#'
#' # Run some methods
#' # lapply -> Subset1 and Subset2
#' Plaque_16S_splitsDA <- lapply(X = split_list, FUN = function(subset){
#'     # apply -> Comparison1, Comparison2, ..., Comparison10
#'     apply(X = subset, MARGIN = 1, FUN = function(splits) {
#'         # Splitting
#'         ps <- phyloseq::prune_samples(phyloseq::sample_names(ps_plaque_16S)
#'             %in% splits, ps_plaque_16S)
#'         # Keep only present taxa
#'         ps <- phyloseq::filter_taxa(ps, function(x) sum(x>0)>0, 1)
#'         # Adding scaling and normalization factors
#'         ps <- norm_edgeR(object = ps, method = "TMM")
#'         ps <- norm_CSS(object = ps, method = "median")
#'         # DA analysis
#'         returnList = list()
#'         returnList = within(returnList, {
#'             da.limma <- DA_limma(object = ps, design = ~ 1 +
#'                 HMP_BODY_SUBSITE, coef = 2, norm = "TMM")
#'             da.limma.css <- DA_limma(object = ps, design = ~ 1 +
#'                 HMP_BODY_SUBSITE, coef = 2, norm = "CSSmedian")
#'         })
#'         return(returnList)
#'     })
#' })
#'
#' # Concordance for p-values
#' concordance_pvalues = createConcordance(object = Plaque_16S_splitsDA, slot =
#'     "pValMat", colName = "rawP", type = "pvalue")
#'
#' # Add area over the concordance curve
#' concordance_area <- areaCAT(concordance = concordance_pvalues)

areaCAT <- function(concordance, plotIt = FALSE) {
    plyr::ddply(.data = concordance, .variables = ~ comparison + n_features +
        method1 + method2, .fun = function(conc){
        MaxArea <- unique(conc[, "n_features"])
        estimated <- conc[, "concordance"]
        # y = x values -> bisector
        theoretical <- seq_along(estimated)/unique(conc[, "n_features"])
        if (plotIt) {
            graphics::plot(x = theoretical, y = estimated, type = "l", xlim = c(
                0,1), ylim = c(0,1), main = paste0("CAT for ", unique(conc[,
                "method1"]), " & ", unique(conc[, "method2"])))
            graphics::abline(0, 1)
        }
        # Consider the y = x line
        difference = estimated - theoretical
        HeightOver <- (estimated - theoretical)/MaxArea # rescaling
        HeightOver[difference <= 0] <- 0 # removing negative values
        return(data.frame("rank" = conc[, "rank"],
            "concordance" = conc[, "concordance"], "heightOver" = HeightOver,
            "areaOver" = cumsum(HeightOver)))
    })
}

#' @title plotConcordanceHeatmap
#'
#' @keywords internal
#' @importFrom ggplot2 ggplot aes facet_grid geom_ribbon geom_line geom_segment
#' @importFrom ggplot2 theme element_blank unit coord_cartesian geom_rect
#' @importFrom ggplot2 unit scale_fill_distiller scale_color_manual
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous xlab ylab
#' @importFrom ggplot2 guides guide_colorbar guide_legend
#' @description
#' Plots the heatmap of concordances.
#'
#' @param c_df A simplified concordance \code{data.frame} produced in
#' \link{plotConcordance} function.
#' @param threshold The threshold for rank (x-axis upper limit if all methods
#' have a higher number of computed statistics).
#' @param cols A named vector containing the color hex codes.
#'
#' @return a \code{ggplot2} object
#'
#' @seealso \code{\link{createConcordance}} and \code{\link{plotConcordance}}

plotConcordanceHeatmap <- function(c_df, threshold, cols){
    concordance <- max_areaOver <- method1 <- method2 <- NULL
    # If we consider the trapezoid with bases B = 1, b = 1-threshold/n_features
    # and height h = threshold, its area will be (B+b)*h/2.
    # Then rescale the area by n_features in order to make it from 0 to 0.5
    n_features <- max(c_df[, "n_features"])
    rescaled_area <- (2 - threshold/n_features) * threshold/2 / n_features
    ggplot2::ggplot(data = c_df, mapping = ggplot2::aes(x = rank, y =
        concordance)) +
        ggplot2::facet_grid(method1 ~ method2, scales = "free_x", switch = "y"
            ) +
        ggplot2::geom_ribbon(aes(ymin = rank/n_features, ymax = concordance,
            fill = max_areaOver)) +
        ggplot2::geom_line() +
        ggplot2::geom_line(mapping = ggplot2::aes(x = rank, y = rank/n_features)
            , lty = 2) +
        ggplot2::theme(legend.position = "bottom",
            axis.line.x.bottom = ggplot2::element_blank(),
            axis.line.y.right = ggplot2::element_blank(),
            strip.text = ggplot2::element_blank(),
            strip.background = ggplot2::element_blank(),
            panel.spacing = ggplot2::unit(0.1,"cm"),
            plot.margin = ggplot2::unit(c(0.1,0.1,0.1,0.1), "cm")) +
        ggplot2::coord_cartesian(xlim = c(1, threshold), ylim = c(0,1),
            expand = FALSE) +
        geom_segment(data = c_df[c_df[, "method1"] == levels(c_df[, "method1"])[
            1],], mapping = ggplot2::aes(x = 1, xend = threshold, y = 1,
            yend = 1, size = 1, color = method2), show.legend = FALSE) +
        geom_segment(data = c_df[c_df[, "method2"] == levels(c_df[, "method1"])[
            1],], mapping = ggplot2::aes(x = 1, xend = 1, y = 0,
            yend = 1, size = 1, color = method1), show.legend = FALSE) +
        ggplot2::geom_rect(data = c_df[c_df[, "method1"] == c_df[, "method2"],],
            mapping = ggplot2::aes(xmin = 1, xmax = threshold, ymin = 0, ymax =
            1), color = "red", fill = NA, show.legend = FALSE) +
        ggplot2::scale_fill_distiller(palette = "RdYlBu", limits = c(0,
            rescaled_area), direction = 1) +
        ggplot2::scale_color_manual(values = cols) +
        ggplot2::scale_y_continuous(breaks = c(0.5, 1), position = "right") +
        ggplot2::scale_x_continuous(limits = c(1, threshold)) +
        ggplot2::xlab("Rank") + ggplot2::ylab("Average concordance") +
        ggplot2::guides(fill = ggplot2::guide_colorbar(title =
            "Average area over the bisector:", direction = "horizontal",
            ticks.colour = "black", frame.colour = "black"),
            colour = ggplot2::guide_legend(NULL))
}

#' @title plotConcordanceDendrogram
#'
#' @keywords internal
#' @importFrom ggdendro dendro_data label
#' @importFrom ggplot2 ggplot aes geom_segment theme element_blank element_rect
#' @importFrom ggplot2 unit scale_fill_manual scale_color_manual geom_point
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous xlab ylab
#' @importFrom ggplot2 geom_label coord_flip scale_y_reverse scale_x_reverse
#' @description
#' Plots the method's dendrogram of concordances.
#'
#' @param hc Hierarchical clustering results produced in
#' \link{plotConcordance} function.
#' @param direction vertical (default \code{direction = "v"}) or horizontal
#' (\code{direction = "h"}).
#' @param cols A named vector containing the color hex codes.
#'
#' @return a \code{ggplot2} object
#'
#' @seealso \code{\link{createConcordance}} and \code{\link{plotConcordance}}

plotConcordanceDendrogram <- function(hc, direction = "v", cols){
    x <- y <- xend <- yend <- NULL
    g_dendro <- ggplot2::ggplot() +
        ggplot2::geom_segment(data = ggdendro::dendro_data(hc)[["segments"]],
            mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
        ggplot2::theme(axis.line = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(fill = "white"),
            panel.grid = ggplot2::element_blank(),
            legend.position = "none",
            panel.spacing = ggplot2::unit(0, "lines"),
            plot.margin = ggplot2::unit(c(0.1,0.1,0.1,0.1), "cm"))
    if(direction == "v"){
        g_dendro <- g_dendro +
            ggplot2::geom_label(data = ggdendro::dendro_data(hc)[["labels"]],
                mapping = ggplot2::aes(x = x, y = y, label = label, hjust = 1,
                fill = label), nudge_y = 0) +
            ggplot2::coord_flip() +
            ggplot2::scale_y_reverse(expand = c(0,0,0,0)) +
            ggplot2::scale_x_reverse() +
            ggplot2::scale_fill_manual(values = cols)
    } else if(direction == "h"){
        g_dendro <- g_dendro +
            ggplot2::geom_point(data = ggdendro::dendro_data(hc)[["labels"]],
                mapping = ggplot2::aes(x = x, y = y, color = label), nudge_y =
                0.1, size = 5) +
            ggplot2::scale_y_continuous(expand = c(0,0,0,0)) +
            ggplot2::scale_x_continuous(expand = c(0,0,0,0)) +
            ggplot2::scale_color_manual(values = cols)
    } else stop("Please supply a direction value between 'h' for horizontal or
                'v' for vertical.")
    return(g_dendro)
}

#' @title plotConcordance
#'
#' @export
#' @importFrom plyr ddply
#' @importFrom reshape2 dcast
#' @importFrom stats hclust as.dist
#' @description
#' Produce a list of graphical outputs summarizing the between and
#' within method concordance.
#'
#' @param concordance A long format \code{data.frame} produced by
#' \link{createConcordance} function.
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
#' data("ps_plaque_16S")
#'
#' set.seed(123)
#' phyloseq::sample_data(ps_plaque_16S)$HMP_BODY_SUBSITE <- factor(
#'     phyloseq::sample_data(ps_plaque_16S)$HMP_BODY_SUBSITE)
#' split_list <- createSplits(object = ps_plaque_16S,
#'                            varName = "HMP_BODY_SUBSITE",
#'                            paired = "RSID",
#'                            balanced = TRUE,
#'                            N = 10) # At least 100 is suggested
#'
#' # Run some methods
#' # lapply -> Subset1 and Subset2
#' Plaque_16S_splitsDA <- lapply(X = split_list, FUN = function(subset){
#'     # apply -> Comparison1, Comparison2, ..., Comparison10
#'     apply(X = subset, MARGIN = 1, FUN = function(splits) {
#'         # Splitting
#'         ps <- phyloseq::prune_samples(phyloseq::sample_names(ps_plaque_16S)
#'             %in% splits, ps_plaque_16S)
#'         # Keep only present taxa
#'         ps <- phyloseq::filter_taxa(ps, function(x) sum(x>0)>0, 1)
#'         # Adding scaling and normalization factors
#'         ps <- norm_edgeR(object = ps, method = "TMM")
#'         ps <- norm_CSS(object = ps, method = "median")
#'         # DA analysis
#'         returnList = list()
#'         returnList = within(returnList, {
#'             da.limma <- DA_limma(object = ps, design = ~ 1 +
#'                 HMP_BODY_SUBSITE, coef = 2, norm = "TMM")
#'             da.limma.css <- DA_limma(object = ps, design = ~ 1 +
#'                 HMP_BODY_SUBSITE, coef = 2, norm = "CSSmedian")
#'         })
#'         return(returnList)
#'     })
#' })
#'
#' # Concordance for p-values
#' concordance_pvalues = createConcordance(object = Plaque_16S_splitsDA, slot =
#'     "pValMat", colName = "rawP", type = "pvalue")
#'
#' # plot concordances from rank 1 to 50.
#' concordance_plots <- plotConcordance(concordance = concordance_pvalues,
#'     threshold = 50)

plotConcordance <- function(concordance, threshold = NULL, cols = NULL){
    area_df <- areaCAT(concordance = concordance, plotIt = FALSE)
    # Maximum rank for each method and comparison
    maxRank <- plyr::ddply(.data = area_df, .variables = ~ comparison +
        method1 + method2, .fun = function(x){
            c(maxRank = max(x[, "rank"]))
        })
    # The maximum threshold allowed for each pair of methods is the minimum
    # maxRank
    threshold_df <- plyr::ddply(.data = maxRank, .variables = ~ method1 +
        method2, .fun = function(x){
            c(threshold = min(x[, "maxRank"]))
        })
    # Common threshold
    common_threshold <- min(threshold_df[, "threshold"])
    # cat(common_threshold,"\n")
    if(is.null(threshold)){
        threshold <- common_threshold
    } else {
        if(threshold > common_threshold){
            cat(threshold,"\n")
            threshold <- common_threshold
            cat(threshold,"\n")
        }
    }
    # Keep only ranks that are lower than threshold
    concordance_df <- area_df[area_df[, "rank"] <= threshold, ]
    # Average values for plotting
    c_df <- plyr::ddply(.data = concordance_df, .variables = ~ method1 +
        method2 + rank, .fun = function(x){
            colMeans(x[, c("n_features", "concordance", "areaOver")])
        })
    # Extract the area at the maximum rank (= threshold)
    max_areaOver <- plyr::ddply(c_df, .variables = ~ method1 + method2, .fun =
        function(x) c(areaOver = max(x[, "areaOver"])))
    # Repeat the max_areaOver for all ranks in pair of method1 and method2
    c_df[, "max_areaOver"] <- rep(max_areaOver[, "areaOver"], each = threshold)
    # Clustering of methods
    # Organize area values in a n_methods x n_methods matrix
    dist_df <- reshape2::dcast(data = max_areaOver, formula = method1 ~ method2,
        value.var = "areaOver")
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
    if(is.null(cols))
        cols <- createColors(hc[["labels"]])
    g_list <- list(
        "concordanceDendrogram" = plotConcordanceDendrogram(hc = hc,
            direction = "v", cols = cols),
        "concordanceHeatmap" = plotConcordanceHeatmap(c_df = c_df, threshold =
            threshold, cols = cols))
    return(g_list)
}


