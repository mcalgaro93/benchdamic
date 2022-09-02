#' @title plotEnrichment
#'
#' @export
#' @importFrom plyr ldply ddply
#' @import ggplot2
#' @description Summary plot for the number of differentially abundant (DA)
#' features and their association with enrichment variable. If some features are
#' UP-Abundant or DOWN-Abundant (or just DA), several bars will be represented
#' in the corresponding direction. Significance thresholds are reported
#' over/under each bar, according to the Fisher exact tests.
#'
#' @inheritParams plotContingency
#' @inheritParams createEnrichment
#'
#' @return a \code{ggplot2} object.
#'
#' @seealso \code{\link{createEnrichment}}, \code{\link{plotContingency}}, and
#' \code{\link{plotMutualFindings}}.
#'
#' @examples 
#' data("ps_plaque_16S")
#' data("microbial_metabolism")
#'
#' # Extract genera from the phyloseq tax_table slot
#' genera <- phyloseq::tax_table(ps_plaque_16S)[, "GENUS"]
#' # Genera as rownames of microbial_metabolism data.frame
#' rownames(microbial_metabolism) <- microbial_metabolism$Genus
#' # Match OTUs to their metabolism
#' priorInfo <- data.frame(genera,
#'     "Type" =  microbial_metabolism[genera, "Type"])
#' # Unmatched genera becomes "Unknown"
#' unknown_metabolism <- is.na(priorInfo$Type)
#' priorInfo[unknown_metabolism, "Type"] <- "Unknown"
#' priorInfo$Type <- factor(priorInfo$Type)
#' # Add a more informative names column
#' priorInfo[, "newNames"] <- paste0(rownames(priorInfo), priorInfo[, "GENUS"])
#'
#' # Add some normalization/scaling factors to the phyloseq object
#' my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
#'     method = c("TMM", "CSS"))
#' ps_plaque_16S <- runNormalizations(normalization_list = my_norm,
#'     object = ps_plaque_16S)
#'
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ 1 + RSID + HMP_BODY_SUBSITE, 
#'     coef = "HMP_BODY_SUBSITESupragingival Plaque",
#'     norm = c("TMM", "CSS"))
#'
#' # Make sure the subject ID variable is a factor
#' phyloseq::sample_data(ps_plaque_16S)[, "RSID"] <- as.factor(
#'     phyloseq::sample_data(ps_plaque_16S)[["RSID"]])
#'     
#' # Perform DA analysis
#' Plaque_16S_DA <- runDA(method_list = my_limma, object = ps_plaque_16S)
#'
#' # Enrichment analysis
#' enrichment <- createEnrichment(object = Plaque_16S_DA,
#'     priorKnowledge = priorInfo, enrichmentCol = "Type", namesCol = "GENUS",
#'     slot = "pValMat", colName = "adjP", type = "pvalue", direction = "logFC",
#'     threshold_pvalue = 0.1, threshold_logfc = 1, top = 10, verbose = TRUE)
#'     
#' # Contingency tables
#' plotContingency(enrichment = enrichment, method = "limma.TMM")
#' # Barplots
#' plotEnrichment(enrichment, enrichmentCol = "Type")
#' # Mutual findings
#' plotMutualFindings(
#'     enrichment = enrichment, enrichmentCol = "Type",
#'     n_methods = 1
#' )
plotEnrichment <- function(enrichment, enrichmentCol, levels_to_plot) {
    # Extract freq and p-values
    df_to_plot <- plyr::ldply(.data = enrichment, .fun = function(x) {
        plyr::ldply(x[["summaries"]], .fun = function(tab) {
            out <- cbind(tab,
                "direction" = rownames(tab),
                "variable" = colnames(tab)[1]
            )
            colnames(out) <- c("value", "pvalue", "direction", "variable")
            return(out)
        }, .id = NULL)
    }, .id = "method")
    if (!missing(levels_to_plot)) { # Check which levels to plot
        keep <- is.element(
            el = df_to_plot[, "variable"],
            set = levels_to_plot
        )
        df_to_plot <- df_to_plot[keep, ]
    }
    # Order methods
    method <- direction <- value <- variable <- pvalue <- NULL
    method_DA <- plyr::ddply(
        .data = df_to_plot, .variables = c("method"),
        function(x) c("freq" = sum(x[, "value"]))
    )
    ord <- method_DA[order(method_DA[, "freq"], decreasing = TRUE), "method"]
    max_DA <- max(abs(df_to_plot[, "value"])) # Max freq available
    g_enrichment <- ggplot(
        data = df_to_plot,
        mapping = aes(
            x = method,
            y = ifelse(test = (direction == "UP Abundant" | direction == "DA"),
                yes = value, no = -value
            ), fill = variable
        )
    ) +
        geom_col(position = "dodge") +
        geom_text(mapping = aes(
            x = method,
            y = ifelse(test = (direction == "UP Abundant" | direction == "DA"),
                yes = value + 5, no = -value - 5
            ),
            label = ifelse(test = (pvalue > 0.1), yes = "", no = ifelse(
                test = (pvalue > 0.05), yes = ".", no = ifelse(
                    test = (pvalue > 0.01), yes = "*", no = ifelse(
                        test = (pvalue > 0.001), yes = "**", no = "***"
                    )
                )
            )),
            color = variable
        ), position = position_dodge(width = 1)) +
        geom_hline(
            mapping = aes(yintercept = 0),
            color = "black"
        ) +
        scale_x_discrete(limits = ord) +
        ylab("Number of DA taxa") +
        theme(
            axis.text.x = element_text(
                angle = 90,
                hjust = 1, vjust = 0.5
            ),
            axis.text.y.right = element_text(angle = -90, hjust = 0.5),
            axis.ticks.y.right = element_blank(),
            legend.position = "bottom"
        ) +
        guides(fill = guide_legend(
            title = enrichmentCol,
            override.aes = aes(label = "")
        ), color = "none") +
        ggtitle("Enrichment analysis",
            subtitle =
                "Fisher exact test pvalue: 0 *** 0.001 ** 0.01 * 0.05 . 0.1"
        )
    if (!is.element(el = "DA", set = df_to_plot[, "direction"])) {
        g_enrichment <- g_enrichment +
            scale_y_continuous(sec.axis = sec_axis(
                trans = ~ . / max_DA, name = "DA direction",
                    breaks = c(-0.5, 0.5),
                labels = c("DOWN Abundant", "UP Abundant")
            ))
    }
    return(g_enrichment)
}

#' @title plotContingency
#'
#' @export
#' @importFrom plyr ldply
#' @importFrom reshape2 melt
#' @import ggplot2
#' @description Plot of the contingency tables for the chosen method. The
#' top-left cells are colored, according to Fisher exact tests' p-values, if the
#' number of features in those cells are enriched.
#'
#' @param enrichment enrichment object produced by \link{createEnrichment}
#' function.
#' @param method name of the method to plot.
#' @param levels_to_plot A character vector containing the levels of the
#' enrichment variable to plot.
#'
#' @return a \code{ggplot2} object.
#'
#' @seealso \code{\link{createEnrichment}}, \code{\link{plotEnrichment}}, and
#' \code{\link{plotMutualFindings}}.
#'
#' @examples
#' data("ps_plaque_16S")
#' data("microbial_metabolism")
#'
#' # Extract genera from the phyloseq tax_table slot
#' genera <- phyloseq::tax_table(ps_plaque_16S)[, "GENUS"]
#' # Genera as rownames of microbial_metabolism data.frame
#' rownames(microbial_metabolism) <- microbial_metabolism$Genus
#' # Match OTUs to their metabolism
#' priorInfo <- data.frame(genera,
#'     "Type" =  microbial_metabolism[genera, "Type"])
#' # Unmatched genera becomes "Unknown"
#' unknown_metabolism <- is.na(priorInfo$Type)
#' priorInfo[unknown_metabolism, "Type"] <- "Unknown"
#' priorInfo$Type <- factor(priorInfo$Type)
#' # Add a more informative names column
#' priorInfo[, "newNames"] <- paste0(rownames(priorInfo), priorInfo[, "GENUS"])
#'
#' # Add some normalization/scaling factors to the phyloseq object
#' my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
#'     method = c("TMM", "CSS"))
#' ps_plaque_16S <- runNormalizations(normalization_list = my_norm,
#'     object = ps_plaque_16S)
#'
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ 1 + RSID + HMP_BODY_SUBSITE, 
#'     coef = "HMP_BODY_SUBSITESupragingival Plaque",
#'     norm = c("TMM", "CSS"))
#'
#' # Make sure the subject ID variable is a factor
#' phyloseq::sample_data(ps_plaque_16S)[, "RSID"] <- as.factor(
#'     phyloseq::sample_data(ps_plaque_16S)[["RSID"]])
#'     
#' # Perform DA analysis
#' Plaque_16S_DA <- runDA(method_list = my_limma, object = ps_plaque_16S)
#'
#' # Enrichment analysis
#' enrichment <- createEnrichment(object = Plaque_16S_DA,
#'     priorKnowledge = priorInfo, enrichmentCol = "Type", namesCol = "GENUS",
#'     slot = "pValMat", colName = "adjP", type = "pvalue", direction = "logFC",
#'     threshold_pvalue = 0.1, threshold_logfc = 1, top = 10, verbose = TRUE)
#'     
#' # Contingency tables
#' plotContingency(enrichment = enrichment, method = "limma.TMM")
#' # Barplots
#' plotEnrichment(enrichment, enrichmentCol = "Type")
#' # Mutual findings
#' plotMutualFindings(
#'     enrichment = enrichment, enrichmentCol = "Type",
#'     n_methods = 1
#' )
plotContingency <- function(enrichment, method, levels_to_plot) {
    # Extract enrichment tables
    enrichment_table <- enrichment[[method]][["tables"]]
    if (length(enrichment_table) == 0) {
        stop("No DA features for ", method)
    }
    df_to_plot <- reshape2::melt(plyr::ldply(enrichment_table,
        .fun = function(tab) {
            out <- cbind(
                "direction" = rownames(tab), tab,
                "level" = colnames(tab)[1]
            )
            return(out)
        }
    ), na.rm = TRUE)
    df_to_plot$dir <- gsub(
        pattern = "non-", replacement = "",
        x = df_to_plot[, "direction"]
    )
    if (!missing(levels_to_plot)) { # Check which levels to plot
        keep <- is.element(el = df_to_plot[, "level"], set = levels_to_plot)
        df_to_plot <- df_to_plot[keep, ]
    }
    df_to_plot[, "direction"] <- factor(df_to_plot[, "direction"],
        levels =
            c(
                "non-DOWN Abundant", "DOWN Abundant", "non-UP Abundant",
                "UP Abundant", "non-DA", "DA"
            ), ordered = TRUE
    )
    # Extract p-values
    pvalue_table <- plyr::ldply(enrichment[[method]][["tests"]])
    colnames(pvalue_table) <- c("direction_variable", "pvalue")
    pvalue_table[, c("direction", "variable")] <- plyr::ldply(strsplit(
        pvalue_table[, "direction_variable"],
        split = "-"
    ))
    df_to_plot <- merge(
        x = df_to_plot, y = pvalue_table,
        by.x = c("variable", "direction"), by.y = c("variable", "direction"),
        all.x = TRUE
    )
    # Plot
    variable <- direction <- pvalue <- value <- NULL
    ggplot(data = df_to_plot, mapping = aes(
        x = variable,
        y = direction, fill = pvalue
    )) +
        facet_grid(dir ~ level, scales = "free") +
        geom_tile(width = 0.8, height = 0.8) +
        geom_text(mapping = aes(label = round(
            x = value,
            digits = 2
        )), color = "black") +
        theme(axis.text.x = element_text(
            angle = 90,
            hjust = 1, vjust = 0.5
        )) +
        ggtitle(
            label = paste0("Enrichment tables - ", method),
            subtitle = "Colored by Fisher p-value"
        ) +
        scale_fill_distiller(
            type = "seq", palette = "RdYlBu",
            breaks = seq.int(from = 0, to = 0.2, by = 0.02), limits = c(0, 0.1),
            direction = 1
        ) +
        guides(fill = guide_colorbar(
            ticks.colour = "black",
            frame.colour = "black"
        ))
}

#' @title plotMutualFindings
#'
#' @export
#' @importFrom plyr ldply ddply
#' @import ggplot2
#' @description Plot and filter the features which are considered differentially
#' abundant, simultaneously, by a specified number of methods.
#'
#' @inheritParams plotEnrichment
#' @param n_methods minimum number of method that mutually find the features.
#'
#' @return a \code{ggplot2} object.
#'
#' @seealso \code{\link{createEnrichment}}, \code{\link{plotEnrichment}}, and
#' \code{\link{plotContingency}}.
#'
#' @examples
#' data("ps_plaque_16S")
#' data("microbial_metabolism")
#' 
#' # Extract genera from the phyloseq tax_table slot
#' genera <- phyloseq::tax_table(ps_plaque_16S)[, "GENUS"]
#' # Genera as rownames of microbial_metabolism data.frame
#' rownames(microbial_metabolism) <- microbial_metabolism$Genus
#' # Match OTUs to their metabolism
#' priorInfo <- data.frame(genera,
#'     "Type" =  microbial_metabolism[genera, "Type"])
#' # Unmatched genera becomes "Unknown"
#' unknown_metabolism <- is.na(priorInfo$Type)
#' priorInfo[unknown_metabolism, "Type"] <- "Unknown"
#' priorInfo$Type <- factor(priorInfo$Type)
#' # Add a more informative names column
#' priorInfo[, "newNames"] <- paste0(rownames(priorInfo), priorInfo[, "GENUS"])
#' 
#' # Add some normalization/scaling factors to the phyloseq object
#' my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
#'    method = c("TMM", "CSS"))
#' ps_plaque_16S <- runNormalizations(normalization_list = my_norm,
#'     object = ps_plaque_16S)
#' 
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ 1 + RSID + HMP_BODY_SUBSITE,
#'     coef = "HMP_BODY_SUBSITESupragingival Plaque",
#'     norm = c("TMM", "CSS"))
#' 
#' # Make sure the subject ID variable is a factor
#' phyloseq::sample_data(ps_plaque_16S)[, "RSID"] <- as.factor(
#'     phyloseq::sample_data(ps_plaque_16S)[["RSID"]])
#' 
#' # Perform DA analysis
#' Plaque_16S_DA <- runDA(method_list = my_limma, object = ps_plaque_16S)
#' 
#' # Enrichment analysis
#' enrichment <- createEnrichment(object = Plaque_16S_DA,
#'     priorKnowledge = priorInfo, enrichmentCol = "Type", namesCol = "GENUS",
#'     slot = "pValMat", colName = "adjP", type = "pvalue", direction = "logFC",
#'     threshold_pvalue = 0.1, threshold_logfc = 1, top = 10, verbose = TRUE)
#' 
#' # Contingency tables
#' plotContingency(enrichment = enrichment, method = "limma.TMM")
#' # Barplots
#' plotEnrichment(enrichment, enrichmentCol = "Type")
#' # Mutual findings
#' plotMutualFindings(
#'     enrichment = enrichment, enrichmentCol = "Type",
#'     n_methods = 1
#' )
plotMutualFindings <- function(enrichment, enrichmentCol, levels_to_plot,
                               n_methods = 1) {
    DA_df <- plyr::ldply(enrichment, function(x) x[["data"]], .id = "method")
    df_to_plot <- DA_df[DA_df[, "DA"] != "non-DA", ]
    if (n_methods > 1) { # Filter only DA mutually found by n_methods methods
        feature_df <- plyr::ddply(df_to_plot, .variables = "feature", nrow)
        feature_to_keep <- feature_df[
            which(feature_df[, "V1"] >= n_methods),
            "feature"
        ]
        keep <- is.element(el = df_to_plot[, "feature"], set = feature_to_keep)
        df_to_plot <- df_to_plot[keep, ]
    }
    if (!missing(levels_to_plot)) { # Check which levels to plot
        keep <- is.element(
            el = df_to_plot[, enrichmentCol],
            set = levels_to_plot
        )
        df_to_plot <- df_to_plot[keep, ]
    }
    df_to_plot <- iterative_ordering(
        df = df_to_plot,
        var_names = c("method", "feature", enrichmentCol),
        decreasing = c(TRUE, FALSE, TRUE)
    )
    method <- feature <- DA <- NULL
    ggplot(data = df_to_plot, mapping = aes(
        x = method,
        y = feature, fill = DA
    )) +
        facet_grid(get(enrichmentCol) ~ .,
            scales = "free",
            space = "free", drop = TRUE
        ) +
        geom_tile(width = 0.8, height = 0.8) +
        theme(
            axis.text.x = element_text(
                angle = 90,
                hjust = 1, vjust = 0.5
            ),
            strip.text.y = element_text(angle = 0)
        ) +
        ggtitle(
            label = "Summary of DA features",
            subtitle = paste0("By ", enrichmentCol)
        )
}

#' @title plotPositives
#'
#' @export
#' @importFrom plyr ldply ddply
#' @import ggplot2
#' @description Plot the difference between the number of true positives (TP)
#' and false positives (FP) for each method and for each 'top' threshold
#' provided by the \code{createPositives()} function.
#'
#' @param positives \code{data.frame} object produced by
#' \code{createPositives()} function.
#' @param cols named vector of cols (default \code{cols = NULL}).
#'
#' @return a \code{ggplot2} object.
#'
#' @seealso \code{\link{getPositives}}, \code{\link{createPositives}}.
#'
#' @examples
#' data("ps_plaque_16S")
#' data("microbial_metabolism")
#'
#' # Extract genera from the phyloseq tax_table slot
#' genera <- phyloseq::tax_table(ps_plaque_16S)[, "GENUS"]
#' # Genera as rownames of microbial_metabolism data.frame
#' rownames(microbial_metabolism) <- microbial_metabolism$Genus
#' # Match OTUs to their metabolism
#' priorInfo <- data.frame(genera,
#'     "Type" =  microbial_metabolism[genera, "Type"])
#' # Unmatched genera becomes "Unknown"
#' unknown_metabolism <- is.na(priorInfo$Type)
#' priorInfo[unknown_metabolism, "Type"] <- "Unknown"
#' priorInfo$Type <- factor(priorInfo$Type)
#' # Add a more informative names column
#' priorInfo[, "newNames"] <- paste0(rownames(priorInfo), priorInfo[, "GENUS"])
#'
#' # Add some normalization/scaling factors to the phyloseq object
#' my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
#'     method = c("TMM", "CSS"))
#' ps_plaque_16S <- runNormalizations(normalization_list = my_norm,
#'     object = ps_plaque_16S)
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ 1 + RSID + HMP_BODY_SUBSITE, 
#'     coef = "HMP_BODY_SUBSITESupragingival Plaque",
#'     norm = c("TMM", "CSS"))
#'
#' # Make sure the subject ID variable is a factor
#' phyloseq::sample_data(ps_plaque_16S)[, "RSID"] <- as.factor(
#'     phyloseq::sample_data(ps_plaque_16S)[["RSID"]])
#'     
#' # Perform DA analysis
#' Plaque_16S_DA <- runDA(method_list = my_limma, object = ps_plaque_16S)
#'
#' # Count TPs and FPs, from the top 1 to the top 20 features.
#' # As direction is supplied, features are ordered by "logFC" absolute values.
#' positives <- createPositives(object = Plaque_16S_DA,
#'     priorKnowledge = priorInfo, enrichmentCol = "Type", 
#'     namesCol = "newNames", slot = "pValMat", colName = "rawP", 
#'     type = "pvalue", direction = "logFC", threshold_pvalue = 1, 
#'     threshold_logfc = 0, top = 1:20, alternative = "greater", 
#'     verbose = FALSE,
#'     TP = list(c("DOWN Abundant", "Anaerobic"), c("UP Abundant", "Aerobic")),
#'     FP = list(c("DOWN Abundant", "Aerobic"), c("UP Abundant", "Anaerobic")))
#'
#' # Plot the TP-FP differences for each threshold
#' plotPositives(positives = positives)
plotPositives <- function(positives, cols = NULL) {
    if (is.null(cols)) {
        cols <- createColors(positives[, "method"])
    }
    top <- TP <- FP <- method <- NULL
    ggplot(data = positives, mapping = aes(
        x = top,
        y = TP - FP, color = method
    )) +
        geom_point() +
        geom_path(mapping = aes(group = method)) +
        scale_color_manual(values = cols) +
        ggtitle(
            label = "Putative TP - Putative FP",
            subtitle = paste0(
                "From ", min(positives[, "top"]), " to ",
                max(positives[, "top"]), " top features."
            )
        )
}
