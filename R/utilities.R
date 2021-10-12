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
#' getStatistics(
#'     method = da.limma, slot = "pValMat", colName = "rawP",
#'     type = "pvalue", direction = NULL
#' )
#' # get negative abs(logFC) values
#' getStatistics(
#'     method = da.limma, slot = "statInfo", colName = "logFC",
#'     type = "logfc", direction = NULL
#' )
#' # get p-values and logFC
#' getStatistics(
#'     method = da.limma, slot = "pValMat", colName = "rawP",
#'     type = "pvalue", direction = "logFC"
#' )
getStatistics <- function(method, slot = "pValMat", colName = "rawP",
                          type = "pvalue", direction = NULL, verbose = FALSE) {
    # Info extraction
    method_name <- method[["name"]]
    if (!is.element(slot, names(method))) {
        stop("'", slot, "' slot not found for ", method_name)
    } else {
        info <- method[[slot]]
    }
    if (!is.element(colName, colnames(info))) {
        stop("'", colName, "' column not found in '", slot, "' slot for ",
            method_name)
    } else {
        info_col <- info[, colName]
    }
    names(info_col) <- rownames(info)
    msg <- paste0("\nMethod: ", method_name)
    # Output vector
    if (type == "pvalue") {
        msg <- paste0(msg, "\n * '", colName, "'")
        out <- info_col # [!is.na(info_col) & info_col < 1]
    } else if (type == "logfc") {
        msg <- paste0(msg, "\n * -|", colName, "|")
        out <- -abs(info_col) # [!is.na(info_col)])
    } else {
        stop("Please choose between type: pvalue or logfc.")
    }
    msg <- paste0(
        msg, " column, as ", type, " type, of ", slot,
        " matrix."
    )
    # Check if direction is not null
    if (!is.null(direction)) {
        msg <- paste0(
            msg, "\n * '", direction,
            "' column, as direction, of statInfo matrix."
        )
        # Extract statInfo
        statInfo <- method[["statInfo"]]
        if (!is.element(direction, colnames(statInfo))) {
            stop(direction, " column not found for ", method_name)
        } else {
            out <- data.frame(out, statInfo[, direction])
            colnames(out) <- c(colName, direction)
        }
    }
    if (verbose) {
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
#' # Add scaling factors
#' my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
#'     method = c("TMM", "median"))
#' ps_plaque_16S <- runNormalizations(normalization_list = my_norm,
#'     object = ps_plaque_16S)
#' # Perform DA analysis
#' my_methods <- set_limma(design = ~ 1 + HMP_BODY_SUBSITE, coef = 2,
#'     norm = c("TMM", "CSSmedian"))
#' Plaque_16S_DA <- runDA(method_list = my_methods, object = ps_plaque_16S)
#' ### Extract statistics for concordance analysis:
#' # Only p-values
#' extracted_pvalues <- extractStatistics(
#'     object = Plaque_16S_DA, slot =
#'         "pValMat", colName = "rawP", type = "pvalue"
#' )
#' # Only transformed log fold changes -abs(logFC)
#' extracted_abslfc <- extractStatistics(
#'     object = Plaque_16S_DA, slot =
#'         "statInfo", colName = "logFC", type = "logfc"
#' )
#' # Only transformed log fold changes for a method and p-values for the other
#' extracted_abslfc_pvalues <- extractStatistics(
#'     object = Plaque_16S_DA,
#'     slot = c("statInfo", "pValMat"), colName = c("logFC", "rawP"), type =
#'         c("logfc", "pvalue")
#' )
#' ### Extract statistics for enrichment analysis:
#' # p-values and log fold changes
#' extracted_pvalues_and_lfc <- extractStatistics(
#'     object = Plaque_16S_DA,
#'     slot = "pValMat", colName = "rawP", type = "pvalue", direction = "logFC"
#' )
#' # transformed log fold changes and untouched log fold changes
#' extracted_abslfc_and_lfc <- extractStatistics(
#'     object = Plaque_16S_DA,
#'     slot = "statInfo", colName = "logFC", type = "logfc", direction =
#'         "logFC"
#' )
#' # Only transformed log fold changes for a method and p-values for the other
#' extracted_mix <- extractStatistics(
#'     object = Plaque_16S_DA,
#'     slot = c("statInfo", "pValMat"), colName = c("logFC", "rawP"), type =
#'         c("logfc", "pvalue"), direction = "logFC"
#' )
extractStatistics <- function(object, slot = "pValMat", colName = "rawP",
                              type = "pvalue", direction = NULL, verbose = FALSE) {
    n_methods <- length(object)
    # Check the dimension of slot, colName, and type.
    if (length(slot) == 1) {
        slot <- rep(slot, n_methods)
    }
    if (length(colName) == 1) {
        colName <- rep(colName, n_methods)
    }
    if (length(type) == 1) {
        type <- rep(type, n_methods)
    }
    # Error if they have unequal lengths
    if (length(slot) != n_methods | length(colName) != n_methods |
        length(type) != n_methods) {
        stop("Unequal lengths for slot, colName, or type arguments.")
    }
    # Check if direction is defined
    if (!is.null(direction)) {
        if (length(direction) == 1) {
            direction <- rep(direction, n_methods)
        }
        if (length(direction) != n_methods) {
            stop("Wrong length for direction argument.")
        }
    }
    # Rename method names
    names(object) <- unlist(lapply(object, function(method) method[["name"]]))
    if (is.null(direction)) {
        out <- mapply(getStatistics,
            method = object, slot = slot,
            colName = colName, type = type, MoreArgs = list(
                direction = NULL,
                verbose = verbose
            ), SIMPLIFY = FALSE
        )
    } else {
        out <- mapply(getStatistics,
            method = object, slot = slot,
            colName = colName, type = type, direction = direction, MoreArgs =
                list(verbose = verbose), SIMPLIFY = FALSE
        )
    }
    return(out)
}

#' @title getDA
#'
#' @export
#' @description
#' Inspect the list of p-values or/and log fold changes from the output of a
#' differential abundance detection method.
#'
#' @inheritParams getStatistics
#' @param threshold_pvalue Threshold value for p-values. If present, features
#' with p-values lower than \code{threshold_pvalue} are considered
#' differentially abundant. Set \code{threshold_pvalue = 1} to not filter by
#' p-values.
#' @param threshold_logfc Threshold value for log fold changes. If present,
#' features with log fold change absolute values higher than
#' \code{threshold_logfc} are considered differentially abundant. Set
#' \code{threshold_logfc = 0} to not filter by log fold change values.
#' @param top If not null, the \code{top} number of features, ordered by
#' p-values or log fold change values, are considered as differentially
#' abundant (default \code{top = NULL}).
#' @param verbose Boolean to display the kind of extracted values
#' (default \code{verbose = FALSE}).
#'
#' @return A \code{data.frame} with several columns:
#' \itemize{
#'     \item{\code{stat}}{ which contains the p-values or the absolute log fold
#'     change values;}
#'     \item{\code{direction}}{ which is present if \code{method} was a
#'     \code{data.frame} object, it contains the information about
#'     directionality of differential abundance (usually log fold changes);}
#'     \item{\code{DA}}{ which can contain several values according to
#'     thresholds and inputs. \code{"DA"} or \code{"non-DA"} if \code{method}
#'     object was a vector, \code{"UP Abundant"}, \code{"DOWN Abundant"}, or
#'     \code{"non-DA"} if \code{method} was a \code{data.frame}.}}
#'
#' @seealso \code{\link{getStatistics}}, \code{\link{extractDA}}
#'
#' @examples
#' data("ps_plaque_16S")
#' # Add scaling factors
#' ps_plaque_16S <- norm_edgeR(object = ps_plaque_16S, method = "TMM")
#' # DA analysis
#' da.limma <- DA_limma(
#'     object = ps_plaque_16S,
#'     design = ~ 1 + HMP_BODY_SUBSITE,
#'     coef = 2,
#'     norm = "TMM"
#' )
#' # features with p-value < 0.1 as DA
#' getDA(
#'     method = da.limma, slot = "pValMat", colName = "rawP", type = "pvalue",
#'     direction = NULL, threshold_pvalue = 0.1, threshold_logfc = 0,
#'     top = NULL
#' )
#' # top 10 feature with highest logFC are DA
#' getDA(
#'     method = da.limma, slot = "pValMat", colName = "rawP", type = "pvalue",
#'     direction = "logFC", threshold_pvalue = 1, threshold_logfc = 0, top = 10
#' )
#' # features with p-value < 0.1 and |logFC| > 1 are DA
#' getDA(
#'     method = da.limma, slot = "pValMat", colName = "rawP", type = "pvalue",
#'     direction = "logFC", threshold_pvalue = 0.1, threshold_logfc = 1, top =
#'         NULL
#' )
#' # top 10 features with |logFC| > 1 are DA
#' getDA(
#'     method = da.limma, slot = "pValMat", colName = "rawP", type = "pvalue",
#'     direction = "logFC", threshold_pvalue = 1, threshold_logfc = 1, top = 10
#' )
getDA <- function(method, slot = "pValMat", colName = "rawP", type = "pvalue",
                  direction = NULL, threshold_pvalue = 1, threshold_logfc = 0, top = NULL,
                  verbose = FALSE) {
    out <- data.frame(getStatistics(
        method = method, slot = slot, colName =
            colName, type = type, direction = direction, verbose = verbose
    ))
    out[, "DA"] <- "non-DA"
    msg <- "\n * DA: "
    if (is.null(top)) {
        top <- nrow(out)
    } else {
        msg <- paste0(msg, "top ", top, " features")
    }
    # Only one between the p-values or -|lfc| were supplied
    if (is.null(direction)) {
        colnames(out) <- c("stat", "DA")
        if (type == "pvalue") { # p-value <= threshold_pvalue are DA
            msg <- paste0(msg, ", ", colName, "<=", threshold_pvalue, "\n")
            index_DA <- which(out[, "stat"] <= threshold_pvalue)
        } else if (type == "logfc") { # -|lfc| <= -threshold_logfc are DA
            msg <- paste0(msg, ", |", colName, "|>=", threshold_logfc)
            index_DA <- which(out[, "stat"] <= -threshold_logfc, "\n")
        } else {
            stop("type should be one between 'pvalue' or 'logfc'")
        }
        DA_taxa <- rownames(out)[index_DA]
        out_DA <- out[DA_taxa, ]
        ordered_out <- out_DA[order(out_DA[, "stat"]), ]
        DA_taxa <- rownames(ordered_out)[seq_len(min(length(DA_taxa), top))]
        if (length(index_DA) != 0) { # Check if all non-DA
            out[DA_taxa, "DA"] <- "DA"
        }
        # p-values and lfc, or -|lfc| and lfc were supplied
    } else {
        colnames(out) <- c("stat", "direction", "DA")
        if (type == "pvalue") { # Check if the first column is a p-value
            msg <- paste0(
                msg, ", ", colName, "<=", threshold_pvalue,
                ", |", direction, "|>=", threshold_logfc, "\n"
            )
            index_DA <- intersect(
                which(out[, "stat"] <= threshold_pvalue),
                which(abs(out[, "direction"]) >= threshold_logfc)
            )
        } else if (type == "logfc") { # If the first column is a -|lfc|
            msg <- paste0(msg, ", |", colName, "|>=", threshold_logfc, "\n")
            index_DA <- which(out[, "stat"] <= -threshold_logfc)
        } else {
            stop("type should be one between 'pvalue' or 'logfc'")
        }
        DA_taxa <- rownames(out)[index_DA]
        out_DA <- out[DA_taxa, ]
        ordered_out <- out_DA[order(-abs(out_DA[, "direction"])), ]
        DA_taxa <- rownames(ordered_out)[seq_len(min(length(DA_taxa), top))]
        if (length(index_DA) != 0) { # Check if all non-DA
            out[DA_taxa, "DA"] <- ifelse(test = sign(out[DA_taxa, "direction"])
            == 1, yes = "UP Abundant", no = "DOWN Abundant")
        }
    }
    if (verbose) {
        message(msg, appendLF = FALSE)
    }
    return(out)
}

#' @title extractDA
#'
#' @export
#' @description
#' Inspect the list of p-values or/and log fold changes from the output of
#' differential abundance detection methods.
#'
#' @inheritParams extractStatistics
#' @param threshold_pvalue A single or a numeric vector of thresholds for
#' p-values. If present, features with p-values lower than
#' \code{threshold_pvalue} are considered differentially abundant. Set
#' \code{threshold_pvalue = 1} to not filter by p-values.
#' @param threshold_logfc A single or a numeric vector of thresholds for log
#' fold changes. If present, features with log fold change absolute values
#' higher than \code{threshold_logfc} are considered differentially abundant.
#' Set \code{threshold_logfc = 0} to not filter by log fold change values.
#' @param top If not null, the \code{top} number of features, ordered by
#' p-values or log fold change values, are considered as differentially
#' abundant (default \code{top = NULL}).
#' @param verbose Boolean to display the kind of extracted values
#' (default \code{verbose = FALSE}).
#'
#' @return A \code{data.frame} with several columns for each method:
#' \itemize{
#'     \item{\code{stat}}{ which contains the p-values or the absolute log fold
#'     change values;}
#'     \item{\code{direction}}{ which is present if \code{direction} was
#'     supplied, it contains the information about directionality of
#'     differential abundance (usually log fold changes);}
#'     \item{\code{DA}}{ which can contain several values according to
#'     thresholds and inputs. \code{"DA"} or \code{"non-DA"} if
#'     \code{direction = NULL}, \code{"UP Abundant"}, \code{"DOWN Abundant"}, or
#'     \code{"non-DA"} otherwise.}}
#'
#' @seealso \code{\link{getDA}}, \code{\link{extractStatistics}}
#'
#' @examples
#' data("ps_plaque_16S")
#' # Add scaling factors
#' my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
#'     method = c("TMM", "median"))
#' ps_plaque_16S <- runNormalizations(normalization_list = my_norm,
#'     object = ps_plaque_16S)
#' # Perform DA analysis
#' my_methods <- set_limma(design = ~ 1 + HMP_BODY_SUBSITE, coef = 2,
#'     norm = c("TMM", "CSSmedian"))
#' Plaque_16S_DA <- runDA(method_list = my_methods, object = ps_plaque_16S)
#' # Top 10 features (ordered by 'direction') are DA
#' DA_1 <- extractDA(
#'     object = Plaque_16S_DA, slot = "pValMat", colName = "adjP",
#'     type = "pvalue", direction = "logFC", threshold_pvalue = 1,
#'     threshold_logfc = 0, top = 10
#' )
#' # Features with p-value < 0.05 and |logFC| > 1 are DA
#' DA_2 <- extractDA(
#'     object = Plaque_16S_DA, slot = "pValMat", colName = "adjP",
#'     type = "pvalue", direction = "logFC", threshold_pvalue = 0.05,
#'     threshold_logfc = 1, top = NULL
#' )
extractDA <- function(object, slot = "pValMat", colName = "adjP", type =
                          "pvalue", direction = NULL, threshold_pvalue = 1, threshold_logfc = 0,
                      top = NULL, verbose = FALSE) {
    if (is.null(direction)) {
        if (is.null(top)) {
            out <- mapply(getDA,
                method = object, slot = slot,
                colName = colName, type = type,
                threshold_pvalue = threshold_pvalue,
                threshold_logfc = threshold_logfc, MoreArgs = list(
                    direction =
                        direction, top = top, verbose = verbose
                ), SIMPLIFY = FALSE
            )
        } else {
            out <- mapply(getDA,
                method = object, slot = slot,
                colName = colName, type = type,
                threshold_pvalue = threshold_pvalue,
                threshold_logfc = threshold_logfc, top = top, MoreArgs = list(
                    direction = direction, verbose = verbose
                ), SIMPLIFY = FALSE
            )
        }
    } else {
        if (is.null(top)) {
            out <- mapply(getDA,
                method = object, slot = slot,
                colName = colName, type = type, direction = direction,
                threshold_pvalue = threshold_pvalue,
                threshold_logfc = threshold_logfc, MoreArgs =
                    list(top = top, verbose = verbose), SIMPLIFY = FALSE
            )
        } else {
            out <- mapply(getDA,
                method = object, slot = slot,
                colName = colName, type = type, direction = direction,
                threshold_pvalue = threshold_pvalue,
                threshold_logfc = threshold_logfc, top = top,
                MoreArgs = list(verbose = verbose), SIMPLIFY = FALSE
            )
        }
    }
    names(out) <- unlist(lapply(object, function(method) {
        method[["name"]]
    }))
    return(out)
}

#' @title getPositives
#'
#' @export
#' @description
#' Inspect the list of p-values or/and log fold changes from the output of a
#' differential abundance detection method and count the True Positives (TP) and
#' the False Positives (FP).
#'
#' @param method Output of differential abundance detection method in which
#' DA information is extracted by the \code{getDA} function, information
#' related to enrichment is appropriately added through the \code{addKnowledge}
#' function and the Fisher exact tests is performed for the contingency tables
#' by the \code{enrichmentTests} function.
#' @inheritParams addKnowledge
#' @param TP A list of length-2 vectors. The entries in the vector are the
#' direction ("UP Abundant", "DOWN Abundant", or "non-DA") in the first
#' position, and the level of the enrichment variable (\code{enrichmentCol})
#' which is expected in that direction, in the second position.
#' @param FP A list of length-2 vectors. The entries in the vector are the
#' direction ("UP Abundant", "DOWN Abundant", or "non-DA") in the first
#' position, and the level of the enrichment variable (\code{enrichmentCol})
#' which is not expected in that direction, in the second position.
#'
#' @return A named vector containing the number of TPs and FPs.
#'
#' @seealso \code{\link{createPositives}}.
#'
#' @examples
#' data("ps_plaque_16S")
#' data("microbial_metabolism")
#' # Extract genera from the phyloseq tax_table slot
#' genera <- phyloseq::tax_table(ps_plaque_16S)[, "GENUS"]
#' # Genera as rownames of microbial_metabolism data.frame
#' rownames(microbial_metabolism) <- microbial_metabolism$Genus
#' # Match OTUs to their metabolism
#' priorInfo <- data.frame(genera,
#'     "Type" = microbial_metabolism[genera, "Type"]
#' )
#' # Unmatched genera becomes "Unknown"
#' unknown_metabolism <- is.na(priorInfo$Type)
#' priorInfo[unknown_metabolism, "Type"] <- "Unknown"
#' priorInfo$Type <- factor(priorInfo$Type)
#' # Add a more informative names column
#' priorInfo[, "newNames"] <- paste0(rownames(priorInfo), priorInfo[, "GENUS"])
#'
#' # DA Analysis
#' # Add scaling factors
#' ps_plaque_16S <- norm_edgeR(object = ps_plaque_16S, method = "TMM")
#' # DA analysis
#' da.limma <- DA_limma(
#'     object = ps_plaque_16S,
#'     design = ~ 1 + HMP_BODY_SUBSITE,
#'     coef = 2,
#'     norm = "TMM"
#' )
#'
#' DA <- getDA(
#'     method = da.limma, slot = "pValMat", colName = "adjP",
#'     type = "pvalue", direction = "logFC", threshold_pvalue = 0.05,
#'     threshold_logfc = 1, top = NULL
#' )
#' # Add a priori information
#' DA_info <- addKnowledge(
#'     method = DA, priorKnowledge = priorInfo,
#'     enrichmentCol = "Type", namesCol = "newNames"
#' )
#' # Create contingency tables and compute F tests
#' DA_info_enriched <- enrichmentTest(
#'     method = DA_info, enrichmentCol = "Type",
#'     alternative = "greater"
#' )
#' # Count True and False Positives
#' DA_TP_FP <- getPositives(
#'     method = DA_info_enriched, enrichmentCol = "Type",
#'     TP = list(c("UP Abundant", "Aerobic"), c("DOWN Abundant", "Anaerobic")),
#'     FP = list(c("UP Abundant", "Anaerobic"), c("DOWN Abundant", "Aerobic"))
#' )
getPositives <- function(method, enrichmentCol, TP, FP) {
    data <- method[["data"]][, c("DA", enrichmentCol)]
    # contingency table
    tab <- table(data)
    # True Positives
    TPs <- unlist(lapply(X = TP, FUN = function(x) {
        if (is.element(el = x[1], set = rownames(tab)) &
            is.element(el = x[2], set = colnames(tab))) {
            return(tab[x[1], x[2]])
        } else {
            return(0)
        }
    }))
    # False Positives
    FPs <- unlist(lapply(X = FP, FUN = function(x) {
        if (is.element(el = x[1], set = rownames(tab)) &
            is.element(el = x[2], set = colnames(tab))) {
            return(tab[x[1], x[2]])
        } else {
            return(0)
        }
    }))
    TP_tot <- sum(TPs)
    FP_tot <- sum(FPs)
    return(c("TP" = TP_tot, "FP" = FP_tot))
}

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
#'
#' @examples
#' # Given qualitative variable
#' cond <- factor(c("A", "A", "B", "B", "C", "D"),
#'     levels = c("A", "B", "C", "D"))
#'
#' # Associate a color to each level (or unique value, if not a factor)
#' cond_colors <- createColors(cond)
createColors <- function(variable) {
    if (is.factor(variable)) {
        levels <- levels(variable)
    } else {
        levels <- unique(variable)
    }
    n_levels <- length(levels)
    if (n_levels > 74) {
          stop("To many levels to color them differently.")
      }
    pal.info <- RColorBrewer::brewer.pal.info
    qual_col_pals <- pal.info[pal.info$category == "qual", ]
    col_vector <- unlist(mapply(
        RColorBrewer::brewer.pal,
        qual_col_pals$maxcolors, rownames(qual_col_pals)
    ))
    cols <- col_vector[seq_len(n_levels)]
    names(cols) <- levels
    return(cols)
} # END - function: createColors

#' @title iterativeOrdering
#'
#' @export
#' @importFrom plyr ddply
#' @description Turn the chosen columns (factor) of the input \code{data.frame}
#' into ordered factors. For each factor, the order is given by the number of
#' elements in each level of that factor.
#'
#' @param df a \code{data.frame} object.
#' @param var_names character vector containing the names of one or more columns
#' of \code{df}.
#' @param i iteration index (default \code{i = 1}).
#' @param decreasing logical value or a vector of them. Each value should be
#' associated to a \code{var_name} vector's element. Should the sort order be
#' increasing or decreasing?
#'
#' @return the input \code{data.frame} with the \code{var_names} variables as
#' ordered factors.
#'
#' @seealso \code{\link{plotMutualFindings}}
#'
#' @examples
#' # From a dataset with some factor columns
#' mpg <- data.frame(ggplot2::mpg)
#' # Order the levels of the desired factors based on the cardinality of each
#' # level.
#' ordered_mpg <- iterative_ordering(df = mpg,
#'    var_names = c("manufacturer", "model"),
#'    decreasing = c(TRUE, TRUE))
#' # Now the levels of the factors are ordered in a decreasing way
#' levels(ordered_mpg$manufacturer)
#' levels(ordered_mpg$model)

iterative_ordering <- function(df, var_names, i = 1, decreasing = TRUE) {
    if (length(decreasing) == 1) {
        decreasing <- rep(decreasing, length(var_names))
    } else {
        if (length(decreasing) != length(var_names)) {
              stop("Wrong length of 'decreasing' vector\n")
          }
    }
    # Count elements
    df_var <- plyr::ddply(df, .variables = var_names[i], nrow)
    # Order levels
    ord <- df_var[
        order(df_var[, "V1"], decreasing = decreasing[i]),
        var_names[i]
    ]
    # Rebuild data.frame with the ordered factor
    df[, var_names[i]] <- factor(
        x = df[, var_names[i]], levels = ord,
        ordered = TRUE
    )
    if (i < length(var_names)) { # Iterative
          iterative_ordering(
              df = df, var_names = var_names, i = i + 1,
              decreasing = decreasing
          )
      } else {
        return(df)
    }
}
