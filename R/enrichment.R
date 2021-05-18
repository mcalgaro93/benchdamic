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
#'
#' # features with p-value < 0.1 as DA
#' getDA(method = da.limma, slot = "pValMat", colName = "rawP", type = "pvalue",
#'     direction = NULL, threshold_pvalue = 0.1, threshold_logfc = 0,
#'     top = NULL)
#' # top 10 feature with highest logFC are DA
#' getDA(method = da.limma, slot = "pValMat", colName = "rawP", type = "pvalue",
#'     direction = "logFC", threshold_pvalue = 1, threshold_logfc = 0, top = 10)
#' # features with p-value < 0.1 and |logFC| > 1 are DA
#' getDA(method = da.limma, slot = "pValMat", colName = "rawP", type = "pvalue",
#'     direction = "logFC", threshold_pvalue = 0.1, threshold_logfc = 1, top =
#'     NULL)
#' # top 10 features with |logFC| > 1 are DA
#' getDA(method = da.limma, slot = "pValMat", colName = "rawP", type = "pvalue",
#'     direction = "logFC", threshold_pvalue = 1, threshold_logfc = 1, top = 10)

getDA <- function(method, slot = "pValMat", colName = "rawP", type = "pvalue",
    direction = NULL, threshold_pvalue = 1, threshold_logfc = 0, top = NULL,
    verbose = FALSE){
    out <- data.frame(getStatistics(method = method, slot = slot, colName =
        colName, type = type, direction = direction, verbose = verbose))
    out[, "DA"] <- "non-DA"
    msg <- "\n * DA: "
    if(is.null(top)){
        top <- nrow(out)
    } else {
        msg <- paste0(msg, "top ", top, " features")
    }
    # Only one between the p-values or -|lfc| were supplied
    if(is.null(direction)){
        colnames(out) <- c("stat", "DA")
        if(type == "pvalue"){ # p-value <= threshold_pvalue are DA
            msg <- paste0(msg, ", ", colName, "<=", threshold_pvalue, "\n")
            index_DA <- which(out[, "stat"] <= threshold_pvalue)
        } else if(type == "logfc") { # -|lfc| <= -threshold_logfc are DA
            msg <- paste0(msg, ", |", colName, "|>=", threshold_logfc)
            index_DA <- which(out[, "stat"] <= -threshold_logfc, "\n")
        } else stop("type should be one between 'pvalue' or 'logfc'")
        DA_taxa <- rownames(out)[index_DA]
        out_DA <- out[DA_taxa, ]
        ordered_out <- out_DA[order(out_DA[, "stat"]), ]
        DA_taxa <- rownames(ordered_out)[seq_len(min(length(DA_taxa), top))]
        if(length(index_DA) != 0){ # Check if all non-DA
            out[DA_taxa, "DA"] <- "DA"
        }
    # p-values and lfc, or -|lfc| and lfc were supplied
    } else {
        colnames(out) <- c("stat", "direction", "DA")
        if(type == "pvalue"){ # Check if the first column is a p-value
            msg <- paste0(msg, ", ", colName, "<=", threshold_pvalue,
                ", |", direction, "|>=", threshold_logfc, "\n")
            index_DA <- intersect(which(out[, "stat"] <= threshold_pvalue),
                which(abs(out[, "direction"]) >= threshold_logfc))
        } else if(type == "logfc") { # If the first column is a -|lfc|
            msg <- paste0(msg, ", |", colName, "|>=", threshold_logfc, "\n")
            index_DA <- which(out[, "stat"] <= -threshold_logfc)
        } else stop("type should be one between 'pvalue' or 'logfc'")
        DA_taxa <- rownames(out)[index_DA]
        out_DA <- out[DA_taxa, ]
        ordered_out <- out_DA[order(-abs(out_DA[, "direction"])), ]
        DA_taxa <- rownames(ordered_out)[seq_len(min(length(DA_taxa), top))]
        if(length(index_DA) != 0){ # Check if all non-DA
            out[DA_taxa, "DA"] <- ifelse(test = sign(out[DA_taxa, "direction"])
                == 1, yes = "UP Abundant", no = "DOWN Abundant")
        }
    }
    if(verbose){
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
#' # Top 10 features (ordered by 'direction') are DA
#' DA_1 = extractDA(object = Plaque_16S_DA, slot = "pValMat", colName = "adjP",
#'     type = "pvalue", direction = "logFC", threshold_pvalue = 1,
#'     threshold_logfc = 0, top = 10)
#' # Features with p-value < 0.05 and |logFC| > 1 are DA
#' DA_2 = extractDA(object = Plaque_16S_DA, slot = "pValMat", colName = "adjP",
#'     type = "pvalue", direction = "logFC", threshold_pvalue = 0.05,
#'     threshold_logfc = 1, top = NULL)

extractDA <- function(object, slot = "pValMat", colName = "adjP", type =
    "pvalue", direction = NULL, threshold_pvalue = 1, threshold_logfc = 0,
    top = NULL, verbose = FALSE){
    if(is.null(direction)){
        if(is.null(top)){
            out <- mapply(getDA, method = object, slot = slot,
                colName = colName, type = type,
                threshold_pvalue = threshold_pvalue,
                threshold_logfc = threshold_logfc, MoreArgs = list(direction =
                direction, top = top, verbose = verbose), SIMPLIFY = FALSE)
        } else {
            out <- mapply(getDA, method = object, slot = slot,
                colName = colName, type = type,
                threshold_pvalue = threshold_pvalue,
                threshold_logfc = threshold_logfc, top = top, MoreArgs = list(
                direction = direction, verbose = verbose), SIMPLIFY = FALSE)
        }
    } else {
        if(is.null(top)){
            out <- mapply(getDA, method = object, slot = slot,
                colName = colName, type = type, direction = direction,
                threshold_pvalue = threshold_pvalue,
                threshold_logfc = threshold_logfc, MoreArgs =
                list(top = top, verbose = verbose), SIMPLIFY = FALSE)
        } else {
            out <- mapply(getDA, method = object, slot = slot,
                colName = colName, type = type, direction = direction,
                threshold_pvalue = threshold_pvalue,
                threshold_logfc = threshold_logfc, top = top,
                MoreArgs = list(verbose = verbose), SIMPLIFY = FALSE)
        }
    }
    names(out) <- unlist(lapply(object, function(method){method[["name"]]}))
    return(out)
}

#' @title addKnowledge
#'
#' @export
#' @description
#' Add a priori knowledge for each feature tested by a method.
#'
#' @param method Output of differential abundance detection method in which
#' DA information is extracted by the \code{getDA} function.
#' @param priorKnowledge \code{data.frame} (with feature names as
#' \code{row.names}) containing feature level metadata.
#' @param enrichmentCol name of the column containing information for enrichment
#' analysis.
#' @param namesCol name of the column containing new names for features (default
#' \code{namesCol = NULL}).
#'
#' @return A \code{data.frame} with a new column containing information for
#' enrichment analysis.
#'
#' @seealso \code{\link{createEnrichment}}.
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
#' DA <- getDA(method = da.limma, slot = "pValMat", colName = "adjP",
#'     type = "pvalue", direction = "logFC", threshold_pvalue = 0.05,
#'     threshold_logfc = 1, top = NULL)
#' # Add a priori information
#' DA_info <- addKnowledge(method = DA, priorKnowledge = priorInfo,
#'     enrichmentCol = "Type", namesCol = "newNames")

addKnowledge <- function(method, priorKnowledge, enrichmentCol, namesCol = NULL)
{
    if(!is.element(enrichmentCol, colnames(priorKnowledge))){
        stop(paste0(enrichmentCol, " column is not present."))
    }
    if(!is.null(namesCol)){
        if(!is.element(namesCol, colnames(priorKnowledge))){
            stop(paste0(namesCol, " column is not present."))
        }
    }
    if(!is.data.frame(priorKnowledge)){
        stop("'priorKnowledge' must be a data.frame object.")
    }
    if(is.null(rownames(priorKnowledge))){
        stop("'priorKnowledge' argument must have the row names.")
    }
    out <- method
    out[, enrichmentCol] <- NA
    out[, enrichmentCol] <- priorKnowledge[rownames(method), enrichmentCol]
    if(is.null(namesCol)){
        out[, "feature"] <- rownames(out)
    } else out[, "feature"] <- priorKnowledge[rownames(method), namesCol]
    return(out)
}

#' @title enrichmentTest
#'
#' @export
#' @importFrom plyr llply
#' @importFrom stats fisher.test
#' @description
#' Perform the Fisher exact test for all the possible 2x2 contingency tables,
#' considering differential abundance direction and enrichment variable.
#'
#' @param method Output of differential abundance detection method in which
#' DA information is extracted by the \code{getDA} function and the information
#' related to enrichment is appropriately added through the \code{addKnowledge}.
#' @inheritParams addKnowledge
#' @inheritParams extractDA
#' @inheritParams stats::fisher.test
#'
#' @return a list of objects:
#' \itemize{
#'     \item{\code{data}}{ a \code{data.frame} object with DA directions,
#'     statistics, and feature names;}
#'     \item{\code{tables}}{ a list of 2x2 contingency tables;}
#'     \item{\code{tests}}{ the list of Fisher exact tests' p-values for each
#'     contingency table;}
#'     \item{\code{summaries}}{ a list with the first element of each
#'     contingency table and its p-value (for graphical purposes);}}
#'
#' @seealso \code{\link{extractDA}}, \code{\link{addKnowledge}}, and
#' \code{\link{createEnrichment}}
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
#' DA <- getDA(method = da.limma, slot = "pValMat", colName = "adjP",
#'     type = "pvalue", direction = "logFC", threshold_pvalue = 0.05,
#'     threshold_logfc = 1, top = NULL)
#' # Add a priori information
#' DA_info <- addKnowledge(method = DA, priorKnowledge = priorInfo,
#'     enrichmentCol = "Type", namesCol = "newNames")
#' # Create contingency tables and compute F tests
#' DA_info_enriched <- enrichmentTest(method = DA_info, enrichmentCol = "Type",
#'     alternative = "greater")

enrichmentTest <- function(method, enrichmentCol, alternative = "greater"){
    # Create enrichment tables and perform Fisher exact tests
    tab <- table(method[, c("DA", enrichmentCol)])
    # Only DA, UP Abundant or DOWN Abundant are of interest
    rows <- rownames(tab)
    row_levels <- rows[rows != "non-DA"]
    # All levels of enrichmentCol are of interest
    col_levels <- colnames(tab)
    out <- list() # output tables
    test <- list() # output F tests
    summary <- list() # output of summary info
    out_tab <- matrix(NA, nrow = 2, ncol = 2)
    for(row in row_levels){
        for(col in col_levels){ # Create enrichment table
            out_tab[1, 1] <- tab[row, col]
            out_tab[1, 2] <- sum(tab[row, which(colnames(tab) != col)])
            out_tab[2, 1] <- sum(tab[which(rownames(tab) != row), col])
            out_tab[2, 2] <- sum(tab[which(rownames(tab) != row), which(
                colnames(tab) != col)])
            rownames(out_tab) <- c(row, paste0("non-", row))
            colnames(out_tab) <- c(col, paste0("non-", col))
            # Fisher exact test
            f_test <- fisher.test(out_tab,
                alternative = alternative)[["p.value"]]
            out[[paste0(row, "-", col)]] <- data.frame(out_tab)
            test[[paste0(row, "-", col)]] <- f_test
            # Only DA number and p-value
            summary_df <- data.frame(out_tab[1,1], f_test)
            colnames(summary_df) <- c(col, "pvalue")
            rownames(summary_df) <- row
            summary[[paste0(row, "-", col)]] <- summary_df
        }
    }
    return(list("data" = method, "tables" = out, "tests" = test,
        "summaries" = summary))
}

#' @title createEnrichment
#'
#' @export
#' @importFrom plyr ldply
#' @description
#' Create a \code{data.frame} object with several information to perform
#' enrichment analysis.
#'
#' @inheritParams addKnowledge
#' @inheritParams extractDA
#' @inheritParams enrichmentTest
#'
#' @return a list of objects for each method. Each list contains:
#' \itemize{
#'     \item{\code{data}}{ a \code{data.frame} object with DA directions,
#'     statistics, and feature names;}
#'     \item{\code{tables}}{ a list of 2x2 contingency tables;}
#'     \item{\code{tests}}{ the list of Fisher exact tests' p-values for each
#'     contingency table;}
#'     \item{\code{summaries}}{ a list with the first element of each
#'     contingency table and its p-value (for graphical purposes);}}
#'
#' @seealso \code{\link{addKnowledge}}, \code{\link{extractDA}}, and
#' \code{\link{enrichmentTest}}.
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
#' # DA analysis
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
#' enrichment <- createEnrichment(object = Plaque_16S_DA,
#'     priorKnowledge = priorInfo, enrichmentCol = "Type", namesCol = "GENUS",
#'     slot = "pValMat", colName = "adjP", type = "pvalue", direction = "logFC",
#'     threshold_pvalue = 0.1, threshold_logfc = 1, top = 10, verbose = TRUE)

createEnrichment <- function(object, priorKnowledge, enrichmentCol,
    namesCol = NULL, slot = "pValMat", colName = "adjP", type = "pvalue",
    direction = NULL, threshold_pvalue = 1, threshold_logfc = 0, top = NULL,
    alternative = "greater", verbose = FALSE){
    # Extract DA information
    object_DA <- extractDA(object = object, slot = slot, colName = colName,
        type = type, direction = direction, threshold_pvalue = threshold_pvalue,
        threshold_logfc = threshold_logfc, top = top, verbose = verbose)
    # Add priorKnowledge to previously produced object
    object_DA_info <- mapply(addKnowledge, method = object_DA, MoreArgs = list(
        priorKnowledge = priorKnowledge, enrichmentCol = enrichmentCol,
        namesCol = namesCol), SIMPLIFY = FALSE)
    # Create contingency tables and perform Fisher exact tests on them
    object_enriched <- lapply(X = object_DA_info, FUN = enrichmentTest,
        enrichmentCol = enrichmentCol, alternative = alternative)
    return(object_enriched)
}

#' @title plotEnrichment
#'
#' @export
#' @importFrom plyr ldply ddply
#' @importFrom ggplot2 ggplot aes geom_col geom_text position_dodge theme
#' @importFrom ggplot2 sec_axis geom_hline element_text scale_y_continuous
#' @importFrom ggplot2 scale_x_discrete guide_legend ggtitle guides
#' @importFrom ggplot2 element_blank
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
#' # DA analysis
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
#' enrichment <- createEnrichment(object = Plaque_16S_DA,
#'     priorKnowledge = priorInfo, enrichmentCol = "Type", namesCol = "GENUS",
#'     slot = "pValMat", colName = "adjP", type = "pvalue", direction = "logFC",
#'     threshold_pvalue = 0.1, threshold_logfc = 1, top = 10, verbose = TRUE)
#' # Contingency tables
#' plotContingency(enrichment = enrichment, method = "limma.TMM")
#' # Barplots
#' plotEnrichment(enrichment, enrichmentCol = "Type")
#' # Mutual findings
#' plotMutualFindings(enrichment = enrichment, enrichmentCol = "Type",
#'     n_methods = 1)

plotEnrichment <- function(enrichment, enrichmentCol, levels_to_plot){
    # Extract freq and p-values
    df_to_plot <- plyr::ldply(.data = enrichment, .fun = function(x){
        plyr::ldply(x[["summaries"]], .fun = function(tab){
            out <- cbind(tab, "direction" = rownames(tab),
                "variable" = colnames(tab)[1])
            colnames(out) <- c("value", "pvalue", "direction", "variable")
            return(out)
        }, .id = NULL)
    },.id = "method")
    if(!missing(levels_to_plot)){ # Check which levels to plot
        keep <- is.element(el = df_to_plot[ , "variable"],
            set = levels_to_plot)
        df_to_plot <- df_to_plot[keep, ]
    }
    # Order methods
    method <- direction <- value <- variable <- pvalue <- NULL
    method_DA <- plyr::ddply(.data = df_to_plot, .variables = c("method"),
        function(x) c("freq" = sum(x[, "value"])))
    ord <- method_DA[order(method_DA[, "freq"], decreasing = TRUE), "method"]
    max_DA <- max(abs(df_to_plot[, "value"])) # Max freq available
    g_enrichment <- ggplot2::ggplot(data = df_to_plot,
        mapping = ggplot2::aes(x = method,
            y = ifelse(test = (direction == "UP Abundant" | direction == "DA"),
            yes = value, no = -value), fill = variable)) +
        ggplot2::geom_col(position = "dodge") +
        ggplot2::geom_text(mapping = ggplot2::aes(x = method,
            y = ifelse(test = (direction == "UP Abundant" | direction == "DA"),
                yes = value + 5, no =  -value - 5),
            label = ifelse(test = (pvalue > 0.1), yes = "", no = ifelse(
                test = (pvalue > 0.05), yes = ".", no = ifelse(
                    test = (pvalue > 0.01), yes = "*", no = ifelse(
                        test = (pvalue > 0.001), yes = "**", no = "***")))),
            color = variable), position = ggplot2::position_dodge(width = 1)) +
        ggplot2::geom_hline(mapping = ggplot2::aes(yintercept = 0),
            color = "black") +
        ggplot2::scale_x_discrete(limits = ord) +
        ggplot2::ylab("Number of DA taxa") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
            hjust = 1, vjust = 0.5),
            axis.text.y.right = ggplot2::element_text(angle = -90, hjust = 0.5),
            axis.ticks.y.right = ggplot2::element_blank(),
            legend.position = "bottom") +
        ggplot2::guides(fill = ggplot2::guide_legend(title = enrichmentCol,
            override.aes = ggplot2::aes(label = "")), color = "none") +
        ggplot2::ggtitle("Enrichment analysis", subtitle =
            "Fisher exact test pvalue: 0 *** 0.001 ** 0.01 * 0.05 . 0.1")
    if(!is.element(el = "DA", set = df_to_plot[, "direction"]))
        g_enrichment <- g_enrichment +
            ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(
            trans = ~./max_DA, name = "DA direction", breaks = c(-0.5, 0.5),
            labels = c("DOWN Abundant","UP Abundant")))
    return(g_enrichment)
}

#' @title plotContingency
#'
#' @export
#' @importFrom plyr ldply
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes facet_grid geom_tile geom_text theme
#' @importFrom ggplot2 element_text ggtitle scale_fill_distiller guides
#' @importFrom ggplot2 guide_colorbar
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
#' # DA analysis
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
#' enrichment <- createEnrichment(object = Plaque_16S_DA,
#'     priorKnowledge = priorInfo, enrichmentCol = "Type", namesCol = "GENUS",
#'     slot = "pValMat", colName = "adjP", type = "pvalue", direction = "logFC",
#'     threshold_pvalue = 0.1, threshold_logfc = 1, top = 10, verbose = TRUE)
#' # Contingency tables
#' plotContingency(enrichment = enrichment, method = "limma.TMM")
#' # Barplots
#' plotEnrichment(enrichment, enrichmentCol = "Type")
#' # Mutual findings
#' plotMutualFindings(enrichment = enrichment, enrichmentCol = "Type",
#'     n_methods = 1)

plotContingency <- function(enrichment, method, levels_to_plot){
    # Extract enrichment tables
    enrichment_table <- enrichment[[method]][["tables"]]
    if(length(enrichment_table) == 0)
        stop(paste0("No DA features for ", method))
    df_to_plot <- reshape2::melt(plyr::ldply(enrichment_table,
        .fun = function(tab){
            out <- cbind("direction" = rownames(tab), tab,
                "level" = colnames(tab)[1])
        return(out)
    }), na.rm = TRUE)
    df_to_plot$dir <- gsub(pattern = "non-", replacement = "",
        x = df_to_plot[,"direction"])
    if(!missing(levels_to_plot)){ # Check which levels to plot
        keep <- is.element(el = df_to_plot[ , "level"], set = levels_to_plot)
        df_to_plot <- df_to_plot[keep, ]
    }
    df_to_plot[, "direction"] <- factor(df_to_plot[, "direction"], levels =
        c("non-DOWN Abundant", "DOWN Abundant", "non-UP Abundant",
        "UP Abundant", "non-DA", "DA"), ordered = TRUE)
    # Extract p-values
    pvalue_table <- plyr::ldply(enrichment[[method]][["tests"]])
    colnames(pvalue_table) <- c("direction_variable", "pvalue")
    pvalue_table[, c("direction","variable")] <- plyr::ldply(strsplit(
        pvalue_table[, "direction_variable"], split = "-"))
    df_to_plot <- merge(x = df_to_plot, y = pvalue_table,
        by.x = c("variable","direction"), by.y = c("variable","direction"),
        all.x = TRUE)
    # Plot
    variable <- direction <- pvalue <- value <- NULL
    ggplot2::ggplot(data = df_to_plot, mapping = ggplot2::aes(x = variable,
        y = direction, fill = pvalue)) +
        ggplot2::facet_grid(dir ~ level, scales = "free") +
        ggplot2::geom_tile(width = 0.8, height = 0.8) +
        ggplot2::geom_text(mapping = ggplot2::aes(label = round(x = value,
            digits = 2)), color = "black") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
            hjust = 1, vjust = 0.5)) +
        ggplot2::ggtitle(label = paste0("Enrichment tables - ", method),
            subtitle = "Colored by Fisher p-value") +
        ggplot2::scale_fill_distiller(type = "seq", palette = "RdYlBu",
            breaks = seq.int(from = 0, to = 0.2, by = 0.02), limits = c(0,0.1),
            direction = 1) +
        ggplot2::guides(fill = ggplot2::guide_colorbar(ticks.colour = "black",
            frame.colour = "black"))
}

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

iterative_ordering <- function(df, var_names, i = 1, decreasing = TRUE){
    if(length(decreasing) == 1){
        decreasing = rep(decreasing, length(var_names))
    } else {
        if(length(decreasing) != length(var_names))
            stop("Wrong length of 'decreasing' vector\n")
    }
    # Count elements
    df_var <- plyr::ddply(df, .variables = var_names[i], nrow)
    # Order levels
    ord <- df_var[order(df_var[, "V1"], decreasing = decreasing[i]),
        var_names[i]]
    # Rebuild data.frame with the ordered factor
    df[, var_names[i]] <- factor(x = df[, var_names[i]], levels = ord,
        ordered = TRUE)
    if(i < length(var_names)) # Iterative
        iterative_ordering(df = df, var_names = var_names, i = i + 1,
            decreasing = decreasing)
    else return(df)
}

#' @title plotMutualFindings
#'
#' @export
#' @importFrom plyr ldply ddply
#' @importFrom ggplot2 ggplot aes facet_grid geom_tile theme
#' @importFrom ggplot2 element_text ggtitle
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
#' # DA analysis
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
#' enrichment <- createEnrichment(object = Plaque_16S_DA,
#'     priorKnowledge = priorInfo, enrichmentCol = "Type", namesCol = "GENUS",
#'     slot = "pValMat", colName = "adjP", type = "pvalue", direction = "logFC",
#'     threshold_pvalue = 0.1, threshold_logfc = 1, top = 10, verbose = TRUE)
#' # Contingency tables
#' plotContingency(enrichment = enrichment, method = "limma.TMM")
#' # Barplots
#' plotEnrichment(enrichment, enrichmentCol = "Type")
#' # Mutual findings
#' plotMutualFindings(enrichment = enrichment, enrichmentCol = "Type",
#'     n_methods = 1)

plotMutualFindings <- function(enrichment, enrichmentCol, levels_to_plot,
    n_methods = 1){
    DA_df <- plyr::ldply(enrichment, function(x) x[["data"]], .id = "method")
    df_to_plot <- DA_df[DA_df[, "DA"] != "non-DA",]
    if(n_methods > 1){ # Filter only DA mutually found by n_methods methods
        feature_df <- plyr::ddply(df_to_plot, .variables = "feature", nrow)
        feature_to_keep <- feature_df[which(feature_df[, "V1"] >= n_methods),
            "feature"]
        keep <- is.element(el = df_to_plot[, "feature"], set = feature_to_keep)
        df_to_plot <- df_to_plot[keep, ]
    }
    if(!missing(levels_to_plot)){ # Check which levels to plot
        keep <- is.element(el = df_to_plot[, enrichmentCol],
            set = levels_to_plot)
        df_to_plot <- df_to_plot[keep, ]
    }
    df_to_plot <- iterative_ordering(df = df_to_plot,
        var_names = c("method", "feature", enrichmentCol),
        decreasing = c(TRUE,FALSE,TRUE))
    method <- feature <- DA <- NULL
    ggplot2::ggplot(data = df_to_plot, mapping = ggplot2::aes(x = method,
        y = feature, fill = DA)) +
        ggplot2::facet_grid(get(enrichmentCol) ~ . , scales = "free",
            space = "free", drop = TRUE) +
        ggplot2::geom_tile(width = 0.8, height = 0.8) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
            hjust = 1, vjust = 0.5),
            strip.text.y = ggplot2::element_text(angle = 0)) +
        ggplot2::ggtitle(label = "Summary of DA features",
            subtitle = paste0("By ", enrichmentCol))
}

#' @title createPositives
#'
#' @export
#' @description
#' Inspect the list of p-values or/and log fold changes from the output of the
#' differential abundance detection methods and count the True Positives (TP)
#' and the False Positives (FP).
#'
#' @inheritParams createEnrichment
#' @inheritParams getPositives
#'
#' @return a \code{data.frame} object which contains the number of TPs and FPs
#' features for each method and for each threshold of the \code{top} argument.
#'
#' @seealso \code{\link{getPositives}}, \code{\link{plotPositives}}.
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
#' # DA analysis
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
#' # Count TPs and FPs, from the top 1 to the top 20 features.
#' # As direction is supplied, features are ordered by "logFC" absolute values.
#' positives <- createPositives(object = Plaque_16S_DA,
#' priorKnowledge = priorInfo, enrichmentCol = "Type", namesCol = "newNames",
#' slot = "pValMat", colName = "rawP", type = "pvalue", direction = "logFC",
#' threshold_pvalue = 1, threshold_logfc = 0, top = 1:20,
#' alternative = "greater", verbose = FALSE,
#' TP = list(c("DOWN Abundant", "Anaerobic"), c("UP Abundant", "Aerobic")),
#' FP = list(c("DOWN Abundant", "Aerobic"), c("UP Abundant", "Anaerobic")))
#' # Plot the TP-FP differences for each threshold
#' plotPositives(positives = positives)

createPositives <- function(object, priorKnowledge, enrichmentCol,
    namesCol = NULL, slot = "pValMat", colName = "adjP", type = "pvalue",
    direction = NULL, threshold_pvalue = 1, threshold_logfc = 0, top = NULL,
    alternative = "greater", verbose = FALSE, TP, FP){
    if(is.null(top)){
        top <- seq.int(from = 0, to = nrow(object[[1]][[slot]]), by = 10)
    }
    top <- top[top != 0 & top <= nrow(object[[1]][[slot]])]
    enrichment_list <- mapply(FUN = createEnrichment, top = top, MoreArgs =
        list(object = object, priorKnowledge = priorKnowledge,
        enrichmentCol = enrichmentCol, namesCol = namesCol, slot = slot,
        colName = colName, type = type, direction = direction,
        threshold_pvalue = threshold_pvalue, threshold_logfc = threshold_logfc,
        alternative, verbose = verbose), SIMPLIFY = FALSE)
    names(enrichment_list) <- top
    out <- ldply(.data = enrichment_list, .fun = function(topN){
        ldply(.data = topN, .fun = function(method){
            getPositives(method = method, enrichmentCol = enrichmentCol,
                TP = TP, FP = FP)
        }, .id = "method")
    }, .id = "top")
    out[, "top"] <- as.numeric(as.character(out[, "top"]))
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
#' DA <- getDA(method = da.limma, slot = "pValMat", colName = "adjP",
#'     type = "pvalue", direction = "logFC", threshold_pvalue = 0.05,
#'     threshold_logfc = 1, top = NULL)
#' # Add a priori information
#' DA_info <- addKnowledge(method = DA, priorKnowledge = priorInfo,
#'     enrichmentCol = "Type", namesCol = "newNames")
#' # Create contingency tables and compute F tests
#' DA_info_enriched <- enrichmentTest(method = DA_info, enrichmentCol = "Type",
#'     alternative = "greater")
#' # Count True and False Positives
#' DA_TP_FP <- getPositives(method = DA_info_enriched, enrichmentCol = "Type",
#'     TP = list(c("UP Abundant", "Aerobic"), c("DOWN Abundant", "Anaerobic")),
#'     FP = list(c("UP Abundant", "Anaerobic"), c("DOWN Abundant", "Aerobic")))

getPositives <- function(method, enrichmentCol, TP, FP){
    data <- method[["data"]][, c("DA", enrichmentCol)]
    # contingency table
    tab <- table(data)
    # True Positives
    TPs <- unlist(lapply(X = TP, FUN = function(x){
        if(is.element(el = x[1], set = rownames(tab)) &
            is.element(el = x[2], set = colnames(tab))){
            return(tab[x[1],x[2]])
        } else return(0)
    }))
    # False Positives
    FPs <- unlist(lapply(X = FP, FUN = function(x){
        if(is.element(el = x[1], set = rownames(tab)) &
            is.element(el = x[2], set = colnames(tab))){
            return(tab[x[1],x[2]])
        } else return(0)
    }))
    TP_tot <- sum(TPs)
    FP_tot <- sum(FPs)
    return(c("TP" = TP_tot, "FP" = FP_tot))
}

#' @title plotPositives
#'
#' @export
#' @importFrom plyr ldply ddply
#' @importFrom ggplot2 ggplot aes geom_path scale_color_manual ggtitle
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
#' # DA analysis
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
#' # Count TPs and FPs, from the top 1 to the top 20 features.
#' # As direction is supplied, features are ordered by "logFC" absolute values.
#' positives <- createPositives(object = Plaque_16S_DA,
#' priorKnowledge = priorInfo, enrichmentCol = "Type", namesCol = "newNames",
#' slot = "pValMat", colName = "rawP", type = "pvalue", direction = "logFC",
#' threshold_pvalue = 1, threshold_logfc = 0, top = 1:20,
#' alternative = "greater", verbose = FALSE,
#' TP = list(c("DOWN Abundant", "Anaerobic"), c("UP Abundant", "Aerobic")),
#' FP = list(c("DOWN Abundant", "Aerobic"), c("UP Abundant", "Anaerobic")))
#' # Plot the TP-FP differences for each threshold
#' plotPositives(positives = positives)

plotPositives <- function(positives, cols = NULL){
    if(is.null(cols))
        cols <- createColors(positives[, "method"])
    top <- TP <- FP <- method <- NULL
    ggplot2::ggplot(data = positives, mapping = ggplot2::aes(x = top,
        y = TP-FP, color = method)) +
        ggplot2::geom_point() +
        ggplot2::geom_path(mapping = ggplot2::aes(group = method)) +
        ggplot2::scale_color_manual(values = cols) +
        ggplot2::ggtitle(label = "Putative TP - Putative FP",
            subtitle = paste0("From ", min(positives[, "top"]), " to ",
            max(positives[, "top"]), " top features."))
}

