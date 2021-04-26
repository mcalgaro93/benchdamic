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
#'
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
#' @param method Extracted statistics and DA information produced by a
#' differential abundance detection method.
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
#' DA <- extractDA(object = Plaque_16S_DA, slot = "pValMat", colName = "adjP",
#'     type = "pvalue", direction = "logFC", threshold_pvalue = 0.05,
#'     threshold_logfc = 1, top = NULL)
#' # Add a priori information only for a method
#' da.limma.css_info <- addKnowledge(method = DA[[1]],
#'     priorKnowledge = priorInfo, enrichmentCol = "Type",
#'     namesCol = "newNames")

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
#'
#' @return A \code{data.frame} object containing statistics, direction (if not
#' \code{NULL}), \code{enrichmentCol}, and feature names for each method.
#'
#' @seealso \code{\link{addKnowledge}}, \code{\link{extractDA}}.
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
#' DA <- extractDA(object = Plaque_16S_DA, slot = "pValMat", colName = "adjP",
#'     type = "pvalue", direction = "logFC", threshold_pvalue = 0.05,
#'     threshold_logfc = 1, top = NULL)
#' # Add a priori information only for a method
#' da.limma.css_info <- addKnowledge(method = DA[[1]],
#'     priorKnowledge = priorInfo, enrichmentCol = "Type",
#'     namesCol = "newNames")
#' # get a data.frame ready for enrichment analysis for all methods
#' enrichment <- createEnrichment(object = Plaque_16S_DA,
#'     priorKnowledge = priorInfo, enrichmentCol = "Type", namesCol = "GENUS",
#'     slot = "pValMat", colName = "adjP", type = "pvalue", direction = "logFC",
#'     threshold_pvalue = 0.1, threshold_logfc = 1, top = 10, verbose = TRUE)

createEnrichment <- function(object, priorKnowledge, enrichmentCol,
    namesCol = NULL, slot = "pValMat", colName = "adjP", type = "pvalue",
    direction = NULL, threshold_pvalue = 1, threshold_logfc = 0, top = NULL,
    verbose = FALSE){
    object_DA <- extractDA(object = object, slot = slot, colName = colName,
        type = type, direction = direction, threshold_pvalue = threshold_pvalue,
        threshold_logfc = threshold_logfc, top = top, verbose = verbose)
    out_list <- mapply(addKnowledge, method = object_DA, MoreArgs = list(
        priorKnowledge = priorKnowledge, enrichmentCol = enrichmentCol,
        namesCol = namesCol), SIMPLIFY = FALSE)
    out_df <- plyr::ldply(.data = out_list, .id = "method")
    return(out_df)
}

