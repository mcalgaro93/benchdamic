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
#' # Add some normalization/scaling factors to the phyloseq object
#' my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
#'     method = c("TMM", "median"))
#' ps_plaque_16S <- runNormalizations(normalization_list = my_norm,
#'     object = ps_plaque_16S)
#'
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ 1 + HMP_BODY_SUBSITE, coef = 2,
#'     norm = c("TMM", "CSSmedian"))
#'
#' # Perform DA analysis
#' Plaque_16S_DA <- runDA(method_list = my_limma, object = ps_plaque_16S)
#'
#' # Enrichment analysis
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
#' # Add some normalization/scaling factors to the phyloseq object
#' my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
#'     method = c("TMM", "median"))
#' ps_plaque_16S <- runNormalizations(normalization_list = my_norm,
#'     object = ps_plaque_16S)
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ 1 + HMP_BODY_SUBSITE, coef = 2,
#'     norm = c("TMM", "CSSmedian"))
#'
#' # Perform DA analysis
#' Plaque_16S_DA <- runDA(method_list = my_limma, object = ps_plaque_16S)
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
#'
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
