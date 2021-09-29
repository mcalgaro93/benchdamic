#' @title DA_Seurat
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data
#' @importFrom Seurat CreateSeuratObject AddMetaData NormalizeData
#' @importFrom Seurat FindVariableFeatures ScaleData FindMarkers
#' @importFrom utils capture.output
#' @export
#' @description
#' Fast run for Seurat differential abundance detection method.
#'
#' @inheritParams DA_DESeq2
#' @inheritParams Seurat::FindMarkers
#'
#' @return A list object containing the matrix of p-values `pValMat`, the matrix
#' of summary statistics for each tag `statInfo`, and a suggested `name` of the
#' final object considering the parameters passed to the function.
#'
#' @seealso \code{\link[Seurat]{CreateSeuratObject}} to create the Seurat
#' object, \code{\link[Seurat]{AddMetaData}} to add metadata information,
#' \code{\link[Seurat]{NormalizeData}} to compute the normalization for the
#' counts, \code{\link[Seurat]{FindVariableFeatures}} to estimate the
#' mean-variance trend, \code{\link[Seurat]{ScaleData}} to scale and center
#' features in the dataset, and \code{\link[Seurat]{FindMarkers}} to perform
#' differential abundance analysis.
#'
#' @examples
#' set.seed(1)
#' # Create a very simple phyloseq object
#' counts <- matrix(rnbinom(n = 60, size = 3, prob = 0.5), nrow = 10, ncol = 6)
#' metadata <- data.frame("Sample" = c("S1", "S2", "S3", "S4", "S5", "S6"),
#'                        "group" = as.factor(c("A", "A", "A", "B", "B", "B")))
#' ps <- phyloseq::phyloseq(phyloseq::otu_table(counts, taxa_are_rows = TRUE),
#'                          phyloseq::sample_data(metadata))
#' # No use of scaling factors
#' ps_NF <- norm_edgeR(object = ps, method = "none")
#' # The phyloseq object now contains the scaling factors:
#' scaleFacts <- phyloseq::sample_data(ps_NF)[, "NF.none"]
#' head(scaleFacts)
#' # Differential abundance
#' DA_Seurat(object = ps_NF, contrast = c("group","B","A"), norm = "none")

DA_Seurat <- function(object, pseudo_count = FALSE, test.use = "wilcox",
    contrast, norm = c("TMM", "TMMwsp", "RLE", "upperquartile",
    "posupperquartile", "none", "ratio", "poscounts", "iterate", "TSS",
    "CSSmedian", "CSSdefault"), verbose = TRUE){
    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)
    # Name building
    name <- "Seurat"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        if(verbose)
            message("Adding a pseudo count...")
        counts <- counts + 1
        name <- paste(name,".pseudo",sep = "")}
    if(length(norm) > 1)
        stop("Please choose one normalization for this istance of differential",
            " abundance analysis.")
    NF.col <- paste("NF", norm, sep = ".")
    # Check if the column with the normalization factors is present
    if(!any(colnames(metadata) == NF.col)){
        stop("Can't find the ", NF.col," column in your object. Make sure to",
            " add the normalization factors column in your object first.")}
    name <- paste(name, ".", norm, sep = "")
    NFs = unlist(metadata[,NF.col])
    # Check if the NFs are scaling factors. If so, make them norm. factors
    if(is.element(norm, c("TMM", "TMMwsp", "RLE", "upperquartile",
        "posupperquartile", "CSSmedian", "CSSdefault", "TSS")))
        NFs <- NFs * colSums(counts)
    NFs <- NFs/exp(mean(log(NFs)))
    norm_counts <- round(counts %*% diag(1/NFs), digits = 0)
    colnames(norm_counts) <- colnames(counts)
    if(verbose)
        message("Differential abundance on ", norm," normalized data")
    # Initialize the Seurat object with the raw (non-normalized data).
    # If the chosen normalization is 'none', then this is true.
    # Keep all features expressed in >= 1 sample
    # Keep all samples with at least 1 detected feature.
    sobj <- Seurat::CreateSeuratObject(counts = norm_counts, min.cells = 1,
        min.features = 1)
    sobj <- Seurat::AddMetaData(object = sobj, metadata = data.frame(metadata),
        col.name = colnames(metadata))
    if(missing(contrast) | (!is.character(contrast) & length(contrast) != 3)){
        stop("Please supply a character vector with exactly three elements:",
            " the name of a factor in the design formula, the name of the",
            " numerator level for the fold change, and the name of the",
            " denominator level for the fold change.")
    } else {
        if(!is(unlist(sobj[[contrast[1]]]), "factor")){
            stop(contrast[1],
                 " variable is not a factor. Please supply a factor.")
        } else {
            if(!is.element(contrast[2], levels(unlist(sobj[[contrast[1]]]))))
                stop(contrast[2], "is not a level of the", contrast[1],
                           "variable. Please supply a present category.")
            if(!is.element(contrast[3], levels(unlist(sobj[[contrast[1]]]))))
                stop(contrast[3], "is not a level of the", contrast[1],
                           "variable. Please supply a present category.")
        }
    }
    if(verbose)
        message("Extracting results for ", contrast[1]," variable, ",
            contrast[2], " / ", contrast[3])
    sobj <- Seurat::NormalizeData(object = sobj, normalization.method =
        "LogNormalize", scale.factor = 10000, verbose = verbose)
    sobj <- Seurat::FindVariableFeatures(object = sobj, nfeatures = round(nrow(
        counts)*0.1, digits = 0), verbose = verbose) # loess.span = 1
    sobj <- Seurat::ScaleData(object = sobj, vars.to.regress = c("nCount_RNA"),
        verbose = verbose)
    if(verbose){ # FindMarkers: If we set verbose = FALSE we'll get an error
        statInfo_ <- Seurat::FindMarkers(sobj, test.use = test.use, group.by =
            contrast[1], ident.1 = contrast[2], ident.2 = contrast[3],
            logfc.threshold = 0, min.cpt = 0)
    } else {
        invisible(utils::capture.output(statInfo_ <- Seurat::FindMarkers(sobj,
            test.use = test.use, group.by = contrast[1], ident.1 = contrast[2],
            ident.2 = contrast[3], logfc.threshold = 0, min.cpt = 0)))
    }
    computed_features <- match(gsub(pattern = "_", x = rownames(counts),
        replacement = "-"),rownames(statInfo_))
    statInfo <- data.frame(matrix(NA, ncol = ncol(statInfo_), nrow = nrow(
        counts)))
    statInfo <- statInfo_[computed_features,]
    name <- paste(name, ".", test.use, sep = "")
    pValMat <- statInfo[, c("p_val", "p_val_adj")]
    colnames(pValMat) <- c("rawP","adjP")
    rownames(pValMat) <- rownames(statInfo) <- rownames(counts)
    list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name)
}# END - function: DA_Seurat

#' @title set_Seurat
#'
#' @export
#' @description
#' Set the parameters for Seurat differential abundance detection method.
#'
#' @inheritParams DA_Seurat
#' @param expand logical, if TRUE create all combinations of input parameters
#' (default \code{expand = TRUE})
#'
#' @return A named list containing the set of parameters for \code{DA_Seurat}
#' method.
#'
#' @seealso \code{\link{DA_Seurat}}
#'
#' @examples
#' # Set some basic combinations of parameters for Seurat
#' base_Seurat <- set_Seurat(contrast = c("group", "B", "A"))
#' # Set a specific set of normalization for Seurat (even of other packages!)
#' setNorm_Seurat <- set_Seurat(contrast = c("group", "B", "A"), norm = c("TSS",
#' "TMM", "poscounts"))
#' # Set many possible combinations of parameters for Seurat
#' all_Seurat <- set_Seurat(pseudo_count = c(TRUE, FALSE),
#'     test.use = c("wilcox", "t", "negbinom", "poisson"),
#'     contrast = c("group", "B", "A"), norm = c("TSS", "TMM"))

set_Seurat <- function(pseudo_count = FALSE, test.use = c("wilcox", "t"),
                       contrast = NULL, norm = "TSS", expand = TRUE) {
    method <- "DA_Seurat"
    if (!is.logical(pseudo_count)) {
        stop("'pseudo_count' and 'boot' must be logical.")
    }
    if (is.null(contrast)) {
        stop("'contrast' is required.")
    }
    if (sum(!is.element(norm, c("TSS", "none"))) > 0) {
        warning("One or more elements into 'norm' are not native to Seurat.")
    }
    if(sum(!is.element(test.use, c("wilcox", "t", "bimod", "roc", "negbinom",
        "poisson", "LR", "DESeq2"))) > 0){
        stop("One or more of the test.use are wrong.")
    }
    if (expand) {
        parameters <- expand.grid(method = method, pseudo_count = pseudo_count,
                                  test.use = test.use, norm = norm,
                                  stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, pseudo_count = pseudo_count,
                                 test.use = test.use, norm = norm)
    }
    # data.frame to list
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        x <- append(x = x, values = list("contrast" = contrast), after = 2)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}
