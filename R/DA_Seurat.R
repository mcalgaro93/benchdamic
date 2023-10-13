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
#' @inheritParams Seurat::NormalizeData
#' @param norm Method for normalization. 
#' \itemize{\item{\code{LogNormalize}}{ Feature counts for each sample are 
#' divided by the total counts of that sample and multiplied by the 
#' scale.factor. This is then natural-log transformed using log1p;}
#' \item{\code{CLR}}{ Applies a centered log ratio transformation;}
#' \item{\code{RC}}{ Relative counts. Feature counts for each sample are 
#' divided by the total counts of that sample and multiplied by the 
#' scale.factor. No log-transformation is applied. For counts per million 
#' (CPM) set scale.factor = 1e6;}
#' \item{\code{none}}{ No normalization}}
#' @param test Denotes which test to use. Available options are:
#' \itemize{\item{\code{"wilcox"}}{ Identifies differentially abundant 
#' features between two groups of samples using a Wilcoxon Rank Sum test 
#' (default).}
#' \item{\code{"bimod"}}{ Likelihood-ratio test for the feature abundances, 
#' (McDavid et al., Bioinformatics, 2013).}
#' \item{\code{"roc"}}{ Identifies 'markers' of feature abundance using ROC 
#' analysis. For each feature, evaluates (using AUC) a classifier built on that 
#' feature alone, to classify between two groups of cells. An AUC value of 1 
#' means that abundance values for this feature alone can perfectly classify the
#' two groupings (i.e. Each of the samples in group.1 exhibit a higher level 
#' than each of the samples in group.2). An AUC value of 0 also means there is 
#' perfect classification, but in the other direction. A value of 0.5 implies 
#' that the feature has no predictive power to classify the two groups. Returns 
#' a 'predictive power' (abs(AUC-0.5) * 2) ranked matrix of putative 
#' differentially expressed genes.}
#' \item{\code{"t"}}{ Identify differentially abundant features between two 
#' groups of samples using the Student's t-test.}
#' \item{\code{"negbinom"}}{ Identifies differentially abundant features 
#' between two groups of samples using a negative binomial generalized linear 
#' model.}
#' \item{\code{"poisson"}}{ Identifies differentially abundant features between 
#' two groups of samples using a poisson generalized linear model.}
#' \item{\code{"LR"}}{ Uses a logistic regression framework to determine 
#' differentially abundant features. Constructs a logistic regression model 
#' predicting group membership based on each feature individually and compares 
#' this to a null model with a likelihood ratio test.}
#' \item{\code{"MAST"}}{ Identifies differentially expressed genes between two 
#' groups of cells using a hurdle model tailored to scRNA-seq data. Utilizes 
#' the MAST package to run the DE testing.}
#' \item{\code{"DESeq2"}}{ Identifies differentially abundant features between 
#' two groups of samples based on a model using DESeq2 which uses a negative 
#' binomial distribution (Love et al, Genome Biology, 2014).}}
#'
#' @return A list object containing the matrix of p-values `pValMat`, the 
#' matrix of summary statistics for each tag `statInfo`, and a suggested `name` 
#' of the final object considering the parameters passed to the function.
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
#'                          
#' # Differential abundance
#' DA_Seurat(object = ps, contrast = c("group","B","A"))
#' 
#' # Perform a simple Wilcoxon test using Seurat on raw data
#' DA_Seurat(object = ps, contrast = c("group","B","A"), norm = "none", 
#'     test = "wilcox")

DA_Seurat <- function(object, assay_name = "counts", pseudo_count = FALSE, 
    norm = "LogNormalize", scale.factor = 10000, test = "wilcox", 
    contrast = NULL, verbose = TRUE){
    counts_and_metadata <- get_counts_metadata(object, assay_name = assay_name)
    counts <- counts_and_metadata[[1]]
    metadata <- counts_and_metadata[[2]]
    is_phyloseq <- counts_and_metadata[[3]]
    # Name building
    name <- "Seurat"
    method <- "DA_Seurat"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        if(verbose)
            message("Adding a pseudo count...")
        counts <- counts + 1
        name <- paste(name,".pseudo",sep = "")}
    # Check the assay
    if (!is_phyloseq){
        if(verbose)
            message("Using the ", assay_name, " assay.")
        name <- paste(name, ".", assay_name, sep = "")
    } 
    # Initialize the Seurat object with the raw (non-normalized data).
    # Keep all features expressed in >= 1 sample
    # Keep all samples with at least 1 detected feature.
    if(verbose){
        sobj <- Seurat::CreateSeuratObject(counts = counts, min.cells = 1, 
            min.features = 1)
    } else {
        sobj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts, 
            min.cells = 1, min.features = 1))
    }  
    if(!is.character(contrast) | length(contrast) != 3)
        stop(method, "\n", 
             "contrast: please supply a character vector with exactly", 
             " three elements: the name of a variable used in",  
             " 'design', the name of the level of interest, and the", 
             " name of the reference level.")
    if(is.element(contrast[1], colnames(metadata))){
        if(!is.factor(metadata[, contrast[1]])){
            if(verbose){
                message("Converting variable ", contrast[1], " to factor.")
            }
            metadata[, contrast[1]] <- as.factor(metadata[, contrast[1]])
        }
        if(!is.element(contrast[2], levels(metadata[, contrast[1]])) | 
           !is.element(contrast[3], levels(metadata[, contrast[1]]))){
            stop(method, "\n", 
                 "contrast: ", contrast[2], " and/or ", contrast[3], 
                 " are not levels of ", contrast[1], " variable.")
        }
        if(verbose){
            message("Setting ", contrast[3], " the reference level for ", 
                    contrast[1], " variable.")
        }
        metadata[, contrast[1]] <- stats::relevel(metadata[, contrast[1]], 
                                                  ref = contrast[3])
    }
    sobj <- Seurat::AddMetaData(object = sobj, metadata = metadata,
        col.name = colnames(metadata))
    if(verbose)
        message("Extracting results for ", contrast[1]," variable, ",
            contrast[2], " / ", contrast[3])
    # Check test
    if(!is.element(test, c("wilcox", "bimod", "roc", "t", "negbinom", 
        "poisson", "LR", "MAST", "DESeq2"))){
        stop(method, "\n", 
            "test: ", test, " is not one of the allowed tests.")
    }
    slot_name <- "data" # Setting the default slot of Seurat object
    # If these tests are used, normalization
    if(is.element(test, c("wilcox", "bimod", "roc", "t", "LR", "MAST"))){
        # Check norm and scale.factor
        if(length(norm) > 1 | length(scale.factor) > 1)
            stop(method, "\n", 
                "Please choose one 'norm' and one 'scale.factor' value for",
                " this istance of differential abundance analysis.")
        if(!is.element(norm, c("LogNormalize", "CLR", "RC", 
            "none"))){
            stop(method, "\n", "norm: ", norm, " normalization is not an", 
                " available option. Please choose between 'LogNormalize',",
                " 'CLR', 'RC', or 'none'.")
        }
        if(norm != "none"){
            sobj <- Seurat::NormalizeData(object = sobj, normalization.method =
                norm, scale.factor = scale.factor, margin = 2,
                verbose = verbose)
            name <- paste(name, ".", norm, sep = "")
            name <- paste(name, ".SF", scale.factor, sep = "")
        } else slot_name <- "counts" # Change the default slot into "counts"
    } else {
        if(verbose & (!is.null(norm) | !is.null(scale.factor))){
            warning(method, "\n", 
                "test = ", test, ": 'norm' and 'scale.factor' won't",
                " be used (raw counts are used instead).")
        }
    }
    if(verbose){ # FindMarkers: If we set verbose = FALSE we'll get an error
        statInfo_ <- Seurat::FindMarkers(sobj, slot = slot_name, 
            test.use = test, group.by = contrast[1], ident.1 = contrast[2], 
            ident.2 = contrast[3], logfc.threshold = 0, min.pct = 0)
    } else {
        invisible(utils::capture.output(statInfo_ <- Seurat::FindMarkers(sobj,
            slot = slot_name, test.use = test, group.by = contrast[1], 
            ident.1 = contrast[2], ident.2 = contrast[3], logfc.threshold = 0, 
            min.pct = 0)))
    }
    computed_features <- match(gsub(pattern = "_", x = rownames(counts),
        replacement = "-"), rownames(statInfo_))
    statInfo <- data.frame(matrix(NA, ncol = ncol(statInfo_), nrow = nrow(
        counts)))
    statInfo <- statInfo_[computed_features,]
    name <- paste(name, ".", test, sep = "")
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
#' # Set many possible combinations of parameters for Seurat
#' all_Seurat <- set_Seurat(test = c("wilcox", "t", "negbinom", "poisson"),
#'     norm = c("LogNormalize", "CLR", "RC", "none"), 
#'     scale.factor = c(1000, 10000), contrast = c("group", "B", "A"))

set_Seurat <- function(assay_name = "counts", pseudo_count = FALSE, 
    test = "wilcox", contrast = NULL, norm = "LogNormalize", 
    scale.factor = 10000, expand = TRUE) {
    method <- "DA_Seurat"
    if (is.null(assay_name)) {
        stop(method, "\n", "'assay_name' is required (default = 'counts').")
    }
    if (!is.logical(pseudo_count)) {
        stop(method, "\n", "'pseudo_count' and 'boot' must be logical.")
    }
    if (is.null(contrast)) {
        stop(method, "\n", "'contrast' is required.")
    }
    if (is.null(norm)) {
        stop(method, "\n", "'norm' is required.")
    }
    if (is.null(scale.factor)) {
        stop(method, "\n", "'scale.factor' is required.")
    }
    if(sum(!is.element(test, c("wilcox", "t", "bimod", "roc", "negbinom",
        "poisson", "LR", "DESeq2"))) > 0){
        stop(method, "\n", "One or more elements of the 'test' are wrong.")
    }
    if(sum(!is.element(norm, c("LogNormalize", "CLR", "RC", "none"))) > 0){
        stop(method, "\n", 
            "One or more elements of the 'norm' are wrong.")
    }
    if (expand) {
        parameters <- expand.grid(method = method, assay_name = assay_name,
            pseudo_count = pseudo_count, 
            norm = norm, 
            scale.factor = scale.factor, test = test, 
            stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, assay_name = assay_name,
            pseudo_count = pseudo_count,
            norm = norm, 
            scale.factor = scale.factor, test = test)
    }
    # Remove cases of normalization with tests that don't need normalization
    count_test <- is.element(parameters[, "test"], c("negbinom", 
        "poisson", "DESeq2"))
    if(sum(count_test) > 0){
        message("Replacing 'norm' and 'scale.factor' with",
            " default values where 'test' is equal to 'negbinom',", 
            " 'poisson', or 'DESeq2'. In those cases, normalization and", 
            " scaling are not performed.")
        parameters[count_test, "norm"] <- "LogNormalize"
        parameters[count_test, "scale.factor"] <- 10000
    }
    # Remove duplicates
    dup <- duplicated(parameters)
    if(sum(dup) > 0){
        message("Removing duplicated set of parameters.")
        parameters <- parameters[-which(dup), ]
    }
    # data.frame to list
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        x <- append(x = x, values = list("contrast" = contrast), after = 3)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}
