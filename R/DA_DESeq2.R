#' @title DA_DESeq2
#'
#' @importFrom DESeq2 sizeFactors DESeq dispersions results 
#' DESeqDataSetFromMatrix
#' @importFrom SummarizedExperiment assays
#' @export
#' @description
#' Fast run for DESeq2 differential abundance detection method.
#'
#' @inheritParams DA_edgeR
#' @param contrast character vector with exactly three elements: the name of a
#' factor in the design formula, the name of the numerator level for the fold
#' change, and the name of the denominator level for the fold change.
#' @param alpha the significance cutoff used for optimizing the independent
#' filtering (by default 0.05). If the adjusted p-value cutoff (FDR) will be a
#' value other than 0.05, alpha should be set to that value.
#' @param norm name of the normalization method to use in the differential 
#' abundance analysis. Choose between the native DESeq2 normalization methods, 
#' such as \code{ratio}, \code{poscounts}, or \code{iterate}. Alternatively 
#' (only for advanced users), if \code{norm} is equal to "TMM", "TMMwsp", 
#' "RLE", "upperquartile", "posupperquartile", or "none" from 
#' \code{\link{norm_edgeR}}, "CSS" from \code{\link{norm_CSS}}, or "TSS" from 
#' \code{\link{norm_TSS}}, the normalization factors are automatically 
#' transformed into size factors. If custom factors are supplied, make sure 
#' they are compatible with DESeq2 size factors.
#'
#' @return A list object containing the matrix of p-values `pValMat`,
#' the dispersion estimates `dispEsts`, the matrix of summary statistics for
#' each tag `statInfo`, and a suggested `name` of the final object considering
#' the parameters passed to the function.
#'
#' @seealso \code{\link[phyloseq]{phyloseq_to_deseq2}} for phyloseq to DESeq2
#' object conversion, \code{\link[DESeq2]{DESeq}} and
#' \code{\link[DESeq2]{results}} for the differential abundance method.
#'
#' @examples
#' set.seed(1)
#' # Create a very simple phyloseq object
#' counts <- matrix(rnbinom(n = 60, size = 3, prob = 0.5), nrow = 10, ncol = 6)
#' metadata <- data.frame("Sample" = c("S1", "S2", "S3", "S4", "S5", "S6"),
#'                        "group" = as.factor(c("A", "A", "A", "B", "B", "B")))
#' ps <- phyloseq::phyloseq(phyloseq::otu_table(counts, taxa_are_rows = TRUE),
#'                          phyloseq::sample_data(metadata))
#' # Calculate the poscounts size factors
#' ps_NF <- norm_DESeq2(object = ps, method = "poscounts")
#' # The phyloseq object now contains the size factors:
#' sizeFacts <- phyloseq::sample_data(ps_NF)[, "NF.poscounts"]
#' head(sizeFacts)
#' # Differential abundance
#' DA_DESeq2(object = ps_NF, pseudo_count = FALSE, design = ~ group, contrast =
#'     c("group", "B", "A"), norm = "poscounts")

DA_DESeq2 <- function(object, assay_name = "counts", pseudo_count = FALSE, 
    design = NULL, contrast = NULL, alpha = 0.05, norm = c("ratio", 
    "poscounts", "iterate"), weights, verbose = TRUE){
    counts_and_metadata <- get_counts_metadata(object, assay_name = assay_name)
    counts <- counts_and_metadata[[1]]
    libSizes <- colSums(counts)
    metadata <- counts_and_metadata[[2]]
    is_phyloseq <- counts_and_metadata[[3]]
    # Name building
    name <- "DESeq2"
    method <- "DA_DESeq2"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        message("Adding a pseudo count... \n")
        counts <- counts + 1
        name <- paste(name, ".pseudo", sep = "")
    }
    # Check the assay
    if (!is_phyloseq){
        if(verbose)
            message("Using the ", assay_name, " assay.")
        name <- paste(name, ".", assay_name, sep = "")
    } 
    if (is.character(design)){
        design <- as.formula(design)
    }
    if(verbose){
        dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, 
            colData = metadata, design = design)
    } else {
        dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(
            countData = counts, colData = metadata, design = design))
    }
    if(length(norm) > 1)
        stop(method, "\n", 
             "norm: please choose one normalization for this istance of", 
            " differential abundance analysis.")
    # Check if the column with the normalization factors is present
    NF.col <- paste("NF", norm, sep = ".")
    if (!any(colnames(metadata) == NF.col)) {
        stop(method, "\n", 
             "norm: can't find the ", NF.col," column in your object.",
             " Make sure to add the normalization factors column in your",
             " object first.")
    }
    NFs <- unlist(metadata[, NF.col])
    # Check for implemented normalizations
    this_norms <- c("ratio", "iterate", "poscounts")
    other_norms <- c("TMM", "TMMwsp", "RLE", "upperquartile", 
        "posupperquartile", "none", "TSS", "CSS")
    if(!is.element(norm, this_norms)) {
        if(verbose){
            warning(method, "\n", 
                norm, " normalization is not a native DESeq2",
                " normalization. Make sure you know what you are doing,", 
                " otherwise choose between 'ratio', 'poscounts', or 'iterate'.")
        } 
        if(is.element(norm, other_norms)){
            if(verbose)
                message("Automatically converting NF.", norm, " factors to", 
                " size factors.")
            if(norm == "TSS"){
                NFs <- 1
                if(verbose)
                    message("'TSS' normalization: automatically setting NFs",
                        " to ones.")
            }
            if(!is.element(norm, c("CSS", "TSS"))){ 
                # CSS scaling factors are ready for DESeq2
                NFs <- (NFs*libSizes)/exp(mean(log(NFs*libSizes)))
            } 
        } else {
            if(verbose){
                warning(method, "\n", 
                    "Make sure that the provided NF.", norm, " factors", 
                    " are compatible with DESeq2 size factors.")
            } 
        }
    }
    suppressMessages(expr = {DESeq2::sizeFactors(dds) <- NFs})
    name <- paste(name, ".", norm, sep = "")
    if(verbose)
        message("Differential abundance on ", norm, " normalized data")
    if(missing(weights)){
        if(verbose)
            message("Estimating Differential Abundance without weighting.",
                " No weights supplied.")
    } else {
        if(is.null(weights)){
            if(verbose)
                message("Estimating Differential Abundance without weighting")
        } else {
            if(verbose)
                message("Estimating Differential Abundance with weights")
            weights[which(weights < 1e-6)] <- 1e-06
            SummarizedExperiment::assays(dds,
                withDimnames = FALSE)[["weights"]] <- weights
            name <- paste(name, ".weighted", sep = "")
        }
    }
    ### Run DESeq
    if(verbose){
        ddsRes <- DESeq2::DESeq(object = dds, test = "LRT", reduced = ~ 1,
            parallel = FALSE, quiet = verbose)
        dispEsts <- DESeq2::dispersions(ddsRes)
    } else {
        suppressMessages(expr = {
            ddsRes <- DESeq2::DESeq(object = dds, test = "LRT", reduced = ~ 1,
                parallel = FALSE)
            dispEsts <- DESeq2::dispersions(ddsRes)})
    }
    if(missing(contrast) | (!is.character(contrast) & length(contrast) != 3))
        stop(method, "\n", 
            "contrast: please supply a character vector with exactly three",
            " elements: the name of a factor in the design formula, the name",
            " of the numerator level for the fold change, and the name of the",
            " denominator level for the fold change.")
    else {
        if(verbose)
            message("Extracting results for ", contrast[1], " variable, ",
                contrast[2], " / ", contrast[3])
    }
    if(verbose){
        res <- DESeq2::results(ddsRes, alpha = alpha, contrast = contrast)
    } else {
        res <- suppressMessages(
            DESeq2::results(ddsRes, alpha = alpha, contrast = contrast))
    }
    statInfo <- as(res, "data.frame")
    pValMat <- statInfo[,c("pvalue", "padj")]
    colnames(pValMat) <- c("rawP", "adjP")
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "dispEsts" =
        dispEsts, "name" = name))
}# END - function: DA_DESeq2

#' @title set_DESeq2
#'
#' @export
#' @description
#' Set the parameters for DESeq2 differential abundance detection method.
#'
#' @inheritParams DA_DESeq2
#' @param weights_logical logical vector, if TRUE a matrix of observational
#' weights will be used for differential abundance analysis (default
#' \code{weights_logical = FALSE}).
#' @param expand logical, if TRUE create all combinations of input parameters
#' (default \code{expand = TRUE}).
#'
#' @return A named list containing the set of parameters for \code{DA_DESeq2}
#' method.
#'
#' @seealso \code{\link{DA_DESeq2}}
#'
#' @examples
#' # Set some basic combinations of parameters for DESeq2
#' base_DESeq2 <- set_DESeq2(design = ~ group, contrast = c("group", "B", "A"))
#' # Set a specific set of normalization for DESeq2
#' setNorm_DESeq2 <- set_DESeq2(design = ~ group, contrast =
#'     c("group", "B", "A"), norm = c("ratio", "poscounts"))
#' # Set many possible combinations of parameters for DESeq2
#' all_DESeq2 <- set_DESeq2(pseudo_count = c(TRUE, FALSE), design = ~ group,
#'     contrast = c("group", "B", "A"), weights_logical = c(TRUE,FALSE))
set_DESeq2 <- function(assay_name = "counts", pseudo_count = FALSE, 
    design = NULL, contrast = NULL, alpha = 0.05, norm = c("ratio", 
    "poscounts", "iterate"), weights_logical = FALSE, expand = TRUE) {
    method <- "DA_DESeq2"
    if (is.null(assay_name)) {
        stop(method, "\n", "'assay_name' is required (default = 'counts').")
    }
    if (!is.logical(pseudo_count) | !is.logical(weights_logical)) {
        stop(method, "\n", 
            "'pseudo_count' and 'weights_logical' must be logical.")
    }
    if (is.null(design) | is.null(contrast)) {
        stop(method, "\n", "'design' and 'contrast' must be specified.")
    }
    if (!is.character(design) & !is(design, "formula")){
        stop(method, "\n", "'design' should be a character or a formula.")
    } else design <- as.formula(design)
    if (!is.character(contrast) & length(contrast) != 3){
        stop(method, "\n", 
            "contrast: please supply a character vector with exactly three", 
            " elements: the name of a factor in the design formula, the name",
            " of the numerator level for the fold change, and the name of the",
            " denominator level for the fold change.")
    }
    if (sum(!is.element(norm, c("ratio", "poscounts", "iterate"))) > 0) {
        warning(method, "\n", 
            "One or more elements into 'norm' are not native to DESeq2")
    }
    if (expand) {
        parameters <- expand.grid(method = method, assay_name = assay_name, 
            pseudo_count = pseudo_count, alpha = alpha, norm = norm, 
            weights = weights_logical, stringsAsFactors = FALSE)
    } else {
        message("Some parameters may be duplicated to fill the matrix.")
        parameters <- data.frame(method = method, assay_name = assay_name, 
            pseudo_count = pseudo_count, alpha = alpha, norm = norm, 
            weights = weights_logical)
    }
    # data.frame to list
    out <- plyr::dlply(.data = parameters, .variables = colnames(parameters))
    out <- lapply(X = out, FUN = function(x){
        x <- append(x = x, values = list("design" = deparse(design),
                                         "contrast" = contrast), after = 3)
    })
    names(out) <- paste0(method, ".", seq_along(out))
    return(out)
}
