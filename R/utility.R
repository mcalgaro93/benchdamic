#' @title createMocks
#'
#' @export
#' @description
#' Given the number of samples of the dataset from which the mocks should be
#' created this function produces a `data.frame` object with as many rows as the
#' number of mocks and as many columns as the number of samples. If an odd
#' number of samples is given, the lower even integer will be considered in
#' order to obtain a balanced design for the mocks.
#' @param nsamples an integer representing the total number of samples.
#' @param N number of mock comparison to generate.
#' @param seed seed to use for random mock generation.
#'
#' @return A `data.frame` containing `N` rows and `nsamples` columns (if even).
#' Each cell of the data frame contains the "grp1" or "grp2" characters which
#' represent the mock groups pattern.
#' @examples
#' # Let's generate the patterns for 100 mock comparison for an experiment with
#' # 15 grp1 and 15 grp2 samples.
#' mocks <- createMocks(nsamples = 30, N = 100)
#' head(mocks)

# Create data.frame with random labels
createMocks <- function(nsamples, N = 1000, seed = 123){
    sample_num <- nsamples %/% 2 * 2 # Balanced design sample numerosity
    mock_df <- matrix(NA, nrow = N, ncol = sample_num)
    set.seed(seed)
    for(i in 1:N){ # N random balanced relabellings
        grps <- rep("grp1", sample_num)
        grps[sample(1:sample_num,size = sample_num/2)] <- "grp2"
        mock_df[i,] <- grps
    }
    rownames(mock_df) <- paste0("Comparison",1:N)
    return(mock_df)
}# END - function: createMocks

### Normalizations ###

#' @title norm_edgeR
#'
#' @importFrom edgeR calcNormFactors
#' @importFrom phyloseq sample_sums taxa_are_rows otu_table
#' @importFrom stats quantile
#' @export
#' @description
#' Calculate normalization factors from a phyloseq object to scale the raw
#' library sizes. Inherited from edgeR `calcNormFactors` function.
#'
#' @param object phyloseq object containing the counts to be normalized.
#' @param method normalization method to be used. Choose between "TMM",
#' "TMMwsp", "RLE", "upperquartile", "posupperquartile" or "none".
#' @inheritParams edgeR::calcNormFactors
#'
#' @return A new column containing the chosen edgeR-based normalization factors
#' is added to the phyloseq `sample_data` slot. The effective library sizes for
#' use in downstream analysis must be multiplied by the normalization factors.
#' @seealso [edgeR::calcNormFactors()] for details.

norm_edgeR <- function(object,
                       method=c("TMM", "TMMwsp",
                                "RLE", "upperquartile",
                                "posupperquartile", "none"),
                       refColumn = NULL, logratioTrim = .3, sumTrim = 0.05,
                       doWeighting = TRUE, Acutoff = -1e10, p = 0.75, ...)
{
    counts <- as(otu_table(object), "matrix")
    if (!phyloseq::taxa_are_rows(object))
    {
        counts <- t(counts)
    } else {}

    if (is.na(method)){
        stop("Please, supply a valid method between 'TMM', 'TMMwsp', 'RLE',
        'upperquartile', 'posupperquartile' or 'none'")
    }
    else if (method == "posupperquartile")
    {
        scaledCounts <- t(counts) / colSums(counts)
        tmpNF <- apply(scaledCounts, MARGIN = 1L, FUN = function(x)
            quantile(x[x != 0], probs = .75))
        normFacts <- tmpNF/exp(mean(log(tmpNF)))
    } else {
        normFacts <- edgeR::calcNormFactors(counts, method=method,
                                            refColumn = refColumn,
                                            logratioTrim = logratioTrim,
                                            sumTrim = sumTrim,
                                            doWeighting= doWeighting,
                                            Acutoff = Acutoff,
                                            p = p, ...)
    }# END - ifelse: posupperquartile = upperquartile only of non-zero counts

    # VERY IMPORTANT: multiply by library sizes and renormalize.
    # edgeR calculates scaling factors, which still have to be multiplied by
    # library sizes to get to the size factors of effective sequencing depth,
    # i.e. robust estimates of the library sizes.
    # normFacts = normFacts*sample_sums(object)
    # normFacts = normFacts/exp(mean(log(normFacts)))

    if (all(is.na(normFacts))) #Resort to proportion normalization in case of
        # failure for all samples
    {
        normFacts = phyloseq::sample_sums(object)
    }
    object@sam_data@.Data <- c(object@sam_data@.Data, list(normFacts))
    aux <- object@sam_data@names
    aux[length(aux)] <- paste("NF", method, sep = ".")
    object@sam_data@names <- aux
    return(object)
}# END - function: norm_edgeR

#' @title norm_DESeq2
#'
#' @importFrom DESeq2 estimateSizeFactors sizeFactors
#' @importFrom phyloseq sample_sums taxa_are_rows otu_table
#' @importFrom stats quantile median
#' @export
#' @description
#' Calculate normalization factors from a phyloseq object to scale the raw
#' library sizes. Inherited from DESeq2 `estimateSizeFactors` function.
#'
#' @param object phyloseq object containing the counts to be normalized.
#' @param method Method for estimation: either "ratio", "poscounts", or
#' "iterate". "ratio" uses the standard median ratio method introduced in DESeq.
#' The size factor is the median ratio of the sample over a "pseudosample": for
#' each gene, the geometric mean of all samples. "poscounts" and "iterate" offer
#' alternative estimators, which can be used even when all features contain a
#' sample with a zero (a problem for the default method, as the geometric mean
#' becomes zero, and the ratio undefined). The "poscounts" estimator deals with
#' a feature with some zeros, by calculating a modified geometric mean by taking
#' the n-th root of the product of the non-zero counts. This evolved out of use
#' cases with Paul McMurdie's phyloseq package for metagenomic samples. The
#' "iterate" estimator iterates between estimating the dispersion with a design
#' of ~1, and finding a size factor vector by numerically optimizing the
#' likelihood of the ~1 model.
#' @param ... other parameters for DESeq2 [DESeq2::estimateSizeFactors()]
#' function.
#'
#' @return A new column containing the chosen DESeq2-based normalization factors
#' is added to the phyloseq `sample_data` slot.
#'
#' @seealso [DESeq2::estimateSizeFactors()] for details.

norm_DESeq2 <- function(object,
                        method = c("ratio", "poscounts", "iterate"),
                        ...)
{
    if (!phyloseq::taxa_are_rows(object))
    {
        object <- t(object)
    } else {}

    ## Calculate size factors
    obj <- phyloseq::phyloseq_to_deseq2(object, design = ~ 1)

    if(missing(method))
        stop("Please supply a normalization method between 'ratio', 'poscounts'
             or 'iterate'.")

    normFacts <- DESeq2::sizeFactors(DESeq2::estimateSizeFactors(obj,
                                                                 type = method,
                                                                 ...))

    object@sam_data@.Data <- c(object@sam_data@.Data, list(normFacts))
    aux <- object@sam_data@names
    aux[length(aux)] <- paste("NF", method, sep = ".")
    object@sam_data@names <- aux
    return(object)
}# END - function: norm_DESeq2

#' @title norm_CSS
#'
#' @importFrom metagenomeSeq newMRexperiment
#' @importFrom phyloseq taxa_are_rows otu_table
#' @importFrom stats median
#' @export
#' @description
#' Calculate normalization factors from a phyloseq object to scale the raw
#' library sizes. Inherited from metagenomeSeq `calcNormFactors` function
#' which performs the Cumulative Sum Scaling normalization.
#'
#' @param object phyloseq object containing the counts to be normalized.
#' @param method  Normalization scaling parameter (default = "1000"). If
#' "median", the median of the normalization factors is used as scaling (Paulson
#' et al. 2013).
#'
#' @return A new column containing the CSS normalization factors is added to the
#' phyloseq `sample_data` slot.
#'
#' @seealso [metagenomeSeq::calcNormFactors()] for details.

norm_CSS <- function(object, method = "default")
{
    if (!phyloseq::taxa_are_rows(object))
    {
        object <- t(object)
    } else {}

    counts <- as(phyloseq::otu_table(object), "matrix")

    obj <- metagenomeSeq::newMRexperiment(counts = counts)
    normFacts <- metagenomeSeq::calcNormFactors(obj = obj)
    normFacts <- drop(as.matrix(normFacts))

    # Default: log2(normFacts/1000 + 1)
    # Original metagenomeSeq paper: log2(normFacts/median(libsize) +1)
    if(method == "default"){
        normFacts <- log2(normFacts/1000 + 1)
    } else if (method == "median"){
        normFacts <- log2(normFacts/stats::median(normFacts) + 1)
    } else stop("Please choose a scaling method between 'default' or 'median'.")
    # Remember to useCSSoffset = FALSE in fitZig function

    object@sam_data@.Data <- c(object@sam_data@.Data, list(normFacts))
    aux <- object@sam_data@names
    aux[length(aux)] <- paste("NF.CSS", method, sep = ".")
    object@sam_data@names <- aux
    return(object)
}# END - function: norm_CSS

#' @title norm_TSS
#'
#' @importFrom phyloseq taxa_are_rows sample_sums
#' @export
#' @description
#' Calculate the raw library sizes from a phyloseq object. If used to divide
#' counts, known as Total Sum Scaling normalization (TSS).
#'
#' @param object phyloseq object containing the counts to be normalized.
#'
#' @return A new column containing the TSS normalization factors is added to the
#' phyloseq `sample_data` slot.

norm_TSS <- function(object)
{
    if (!phyloseq::taxa_are_rows(object))
    {
        object <- t(object)
    } else {}

    normFacts <- 1/phyloseq::sample_sums(object)
    object@sam_data@.Data <- c(object@sam_data@.Data, list(normFacts))
    aux <- object@sam_data@names
    aux[length(aux)] <- "NF.TSS"
    object@sam_data@names <- aux
    return(object)
}# END - function: norm_TSS

### Differential Abundance Methods  ###

#' @title DA_edgeR
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data
#' @importFrom edgeR DGEList estimateDisp estimateGLMRobustDisp glmQLFit
#' glmQLFTest
#' @importFrom stats model.matrix p.adjust
#' @export
#' @description
#' Fast run for edgeR differential abundance detection method.
#'
#' @param object phyloseq object.
#' @param pseudo_count Add 1 to all counts if TRUE (default = FALSE).
#' @inheritParams edgeR::DGEList
#' @inheritParams edgeR::estimateDisp
#' @inheritParams edgeR::glmQLFit
#' @inheritParams edgeR::glmQLFTest
#' @inheritParams edgeR::estimateGLMRobustDisp
#' @param norm name of the normalization method used to compute the
#' normalization factors to use in the differential abundace analysis.
#'
#' @return A list object containing the matrix of p-values, the dispersion
#' estimates, the matrix of summary statistics for each tag, and a suggested
#' name of the final object considering the parameters passed to the function.
#'
#' @seealso [edgeR::DGEList()] for the edgeR DEG object creation,
#' [edgeR::estimateDisp()] and [edgeR::estimateGLMRobustDisp()] for dispersion
#' estimation, and [edgeR::glmQLFit()] and [edgeR::glmQLFTest()] for the
#' quasi-likelihood negative binomial model fit.

DA_edgeR <- function(object,
                     pseudo_count = FALSE,
                     group = NULL,
                     design = NULL,
                     contrast = NULL,
                     robust = FALSE,
                     coef = 2,
                     norm = c("TMM", "TMMwsp", "RLE", "upperquartile", # edgeR
                              "posupperquartile", "none", # edgeR
                              "ratio", "poscounts", "iterate", # DESeq2
                              "TSS", "CSS"), # others
                     weights){

    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
    {
        object <- t(object)
    } else {}

    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)

    # Name building
    name <- "edgeR"

    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count)
    {
        message("Adding a pseudo count... \n")
        counts <- counts + 1
        name <- paste(name,".pseudo",sep = "")
    } else {}

    if(length(norm) > 1)
        stop("Please choose one normalization for this istance of differential
             abundance analysis.")

    if(norm == "TSS"){
        NFs = 1
    } else {
        # Check if the column with the normalization factors is present
        NF.col <- paste("NF", norm, sep = ".")
        if(!any(colnames(metadata) == NF.col))
            stop(paste0("Can't find the ", NF.col," column in your object. Be
            sure to add the normalization factors column in your object
                        first."))

        NFs = unlist(metadata[,NF.col])
        # DESeq2 NFs are supplied -> make them scaling factors!
        if(is.element(norm, c("ratio", "poscounts", "iterate"))){
            NFs <- NFs*colSums(counts)
        }
        NFs = NFs/exp(mean(log(NFs))) # Make NFs multiply to 1
    }

    name <- paste(name, ".", norm, sep = "")
    message(paste0("Differential abundance on ", norm," normalized data"))

    dge <- edgeR::DGEList(counts = counts,
                          norm.factors = NFs,
                          group = group)

    if(missing(weights)){
        message("Estimating Differential Abundance without weighting")
    } else {
        message("Estimating Differential Abundance with weights")
        dge$weights <- weights
        name <- paste(name,".weighted",sep = "")
    }

    if(class(design) == "formula")
        design <- stats::model.matrix(object = design,
                                      data = data.frame(metadata))

    if(!robust){
        dge <- edgeR::estimateDisp(y = dge,
                                   design = design)
    } else {
        message(paste0("Estimating robust dispersions"))
        dge <- edgeR::estimateGLMRobustDisp(y = dge,
                                            design = design)
        name <- paste(name,".robust",sep = "")
    }

    message(paste0("Extracting results"))
    glmFit <- edgeR::glmQLFit(y = dge,
                              dispersion = dge$tagwise.dispersion,
                              robust = robust,
                              design = design)
    glmRes <- edgeR::glmQLFTest(glmFit,
                                coef = coef,
                                contrast = contrast)
    if(missing(contrast)){
        message(paste0("Extracting results for ",
                       colnames(glmRes$coefficients)[coef],
                       " coefficient"))
    } else {
        message(paste0("Extracting results for ",
                       contrast,
                       " contrasts"))
    }


    pval <- glmRes$table$PValue
    padj <- stats::p.adjust(pval, "BH")
    pValMat <- cbind("rawP" = pval, "adjP" = padj)
    rownames(pValMat) = rownames(glmRes$table)

    statInfo <- glmRes$table

    return(list("pValMat" = pValMat,
                "statInfo" = statInfo,
                "dispEsts" = dge$tagwise.dispersion,
                "name" = name))
}# END - function: DA_edgeR

#' @title DA_DESeq2
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data phyloseq_to_deseq2
#' @importFrom DESeq2 sizeFactors DESeq dispersions results
#' @importFrom SummarizedExperiment assays
#' @export
#' @description
#' Fast run for DESeq2 differential abundance detection method.
#'
#' @param object phyloseq object.
#' @param pseudo_count Add 1 to all counts if TRUE (default = FALSE).
#' @inheritParams phyloseq::phyloseq_to_deseq2
#' @param contrast character vector with exactly three elements: the name of a
#' factor in the design formula, the name of the numerator level for the fold
#' change, and the name of the denominator level for the fold change.
#' @inheritParams DESeq2::results
#' @param norm name of the normalization method used to compute the
#' normalization factors to use in the differential abundace analysis.
#' @param weights optional numeric matrix giving observation weights.
#'
#' @return A list object containing the matrix of p-values `pValMat`,
#' the dispersion estimates `dispEsts`, the matrix of summary statistics for
#' each tag `statInfo`, and a suggested `name` of the final object considering
#' the parameters passed to the function.
#'
#' @seealso [phyloseq::phyloseq_to_deseq2()] for phyloseq to DESeq2 object
#' conversion, [DESeq2::DESeq()] and [DESeq2::results()] for the differential
#' abundance method.

DA_DESeq2 <- function(object,
                      pseudo_count = FALSE,
                      design = NULL,
                      contrast = NULL,
                      alpha = 0.05,
                      norm = c("TMM", "TMMwsp", "RLE", "upperquartile", # edgeR
                               "posupperquartile", "none", # edgeR
                               "ratio", "poscounts", "iterate", # DESeq2
                               "TSS", "CSS"), # others
                      weights){

    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
    {
        object <- t(object)
    } else {}

    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)

    # Name building
    name <- "DESeq2"

    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count)
    {
        message("Adding a pseudo count... \n")
        counts <- counts + 1
        phyloseq::otu_table(object) <- counts
        name <- paste(name,".pseudo",sep = "")
    } else {}

    dds <- phyloseq::phyloseq_to_deseq2(object, design = design)

    if(length(norm) > 1)
        stop("Please choose one normalization for this istance of differential
             abundance analysis.")

    # Check if the column with the normalization factors is present
    NF.col <- paste("NF", norm, sep = ".")
    if(!any(colnames(metadata) == NF.col))
        stop(paste0("Can't find the ", NF.col," column in your object. Be sure
        to add the normalization factors column in your object first."))

    NFs = unlist(metadata[,NF.col])
    # edgeR, TSS, and CSS NFs are supplied -> make them normalization factors!
    if(is.element(norm, c("TMM", "TMMwsp", "RLE", "upperquartile",
                          "posupperquartile", "CSS", "TSS"))){
        NFs <- NFs*colSums(counts)
    }
    DESeq2::sizeFactors(dds) = NFs/exp(mean(log(NFs))) # Make NFs multiply to 1

    name <- paste(name, ".", norm, sep = "")
    message(paste0("Differential abundance on ", norm," normalized data"))

    if(missing(weights)){
        message("Estimating Differential Abundance without weighting")
    } else {
        message("Estimating Differential Abundance with weights")
        weights[which(weights < 1e-6)] <- 1e-06
        SummarizedExperiment::assays(dds)[["weights"]] <- weights
        name <- paste(name,".weighted",sep = "")
    }

    ### Run DESeq
    ddsRes <- DESeq2::DESeq(object = dds,
                            test = "LRT",
                            reduced = ~ 1,
                            parallel = FALSE)
    dispEsts <- DESeq2::dispersions(ddsRes)

    if(missing(contrast) | (!is.character(contrast) & length(contrast) != 3)){
        stop(paste0("Please supply a character vector with exactly three
                    elements: the name of a factor in the design formula,
                    the name of the numerator level for the fold change,
                    and the name of the denominator level for the fold
                    change."))
    } else {
        message(paste0("Extracting results for ",
                       contrast[1]," variable, ",
                       contrast[2], " vs ",
                       contrast[3]))
    }

    Res <- DESeq2::results(ddsRes,
                           alpha = alpha,
                           contrast = contrast)

    pValMat <- as.matrix(Res[, c("pvalue", "padj")])
    colnames(pValMat) <- c("rawP", "adjP")
    statInfo <- data.frame(Res@listData)

    return(list("pValMat" = pValMat,
                "statInfo" = statInfo,
                "dispEsts" = dispEsts,
                "name" = name))
}# END - function: DA_DESeq2

