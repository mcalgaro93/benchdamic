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
    # Renormalize: multiply to 1
    # normFacts = normFacts/exp(mean(log(normFacts)))

    if (all(is.na(normFacts)))
    {
        stop("Failed to compute normalization factors!")
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

#' @title weights_ZINB
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data
#' @importFrom stats model.matrix
#' @importFrom zinbwave zinbFit computeObservationalWeights
#' @importFrom BiocParallel SerialParam
#' @export
#' @description
#' Computes the observational weights of the counts under a zero-inflated
#' negative binomial (ZINB) model. For each count, the ZINB distribution is
#' parametrized by three parameters: the mean value and the dispersion of the
#' negative binomial distribution, and the probability of the zero component.
#'
#' @param object phyloseq object containing the counts and the sample data.
#' @param design character name of the metadata columns, formula, or design
#' matrix with rows corresponding to samples and columns to coefficients to be
#' estimated (the user needs to explicitly include the intercept in the design).
#' @inheritParams zinbwave::zinbFit
#'
#' @return A matrix of weights.
#'
#' @seealso [zinbwave::zinbFit()] for zero-inflated negative binomial
#' parameters' estimation and [zinbwave::computeObservationalWeights()] for
#' weights extraction.

weights_ZINB <- function(object,
                         design,
                         K = 0,
                         commondispersion = TRUE,
                         zeroinflation = TRUE,
                         verbose = FALSE,
                         ...){
    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
    {
        object <- t(object)
    } else {}

    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)

    if(is.character(design)){
        if(grepl("~", design))
            design <- as.formula(design)
        else design <- as.formula(paste0("~", design))
    }

    if(class(design) == "formula")
        design <- stats::model.matrix(object = design,
                                      data = data.frame(metadata))

    zinbmodel <- zinbwave::zinbFit(Y = counts,
                                   X = design,
                                   K = K,
                                   commondispersion = commondispersion,
                                   zeroinflation = TRUE,
                                   verbose = verbose,
                                   BPPARAM = BiocParallel::SerialParam(),
                                   ... = ...)

    w <- zinbwave::computeObservationalWeights(model = zinbmodel,
                                               x = counts)
    return(w)
}# END - function: weights_ZINB

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
#' normalization factors to use in the differential abundance analysis.
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
                              "TSS", "CSSmedian", "CSSdefault"), # others
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
        # DESeq2 NFs are size factors. To obtain normalized counts
        # we need to divide the raw counts by the NFs
        # DESeq2 NFs are supplied -> make them scaling factors!
        if(is.element(norm, c("ratio", "poscounts", "iterate"))){
            NFs <- NFs/colSums(counts)
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

    if(is.character(design)){
        if(grepl("~", design))
            design <- as.formula(design)
        else design <- as.formula(paste0("~", design))
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
#' normalization factors to use in the differential abundance analysis.
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
                               "TSS", "CSSmedian", "CSSdefault"), # others
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
        stop(paste0("Can't find the ", NF.col," column in your object. Make sure
        to add the normalization factors column in your object first."))

    NFs = unlist(metadata[,NF.col])
    # edgeR, TSS, and CSS NFs supplied -> make them normalization factors!
    if(is.element(norm, c("TMM", "TMMwsp", "RLE", "upperquartile",
                          "posupperquartile", "CSSmedian", "CSSdefault",
                          "TSS"))){
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
        SummarizedExperiment::assays(dds,withDimnames = FALSE)$weights <- weights
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

#' @title DA_limma
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom stats model.matrix
#' @export
#' @description
#' Fast run for limma voom differential abundance detection method.
#'
#' @param object phyloseq object.
#' @param pseudo_count Add 1 to all counts if TRUE (default = FALSE).
#' @param norm name of the normalization method used to compute the
#' normalization factors to use in the differential abundance analysis.
#' @param design character name of the metadata columns, formula, or design
#' matrix with rows corresponding to samples and columns to coefficients to be
#' estimated.
#' @inheritParams limma::topTable
#' @param weights non-negative precision weights. Can be a numeric matrix of
#' individual weights of same size as the object expression matrix, or a numeric
#' vector of array weights with length equal to `ncol` of the expression matrix,
#' or a numeric vector of feature weights with length equal to `nrow `of the
#' expression matrix.
#'
#' @return A list object containing the matrix of p-values, the matrix of
#' summary statistics for each tag, and a suggested name of the final object
#' considering the parameters passed to the function.
#'
#' @seealso [limma::voom()] for the mean-variance relationship estimation,
#' [limma::lmFit()] for the linear model framework.

DA_limma <- function(object,
                     pseudo_count = FALSE,
                     design = NULL,
                     coef = 2,
                     norm = c("TMM", "TMMwsp", "RLE", "upperquartile", # edgeR
                              "posupperquartile", "none", # edgeR
                              "ratio", "poscounts", "iterate", # DESeq2
                              "TSS", "CSSmedian", "CSSdefault"), # others
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
    name <- "limma"

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
            stop(paste0("Can't find the ", NF.col," column in your object. Make
            sure to add the normalization factors column in your object first."
                        ))

        NFs = unlist(metadata[,NF.col])
        # DESeq2 NFs are supplied -> make them scaling factors!
        if(is.element(norm, c("ratio", "poscounts", "iterate"))){
            NFs <- NFs/colSums(counts)
        }
        NFs = NFs/exp(mean(log(NFs))) # Make NFs multiply to 1
    }

    name <- paste(name, ".", norm, sep = "")
    message(paste0("Differential abundance on ", norm," normalized data"))

    if(is.character(design)){
        if(grepl("~", design))
            design <- as.formula(design)
        else design <- as.formula(paste0("~", design))
    }

    if(class(design) == "formula")
        design <- stats::model.matrix(object = design,
                                      data = data.frame(metadata))

    v <- limma::voom(counts = counts,
                     design = design,
                     lib.size = colSums(counts) * NFs,
                     plot = FALSE)

    if(missing(weights)){
        message("Estimating Differential Abundance without weighting")
        fit <- limma::lmFit(object = v,
                            design = design)
    } else {
        message("Estimating Differential Abundance with weights")
        name <- paste(name,".weighted",sep = "")
        fit <- limma::lmFit(object = v,
                            design = design,
                            weights = v$weights * weights)
    }

    fit <- limma::eBayes(fit)

    message(paste0("Extracting results for ",
                   colnames(fit$coefficients)[coef],
                   " coefficient"))

    statInfo <- limma::topTable(fit,
                                coef = 2,
                                n = nrow(counts),
                                sort.by="none")

    pValMat <- cbind("rawP" = statInfo$P.Value,
                     "adjP" = statInfo$adj.P.Val)
    rownames(pValMat) = rownames(statInfo)

    return(list("pValMat" = pValMat,
                "statInfo" = statInfo,
                "name" = name))

}# END - function: DA_limma

#' @title DA_metagenomeSeq
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data
#' @importFrom phyloseq phyloseq_to_metagenomeSeq
#' @importFrom metagenomeSeq fitZig MRfulltable
#' @importFrom stats model.matrix
#' @export
#' @description
#' Fast run for the metagenomeSeq's differential abundance detection method.
#'
#' @param object phyloseq object.
#' @param pseudo_count Add 1 to all counts if TRUE (default = FALSE).
#' @param design The model for the count distribution. Can be the variable name,
#' or a character similar to "~ 1 + group", or a formula, or a `model.matrix`
#' object.
#' @param norm name of the normalization method used to compute the
#' normalization factors to use in the differential abundance analysis.
#' @inheritParams metagenomeSeq::fitZig
#' @inheritParams metagenomeSeq::MRfulltable
#'
#' @return A list object containing the matrix of p-values, the matrix of
#' summary statistics for each tag, and a suggested name of the final object
#' considering the parameters passed to the function.
#'
#' @seealso [metagenomeSeq::fitZig()] for the Zero-Inflated Gaussian regression
#' model estimation and [metagenomeSeq::MRfulltable()] for results extraction.

DA_metagenomeSeq <- function(object,
                             pseudo_count = FALSE,
                             design = NULL,
                             coef = 2,
                             norm = c("TMM", "TMMwsp", "RLE", "upperquartile",
                                      "posupperquartile", "none", # edgeR
                                      "ratio", "poscounts", "iterate", # DESeq2
                                      "TSS", "CSSmedian",
                                      "CSSdefault")) # others
    {


    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
    {
        object <- t(object)
    } else {}

    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)
    obj <- phyloseq::phyloseq_to_metagenomeSeq(object)

    # Name building
    name <- "metagenomeSeq"

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

    NF.col <- paste("NF", norm, sep = ".")
    if(norm == "TSS"){
        phyloseq::sample_data(obj)[NF.col] <- NFs <- 1L
    } else {
        # Check if the column with the normalization factors is present
        if(!any(colnames(metadata) == NF.col))
            stop(paste0("Can't find the ", NF.col," column in your object. Make
            sure to add the normalization factors column in your object first."
            ))

        NFs = unlist(metadata[,NF.col])
        # DESeq2 NFs are supplied -> make them scaling factors!
        if(is.element(norm, c("ratio", "poscounts", "iterate"))){
            NFs <- NFs/colSums(counts)
        }
        NFs = NFs/exp(mean(log(NFs))) # Make NFs multiply to 1
    }

    name <- paste(name, ".", norm, sep = "")
    message(paste0("Differential abundance on ", norm," normalized data"))

    metagenomeSeq::normFactors(object = obj) <- NFs

    if(is.character(design)){
        if(grepl("~", design))
            design <- as.formula(design)
        else design <- as.formula(paste0("~", design))
    }

    if(class(design) == "formula"){
        design <- as.formula(paste0(paste0(design,collapse = " ")," + ",NF.col))
        design <- stats::model.matrix(object = design,
                                      data = data.frame(metadata))
    }

    suppressWarnings(fit <- try(metagenomeSeq::fitZig(obj = obj,
                                          mod = design,
                                          verbose = FALSE,
                                          useCSSoffset = FALSE,
                                          control = zigControl(maxit = 1000)),
                                silent = TRUE))

    if(class(fit) == "try-error"){
        res = matrix(NA, ncol = 2, nrow = nrow(counts))
        stop("Error! Something went wrong during fitZig estimation.")
    } else {
        message(paste0("Extracting results for ",
                       colnames(fit@eb$coefficients)[coef],
                       " coefficient"))
        # You need to specify all OTUs to get the full table from MRfulltable.
        res <- metagenomeSeq::MRfulltable(obj = fit,
                                          number = nrow(counts),
                                          coef = 2)
    }

    pValMat <- cbind("rawP" = res$pvalues,
                     "adjP" = res$adjPvalues)
    rownames(pValMat) = rownames(res)
    colnames(res) <- c("rawP", "adjP")
    pValMat <- as.matrix(res)

    lods <- fit@eb$lods
    statInfo <- cbind(res,lods)
    return(list("pValMat" = pValMat,
                "statInfo" = statInfo,
                "name" = name))

}# END - function: DA_metagenomeSeq

#' @title DA_ALDEx2
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data
#' @importFrom ALDEx2 aldex
#' @export
#' @description
#' Fast run for the ALDEx2's differential abundance detection method.
#'
#' @param object phyloseq object.
#' @param pseudo_count Add 1 to all counts if TRUE (default = FALSE).
#' @param seed seed to use for the MC sampling.
#' @param norm name of the normalization method used to compute the
#' normalization factors to use in the differential abundance analysis.
#' @inheritParams ALDEx2::aldex
#'
#' @return A list object containing the matrix of p-values, the matrix of
#' summary statistics for each tag, and a suggested name of the final object
#' considering the parameters passed to the function.
#'
#' @seealso [ALDEx2::aldex] for the Dirichlet-Multinomial model estimation.
#' Several and more complex tests are present in the ALDEx2 framework.

DA_ALDEx2 <- function(object,
                      pseudo_count = FALSE,
                      seed = 1,
                      conditions = NULL,
                      mc.samples = 128,
                      test = c("t","wilcox"),
                      denom = "iqlr",
                      norm = c("TMM", "TMMwsp", "RLE", "upperquartile",
                               "posupperquartile", "none", # edgeR
                               "ratio", "poscounts", "iterate", # DESeq2
                               "TSS", "CSSmedian", "CSSdefault")) # others
{
    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
    {
        object <- t(object)
    } else {}

    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)

    # Name building
    name <- "ALDEx2"

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

    NF.col <- paste("NF", norm, sep = ".")

    # Check if the column with the normalization factors is present
    if(!any(colnames(metadata) == NF.col)){
        stop(paste0("Can't find the ", NF.col," column in your object. Make
        sure to add the normalization factors column in your object first."
        ))
    }

    name <- paste(name, ".", norm, sep = "")

    NFs = unlist(metadata[,NF.col])
    # Check if the NFs are scaling factors. If so, make them norm. factors
    if(is.element(norm, c("TMM", "TMMwsp", "RLE", "upperquartile",
                          "posupperquartile", "CSSmedian", "CSSdefault",
                          "TSS"))){
        NFs <- NFs * colSums(counts)
    }

    NFs <- NFs/exp(mean(log(NFs)))
    norm_counts <- round(counts %*% diag(1/NFs),
                         digits = 0)
    colnames(norm_counts) <- colnames(counts)

    if(is.null(conditions))
        stop("Please supply the name of the variable of interest or the
             entire character vector.")
    else if(length(conditions) == 1)
        conditions = unlist(metadata[,conditions])

    name <- paste(name, ".", denom, sep = "")

    if(!is.element(test, c("t","wilcox")) | length(test) != 1)
        stop("Please choose between p-values produced by Welch t-test (t) or by
             the Wilcoxon test (wilcox).")

    name <- paste(name, ".", test, sep = "")

    set.seed(seed = seed)
    statInfo <- ALDEx2::aldex(reads = norm_counts,
                              conditions = conditions,
                              mc.samples = mc.samples,
                              test = test,
                              effect = TRUE,
                              include.sample.summary = FALSE,
                              denom = denom,
                              verbose = TRUE)

    if(test == "t"){
        pValMat <- cbind("rawP" = statInfo$we.ep,
                         "adjP" = statInfo$we.eBH)
    } else
        pValMat <- cbind("rawP" = statInfo$wi.ep,
                         "adjP" = statInfo$wi.eBH)

    return(list("pValMat" = pValMat,
                "statInfo" = statInfo,
                "name" = name))

}# END - function: DA_ALDEx2

#' @title DA_MAST
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data
#' @importFrom corncob differentialTest
#' @importFrom stats coef
#' @export
#' @description
#' Fast run for corncob differential abundance detection method.
#'
#' @param object phyloseq object.
#' @param pseudo_count Add 1 to all counts if TRUE (default = FALSE).
#' @param coefficient The coefficient of interest as a single word formed by the
#' variable name and the non reference level. (e.g.: 'ConditionDisease' if the
#' reference level for the variable 'Condition' is 'control').
#' @param norm name of the normalization method used to compute the
#' normalization factors to use in the differential abundance analysis.
#' @inheritParams corncob::differentialTest
#'
#' @return A list object containing the matrix of p-values, the matrix of
#' summary statistics for each tag, and a suggested name of the final object
#' considering the parameters passed to the function.
#'
#' @seealso [corncob::bbdml()] and [corncob::differentialTest] for differential
#' abundance and differential variance evaluation.

DA_corncob <- function(object,
                       pseudo_count = FALSE,
                       formula,
                       phi.formula,
                       formula_null,
                       phi.formula_null,
                       test,
                       boot = FALSE,
                       coefficient = NULL,
                       norm = c("TMM", "TMMwsp", "RLE", "upperquartile",
                                "posupperquartile", "none", # edgeR
                                "ratio", "poscounts", "iterate", # DESeq2
                                "TSS", "CSSmedian", "CSSdefault")) # others)
{

    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
    {
        object <- t(object)
    } else {}

    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)

    # Name building
    name <- "corncob"

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

    NF.col <- paste("NF", norm, sep = ".")
    # Check if the column with the normalization factors is present
    if(!any(colnames(metadata) == NF.col)){
        stop(paste0("Can't find the ", NF.col," column in your object. Make
        sure to add the normalization factors column in your object first."
        ))
    }

    name <- paste(name, ".", norm, sep = "")

    NFs = unlist(metadata[,NF.col])
    # Check if the NFs are scaling factors. If so, make them norm. factors
    if(is.element(norm, c("TMM", "TMMwsp", "RLE", "upperquartile",
                          "posupperquartile", "CSSmedian", "CSSdefault",
                          "TSS"))){
        NFs <- NFs * colSums(counts)
    }

    NFs <- NFs/exp(mean(log(NFs)))
    norm_counts <- round(counts %*% diag(1/NFs),
                         digits = 0)
    colnames(norm_counts) <- colnames(counts)

    message(paste0("Differential abundance on ", norm," normalized data"))

    if(missing(test)){
        stop("Please supply the test to perform, 'Wald' or 'LRT'.")
    } else name <- paste(name, ".", test, sep = "")

    if(boot)
        name <- paste(name, ".", "boot", sep = "")

    ## differential expression
    requireNamespace("phyloseq")
    fit <- corncob::differentialTest(formula = formula,
                                     phi.formula = phi.formula,
                                     formula_null = formula_null,
                                     phi.formula_null = phi.formula_null,
                                     data = norm_counts,
                                     sample_data = metadata,
                                     test = test,
                                     boot = boot)

    # differentialTest's output has a complex structure:
    # extraction of summary table for each estimated model
    # mu.(Intercept), mu.condition, and phi.(Intercept), phi.condition
    # only the mu.condition is of interest.

    if(is.null(coefficient) | !is.element(coefficient,fit$restrictions_DA)){
        stop("Please supply the coefficient of interest as a single word formed
        by the variable name and the non reference level. (e.g.:
        'ConditionDisease' if the reference level for the variable 'Condition'
        is 'control')")
    }

    statInfo <- plyr::ldply(.data = fit$all_models,
                            .fun = function(model) {
                                stats::coef(model)[paste0("mu.",coefficient),]
                            })

    if(length(fit$restrictions_DA)>0){
        message(paste("Differential abundance across",
                      paste0(fit$restrictions_DA, collapse = " and ")))
    }
    if(length(fit$restrictions_DV)>0){
        message(paste("Differential variability across",
                      paste0(fit$restrictions_DV, collapse = " and ")))
    }

    pValMat <- data.frame("pval" = fit$p, "adjp" = fit$p_fdr)
    rownames(statInfo) <- rownames(pValMat) <- rownames(counts)
    list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name)

}# END - function: corncob

#' @title DA_MAST
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data
#' @importFrom MAST zlm summary FromMatrix
#' @importFrom SummarizedExperiment colData assay
#' @export
#' @description
#' Fast run for MAST differential abundance detection method.
#'
#' @param object phyloseq object.
#' @param pseudo_count Add 1 to all counts if TRUE (default = FALSE).
#' @param rescale Rescale count data, per million if 'default', or per median
#' library size if 'median' ('median' is suggested for metagenomics data).
#' @param design The model for the count distribution. Can be the variable name,
#' or a character similar to "~ 1 + group", or a formula, or a `model.matrix`
#' object.
#' @param coefficient The coefficient of interest as a single word formed by the
#' variable name and the non reference level. (e.g.: 'ConditionDisease' if the
#' reference level for the variable 'Condition' is 'control')
#' @param norm name of the normalization method used to compute the
#' normalization factors to use in the differential abundance analysis.
#'
#' @return A list object containing the matrix of p-values, the matrix of
#' summary statistics for each tag, and a suggested name of the final object
#' considering the parameters passed to the function.
#'
#' @seealso [MAST::zlm] for the Truncated Gaussian Hurdle model estimation.

DA_MAST <- function(object,
                    pseudo_count = FALSE,
                    rescale = c("median","default"),
                    design,
                    coefficient = NULL,
                    norm = c("TMM", "TMMwsp", "RLE", "upperquartile",
                             "posupperquartile", "none", # edgeR
                             "ratio", "poscounts", "iterate", # DESeq2
                             "TSS", "CSSmedian", "CSSdefault")) # others
{

    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
    {
        object <- t(object)
    } else {}

    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)

    # Name building
    name <- "MAST"

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

    NF.col <- paste("NF", norm, sep = ".")
    # Check if the column with the normalization factors is present
    if(!any(colnames(metadata) == NF.col)){
        stop(paste0("Can't find the ", NF.col," column in your object. Make
        sure to add the normalization factors column in your object first."
        ))
    }

    name <- paste(name, ".", norm, sep = "")

    NFs = unlist(metadata[,NF.col])
    # Check if the NFs are scaling factors. If so, make them norm. factors
    if(is.element(norm, c("TMM", "TMMwsp", "RLE", "upperquartile",
                          "posupperquartile", "CSSmedian", "CSSdefault",
                          "TSS"))){
        NFs <- NFs * colSums(counts)
    }

    NFs <- NFs/exp(mean(log(NFs)))
    norm_counts <- round(counts %*% diag(1/NFs),
                         digits = 0)
    colnames(norm_counts) <- colnames(counts)

    message(paste0("Differential abundance on ", norm," normalized data"))

    if(length(rescale) > 1 | !is.element(rescale,c("default","median"))){
        stop("Please choose between 'default' or 'median' for the rescale
             parameter. 'median' is suggested for metagenomics data.")
    } else if(rescale == "median"){
        message(paste0("per ",rescale,"-lib.size rescaled data"))
        tpm <- norm_counts * median(colSums(norm_counts)) / colSums(norm_counts)
    } else {
        message(paste0(rescale," (per million) rescaled data"))
        tpm <- norm_counts * 10^6 / colSums(norm_counts)
    }

    name <- paste(name, ".", rescale, sep = "")
    tpm <- log2(tpm + 1)
    sca <- MAST::FromMatrix(exprsArray = tpm,
                            cData = metadata)
    # here, we keep all OTUs so that we can fairly compare MAST and the other
    # methods. So, no adaptive thresholding or filtering by gene expression.
    SummarizedExperiment::assays(sca) <-
        list(tpm = SummarizedExperiment::assay(sca))
    ngeneson <- apply(norm_counts, 2, function(x) mean(x>0))
    metadata$cngeneson <- ngeneson - mean(ngeneson)
    SummarizedExperiment::colData(sca)$cngeneson <- metadata$cngeneson

    if(is.character(design)){
        if(grepl("~", design))
            design <- as.formula(design)
        else design <- as.formula(paste0("~", design))
    }

    if(class(design) == "formula"){
        design <- as.formula(paste0(paste0(design,
                                           collapse = " "), " + cngeneson"))
    }

    ## differential expression
    fit <- MAST::zlm(formula = design,
                     sca = sca,
                     method = "bayesglm",
                     ebayes = TRUE)
    if(!is.element(coefficient, colnames(fit@coefC)))
        stop("Please supply the coefficient of interest as a single word formed
        by the variable name and the non reference level. (e.g.:
        'ConditionDisease' if the reference level for the variable 'Condition'
        is 'control')")
    summaryDt = MAST::summary(fit, doLRT = coefficient)$datatable

    contrast <- component <- NULL
    fcHurdle <- merge(x = summaryDt[contrast == coefficient &
                                        component == 'H',
                                c("primerid", "Pr(>Chisq)")],
                      y = summaryDt[contrast == coefficient &
                                        component == "logFC",
                                c("primerid", "coef", "ci.hi", "ci.lo")],
                      by = "primerid")
    statInfo = data.frame(logFC = fcHurdle$coef,
                          logFC.lo = fcHurdle$ci.lo,
                          logFC.hi = fcHurdle$ci.hi,
                          pval = fcHurdle$`Pr(>Chisq)`,
                          padj = p.adjust(fcHurdle$`Pr(>Chisq)`, 'BH'))
    rownames(statInfo) <- fcHurdle$primerid

    pValMat <- statInfo[, c("pval", "padj")]
    list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name)
}# END - function: MAST

#' @title DA_Seurat
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data
#' @importFrom Seurat CreateSeuratObject AddMetaData NormalizeData
#' @importFrom Seurat FindVariableFeatures ScaleData FindMarkers
#' @export
#' @description
#' Fast run for Seurat differential abundance detection method.
#'
#' @param object phyloseq object.
#' @param pseudo_count Add 1 to all counts if TRUE (default = FALSE).
#' @inheritParams Seurat::FindMarkers
#' @param group.by the name of a factor variable to group samples.
#' @param reference.level the level of the factor variable to take as reference
#' to compute differential abundance.
#' @param norm name of the normalization method used to compute the
#' normalization factors to use in the differential abundance analysis.
#'
#' @return A list object containing the matrix of p-values, the matrix of
#' summary statistics for each tag, and a suggested name of the final object
#' considering the parameters passed to the function.
#'
#' @seealso [Seurat::CreateSeuratObject] to create the Seurat object,
#' [Seurat::AddMetaData()] to add metadata information,
#' [Seurat::NormalizeData()] to compute the normalization for the counts,
#' [Seurat::FindVariableFeatures()] to estimate the mean-variance trend,
#' [Seurat::ScaleData()] to scale and center features in the dataset, and
#' [Seurat::FindMarkers()] to perform differential abundance analysis.

DA_Seurat <- function(object,
                      pseudo_count = FALSE,
                      test.use = "wilcox",
                      group.by = NULL,
                      reference.level = NULL,
                      norm = c("TMM", "TMMwsp", "RLE", "upperquartile",
                               "posupperquartile", "none", # edgeR
                               "ratio", "poscounts", "iterate", # DESeq2
                               "TSS", "CSSmedian", "CSSdefault")) # others
{

    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
    {
        object <- t(object)
    } else {}

    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)

    # Name building
    name <- "Seurat"

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

    NF.col <- paste("NF", norm, sep = ".")
    # Check if the column with the normalization factors is present
    if(!any(colnames(metadata) == NF.col)){
        stop(paste0("Can't find the ", NF.col," column in your object. Make
        sure to add the normalization factors column in your object first."
        ))
    }

    name <- paste(name, ".", norm, sep = "")

    NFs = unlist(metadata[,NF.col])
    # Check if the NFs are scaling factors. If so, make them norm. factors
    if(is.element(norm, c("TMM", "TMMwsp", "RLE", "upperquartile",
                          "posupperquartile", "CSSmedian", "CSSdefault",
                          "TSS"))){
        NFs <- NFs * colSums(counts)
    }

    NFs <- NFs/exp(mean(log(NFs)))
    norm_counts <- round(counts %*% diag(1/NFs),
                         digits = 0)
    colnames(norm_counts) <- colnames(counts)

    message(paste0("Differential abundance on ", norm," normalized data"))

    # Initialize the Seurat object with the raw (non-normalized data).
    # If the chosen normalization is 'none', then this is true.
    # Keep all features expressed in >= 1 sample
    # Keep all samples with at least 1 detected feature.
    sobj <- Seurat::CreateSeuratObject(counts = norm_counts,
                                       min.cells = 1,
                                       min.features = 1)

    sobj <- Seurat::AddMetaData(object = sobj,
                                metadata = data.frame(metadata),
                                col.name = colnames(metadata))

    if(is.null(group.by) | is.null(reference.level)){
        stop("Please supply the variable to study differential abundance for and
             its reference level.")
    } else{
        if(!is(unlist(sobj[[group.by]]),"factor")){
            stop(paste(group.by,"variable is not a factor. Please supply a
                       factor."))
        } else{
            if(!is.element(reference.level,levels(unlist(sobj[[group.by]])))){
                stop(paste(reference.level, "is not a level of the",  group.by,
                           "variable. Please supply a present category."))
            }
        }
    }

    message("Differential abundance analysis on ", group.by, " variable.")

    sobj <- Seurat::NormalizeData(object = sobj,
                                  normalization.method = "LogNormalize",
                                  scale.factor = 10000)
    sobj <- Seurat::FindVariableFeatures(object = sobj,
                                         nfeatures = round(nrow(counts)*0.1,
                                                           digits = 0))
    sobj <- Seurat::ScaleData(object = sobj,
                              vars.to.regress = c("nCount_RNA"))

    statInfo <- Seurat::FindMarkers(sobj,
                                    test.use = test.use,
                                    group.by = group.by,
                                    ident.1 = reference.level,
                                    logfc.threshold = 0)

    name <- paste(name, ".", test.use, sep = "")

    pValMat <- data.frame(pval = statInfo$p_val,
                          adjp = statInfo$p_val_adj)
    rownames(pValMat) <- rownames(statInfo)

    list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name)

}# END - function: DA_Seurat
