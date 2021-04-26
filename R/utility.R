#' @title createMocks
#'
#' @export
#' @description
#' Given the number of samples of the dataset from which the mocks should be
#' created, this function produces a `data.frame` object with as many rows as
#' the number of mocks and as many columns as the number of samples. If an odd
#' number of samples is given, the lower even integer will be considered in
#' order to obtain a balanced design for the mocks.
#' @param nsamples an integer representing the total number of samples.
#' @param N number of mock comparison to generate.
#'
#' @return A `data.frame` containing `N` rows and `nsamples` columns (if even).
#' Each cell of the data frame contains the "grp1" or "grp2" characters which
#' represent the mock groups pattern.
#'
#' @examples
#' data(ps_stool_16S)
#' # Generate the patterns for 100 mock comparison for an experiment
#' mocks <- createMocks(nsamples = phyloseq::nsamples(ps_stool_16S), N = 100)
#' head(mocks)

# Create data.frame with random labels
createMocks <- function(nsamples, N = 1000){
    sample_num <- nsamples %/% 2 * 2 # Balanced design sample numerosity
    mock_df <- matrix(NA, nrow = N, ncol = sample_num)
    for(i in seq_len(N)){ # N random balanced relabellings
        grps <- rep("grp1", sample_num)
        grps[sample(seq_len(sample_num),size = sample_num/2)] <- "grp2"
        mock_df[i,] <- grps}
    rownames(mock_df) <- paste0("Comparison",seq_len(N))
    return(mock_df)
}# END - function: createMocks

#' @title createSplits
#'
#' @export
#' @description
#' Given the phyloseq object from which the random splits should be created,
#' this function produces a list of 2 `data.frame` objects: \code{Subset1} and
#' \code{Subset2} with as many rows as the number of splits and as many columns
#' as the half of the number of samples.
#' @param object a \code{phyloseq} object.
#' @param varName name of a factor variable with 2 levels.
#' @param paired name of the unique subject identifier variable. If specified,
#' paired samples will remain in the same split. (default = NULL).
#' @param balanced If \code{TRUE} a balanced design will be created for the
#' splits. (Ignored if paired is supplied).
#' @param N number of splits to generate.
#'
#' @return A list of 2 `data.frame` objects: \code{Subset1} and
#' \code{Subset2} containing `N` rows and half of the total number of samples
#' columns. Each cell contains a unique sample identifier.
#' @examples
#' data(ps_plaque_16S)
#' set.seed(123)
#' # Balanced design for repeated measures
#' splits_df <- createSplits(object = ps_plaque_16S, varName =
#'     "HMP_BODY_SUBSITE", paired = "RSID", balanced = TRUE, N = 100)
#' # Balanced design for independent samples
#' splits_df <- createSplits(object = ps_plaque_16S, varName =
#'     "HMP_BODY_SUBSITE", balanced = TRUE, N = 100)
#' # Unbalanced design
#' splits_df <- createSplits(object = ps_plaque_16S, varName =
#'     "HMP_BODY_SUBSITE", balanced = FALSE, N = 100)

createSplits <- function(object, varName = NULL, paired = NULL, balanced = TRUE,
    N = 1000){
    metadata <- data.frame(phyloseq::sample_data(object))
    # Take variable names and levels
    if(is.null(varName)){
        stop("Please supply the name of the variable to perform the splitting")}
    else if(!is.element(varName,colnames(metadata))){
        stop("Variable not found")}
    else {
        variable <- metadata[,varName]
        if(!is.factor(variable)){
            warning(paste("The variable", varName, "is not a factor.",
                "Coercing to factor."))
            variable <- as.factor(variable)}
        var_levels <- levels(variable)
        if(length(var_levels) != 2)
            stop(paste("The variable", varName, "has not 2 levels."))
        else {
            grp1_name = var_levels[1]
            grp2_name = var_levels[2]}}
    if(!is.null(paired)){
        # 1 ID is enough to identify the two paired samples
        ID <- unique(metadata[,paired])
        num_ID <- length(ID)
        half = num_ID/2
        Subset1 <- matrix(NA, nrow = N, ncol = 2 * half)
        Subset2 <- matrix(NA, nrow = N, ncol = 2 * half)
        for(i in seq_len(N)){
            chosenID <- sample(x = ID, size = half, replace = FALSE)
            Subset1[i,] <- rownames(metadata[metadata[,paired] %in% chosenID,])
            Subset2[i,] <- rownames(metadata[metadata[,paired] %in% setdiff(ID,
                chosenID),])}}
    else {
        # find indexes for the 2 levels
        ID1 <- rownames(metadata[metadata[,varName] == grp1_name,])
        ID2 <- rownames(metadata[metadata[,varName] == grp2_name,])
        num_ID1 <- length(ID1)
        num_ID2 <- length(ID2)
        if(balanced){
            min_length <- min(num_ID1,num_ID2)
            ID1 <- ID1[seq_len(min_length)]
            ID2 <- ID2[seq_len(min_length)]
            num_ID1 <- num_ID2 <- min_length
            half <- rep(min_length %/% 2,2)
        } else half <- c(num_ID1,num_ID2) %/% 2
        Subset1 <- matrix(NA, nrow = N, ncol = sum(half))
        Subset2 <- matrix(NA, nrow = N, ncol = sum(c(num_ID1,num_ID2))-sum(
            half))
        for(i in seq_len(N)){
            chosenID1 <- sample(x = ID1, size = half[1], replace = FALSE)
            chosenID2 <- sample(x = ID2, size = half[2], replace = FALSE)
            Subset1[i,] <- c(chosenID1, chosenID2)
            Subset2[i,] <- c(setdiff(ID1,chosenID1),setdiff(ID2,chosenID2))}}
    rownames(Subset1) <- rownames(Subset2) <- paste0("Comparison",seq_len(N))
    return(list("Subset1" = Subset1, "Subset2" = Subset2))
}# END - function: createSplits

### Normalizations ###

#' @title norm_edgeR
#'
#' @importFrom edgeR calcNormFactors
#' @importFrom phyloseq sample_sums taxa_are_rows otu_table
#' @importFrom stats quantile
#' @export
#' @description
#' Calculate scaling factors from a phyloseq object to scale the raw library
#' sizes. Inherited from edgeR `calcNormFactors` function.
#'
#' @param object phyloseq object containing the counts to be normalized.
#' @param method normalization method to be used. Choose between "TMM",
#' "TMMwsp", "RLE", "upperquartile", "posupperquartile" or "none".
#' @inheritParams edgeR::calcNormFactors
#'
#' @return A new column containing the chosen edgeR-based scaling factors is
#' added to the phyloseq `sample_data` slot. The effective library sizes for
#' use in downstream analysis must be multiplied by the normalization factors.
#' @seealso \code{\link[edgeR]{calcNormFactors}} for details.
#'
#' @examples
#' data(ps_stool_16S)
#' # Calculate the scaling factors
#' ps_stool_16S <- norm_edgeR(object = ps_stool_16S, method = "TMM")
#'
#' # The phyloseq object now contains the scaling factors:
#' scaleFacts <- phyloseq::sample_data(ps_stool_16S)[,"NF.TMM"]
#' head(scaleFacts)
#'
#' # VERY IMPORTANT: to convert scaling factors to normalization factors
#' # multiply them by the library sizes and renormalize.
#' normFacts = scaleFacts * phyloseq::sample_sums(ps_stool_16S)
#' # Renormalize: multiply to 1
#' normFacts = normFacts/exp(mean(log(normFacts)))

norm_edgeR <- function(object, method = c("TMM", "TMMwsp", "RLE",
    "upperquartile", "posupperquartile", "none"), refColumn = NULL,
    logratioTrim = .3, sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e10,
    p = 0.75, ...){
    counts <- as(phyloseq::otu_table(object), "matrix")
    if (!phyloseq::taxa_are_rows(object))
        counts <- t(counts)
    if (is.na(method))
        stop("Please, supply a valid method between 'TMM', 'TMMwsp', 'RLE',
            'upperquartile', 'posupperquartile' or 'none'")
    else if (method == "posupperquartile"){
        scaledCounts <- t(counts) / colSums(counts)
        tmpNF <- apply(scaledCounts, MARGIN = 1L, FUN = function(x)
            quantile(x[x != 0], probs = .75))
        normFacts <- tmpNF/exp(mean(log(tmpNF)))
    } else {
        normFacts <- edgeR::calcNormFactors(counts, method = method, refColumn =
            refColumn, logratioTrim = logratioTrim, sumTrim = sumTrim,
            doWeighting = doWeighting, Acutoff = Acutoff, p = p, ...) }
    # END - ifelse: posupperquartile = upperquartile only of non-zero counts
    if (all(is.na(normFacts)))
        stop("Failed to compute normalization factors!")
    phyloseq::sample_data(object)[,paste("NF", method, sep = ".")] <- normFacts
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
#' @param ... other parameters for DESeq2
#' \code{\link[DESeq2]{estimateSizeFactors}} function.
#'
#' @return A new column containing the chosen DESeq2-based normalization factors
#' is added to the phyloseq `sample_data` slot.
#'
#' @seealso \code{\link[DESeq2]{estimateSizeFactors}} for details.
#'
#' @examples
#' data(ps_stool_16S)
#' # Calculate the normalization factors
#' ps_stool_16S <- norm_DESeq2(object = ps_stool_16S, method = "poscounts")
#'
#' # The phyloseq object now contains the normalization factors:
#' normFacts <- phyloseq::sample_data(ps_stool_16S)[,"NF.poscounts"]
#' head(normFacts)
#'
#' # VERY IMPORTANT: to convert normalization factors to scaling factors divide
#' # them by the library sizes and renormalize.
#' scaleFacts = normFacts / phyloseq::sample_sums(ps_stool_16S)
#' # Renormalize: multiply to 1
#' scaleFacts = scaleFacts/exp(mean(log(scaleFacts)))

norm_DESeq2 <- function(object, method = c("ratio", "poscounts", "iterate"),
    ...){
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    ## Calculate size factors
    obj <- phyloseq::phyloseq_to_deseq2(object, design = ~ 1)
    if(missing(method))
        stop("Please supply a normalization method between 'ratio', 'poscounts'
            or 'iterate'.")
    normFacts <- DESeq2::sizeFactors(DESeq2::estimateSizeFactors(obj,
        type = method, ...))
    phyloseq::sample_data(object)[,paste("NF", method, sep = ".")] <- normFacts
    return(object)
}# END - function: norm_DESeq2

#' @title norm_CSS
#'
#' @importFrom metagenomeSeq newMRexperiment
#' @importFrom phyloseq taxa_are_rows otu_table
#' @importFrom stats median
#' @export
#' @description
#' Calculate scaling factors from a phyloseq object to scale the raw library
#' sizes. Inherited from metagenomeSeq `calcNormFactors` function which performs
#' the Cumulative Sum Scaling normalization.
#'
#' @param object phyloseq object containing the counts to be normalized.
#' @param method  Normalization scaling parameter (default = "1000"). If
#' "median", the median of the normalization factors is used as scaling (Paulson
#' et al. 2013).
#'
#' @return A new column containing the CSS scaling factors is added to the
#' phyloseq `sample_data` slot.
#'
#' @seealso \code{\link[metagenomeSeq]{calcNormFactors}} for details.
#'
#' @examples
#' data(ps_stool_16S)
#' # Calculate the scaling factors
#' ps_stool_16S <- norm_CSS(object = ps_stool_16S, method = "median")
#'
#' # The phyloseq object now contains the scaling factors:
#' scaleFacts <- phyloseq::sample_data(ps_stool_16S)[, "NF.CSSmedian"]
#' head(scaleFacts)
#'
#' # VERY IMPORTANT: to convert scaling factors to normalization factors
#' # multiply them by the library sizes and renormalize.
#' normFacts = scaleFacts * phyloseq::sample_sums(ps_stool_16S)
#' # Renormalize: multiply to 1
#' normFacts = normFacts/exp(mean(log(normFacts)))

norm_CSS <- function(object, method = "default")
{
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    counts <- as(phyloseq::otu_table(object), "matrix")
    obj <- metagenomeSeq::newMRexperiment(counts = counts)
    normFacts <- metagenomeSeq::calcNormFactors(obj = obj)
    normFacts <- drop(as(normFacts, "matrix"))
    # Default: log2(normFacts/1000 + 1)
    # Original metagenomeSeq paper: log2(normFacts/median(libsize) +1)
    if(method == "default")
        normFacts <- log2(normFacts/1000 + 1)
    else if (method == "median")
        normFacts <- log2(normFacts/stats::median(normFacts) + 1)
    else stop("Please choose a scaling method between 'default' or 'median'.")
    # Remember to useCSSoffset = FALSE in fitZig function
    phyloseq::sample_data(object)[,paste0("NF.CSS", method)] <-
        normFacts
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
#' @return A new column containing the TSS scaling factors is added to the
#' phyloseq `sample_data` slot.
#'
#' @examples
#' data(ps_stool_16S)
#' # Calculate the scaling factors
#' ps_stool_16S <- norm_TSS(object = ps_stool_16S)
#'
#' # The phyloseq object now contains the scaling factors:
#' scaleFacts <- phyloseq::sample_data(ps_stool_16S)[, "NF.TSS"]
#' head(scaleFacts)
#'
#' # VERY IMPORTANT: to convert scaling factors to normalization factors
#' # multiply them by the library sizes and renormalize.
#' normFacts = scaleFacts * phyloseq::sample_sums(ps_stool_16S)
#' # Renormalize: multiply to 1
#' normFacts = normFacts/exp(mean(log(normFacts)))
#' # In this case they will be ones.

norm_TSS <- function(object)
{
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    normFacts <- 1/phyloseq::sample_sums(object)
    phyloseq::sample_data(object)[,"NF.TSS"] <- normFacts
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
#' @seealso \code{\link[zinbwave]{zinbFit}} for zero-inflated negative binomial
#' parameters' estimation and
#' \code{\link[zinbwave]{computeObservationalWeights}} for weights extraction.
#'
#' @examples
#' data(ps_stool_16S)
#' # Calculate the ZINB weights
#' zinbweights <- weights_ZINB(object = ps_stool_16S, K = 0, design = "~ 1")
#' head(zinbweights)

weights_ZINB <- function(object, design, K = 0, commondispersion = TRUE,
    zeroinflation = TRUE, verbose = FALSE, ...){
    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)
    if(is.character(design)){
        if(grepl("~", design))
            design <- as.formula(design)
        else design <- as.formula(paste0("~", design))}
    if(is(design, "formula"))
        design <- stats::model.matrix(object = design,
            data = data.frame(metadata))
    zinbmodel <- zinbwave::zinbFit(Y = counts, X = design, K = K,
        commondispersion = commondispersion, zeroinflation = TRUE, verbose =
        verbose, BPPARAM = BiocParallel::SerialParam(), ... = ...)
    w <- zinbwave::computeObservationalWeights(model = zinbmodel, x = counts)
    return(w)
}# END - function: weights_ZINB

### Differential Abundance Methods  ###

#' @title DA_edgeR
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data
#' @importFrom edgeR DGEList estimateDisp estimateGLMRobustDisp glmQLFit
#' glmQLFTest getDispersion
#' @importFrom stats model.matrix p.adjust
#' @export
#' @description
#' Fast run for edgeR differential abundance detection method.
#'
#' @param object phyloseq object.
#' @param pseudo_count Add 1 to all counts if TRUE (default = FALSE).
#' @param group vector or factor giving the experimental group/condition for
#' each sample/library.
#' @param design numeric design matrix. Defaults to \code{model.matrix(~group)}
#' if \code{group} is specified and otherwise to a single column of ones.
#' @param contrast numeric vector or matrix specifying one or more contrasts of
#' the linear model coefficients to be tested equal to zero.
#' @param robust logical, should the estimation of \code{prior.df} be
#' robustified against outliers?
#' @param coef integer or character index vector indicating which coefficients
#' of the linear model are to be tested equal to zero. Ignored if
#' \code{contrast} is not \code{NULL}.
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
#' @seealso \code{\link[edgeR]{DGEList}} for the edgeR DEG object creation,
#' \code{\link[edgeR]{estimateDisp}} and
#' \code{\link[edgeR]{estimateGLMRobustDisp}} for dispersion estimation, and
#' \code{\link[edgeR]{glmQLFit}} and \code{\link[edgeR]{glmQLFTest}} for the
#' quasi-likelihood negative binomial model fit.
#'
#' @examples
#' data(ps_stool_16S)
#' # Calculate the scaling factors
#' ps_stool_16S <- norm_edgeR(object = ps_stool_16S, method = "TMM")
#'
#' # The phyloseq object now contains the scaling factors:
#' scaleFacts <- phyloseq::sample_data(ps_stool_16S)[, "NF.TMM"]
#' head(scaleFacts)
#'
#' # Differential abundance
#' group <- sample(x = c("grp1","grp2"), size = phyloseq::nsamples(
#'     ps_stool_16S), replace = TRUE)
#' phyloseq::sample_data(ps_stool_16S)$group <- group
#' DA_edgeR(ps_stool_16S, group = group, design = ~ group, coef = 2,
#'     norm = "TMM")

DA_edgeR <- function(object, pseudo_count = FALSE, group = NULL, design = NULL,
    contrast = NULL, robust = FALSE, coef = 2, norm = c("TMM", "TMMwsp", "RLE",
    "upperquartile", "posupperquartile", "none", "ratio", "poscounts",
    "iterate", "TSS", "CSSmedian", "CSSdefault"), weights){
    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)
    # Name building
    name <- "edgeR"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        message("Adding a pseudo count... \n")
        counts <- counts + 1
        name <- paste(name,".pseudo",sep = "")}
    if(length(norm) > 1)
        stop("Please choose one normalization for this istance of differential
            abundance analysis.")
    if(norm == "TSS")
        NFs = 1
    else {
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
        if(is.element(norm, c("ratio", "poscounts", "iterate")))
            NFs <- NFs/colSums(counts)
        NFs = NFs/exp(mean(log(NFs)))} # Make NFs multiply to 1
    name <- paste(name, ".", norm, sep = "")
    message(paste0("Differential abundance on ", norm," normalized data"))
    dge <- edgeR::DGEList(counts = counts, norm.factors = NFs, group = group)
    if(missing(weights))
        message("Estimating Differential Abundance without weighting")
    else {
        message("Estimating Differential Abundance with weights")
        dge$weights <- weights
        name <- paste(name,".weighted",sep = "")}
    if(is.character(design)){
        if(grepl("~", design))
            design <- as.formula(design)
        else design <- as.formula(paste0("~", design))}
    if(is(design, "formula"))
        design <- stats::model.matrix(object = design, data = data.frame(
            metadata))
    if(!robust)
        dge <- edgeR::estimateDisp(y = dge, design = design)
    else {
        message(paste0("Estimating robust dispersions"))
        dge <- edgeR::estimateGLMRobustDisp(y = dge, design = design)
        name <- paste(name,".robust",sep = "")}
    message(paste0("Extracting results"))
    dispEsts <- edgeR::getDispersion(dge)
    glmFit <- edgeR::glmQLFit(y = dge, dispersion = dispEsts,
        robust = robust, design = design)
    glmRes <- edgeR::glmQLFTest(glmFit, coef = coef, contrast = contrast)
    if(missing(contrast))
        message(paste0("Extracting results for ", colnames(coef(glmRes))[coef],
            " coefficient"))
    else message(paste0("Extracting results for ", contrast, " contrasts"))
    statInfo <- glmRes[["table"]]
    pval <- statInfo[, "PValue"]
    pValMat <- data.frame("rawP" = pval, "adjP" = stats::p.adjust(pval, "BH"))
    rownames(pValMat) <- rownames(statInfo)
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "dispEsts" =
        dispEsts, "name" = name))
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
#' @param design (Required). A \code{\link{formula}} which specifies the design
#' of the experiment, taking the form \code{formula(~ x + y + z)}. That is, a
#' formula with right-hand side only. By default, the functions in this package
#' and DESeq2 will use the last variable in the formula (e.g. \code{z}) for
#' presenting results (fold changes, etc.) and plotting. When considering your
#' specification of experimental design, you will want to  re-order the levels
#' so that the \code{NULL} set is first. For example, the following line of code
#' would ensure that Enterotype 1 is used as the reference sample class in tests
#' by setting it to the first of the factor levels using the
#' \code{\link{relevel}} function:
#' \code{sample_data(entill)$Enterotype <-
#' relevel(sample_data(entill)$Enterotype, "1")}
#' @param contrast character vector with exactly three elements: the name of a
#' factor in the design formula, the name of the numerator level for the fold
#' change, and the name of the denominator level for the fold change.
#' @param alpha the significance cutoff used for optimizing the independent
#' filtering (by default 0.05). If the adjusted p-value cutoff (FDR) will be a
#' value other than 0.05, alpha should be set to that value.
#' @param norm name of the normalization method used to compute the
#' normalization factors to use in the differential abundance analysis.
#' @param weights optional numeric matrix giving observation weights.
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
#' data(ps_stool_16S)
#' # Calculate the scaling factors
#' ps_stool_16S <- norm_DESeq2(object = ps_stool_16S, method = "poscounts")
#'
#' # The phyloseq object now contains the normalization factors:
#' normFacts <- phyloseq::sample_data(ps_stool_16S)[, "NF.poscounts"]
#' head(normFacts)
#'
#' # Differential abundance
#' group <- sample(x = c("grp1","grp2"), size = phyloseq::nsamples(
#'     ps_stool_16S), replace = TRUE)
#' phyloseq::sample_data(ps_stool_16S)$group <- group
#' DA_DESeq2(ps_stool_16S, design = ~ group, contrast = c("group", "grp2",
#'     "grp1"), norm = "poscounts")

DA_DESeq2 <- function(object, pseudo_count = FALSE, design = NULL, contrast =
    NULL, alpha = 0.05, norm = c("TMM", "TMMwsp", "RLE", "upperquartile",
    "posupperquartile", "none", "ratio", "poscounts", "iterate", "TSS",
    "CSSmedian", "CSSdefault"), weights){
    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)
    # Name building
    name <- "DESeq2"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        message("Adding a pseudo count... \n")
        counts <- counts + 1
        phyloseq::otu_table(object) <- counts
        name <- paste(name,".pseudo",sep = "")}
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
        "posupperquartile", "CSSmedian", "CSSdefault", "TSS"))){
        NFs <- NFs*colSums(counts)}
    DESeq2::sizeFactors(dds) = NFs/exp(mean(log(NFs))) # Make NFs multiply to 1
    name <- paste(name, ".", norm, sep = "")
    message(paste0("Differential abundance on ", norm," normalized data"))
    if(missing(weights))
        message("Estimating Differential Abundance without weighting")
    else {
        message("Estimating Differential Abundance with weights")
        weights[which(weights < 1e-6)] <- 1e-06
        SummarizedExperiment::assays(dds, withDimnames = FALSE)$weights <-
            weights
        name <- paste(name,".weighted",sep = "")}
    ### Run DESeq
    ddsRes <- DESeq2::DESeq(object = dds, test = "LRT", reduced = ~ 1,
        parallel = FALSE)
    dispEsts <- DESeq2::dispersions(ddsRes)
    if(missing(contrast) | (!is.character(contrast) & length(contrast) != 3))
        stop(paste0("Please supply a character vector with exactly three
            elements: the name of a factor in the design formula, the name of
            the numerator level for the fold change, and the name of the
            denominator level for the fold change."))
    else message(paste0("Extracting results for ", contrast[1]," variable, ",
        contrast[2], " / ", contrast[3]))
    res <- DESeq2::results(ddsRes, alpha = alpha, contrast = contrast)
    statInfo <- as(res, "data.frame")
    pValMat <- statInfo[,c("pvalue", "padj")]
    colnames(pValMat) <- c("rawP", "adjP")
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "dispEsts" =
        dispEsts, "name" = name))
}# END - function: DA_DESeq2

#' @title DA_limma
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom stats model.matrix weights
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
#' @param coef column number or column name specifying which coefficient or
#' contrast of the linear model is of interest. For \code{topTable}, can also be
#' a vector of column subscripts, in which case the gene ranking is by
#' F-statistic for that set of contrasts.
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
#' @seealso \code{\link[limma]{voom}} for the mean-variance relationship
#' estimation, \code{\link[limma]{lmFit}} for the linear model framework.
#'
#' @examples
#' data(ps_stool_16S)
#' # Calculate the scaling factors
#' ps_stool_16S <- norm_edgeR(object = ps_stool_16S, method = "TMM")
#'
#' # The phyloseq object now contains the scaling factors:
#' scaleFacts <- phyloseq::sample_data(ps_stool_16S)[, "NF.TMM"]
#' head(scaleFacts)
#'
#' # Differential abundance
#' group <- sample(x = c("grp1","grp2"), size = phyloseq::nsamples(
#'     ps_stool_16S), replace = TRUE)
#' phyloseq::sample_data(ps_stool_16S)$group <- group
#' DA_limma(ps_stool_16S, design = ~ group, coef = 2, norm = "TMM")

DA_limma <- function(object, pseudo_count = FALSE, design = NULL, coef = 2,
    norm = c("TMM", "TMMwsp", "RLE", "upperquartile", "posupperquartile",
    "none", "ratio", "poscounts", "iterate", "TSS", "CSSmedian", "CSSdefault"),
    weights){
    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)
    # Name building
    name <- "limma"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        message("Adding a pseudo count... \n")
        counts <- counts + 1
        name <- paste(name,".pseudo",sep = "")}
    if(length(norm) > 1)
        stop("Please choose one normalization for this istance of differential
        abundance analysis.")
    if(norm == "TSS")
        NFs = 1
    else {
        # Check if the column with the normalization factors is present
        NF.col <- paste("NF", norm, sep = ".")
        if(!any(colnames(metadata) == NF.col))
            stop(paste0("Can't find the ", NF.col," column in your object. Make
                sure to add the normalization factors column in your object
                first."))
        NFs = unlist(metadata[,NF.col])
        # DESeq2 NFs are supplied -> make them scaling factors!
        if(is.element(norm, c("ratio", "poscounts", "iterate")))
            NFs <- NFs/colSums(counts)
        NFs = NFs/exp(mean(log(NFs)))} # Make NFs multiply to 1
    name <- paste(name, ".", norm, sep = "")
    message(paste0("Differential abundance on ", norm," normalized data"))
    if(is.character(design)){
        if(grepl("~", design))
            design <- as.formula(design)
        else design <- as.formula(paste0("~", design))}
    if(is(design, "formula"))
        design <- stats::model.matrix(object = design, data = data.frame(
            metadata))
    v <- limma::voom(counts = counts, design = design, lib.size = colSums(
        counts) * NFs, plot = FALSE)
    if(missing(weights)){
        message("Estimating Differential Abundance without weighting")
        fit <- limma::lmFit(object = v, design = design)
    } else {
        message("Estimating Differential Abundance with weights")
        name <- paste(name,".weighted",sep = "")
        fit <- limma::lmFit(object = v, design = design, weights =
            stats::weights(v) * weights)}
    fit <- limma::eBayes(fit)
    message(paste0("Extracting results for ", colnames(coef(fit[, coef])),
        " coefficient"))
    statInfo <- limma::topTable(fit, coef = coef, n = nrow(counts), sort.by =
        "none")
    pValMat <- statInfo[, c("P.Value", "adj.P.Val")]
    colnames(pValMat) <- c("rawP", "adjP")
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
}# END - function: DA_limma

#' @title DA_metagenomeSeq
#'
#' @importFrom phyloseq taxa_are_rows otu_table sample_data taxa_names
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
#' @inheritParams metagenomeSeq::MRcoefs
#'
#' @return A list object containing the matrix of p-values, the matrix of
#' summary statistics for each tag, and a suggested name of the final object
#' considering the parameters passed to the function.
#'
#' @seealso \code{\link[metagenomeSeq]{fitZig}} for the Zero-Inflated Gaussian
#' regression model estimation and \code{\link[metagenomeSeq]{MRfulltable}}
#' for results extraction.
#'
#' @examples
#' data(ps_stool_16S)
#' # Calculate the scaling factors
#' ps_stool_16S <- norm_CSS(object = ps_stool_16S, method = "median")
#'
#' # The phyloseq object now contains the scaling factors:
#' scaleFacts <- phyloseq::sample_data(ps_stool_16S)[, "NF.CSSmedian"]
#' head(scaleFacts)
#'
#' # Differential abundance
#' group <- sample(x = c("grp1","grp2"), size = phyloseq::nsamples(
#'     ps_stool_16S), replace = TRUE)
#' phyloseq::sample_data(ps_stool_16S)$group <- group
#' DA_metagenomeSeq(ps_stool_16S, design = ~ group, coef = 2, norm =
#'     "CSSmedian")

DA_metagenomeSeq <- function(object, pseudo_count = FALSE, design = NULL, coef =
    2, norm = c("TMM", "TMMwsp", "RLE", "upperquartile", "posupperquartile",
    "none", "ratio", "poscounts", "iterate", "TSS", "CSSmedian", "CSSdefault")){
    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)
    obj <- phyloseq::phyloseq_to_metagenomeSeq(object)
    # Name building
    name <- "metagenomeSeq"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        message("Adding a pseudo count... \n")
        counts <- counts + 1
        name <- paste(name,".pseudo",sep = "")}
    if(length(norm) > 1)
        stop("Please choose one normalization for this istance of differential
        abundance analysis.")
    NF.col <- paste("NF", norm, sep = ".")
    if(norm == "TSS")
        phyloseq::sample_data(obj)[NF.col] <- NFs <- 1L
    else {
        # Check if the column with the normalization factors is present
        if(!any(colnames(metadata) == NF.col))
            stop(paste0("Can't find the ", NF.col," column in your object. Make
            sure to add the normalization factors column in your object first."
            ))
        NFs = unlist(metadata[,NF.col])
        # DESeq2 NFs are supplied -> make them scaling factors!
        if(is.element(norm, c("ratio", "poscounts", "iterate")))
            NFs <- NFs/colSums(counts)
        NFs = NFs/exp(mean(log(NFs)))} # Make NFs multiply to 1
    name <- paste(name, ".", norm, sep = "")
    message(paste0("Differential abundance on ", norm," normalized data"))
    metagenomeSeq::normFactors(object = obj) <- NFs
    if(is.character(design)){
        if(grepl("~", design))
            design <- as.formula(design)
        else design <- as.formula(paste0("~", design))}
    if(is(design, "formula")){
        design <- as.formula(paste0(paste0(design,collapse = " "), " + ",
            NF.col))
        design <- stats::model.matrix(object = design, data = data.frame(
            metadata))}
    suppressWarnings(fit <- try(metagenomeSeq::fitZig(obj = obj, mod = design,
        verbose = FALSE, useCSSoffset = FALSE, control = zigControl(maxit =
        1000)), silent = TRUE))
    if(is(fit, "try-error")){
        res = matrix(NA, ncol = 2, nrow = nrow(counts))
        stop("Error! Something went wrong during fitZig estimation.")
    } else {
        statInfo <- metagenomeSeq::MRcoefs(obj = fit, by = coef, number = nrow(
            counts))
        statInfo <- statInfo[phyloseq::taxa_names(object),]
        message(paste0("Extracting results for ", colnames(statInfo[coef]),
            " coefficient"))}
    pValMat <- statInfo[, c("pvalues", "adjPvalues")]
    colnames(pValMat) = c("rawP", "adjP")
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
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
#' @param norm name of the normalization method used to compute the
#' normalization factors to use in the differential abundance analysis.
#' @param mc.samples An integer. The number of Monte Carlo samples to use when
#' estimating the underlying distributions. Since we are estimating central
#' tendencies, 128 is usually sufficient.
#' @param denom A character string. Indicates which features to retain as the
#' denominator for the Geometric Mean calculation. Using "iqlr" accounts for
#' data with systematic variation and centers the features on the set features
#' that have variance that is between the lower and upper quartile of variance.
#' Using "zero" is a more extreme case where there are many non-zero features in
#' one condition but many zeros in another. In this case the geometric mean of
#' each group is calculated using the set of per-group non-zero features.
#' @inheritParams ALDEx2::aldex
#'
#' @return A list object containing the matrix of p-values, the matrix of
#' summary statistics for each tag, and a suggested name of the final object
#' considering the parameters passed to the function.
#'
#' @seealso \code{\link[ALDEx2]{aldex}} for the Dirichlet-Multinomial model
#' estimation. Several and more complex tests are present in the ALDEx2
#' framework.
#'
#' @examples
#' data(ps_stool_16S)
#' # Calculate the scaling factors
#' ps_stool_16S <- norm_edgeR(object = ps_stool_16S, method = "none")
#'
#' # The phyloseq object now contains the scaling factors:
#' scaleFacts <- phyloseq::sample_data(ps_stool_16S)[,"NF.none"]
#' head(scaleFacts)
#'
#' # Differential abundance
#' group <- sample(x = c("grp1","grp2"), size = phyloseq::nsamples(
#'     ps_stool_16S), replace = TRUE)
#' phyloseq::sample_data(ps_stool_16S)$group <- group
#' DA_ALDEx2(ps_stool_16S, conditions = group, test = "t", denom = "iqlr",
#'     norm = "none")

DA_ALDEx2 <- function(object, pseudo_count = FALSE, conditions = NULL,
    mc.samples = 128, test = c("t","wilcox"), denom = "iqlr", norm = c("TMM",
    "TMMwsp", "RLE", "upperquartile", "posupperquartile", "none", "ratio",
    "poscounts", "iterate", "TSS", "CSSmedian", "CSSdefault")){
    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)
    # Name building
    name <- "ALDEx2"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        message("Adding a pseudo count... \n")
        counts <- counts + 1
        name <- paste(name,".pseudo",sep = "")}
    if(length(norm) > 1)
        stop("Please choose one normalization for this istance of differential
            abundance analysis.")
    NF.col <- paste("NF", norm, sep = ".")
    # Check if the column with the normalization factors is present
    if(!any(colnames(metadata) == NF.col)){
        stop(paste0("Can't find the ", NF.col," column in your object. Make
        sure to add the normalization factors column in your object first."
        ))}
    name <- paste(name, ".", norm, sep = "")
    NFs = unlist(metadata[, NF.col])
    # Check if the NFs are scaling factors. If so, make them norm. factors
    if(is.element(norm, c("TMM", "TMMwsp", "RLE", "upperquartile",
        "posupperquartile", "CSSmedian", "CSSdefault", "TSS")))
        NFs <- NFs * colSums(counts)
    NFs <- NFs/exp(mean(log(NFs)))
    norm_counts <- round(counts %*% diag(1/NFs), digits = 0)
    colnames(norm_counts) <- colnames(counts)
    if(is.null(conditions))
        stop("Please supply the name of the variable of interest or the
            entire character vector.")
    else if(length(conditions) == 1)
        conditions = unlist(metadata[, conditions])
    name <- paste(name, ".", denom, sep = "")
    if(!is.element(test, c("t","wilcox")) | length(test) != 1)
        stop("Please choose between p-values produced by Welch t-test (t) or by
            the Wilcoxon test (wilcox).")
    name <- paste(name, ".", test, sep = "")
    statInfo <- ALDEx2::aldex(reads = norm_counts, conditions = conditions,
        mc.samples = mc.samples, test = test, effect = TRUE,
        include.sample.summary = FALSE, denom = denom, verbose = TRUE)
    if(test == "t")
        pValMat <- data.frame(statInfo[, c("we.ep", "we.eBH")])
    else pValMat <- data.frame(statInfo[, c("wi.ep", "wi.eBH")])
    colnames(pValMat) <- c("rawP", "adjP")
    return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
}# END - function: DA_ALDEx2

#' @title DA_corncob
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
#' @param formula an object of class \code{formula} without the response: a
#' symbolic description of the model to be fitted to the abundance.
#' @param phi.formula an object of class \code{formula} without the response: a
#' symbolic description of the model to be fitted to the dispersion.
#' @param coefficient The coefficient of interest as a single word formed by the
#' variable name and the non reference level. (e.g.: 'ConditionDisease' if the
#' reference level for the variable 'Condition' is 'control').
#' @param norm name of the normalization method used to compute the
#' normalization factors to use in the differential abundance analysis.
#' @param test Character. Hypothesis testing procedure to use. One of
#' \code{"Wald"} or \code{"LRT"} (likelihood ratio test).
#' @param boot Boolean. Defaults to \code{FALSE}. Indicator of whether or not to
#' use parametric bootstrap algorithm. (See \code{\link[corncob]{pbWald}} and
#' \code{\link[corncob]{pbLRT}}).
#' @inheritParams corncob::differentialTest
#'
#' @return A list object containing the matrix of p-values, the matrix of
#' summary statistics for each tag, and a suggested name of the final object
#' considering the parameters passed to the function.
#'
#' @seealso \code{\link[corncob]{bbdml}} and
#' \code{\link[corncob]{differentialTest}} for differential abundance and
#' differential variance evaluation.
#'
#' @examples
#' library(phyloseq)
#' data(ps_stool_16S)
#' # Calculate the scaling factors
#' ps_stool_16S <- norm_edgeR(object = ps_stool_16S, method = "none")
#'
#' # The phyloseq object now contains the scaling factors:
#' scaleFacts <- phyloseq::sample_data(ps_stool_16S)[, "NF.none"]
#' head(scaleFacts)
#'
#' # Differential abundance
#' group <- sample(x = c("grp1","grp2"), size = phyloseq::nsamples(
#'     ps_stool_16S), replace = TRUE)
#' phyloseq::sample_data(ps_stool_16S)$group <- group
#' DA_corncob(ps_stool_16S, formula = ~ group, phi.formula = ~ group,
#'     formula_null = ~ 1, phi.formula_null = ~ 1, coefficient = "groupgrp2",
#'     norm = "none", test = "Wald")

DA_corncob <- function(object, pseudo_count = FALSE, formula, phi.formula,
    formula_null, phi.formula_null, test, boot = FALSE, coefficient = NULL,
    norm = c("TMM", "TMMwsp", "RLE", "upperquartile", "posupperquartile",
    "none", "ratio", "poscounts", "iterate", "TSS", "CSSmedian", "CSSdefault")){
    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)
    # Name building
    name <- "corncob"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        message("Adding a pseudo count... \n")
        counts <- counts + 1
        name <- paste(name,".pseudo",sep = "")}
    if(length(norm) > 1)
        stop("Please choose one normalization for this istance of differential
            abundance analysis.")
    NF.col <- paste("NF", norm, sep = ".")
    # Check if the column with the normalization factors is present
    if(!any(colnames(metadata) == NF.col))
        stop(paste0("Can't find the ", NF.col," column in your object. Make
            sure to add the normalization factors column in your object first."
            ))
    name <- paste(name, ".", norm, sep = "")
    NFs = unlist(metadata[, NF.col])
    # Check if the NFs are scaling factors. If so, make them norm. factors
    if(is.element(norm, c("TMM", "TMMwsp", "RLE", "upperquartile",
        "posupperquartile", "CSSmedian", "CSSdefault", "TSS")))
        NFs <- NFs * colSums(counts)
    NFs <- NFs/exp(mean(log(NFs)))
    norm_counts <- round(counts %*% diag(1/NFs), digits = 0)
    colnames(norm_counts) <- colnames(counts)
    message(paste0("Differential abundance on ", norm," normalized data"))
    if(missing(test))
        stop("Please supply the test to perform, 'Wald' or 'LRT'.")
    else name <- paste(name, ".", test, sep = "")
    if(boot)
        name <- paste(name, ".", "boot", sep = "")
    ## differential expression
    requireNamespace("phyloseq")
    fit <- corncob::differentialTest(formula = formula, phi.formula =
        phi.formula, formula_null = formula_null, phi.formula_null =
        phi.formula_null, data = norm_counts, sample_data = metadata, test =
        test, boot = boot)
    # differentialTest's output has a complex structure:
    # extraction of summary table for each estimated model
    # mu.(Intercept), mu.condition, and phi.(Intercept), phi.condition
    # only the mu.condition is of interest.
    if(is.null(coefficient) | !is.element(coefficient,fit[["restrictions_DA"]]))
        stop("Please supply the coefficient of interest as a single word formed
            by the variable name and the non reference level. (e.g.:
            'ConditionDisease' if the reference level for the variable
            'Condition' is 'control')")
    statInfo <- plyr::ldply(.data = fit[["all_models"]], .fun = function(model){
        stats::coef(model)[paste0("mu.",coefficient),]})
    if(length(fit[["restrictions_DA"]])>0)
        message(paste("Differential abundance across", paste0(
            fit[["restrictions_DA"]], collapse = " and ")))
    if(length(fit[["restrictions_DV"]])>0)
        message(paste("Differential variability across", paste0(
            fit[["restrictions_DV"]], collapse = " and ")))
    pValMat <- data.frame("rawP" = fit[["p"]], "adjP" = fit[["p_fdr"]])
    rownames(statInfo) <- rownames(pValMat) <- names(fit[["p"]])
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
#' @seealso \code{\link[MAST]{zlm}} for the Truncated Gaussian Hurdle model
#' estimation.
#'
#' @examples
#' data(ps_stool_16S)
#' # Calculate the scaling factors
#' ps_stool_16S <- norm_edgeR(object = ps_stool_16S, method = "none")
#'
#' # The phyloseq object now contains the scaling factors:
#' scaleFacts <- phyloseq::sample_data(ps_stool_16S)[, "NF.none"]
#' head(scaleFacts)
#'
#' # Differential abundance
#' group <- sample(x = c("grp1","grp2"), size = phyloseq::nsamples(
#'     ps_stool_16S), replace = TRUE)
#' phyloseq::sample_data(ps_stool_16S)$group <- group
#' DA_MAST(ps_stool_16S, rescale = "median", design = ~ group, norm = "none",
#'     coefficient = "groupgrp2")

DA_MAST <- function(object, pseudo_count = FALSE, rescale = c("median",
    "default"), design, coefficient = NULL, norm = c("TMM", "TMMwsp", "RLE",
    "upperquartile", "posupperquartile", "none", "ratio", "poscounts",
    "iterate", "TSS", "CSSmedian", "CSSdefault")){
    # Check the orientation
    if (!phyloseq::taxa_are_rows(object))
        object <- t(object)
    # Slot extraction of phyloseq object
    counts <- as(phyloseq::otu_table(object), "matrix")
    metadata <- phyloseq::sample_data(object)
    # Name building
    name <- "MAST"
    # add 1 if any zero counts
    if (any(counts == 0) & pseudo_count){
        message("Adding a pseudo count... \n")
        counts <- counts + 1
        name <- paste(name,".pseudo",sep = "")}
    if(length(norm) > 1)
        stop("Please choose one normalization for this istance of differential
            abundance analysis.")
    NF.col <- paste("NF", norm, sep = ".")
    # Check if the column with the normalization factors is present
    if(!any(colnames(metadata) == NF.col))
        stop(paste0("Can't find the ", NF.col," column in your object. Make
            sure to add the normalization factors column in your object first."
            ))
    name <- paste(name, ".", norm, sep = "")
    NFs = unlist(metadata[, NF.col])
    # Check if the NFs are scaling factors. If so, make them norm. factors
    if(is.element(norm, c("TMM", "TMMwsp", "RLE", "upperquartile",
        "posupperquartile", "CSSmedian", "CSSdefault", "TSS")))
        NFs <- NFs * colSums(counts)
    NFs <- NFs/exp(mean(log(NFs)))
    norm_counts <- round(counts %*% diag(1/NFs), digits = 0)
    colnames(norm_counts) <- colnames(counts)
    message(paste0("Differential abundance on ", norm, " normalized data"))
    if(length(rescale) > 1 | !is.element(rescale,c("default","median")))
        stop("Please choose between 'default' or 'median' for the rescale
            parameter. 'median' is suggested for metagenomics data.")
    else if(rescale == "median"){
        message(paste0("per ",rescale,"-lib.size rescaled data"))
        tpm <- norm_counts * median(colSums(norm_counts)) / colSums(norm_counts)
    } else {
        message(paste0(rescale," (per million) rescaled data"))
        tpm <- norm_counts * 10^6 / colSums(norm_counts)}
    name <- paste(name, ".", rescale, sep = "")
    tpm <- log2(tpm + 1)
    sca <- MAST::FromMatrix(exprsArray = tpm, cData = metadata)
    # here, we keep all OTUs so that we can fairly compare MAST and the other
    # methods. So, no adaptive thresholding or filtering by gene expression.
    SummarizedExperiment::assays(sca) <- list(tpm =
        SummarizedExperiment::assay(sca))
    ngeneson <- apply(norm_counts, 2, function(x) mean(x>0))
    metadata[, "cngeneson"] <- ngeneson - mean(ngeneson)
    SummarizedExperiment::colData(sca)[, "cngeneson"] <- metadata[, "cngeneson"]
    if(is.character(design)){
        if(grepl("~", design))
            design <- as.formula(design)
        else design <- as.formula(paste0("~", design))}
    if(is(design, "formula"))
        design <- as.formula(paste0(paste0(design, collapse = " "),
            " + cngeneson"))
    ## differential expression
    fit <- MAST::zlm(formula = design, sca = sca, method = "bayesglm", ebayes =
        TRUE)
    if(!is.element(coefficient, colnames(coef(fit, "C"))))
        stop("Please supply the coefficient of interest as a single word formed
            by the variable name and the non reference level. (e.g.:
            'ConditionDisease' if the reference level for the variable
            'Condition' is 'control')")
    summaryDt = data.frame(MAST::summary(fit, doLRT = coefficient)[[
        "datatable"]])
    contrast <- component <- NULL
    fcHurdle <- merge(x = summaryDt[summaryDt[,"contrast"] == coefficient &
        summaryDt[,"component"] == 'H', c("primerid", "Pr..Chisq.")], y =
        summaryDt[summaryDt[,"contrast"] == coefficient & summaryDt[,
        "component"] == "logFC", c("primerid", "coef", "ci.hi", "ci.lo")],
        by = "primerid")
    statInfo = data.frame(logFC = fcHurdle[, "coef"], logFC.lo = fcHurdle[,
        "ci.lo"], logFC.hi = fcHurdle[, "ci.hi"], rawP = fcHurdle[,
        "Pr..Chisq."], adjP = stats::p.adjust(fcHurdle[, "Pr..Chisq."], 'BH'))
    rownames(statInfo) <- fcHurdle[, "primerid"]
    pValMat <- statInfo[, c("rawP", "adjP")]
    statInfo <- statInfo[rownames(counts),]
    pValMat <- pValMat[rownames(counts),]
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
#' @param contrast character vector with exactly three elements: the name of a
#' factor in the design formula, the name of the numerator level for the fold
#' change, and the name of the denominator level for the fold change.
#' @param norm name of the normalization method used to compute the
#' normalization factors to use in the differential abundance analysis.
#'
#' @return A list object containing the matrix of p-values, the matrix of
#' summary statistics for each tag, and a suggested name of the final object
#' considering the parameters passed to the function.
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
#' data(ps_stool_16S)
#' # Calculate the scaling factors
#' ps_stool_16S <- norm_edgeR(object = ps_stool_16S, method = "none")
#'
#' # The phyloseq object now contains the scaling factors:
#' scaleFacts <- phyloseq::sample_data(ps_stool_16S)[,"NF.none"]
#' head(scaleFacts)
#'
#' # Differential abundance
#' group <- sample(x = c("grp1","grp2"), size = phyloseq::nsamples(
#'     ps_stool_16S), replace = TRUE)
#' phyloseq::sample_data(ps_stool_16S)$group <- as.factor(group)
#' DA_Seurat(ps_stool_16S, contrast = c("group","grp2","grp1"), norm = "none")

DA_Seurat <- function(object, pseudo_count = FALSE, test.use = "wilcox",
    contrast, norm = c("TMM", "TMMwsp", "RLE", "upperquartile",
    "posupperquartile", "none", "ratio", "poscounts", "iterate", "TSS",
    "CSSmedian", "CSSdefault")){
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
        message("Adding a pseudo count... \n")
        counts <- counts + 1
        name <- paste(name,".pseudo",sep = "")}
    if(length(norm) > 1)
        stop("Please choose one normalization for this istance of differential
            abundance analysis.")
    NF.col <- paste("NF", norm, sep = ".")
    # Check if the column with the normalization factors is present
    if(!any(colnames(metadata) == NF.col)){
        stop(paste0("Can't find the ", NF.col," column in your object. Make
            sure to add the normalization factors column in your object first."
            ))}
    name <- paste(name, ".", norm, sep = "")
    NFs = unlist(metadata[,NF.col])
    # Check if the NFs are scaling factors. If so, make them norm. factors
    if(is.element(norm, c("TMM", "TMMwsp", "RLE", "upperquartile",
        "posupperquartile", "CSSmedian", "CSSdefault", "TSS")))
        NFs <- NFs * colSums(counts)
    NFs <- NFs/exp(mean(log(NFs)))
    norm_counts <- round(counts %*% diag(1/NFs), digits = 0)
    colnames(norm_counts) <- colnames(counts)
    message(paste0("Differential abundance on ", norm," normalized data"))
    # Initialize the Seurat object with the raw (non-normalized data).
    # If the chosen normalization is 'none', then this is true.
    # Keep all features expressed in >= 1 sample
    # Keep all samples with at least 1 detected feature.
    sobj <- Seurat::CreateSeuratObject(counts = norm_counts, min.cells = 1,
        min.features = 1)
    sobj <- Seurat::AddMetaData(object = sobj, metadata = data.frame(metadata),
        col.name = colnames(metadata))
    if(missing(contrast) | (!is.character(contrast) & length(contrast) != 3))
        stop(paste0("Please supply a character vector with exactly three
            elements: the name of a factor in the design formula, the name of
            the numerator level for the fold change, and the name of the
            denominator level for the fold change."))
    else {
        if(!is(unlist(sobj[[contrast[1]]]),"factor"))
            stop(paste(contrast[1]," variable is not a factor. Please supply a
                factor."))
        else{
            if(!is.element(contrast[2],levels(unlist(sobj[[contrast[1]]]))))
                stop(paste(contrast[2], "is not a level of the", contrast[1],
                    "variable. Please supply a present category."))
            if(!is.element(contrast[3],levels(unlist(sobj[[contrast[1]]]))))
                stop(paste(contrast[3], "is not a level of the", contrast[1],
                    "variable. Please supply a present category."))
        }
    }
    message(paste0("Extracting results for ", contrast[1]," variable, ",
        contrast[2], " / ", contrast[3]))
    sobj <- Seurat::NormalizeData(object = sobj, normalization.method =
        "LogNormalize", scale.factor = 10000)
    sobj <- Seurat::FindVariableFeatures(object = sobj, nfeatures = round(nrow(
        counts)*0.1, digits = 0))
    sobj <- Seurat::ScaleData(object = sobj, vars.to.regress = c("nCount_RNA"))
    statInfo_ <- Seurat::FindMarkers(sobj, test.use = test.use, group.by =
        contrast[1], ident.1 = contrast[2], ident.2 = contrast[3],
        logfc.threshold = 0, min.cpt = 0)
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
