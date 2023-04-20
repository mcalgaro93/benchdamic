#' @title fitNB
#'
#' @importFrom edgeR calcNormFactors estimateDisp DGEList glmFit getCounts
#' @importFrom edgeR getDispersion
#' @importFrom phyloseq phyloseq otu_table taxa_are_rows
#' @export
#' @description
#' Fit a Negative Binomial (NB) distribution for each taxon of the count data.
#' The NB estimation procedure is performed by edgeR \code{\link{glmFit}}
#' function, using \code{TMM} normalized counts, tag-wise dispersion estimation,
#' and not assuming the presence of any group in the samples (design matrix
#' equal to a column of ones).
#' 
#' @param object a phyloseq object, a TreeSummarizedExperiment object, or a 
#' matrix of counts.
#' @inheritParams get_counts_metadata
#' @param verbose an optional logical value. If \code{TRUE} information on the
#' steps of the algorithm is printed. Default \code{verbose = TRUE}.
#'
#' @return A data frame containing the continuity corrected logarithms of the
#' average fitted values for each row of the `counts` matrix in the `Y` column,
#' and the estimated probability to observe a zero in the `Y0` column.
#' @examples
#' # Generate some random counts
#' counts = matrix(rnbinom(n = 60, size = 3, prob = 0.5), nrow = 10, ncol = 6)
#'
#' # Fit model on the matrix of counts
#' NB <- fitNB(counts)
#' head(NB)
fitNB <- function(object, assay_name = "counts", verbose = TRUE) {
    if(is(object, "phyloseq") | is(object, "TreeSummarizedExperiment")){
        counts_and_metadata <- get_counts_metadata(object, 
            assay_name = assay_name)
        counts <- counts_and_metadata[[1]]
    } else if (is.matrix(object)) {
        NULL
    } else {
        stop("Please supply a phyloseq object, a TreeSummarizedExperiment,", 
            " or a matrix of counts.")
    }
    if(verbose)
        message("Model: Negative Binomial")
    # Default normalization
    normFacts <- edgeR::calcNormFactors(counts)
    # DGEList object creation
    dge <- edgeR::DGEList(counts = counts, norm.factors = normFacts)
    # Dispersion estimate
    if(verbose){
        message("Estimating dispersions")
        disp <- edgeR::estimateDisp(y = dge, tagwise = TRUE)
    } else {
        disp <- suppressMessages(edgeR::estimateDisp(y = dge, tagwise = TRUE))
    }
    # GLM
    disps <- edgeR::getDispersion(disp)
    fit <- edgeR::glmFit(edgeR::getCounts(dge), dispersion = disps)
    # Fitted values extraction
    fitVals <- fit[["fitted.values"]]
    Y <- log1p(rowMeans(fitVals))
    Y0 <- rowMeans((1 + fitVals * disps)^(-1 / disps))
    return(data.frame("Y" = Y, "Y0" = Y0))
}

#' @title fitZINB
#'
#' @importFrom zinbwave zinbFit getMu getPi getPhi
#' @importFrom BiocParallel SerialParam
#' @importFrom phyloseq otu_table taxa_are_rows
#' @export
#' @description
#' Fit a Zero-Inflated Negative Binomial (ZINB) distribution for each taxon of
#' the countdata. The ZINB estimation procedure is performed by zinbwave
#' \code{\link{zinbFit}} function with \code{commondispersion = FALSE},
#' regularization parameter \code{epsilon = 1e10}, and not assuming the presence
#' of any group in the samples (design matrix equal to a column of ones.)
#'
#' @inheritParams fitNB
#'
#' @return A data frame containing the continuity corrected logarithms of the
#' average fitted values for each row of the matrix of counts in the \code{Y}
#' column, and the estimated probability to observe a zero in the \code{Y0}
#' column.
#' @examples
#' # Generate some random counts
#' counts = matrix(rnbinom(n = 60, size = 3, prob = 0.5), nrow = 10, ncol = 6)
#'
#' # Fit model on the counts matrix
#' ZINB <- fitZINB(counts)
#' head(ZINB)
fitZINB <- function(object, assay_name = "counts", verbose = TRUE) {
    if(is(object, "phyloseq") | is(object, "TreeSummarizedExperiment")){
        counts_and_metadata <- get_counts_metadata(object, 
            assay_name = assay_name)
        counts <- counts_and_metadata[[1]]
    } else if (is.matrix(object)) {
        NULL
    } else {
        stop("Please supply a phyloseq object, a TreeSummarizedExperiment,", 
             " or a matrix of counts.")
    }
    if(verbose)
        message("Model: Zero-Inflated Negative Binomial")
    fit <- zinbwave::zinbFit(
        Y = counts, epsilon = 1e10, commondispersion = TRUE,
        BPPARAM = BiocParallel::SerialParam(), verbose = verbose
    )
    mu <- t(zinbwave::getMu(fit))
    pi <- t(zinbwave::getPi(fit))
    phi <- zinbwave::getPhi(fit)
    Y <- log1p(rowMeans((1 - pi) * mu))
    Y0 <- rowMeans(pi + (1 - pi) * (1 + phi * mu)^(-1 / phi))
    return(data.frame("Y" = Y, "Y0" = Y0))
}

#' @title fitHURDLE
#'
#' @importFrom stats model.matrix as.formula median
#' @importFrom MAST FromMatrix zlm invlogit
#' @importFrom SummarizedExperiment assay colData
#' @importFrom phyloseq otu_table taxa_are_rows
#' @export
#' @description
#' Fit a truncated gaussian hurdle model for each taxon of the count data. The
#' hurdle model estimation procedure is performed by MAST \code{\link{zlm}}
#' function without assuming the presence of any group in the samples (design
#' matrix equal to a column of ones.)
#'
#' @inheritParams fitNB
#' @param scale Character vector, either \code{median} or \code{default} to
#' choose between the median of the library size or one million to scale raw
#' counts.
#'
#' @return A data frame containing the continuity corrected logarithms of the
#' average fitted values for each row of the matrix of counts in the \code{Y}
#' column, and the estimated probability to observe a zero in the \code{Y0}
#' column.
#' @examples
#' # Generate some random counts
#' counts = matrix(rnbinom(n = 600, size = 3, prob = 0.5), nrow = 100, ncol = 6)
#'
#' # Fit model on the counts matrix
#' HURDLE <- fitHURDLE(counts, scale = "median")
#' head(HURDLE)
fitHURDLE <- function(object, assay_name = "counts", scale = "default", 
    verbose = TRUE) {
    if(is(object, "phyloseq") | is(object, "TreeSummarizedExperiment")){
        counts_and_metadata <- get_counts_metadata(object, 
            assay_name = assay_name)
        counts <- counts_and_metadata[[1]]
    } else if (is.matrix(object)) {
        NULL
    } else {
        stop("Please supply a phyloseq object, a TreeSummarizedExperiment,", 
             " or a matrix of counts.")
    }
    if(verbose)
        message("Model: Truncated Gaussian Hurdle")
    # tpm scaling
    if (scale == "median") {
        tpm <- counts * stats::median(colSums(counts)) / colSums(counts)
    } else if (scale == "default") {
        tpm <- counts * 1e6 / colSums(counts)
    } else {
        stop("'scale' must be 'median' or 'default'")
    }
    # log2 + 1 transformation
    tpm <- log2(tpm + 1)
    # Single Cell object
    if(verbose){
        message("The counts provided have been rescaled and log2 transformed.")
        message("Making a SingleCellExperiment:")
        sca <- MAST::FromMatrix(exprsArray = tpm)
    } else {
        sca <- suppressMessages(MAST::FromMatrix(tpm))
    }
    # Normalization suggested in MAST vignette
    ngeneson <- colSums(SummarizedExperiment::assay(sca))
    CD <- SummarizedExperiment::colData(sca)
    CD[, "ngeneson"] <- ngeneson
    CD[, "cngeneson"] <- scale(ngeneson)
    SummarizedExperiment::colData(sca) <- CD
    # hurdle model estimation
    hm <- MAST::zlm(~ 1 + cngeneson, sca = sca, onlyCoef = TRUE)
    # C = Continous, positive counts, D = Discrete, 0 or 1
    betaC <- hm[, , "C"]
    fittedC <- tcrossprod(betaC, stats::model.matrix(~ 1 + cngeneson,
        data = SummarizedExperiment::colData(sca)
    ))
    betaD <- hm[, , "D"]
    fittedD <- tcrossprod(betaD, stats::model.matrix(~ 1 + cngeneson,
        data = SummarizedExperiment::colData(sca)
    ))
    # Estimated parameters
    mu <- fittedC
    pi <- MAST::invlogit(fittedD)
    # Back-transformation
    Y <- log(rowMeans(exp((pi * mu) * log(2))))
    Y0 <- rowMeans(1 - pi)
    return(data.frame("Y" = Y, "Y0" = Y0))
}

#' @title fitZIG
#'
#' @importFrom metagenomeSeq newMRexperiment cumNormStat cumNorm normFactors
#'     fitZig zigControl
#' @importFrom utils capture.output
#' @importFrom stats median
#' @importFrom phyloseq otu_table taxa_are_rows
#' @export
#' @description
#' Fit a Zero-Inflated Gaussian (ZIG) distribution for each taxon of the count
#' data. The model estimation procedure is performed by metagenomeSeq
#' \code{\link{fitZig}} function without assuming the presence of any group in
#' the samples (design matrix equal to a column of ones.)
#'
#' @inheritParams fitNB
#'
#' @return A data frame containing the continuity corrected logarithms of the
#' average fitted values for each row of the matrix of counts in the \code{Y}
#' column, and the estimated probability to observe a zero in the \code{Y0}
#' column.
#' @examples
#' # Generate some random counts
#' counts = matrix(rnbinom(n = 60, size = 3, prob = 0.5), nrow = 10, ncol = 6)
#'
#' # Fit model on the counts matrix
#' ZIG <- fitZIG(counts)
#' head(ZIG)
fitZIG <- function(object, assay_name = "counts", verbose = TRUE) {
    if(is(object, "phyloseq") | is(object, "TreeSummarizedExperiment")){
        counts_and_metadata <- get_counts_metadata(object,
            assay_name = assay_name)
        counts <- counts_and_metadata[[1]]
    } else if (is.matrix(object)) {
        NULL
    } else {
        stop("Please supply a phyloseq object, a TreeSummarizedExperiment,", 
             " or a matrix of counts.")
    }
    if(verbose)
        message("Model: Zero-Inflated Gaussian")
    MGS <- metagenomeSeq::newMRexperiment(counts = counts)
    # Normalization
    if(verbose){
        MGSp <- metagenomeSeq::cumNormStat(MGS)
    } else {
        MGSp <- suppressMessages(metagenomeSeq::cumNormStat(MGS))
    }
    MGS <- metagenomeSeq::cumNorm(MGS, MGSp)
    normFactor <- metagenomeSeq::normFactors(MGS)
    # scaling
    normFactor <- log2(normFactor / stats::median(normFactor) + 1)
    # Design matrix
    desMat <- cbind(1, normFactor = normFactor)
    # Estimation
    control <- metagenomeSeq::zigControl(maxit = 1000, verbose = verbose)
    zig <- metagenomeSeq::fitZig(MGS, desMat, control = control,
        useCSSoffset = FALSE)
    # Coefficient extraction (metagenomeSeq::MRcoefs() changes the order, use @)
    mu <- tcrossprod(coef(zig@fit), desMat)
    Y <- rowMeans(mu) * log(2)
    Y0 <- rowMeans(zig@z)
    return(data.frame("Y" = Y, "Y0" = Y0))
}

#' @title fitDM
#'
#' @importFrom MGLM MGLMreg
#' @importFrom stats as.formula model.matrix
#' @importFrom stats4 coef
#' @importFrom phyloseq otu_table taxa_are_rows
#' @import methods
#'
#' @export
#'
#' @description
#' Fit a Dirichlet-Multinomial (DM) distribution for each taxon of the count
#' data. The model estimation procedure is performed by MGLM
#' \code{\link{MGLMreg}} function without assuming the presence of any group in
#' the samples (design matrix equal to a column of ones.)
#' @inheritParams fitNB
#'
#' @return A data frame containing the continuity corrected logarithms of the
#' average fitted values for each row of the matrix of counts in the \code{Y}
#' column, and the estimated probability to observe a zero in the \code{Y0}
#' column.
#' @examples
#' # Generate some random counts
#' counts = matrix(rnbinom(n = 60, size = 3, prob = 0.5), nrow = 10, ncol = 6)
#'
#' # Fit model on the counts matrix
#' DM <- fitDM(counts)
#' head(DM)
fitDM <- function(object, assay_name = "counts", verbose = TRUE) {
    if(is(object, "phyloseq") | is(object, "TreeSummarizedExperiment")){
        counts_and_metadata <- get_counts_metadata(object,
            assay_name = assay_name)
        counts <- counts_and_metadata[[1]]
    } else if (is.matrix(object)) {
        NULL
    } else {
        stop("Please supply a phyloseq object, a TreeSummarizedExperiment,", 
             " or a matrix of counts.")
    }
    if(verbose)
        message("Model: Dirichlet Multinomial")
    # library sizes
    ls <- colSums(counts)
    data <- t(counts)
    # Design only intercept
    X <- stats::model.matrix(data ~ 1)
    # Model fit
    if(verbose){
        dmFit <- MGLM::MGLMreg.fit(Y = data, X = X, dist = "DM", 
            display = verbose)
    } else {
        dmFit <- suppressWarnings(MGLM::MGLMreg.fit(Y = data, X = X, 
            dist = "DM", display = verbose))
    }
    # fitted_values <- dmFit@fitted * ls
    # Coefficent extraction
    alpha_i <- exp(stats4::coef(dmFit))
    alpha_0 <- rowSums(alpha_i)
    Y <- log1p(colMeans(ls %*% alpha_i / alpha_0))
    # The univariate version of a DM is a Beta-Binomial
    Y0 <- rowMeans(vapply(
        X = ls, FUN = function(l) {
            beta(
                a = alpha_i,
                b = l + alpha_0 - alpha_i
            ) / beta(a = alpha_i, b = alpha_0 - alpha_i)
        },
        FUN.VALUE = seq(0, 1, length.out = nrow(counts))
    ))
    return(data.frame("Y" = Y, "Y0" = Y0))
}
