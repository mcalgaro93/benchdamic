#' @title fitNB
#'
#' @importFrom edgeR calcNormFactors estimateDisp DGEList glmFit getCounts
#' @importFrom edgeR getDispersion
#' @export
#' @description
#' Fit a Negative Binomial (NB) distribution for each taxon of the count data.
#' The NB estimation procedure is performed by edgeR \code{\link{glmFit}}
#' function, using \code{TMM} normalized counts, tag-wise dispersion estimation,
#' and not assuming the presence of any group in the samples (design matrix
#' equal to a column of ones.)
#' @param counts a matrix of counts with features (OTUs, ASVs, genes) by row and
#' samples by column.
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
fitNB <- function(counts, verbose = TRUE) {
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
fitZINB <- function(counts, verbose = TRUE) {
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
#' counts = matrix(rnbinom(n = 60, size = 3, prob = 0.5), nrow = 10, ncol = 6)
#'
#' # Fit model on the counts matrix
#' HURDLE <- fitHURDLE(counts, scale = "median")
#' head(HURDLE)
fitHURDLE <- function(counts, scale = "default", verbose = TRUE) {
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
        sca <- MAST::FromMatrix(tpm)
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
#' @export
#' @description
#' Fit a Zero-Inflated Gaussian (ZIG) distribution for each taxon of the count
#' data. The model estimation procedure is performed by metagenomeSeq
#' \code{\link{fitZig}} function without assuming the presence of any group in
#' the samples (design matrix equal to a column of ones.)
#'
#' @inheritParams fitNB
#' @param scale Character vector, either \code{median} or \code{default} to
#' choose between the median of the library size or one thousand to scale
#' normalization factors.
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
#' ZIG <- fitZIG(counts, scale = "median")
#' head(ZIG)
fitZIG <- function(counts, scale = "default", verbose = TRUE) {
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
    if (scale == "median") {
        normFactor <- log2(normFactor / stats::median(normFactor) + 1)
    } else if (scale == "default") {
        normFactor <- log2(normFactor / 1000 + 1)
    } else {
        stop("'scale' must be 'median' or 'default'")
    }
    # Design matrix
    desMat <- cbind(1, normFactor = normFactor)
    # Estimation
    control = metagenomeSeq::zigControl(maxit = 1000, verbose = verbose)
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
#' @importFrom stats as.formula
#' @importFrom stats4 coef
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
fitDM <- function(counts, verbose = TRUE) {
    if(verbose)
        message("Model: Dirichlet Multinomial")
    # library sizes
    ls <- colSums(counts)
    data <- t(counts)
    # Design only intercept
    desFormula <- stats::as.formula("data ~ 1")
    # Model fit
    dmFit <- MGLM::MGLMreg(data ~ 1L, dist = "DM", display = verbose)
    # fitted_values <- dmFit@fitted * ls
    # Coefficent extraction
    alpha_i <- exp(coef(dmFit))
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
