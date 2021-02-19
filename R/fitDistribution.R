#' @title fitNB
#'
#' @importFrom edgeR calcNormFactors estimateDisp DGEList glmFit
#' @export
#' @description
#' Fit a Negative Binomial (NB) distribution for each taxon of the count data.
#' The NB estimation procedure is performed by edgeR `glmFit` function, using
#' "TMM" normalized counts, tag-wise dispersion estimation, and not assuming the
#' presence of any group in the samples (design matrix equal to a column of
#' ones.)
#' @param counts a matrix of counts with features (OTUs, ASVs, genes) by row and
#' samples by column.
#'
#' @return A data frame containing the continuity corrected logarithms of the
#' average fitted values for each row of the `counts` matrix in the `Y` column,
#' and the estimated probability to observe a zero in the `Y0` column.
#' @examples
#' # Let's generate some random count data
#' nlibs <- 10
#' ntaxa <- 100
#' lambda = 2
#' counts <- matrix(rpois(ntaxa*nlibs, lambda),
#'                  nrow = ntaxa,
#'                  ncol = nlibs)
#' # NB estimation
#' f_values <- fitNB(counts)
#' head(f_values)

# Negative Binomial fitting
fitNB <- function(counts){
  cat("Model: Negative Binomial \n")
  # Default normalization
  normFacts <- edgeR::calcNormFactors(counts)
  # DGEList object creation
  dge <- edgeR::DGEList(counts = counts, norm.factors = normFacts)
  # Dispersion estimate
  disp <- edgeR::estimateDisp(y = dge, tagwise = TRUE)
  # GLM
  fit <- edgeR::glmFit(dge$counts,
                       dispersion = disp$tagwise.dispersion)
  # nbloglik <- rowSums(dnbinom(x = counts,
  #                             size=1/disp$tagwise.dispersion,
  #                             mu=rowMeans(fit$fitted.values),
  #                             log = TRUE))
  # Fitted values extraction
  # Return the log(average fitted values + 1) to
  Y = log1p(rowMeans(fit$fitted.values))
  Y0 = rowMeans((1 + fit$fitted.values * disp$tagwise.dispersion)^(-1/disp$tagwise.dispersion))
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
#' `zinbFit` function with `commondispersion = FALSE`, regularization parameter
#' `epsilon = 1e10`, and not assuming the presence of any group in the samples
#' (design matrix equal to a column of ones.)
#'
#' @inheritParams fitNB
#'
#' @return A data frame containing the continuity corrected logarithms of the
#' average fitted values for each row of the `counts` matrix in the `Y` column,
#' and the estimated probability to observe a zero in the `Y0` column.
#' @examples
#' # Let's generate some random count data
#' nlibs <- 10
#' ntaxa <- 100
#' lambda = 2
#' counts <- matrix(rpois(ntaxa*nlibs, lambda),
#'                  nrow = ntaxa,
#'                  ncol = nlibs)
#' # ZINB estimation
#' f_values <- fitZINB(counts)
#' head(f_values)

fitZINB <- function(counts)
{
  cat("Model: Zero-Inflated Negative Binomial \n")
  fit <- zinbwave::zinbFit(Y = counts,
                           epsilon = 1e10, # Regularization parameter fixed
                           commondispersion = TRUE,
                           BPPARAM = BiocParallel::SerialParam())
  mu = t(zinbwave::getMu(fit))
  pi = t(zinbwave::getPi(fit))
  phi = zinbwave::getPhi(fit)


  Y = log1p(rowMeans((1 - pi) * mu))
  Y0 = rowMeans(pi + (1 - pi) * (1 + phi * mu) ^ (-1/phi))
  return(data.frame("Y" = Y, "Y0" = Y0))
}

#' @title fitHURDLE
#'
#' @importFrom stats model.matrix as.formula median coef
#' @importFrom MAST FromMatrix zlm invlogit
#' @export
#' @description
#' Fit a truncated gaussian hurdle model for each taxon of the count data. The
#' hurdle model estimation procedure is performed by MAST `zlm` function without
#' assuming the presence of any group in the samples (design matrix equal to a
#' column of ones.)
#'
#' @inheritParams fitNB
#' @param scale Character vector, either 'median' or 'default' to choose between
#' the median of the library size or one million to scale raw counts.
#'
#' @return A data frame containing the continuity corrected logarithms of the
#' average fitted values for each row of the `counts` matrix in the `Y` column,
#' and the estimated probability to observe a zero in the `Y0` column.
#' @examples
#' # Let's generate some random count data
#' nlibs <- 10
#' ntaxa <- 100
#' lambda = 2
#' counts <- matrix(rpois(ntaxa*nlibs, lambda),
#'                  nrow = ntaxa,
#'                  ncol = nlibs)
#' # HURDLE model estimation
#' f_values <- fitHURDLE(counts,"median")
#' head(f_values)

fitHURDLE <- function(counts, scale = "default"){
  cat("Model: Truncated Gaussian Hurdle \n")
  # tpm scaling
  if(scale == "median"){
    tpm <- counts * stats::median(colSums(counts)) / colSums(counts)
  } else if(scale == "default"){
    tpm <- counts * 1e6 / colSums(counts)
  } else stop("'scale' must be 'median' or 'default'")
  # log2 + 1 transformation
  tpm <- log2(tpm + 1)
  # Single Cell object
  sca <- MAST::FromMatrix(tpm)
  # Normalization suggested in MAST vignette
  ngeneson <- colSums(sca@assays@data$et)
  CD <- sca@colData
  CD$ngeneson <- ngeneson
  CD$cngeneson <- scale(ngeneson)
  sca@colData <- CD
  # hurdle model estimation
  hm <- MAST::zlm(~ 1 + cngeneson, sca = sca)
  # C = Continous, positive counts, D = Discrete, 0 or 1
  betaC <- hm@coefC
  fittedC <- tcrossprod(betaC,
                        stats::model.matrix(~ 1 + cngeneson,
                                            data = sca@colData))
  betaD <- hm@coefD
  fittedD <- tcrossprod(betaD,
                        stats::model.matrix(~ 1 + cngeneson,
                                            data = sca@colData))
  # Estimated parameters
  mu <- fittedC
  pi <- MAST::invlogit(fittedD)
  # Back-transformation
  Y = log(rowMeans(exp((pi * mu) * log(2))))
  Y0 = rowMeans(1 - pi)
  return(data.frame("Y" = Y, "Y0" = Y0))
}

#' @title fitZIG
#'
#' @importFrom metagenomeSeq newMRexperiment cumNormStat cumNorm normFactors
#'     fitZig zigControl
#' @importFrom stats median
#' @export
#' @description
#' Fit a Zero-Inflated Gaussian (ZIG) distribution for each taxon of the count
#' data. The model estimation procedure is performed by metagenomeSeq `fitZig`
#' function without assuming the presence of any group in the samples (design
#' matrix equal to a column of ones.)
#'
#' @inheritParams fitNB
#' @param scale Character vector, either 'median' or 'default' to choose between
#' the median of the library size or one thousand to scale normalization
#' factors.
#'
#' @return A data frame containing the continuity corrected logarithms of the
#' average fitted values for each row of the `counts` matrix in the `Y` column,
#' and the estimated probability to observe a zero in the `Y0` column.
#' @examples
#' # Let's generate some random count data
#' nlibs <- 10
#' ntaxa <- 100
#' lambda = 2
#' counts <- matrix(rpois(ntaxa*nlibs, lambda),
#'                  nrow = ntaxa,
#'                  ncol = nlibs)
#' # ZIG model estimation
#' f_values <- fitZIG(counts,"median")
#' head(f_values)

fitZIG <- function(counts, scale = "default"){
  cat("Model: Zero-Inflated Gaussian \n")
  MGS <- metagenomeSeq::newMRexperiment(counts = counts)
  # Normalization
  MGSp = metagenomeSeq::cumNormStat(MGS)
  MGS <- metagenomeSeq::cumNorm(MGS,MGSp)
  normFactor = metagenomeSeq::normFactors(MGS)
  # scaling
  if(scale == "median"){
    normFactor = log2(normFactor / stats::median(normFactor) + 1)
  } else if(scale == "default"){
    normFactor = log2(normFactor / 1000 + 1)
  } else stop("'scale' must be 'median' or 'default'")
  # Design matrix
  desMat <- cbind(1, normFactor = normFactor)
  # Estimation
  zig <- metagenomeSeq::fitZig(MGS,
                               desMat,
                               control = zigControl(maxit = 1000),
                               useCSSoffset = FALSE) # To allow the "median"
  # Coefficient extraction
  mu <- tcrossprod(zig@fit$coefficients, desMat)
  Y <- rowMeans(mu) * log(2)
  Y0 <- rowMeans(zig@z)
  return(data.frame("Y" = Y, "Y0" = Y0))
}

#' @title fitDM
#'
#' @importFrom MGLM MGLMreg
#' @importFrom stats as.formula
#' @export
#' @description
#' Fit a Dirichlet-Multinomial (DM) distribution for each taxon of the count
#' data. The model estimation procedure is performed by MGLM `MGLMreg` function
#' without assuming the presence of any group in the samples (design matrix
#' equal to a column of ones.)
#'
#' @inheritParams fitNB
#'
#' @return A data frame containing the continuity corrected logarithms of the
#' average fitted values for each row of the `counts` matrix in the `Y` column,
#' and the estimated probability to observe a zero in the `Y0` column.
#' @examples
#' # Let's generate some random count data
#' nlibs <- 10
#' ntaxa <- 100
#' lambda = 2
#' counts <- matrix(rpois(ntaxa*nlibs, lambda),
#'                  nrow = ntaxa,
#'                  ncol = nlibs)
#' # DM model estimation
#' f_values <- fitDM(counts)
#' head(f_values)

fitDM <- function(counts){
  cat("Model: Dirichlet Multinomial \n")
  # library sizes
  ls <- colSums(counts)
  data <- t(counts)
  # Design only intercept
  desFormula = stats::as.formula("data ~ 1")
  # Model fit
  dmFit <- MGLM::MGLMreg(desFormula, dist = "DM", display = TRUE)
  fitted_values <- dmFit@fitted * ls
  # Coefficent extraction
  alpha_i <- t(exp(t(dmFit@coefficients) %*% 1))
  alpha_0 <- rowSums(alpha_i)

  Y <- log1p(colMeans(ls %*% alpha_i/alpha_0))
  # The univariate version of a DM is a Beta-Binomial
  Y0 <- rowMeans(sapply(ls,function(l)
    beta(a = alpha_i,
         b = l + alpha_0 - alpha_i) /
      beta(a = alpha_i,
           b =  alpha_0 - alpha_i)))
  return(data.frame("Y" = Y, "Y0" = Y0))
}

#' @title prepareObserved
#'
#' @importFrom stats median
#' @export
#' @description
#' Average count value continuity corrected logarithms and fraction of zeroes by
#' taxon.
#'
#' @inheritParams fitNB
#' @param scale If specified it refers to the character vector used in
#' `fitHURDLE` function. Either 'median' or 'default' to choose between the
#' median library size or one million as scaling factors for raw counts.
#'
#' @return A data frame containing the continuity corrected logarithm for the
#' raw count mean values for each taxon of the `counts` matrix in the `Y` column
#' and the observed zero rate in the `Y0` column. If `scale` is specified the
#' continuity corrected logarithm for the mean CPM (`scale = "default"`) or the
#' mean counts per median library size (`scale = "median"`) is computed instead.
#'
#' @examples
#' # Let's generate some random count data
#' nlibs <- 10
#' ntaxa <- 100
#' lambda = 2
#' counts <- matrix(rpois(ntaxa*nlibs, lambda),
#'                  nrow = ntaxa,
#'                  ncol = nlibs)
#'
#' # observed values
#' o_values <- prepareObserved(counts)
#' head(o_values)

prepareObserved <- function(counts,
                            scale = NULL){
  if(!is.null(scale)){
    if(scale == "median"){
      counts <- counts * stats::median(colSums(counts)) / colSums(counts)
    } else if(scale == "default"){
      counts <- counts * 1e6 / colSums(counts)
    } else stop("When specified, 'scale' must be 'median' or 'default'")
  }
  Y <- log1p(rowMeans(counts))
  Y0 <- rowMeans(counts == 0)
  return(data.frame("Y" = Y, "Y0" = Y0))
}

#' @title meanDifferences
#'
#' @export
#' @description
#' Compute the differences between the estimated and the observed continuity
#' corrected logarithms of the average count values (MD), and between the
#' estimated average probability to observe a zero and the the observed zero
#' rate (ZPD).
#'
#' @param estimated a two column data.frame, output of `fitNB`, `fitZINB`,
#' `fitDM`, `fitZIG`, or `fitHURDLE` functions. More in general, a data frame
#' containing the continuity corrected logarithm for the average of the fitted
#' values for each row of a matrix of counts in the `Y` column, and the
#' estimated probability to observe a zero in the `Y0` column.
#' @param observed a two column data.frame, output of `prepareObserved`
#' function. More in general, a data frame containing
#' the continuity corrected logarithm for the average of the observed values for
#' each row of a matrix of counts in the `Y` column, and the estimated
#' proportion of zeroes in the `Y0` column.
#'
#' @return `data.frame` containing the differences between the estimated and the
#' observed continuity corrected logarithms of the average count values in the
#' `MD` columns, and between the estimated average probability to observe a zero
#' and the the observed zero rate in the `ZPD` column.

meanDifferences <- function(estimated,
                            observed){
  if(sum(colnames(estimated) != colnames(observed))>0)
    stop("Estimated and Observed data.frames have different colnames")
  if(sum(colnames(estimated) != c("Y","Y0")) > 0 |
     sum(colnames(observed) != c("Y","Y0")) > 0)
    stop("Please rename the colnames of the data.frames as Y and Y0")
  if(nrow(estimated) != nrow(observed))
    stop("Estimated and Observed data.frames have different number of rows")
  MD <- estimated$Y - observed$Y
  ZPD <- estimated$Y0 - observed$Y0
  return(data.frame("MD" = MD, "ZPD" = ZPD))
}

#' @title RMSE
#'
#' @export
#' @description
#' Computes the Root Mean Square Error (RMSE) from a vector of differences.
#'
#' @param differences a vector of differences.
#'
#' @return RMSE value

RMSE <- function(differences){
  sqrt(mean(differences^2,na.rm = TRUE))
}

#' @title fitModels
#'
#' @export
#' @description
#' A wrapper function that fits the specified models for each taxon of the count
#' data and computes the mean difference (MD) and zero probability difference
#' (ZPD) between estimated and observed values.
#'
#' @inheritParams fitNB
#' @param models character vector which assumes the values 'NB', 'ZINB', 'DM',
#' 'ZIG', and 'HURDLE'.
#' @param scale_ZIG character vector, either 'median' or 'default' to choose
#' between the median of the library size or one thousand to scale normalization
#' factors for the zero-inflated gaussian model.
#' @param scale_HURDLE character vector, either 'median' or 'default' to choose
#' between the median of the library size or one million to scale raw counts for
#' the truncated gaussian hurdle model.
#'
#' @return list of `data.frame` objects for each model. The first two columns
#' contain the properly transformed observed values for mean and zero
#' proportion, while the third and the fourth columns contain the estimated
#' values for the mean and the zero rate resplectively.
#'
#' @seealso \code{\link{fitNB}}, \code{\link{fitZINB}}, \code{\link{fitDM}},
#' \code{\link{fitZIG}}, and \code{\link{fitHURDLE}} for the model estimations,
#' \code{\link{prepareObserved}} for raw counts preparation, and
#' \code{\link{meanDifferences}} for the Mean Difference (MD) and Zero
#' Probability Difference (ZPD) computations.

fitModels <- function(counts,
                      models = c("NB","ZINB","DM","ZIG","HURDLE"),
                      scale_ZIG = c("default","median"),
                      scale_HURDLE = c("default","median")){
  fittedModels <- list()
  observed <- prepareObserved(counts)
  if("NB" %in% models){
    fitted <- fitNB(counts)
    MD <- meanDifferences(estimated = fitted,
                          observed = observed)
    fittedModels$NB <- data.frame(observed,MD)
  }
  if("ZINB" %in% models){
    fitted <- fitZINB(counts)
    MD <- meanDifferences(estimated = fitted,
                          observed = observed)
    fittedModels$ZINB <- data.frame(observed,MD)
  }
  if("DM" %in% models){
    fitted <- fitDM(counts)
    MD <- meanDifferences(estimated = fitted,
                          observed = observed)
    fittedModels$DM <- data.frame(observed,MD)
  }
  if("ZIG" %in% models){
    fitted <- fitZIG(counts, scale_ZIG[1])
    MD <- meanDifferences(estimated = fitted,
                          observed = observed)
    name <- paste0("ZIG_",scale_ZIG[1])
    fittedModels[[name]] <- data.frame(observed, MD)
    if(length(scale_ZIG) == 2){
      fitted <- fitZIG(counts, scale_ZIG[2])
      MD <- meanDifferences(estimated = fitted,
                            observed = observed)
      name <- paste0("ZIG_",scale_ZIG[2])
      fittedModels[[name]] <- data.frame(observed, MD)
    }
  }
  if("HURDLE" %in% models){
    fitted <- fitHURDLE(counts, scale_HURDLE[1])
    MD <- meanDifferences(estimated = fitted,
                          observed = observed)
    name <- paste0("HURDLE_",scale_HURDLE[1])
    fittedModels[[name]] <- data.frame(observed, MD)
    if(length(scale_HURDLE) == 2){
      fitted <- fitHURDLE(counts, scale_HURDLE[2])
      MD <- meanDifferences(estimated = fitted,
                            observed = observed)
      name <- paste0("HURDLE_",scale_HURDLE[2])
      fittedModels[[name]] <- data.frame(observed, MD)
    }
  }
  return(fittedModels)
}

#' @title plotMD
#'
#' @importFrom plyr ddply
#' @import ggplot2
#' @export
#' @description
#' A function to plot mean difference (MD) and zero probability difference (ZPD)
#' values between estimated and observed values.
#'
#' @param data a `data.frame` object with a 'Model', 'Y', 'Y0', 'MD', and 'ZPD'
#' columns containing the model name, the observed values for the mean and the
#' zero proportion and the differences between observed and estimated values.
#' @param difference character vector, either 'MD' or 'ZPD' to plot the
#' differences between estimated and observed mean counts or the differences
#' between estimated zero probability and observed zero proportion.
#' @param split Display each model mean differences in different facets (default
#' = TRUE). If FALSE, points are not displayed for a more clear representation.
#'
#' @return a `ggplot` object.
#'
#' @seealso \code{\link{fitModels}} and \code{\link{RMSE}} for the model
#' estimations and the RMSE computations respectively.

plotMD <- function(data, difference = NULL, split = TRUE){
  # To avoid notes like: "no visible binding for global variable..."
  Y <- MD <- Model <- Y0 <- ZPD <- NULL
  if(difference == "MD"){
    if("Model" %in% colnames(data) & "Y" %in% colnames(data) &
       "MD" %in% colnames(data)){

      RMSE_MD <- plyr::ddply(.data = data,
                             .variables = ~ Model,
                             .fun = function(m) cbind("RMSE" = RMSE(m$MD)))

      gobj <- ggplot(data = data, aes(x = Y, y = MD, color = Model)) +
        ggtitle(label = "Mean Differences plot",
                subtitle = paste0("Observed = log(mean(counts*)+1)",
                                  "\n",
                                  "Estimated = log(mean(fitted*)+1)"))

      if(split){
        gobj <- gobj +
          geom_text(data = RMSE_MD, color = "black",
                    aes(x = mean(data$Y),
                        y = max(data$MD,na.rm = TRUE),
                        label = paste0("RMSE:",round(RMSE,2))))
      }

    } else {stop("data should contains 'Model', 'Y', and 'MD' columns for model
              name, observed values and mean difference values respectively.")}
  } else if(difference == "ZPD"){
    if("Model" %in% colnames(data) &
       "Y0" %in% colnames(data) &
       "ZPD" %in% colnames(data)){

      RMSE_ZPD <- plyr::ddply(.data = data,
                              .variables = ~ Model,
                              .fun = function(m) cbind("RMSE" = RMSE(m$ZPD)))

      gobj <- ggplot(data = data, aes(x = Y0,  y = ZPD,  color = Model)) +
        ggtitle(label = "Zero Probability Differences plot", subtitle =
                  "Observed = mean(counts=0)\nEstimated = mean(P(Y=0))")

      if(split){
        gobj <- gobj +
          geom_text(data = RMSE_ZPD, color = "black",
                    aes(x = mean(data$Y0),
                        y = max(data$ZPD,na.rm = TRUE),
                        label = paste0("RMSE:",round(RMSE,4))))
      }

    } else {stop("df should contains 'Model', 'Y0', and 'ZPD' columns for model
              name, zero rate observed values and zero probability difference
              values respectively.")}
  } else stop("Difference must be 'MD' or 'ZPD'")

  gobj <- gobj +
    geom_hline(aes(yintercept = 0), lty = 2, color = "black") +
    theme(legend.position = "bottom") +
    xlab("Observed") +
    ylab("Estimated-Observed")

  if(length(unique(data$Model))>1){
    if(split){
      gobj <- gobj +
        facet_grid(~ Model, labeller = labeller(.cols = label_both)) +
        geom_point(pch = 21) +
        geom_smooth(color = "black")
    } else {
      gobj <- gobj +
        geom_smooth()
    }
  }

  return(gobj)
}

#' @title plotRMSE
#'
#' @importFrom plyr ddply
#' @import ggplot2
#' @export
#' @description
#' A function to plot RMSE values computed for mean difference (MD) and zero
#' probability difference (ZPD) values between estimated and observed values.
#'
#' @param data a `data.frame` object with a 'Model', 'Y', 'Y0', 'MD', and 'ZPD'
#' columns containing the model name, the observed values for the mean and the
#' zero proportion and the differences between observed and estimated values.
#' @param difference character vector, either 'MD' or 'ZPD' to plot the RMSE
#' values for the differences between estimated and observed mean counts or the
#' differences between estimated zero probability and observed zero proportion.
#'
#' @return a `ggplot` object.
#'
#' @seealso \code{\link{fitModels}} and \code{\link{RMSE}} for the model
#' estimations and the RMSE computations respectively.

plotRMSE <- function(data, difference = NULL){
  # To avoid notes like: "no visible binding for global variable Model"
  Model <- NULL
  if(difference == "MD"){
    if("Model" %in% colnames(data) & "Y" %in% colnames(data) &
       "MD" %in% colnames(data)){

      RMSE <- plyr::ddply(.data = data,
                          .variables = ~ Model,
                          .fun = function(m) cbind("RMSE" = RMSE(m$MD)))

      gobj <- ggplot(data = RMSE, aes(x = Model, y = RMSE, fill = Model)) +
        geom_col() +
        geom_label(aes(x = Model, y = RMSE, label = round(RMSE,2)),
                   fill = "white") +
        ggtitle(label = "RMSE",
                subtitle = "Mean differences")

    } else {stop("data should contains 'Model', 'Y', and 'MD' columns for model
              name, observed values and mean difference values respectively.")}
  } else if(difference == "ZPD"){
    if("Model" %in% colnames(data) &
       "Y0" %in% colnames(data) &
       "ZPD" %in% colnames(data)){

      RMSE <- plyr::ddply(.data = data,
                          .variables = ~ Model,
                          .fun = function(m) cbind("RMSE" = RMSE(m$ZPD)))

      gobj <- ggplot(data = RMSE, aes(x = Model, y = RMSE, fill = Model)) +
        geom_col() +
        geom_label(aes(x = Model, y = RMSE, label = round(RMSE,4)),
                   fill = "white") +
        ggtitle(label = "RMSE",
                subtitle = "Zero probability difference")

    } else {stop("df should contains 'Model', 'Y0', and 'ZPD' columns for model
              name, zero rate observed values and zero probability difference
              values respectively.")}

  } else stop("Difference must be 'MD' or 'ZPD'")

  gobj <- gobj +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_x_discrete(limits = RMSE$Model[order(RMSE$RMSE)])

  return(gobj)
}
