#' @title createMocks
#'
#' @export
#' @description
#' Given the number of samples of the dataset from which the mocks should be
#' created, this function produces a \code{data.frame} object with as many rows
#' as the number of mocks and as many columns as the number of samples. If an
#' odd number of samples is given, the lower even integer will be considered in
#' order to obtain a balanced design for the mocks.
#' @param nsamples an integer representing the total number of samples.
#' @param N number of mock comparison to generate.
#'
#' @return a \code{data.frame} containing \code{N} rows and \code{nsamples}
#' columns (if even). Each cell of the data frame contains the "grp1" or "grp2"
#' characters which represent the mock groups pattern.
#'
#' @examples
#' # Generate the pattern for 100 mock comparisons for an experiment with 30
#' # samples
#' mocks <- createMocks(nsamples = 30, N = 100)
#' head(mocks)
createMocks <- function(nsamples, N = 1000) {
    sample_num <- nsamples %/% 2 * 2 # Balanced design sample numerosity
    mock_df <- matrix(NA, nrow = N, ncol = sample_num)
    for (i in seq_len(N)) { # N random balanced relabellings
        grps <- rep("grp1", sample_num)
        grps[sample(seq_len(sample_num), size = sample_num / 2)] <- "grp2"
        mock_df[i, ] <- grps
    }
    rownames(mock_df) <- paste0("Comparison", seq_len(N))
    return(mock_df)
} # END - function: createMocks

#' @title runMocks
#'
#' @export
#' @importFrom phyloseq sample_data
#' @description
#' Run the differential abundance detection methods on mock datasets.
#'
#' @param mocks a \code{data.frame} containing \code{N} rows and \code{nsamples}
#' columns (if even). Each cell of the data frame contains the "grp1" or "grp2"
#' characters which represent the mock groups pattern. Produced by the
#' \code{\link{createMocks}} function.
#' @inheritParams get_counts_metadata
#' @inheritParams runDA
#'
#' @return A named list containing the results for each method.
#'
#' @examples
#' # Load some data
#' data(ps_stool_16S)
#'
#' # Generate the pattern for 10 mock comparisons
#' # (N = 1000 is suggested)
#' mocks <- createMocks(nsamples = phyloseq::nsamples(ps_stool_16S), N = 10)
#' head(mocks)
#'
#' # Add some normalization/scaling factors to the phyloseq object
#' my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
#'     method = c("TMM", "CSS"))
#' ps_stool_16S <- runNormalizations(normalization_list = my_norm,
#'     object = ps_stool_16S)
#'
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ group, coef = 2, norm = c("TMM", "CSS"))
#'
#' # Run methods on mock datasets
#' results <- runMocks(mocks = mocks, method_list = my_limma,
#'     object = ps_stool_16S)
runMocks <- function(mocks, method_list, object, weights = NULL, 
    verbose = TRUE){
    is_phyloseq <- ifelse(is(object, "phyloseq"), TRUE, FALSE)
    index <- seq_len(nrow(mocks))
    out <- apply(X = cbind(index, mocks), MARGIN = 1, FUN = function(x) {
        # Group assignment
        i <- x[1]
        x <- x[-1]
        if(is_phyloseq){
            phyloseq::sample_data(object)[, "group"] <- factor(x)
        } else {
            SummarizedExperiment::colData(object)[, "group"] <- factor(x)
        }
        if(verbose)
            message("  - Comparison", i, "\n")
        runDA(method_list = method_list, object = object,
            weights = weights, verbose = verbose)
    })
    return(out)
}

#' @title createTIEC
#'
#' @export
#' @importFrom plyr ldply ddply
#' @importFrom stats ks.test
#' @importFrom reshape2 melt
#' @description
#' Extract the list of p-values from the outputs of the differential abundance
#' detection methods to compute several statistics to study the ability to
#' control the type I error and the p-values distribution.
#'
#' @param object Output of the differential abundance tests on mock comparisons.
#' Must follow a specific structure with comparison, method, matrix of
#' p-values, and method's name (See vignette for detailed information).
#'
#' @return A \code{list} of \code{data.frame}s:
#' \itemize{
#'     \item{\code{df_pval}}{ 5 columns per number_of_features x methods x
#'     comparisons rows data.frame. The four columns are called Comparison,
#'     Method, variable (containing the feature names), pval, and padj;}
#'     \item{\code{df_FPR}}{ 5 columns per methods x comparisons rows 
#'     data.frame. For each set of method and comparison, the proportion of 
#'     false discoveries, considering 3 threshold (0.01, 0.05, 0.1) are 
#'     reported;}
#'     \item{\code{df_QQ}}{ contains the coordinates to draw the QQ-plot to
#'     compare the mean observed p-value distribution across comparisons, with
#'     the theoretical uniform distribution;}
#'     \item{\code{df_KS}}{ 5 columns and methods x comparisons rows data.frame.
#'     For each set of method and comparison, the Kolmogorov-Smirnov test
#'     statistics and p-values are reported in KS and KS_pval columns
#'     respectively.}}
#'
#' @seealso \code{\link{createMocks}}
#'
#' @examples
#' # Load some data
#' data(ps_stool_16S)
#'
#' # Generate the patterns for 10 mock comparison for an experiment
#' # (N = 1000 is suggested)
#' mocks <- createMocks(nsamples = phyloseq::nsamples(ps_stool_16S), N = 10)
#' head(mocks)
#'
#' # Add some normalization/scaling factors to the phyloseq object
#' my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
#'     method = c("TMM", "CSS"))
#' ps_stool_16S <- runNormalizations(normalization_list = my_norm,
#'     object = ps_stool_16S)
#'
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ group, coef = 2,
#'     norm = c("TMM", "CSS"))
#'
#' # Run methods on mock datasets
#' results <- runMocks(mocks = mocks, method_list = my_limma,
#'     object = ps_stool_16S)
#'
#' # Prepare results for Type I Error Control
#' TIEC_summary <- createTIEC(results)
#'
#' # Plot the results
#' plotFPR(df_FPR = TIEC_summary$df_FPR)
#' plotQQ(df_QQ = TIEC_summary$df_QQ, zoom = c(0, 0.1))
#' plotKS(df_KS = TIEC_summary$df_KS)
#' plotLogP(df_QQ = TIEC_summary$df_QQ)
createTIEC <- function(object) {
    # Create a list of data frames with two columns: pval and method name
    # One data frame for each comparison
    message("1. Extracting statistics")
    df_list_pval <- lapply(X = object, FUN = function(Comparison) {
        pval_df <- reshape2::melt(
            plyr::ldply(
                extractStatistics(object = Comparison,
                    slot = "pValMat", colName = "rawP", type = "pvalue",
                    direction = NULL, verbose = FALSE
                ), .id = "Method"
            ), value.name = "pval"
        )
        padj_df <- reshape2::melt(
            plyr::ldply(
                extractStatistics(object = Comparison,
                    slot = "pValMat", colName = "adjP", type = "pvalue",
                    direction = NULL, verbose = FALSE
                ), .id = "Method"
            ), value.name = "padj"
        )
        return(data.frame(pval_df, "padj" = padj_df[ ,"padj"]))
    })
    # Melt down to a single data frame with the Comparison column added
    df_pval <- plyr::ldply(df_list_pval, .id = "Comparison")
    ### FPR ###
    message("2. Counting p-values lower than some thresholds")
    # Count the p-values which are lower than a selected threshold
    df_pval_FPR <- plyr::ddply(.data = df_pval, .variables = ~ Comparison +
        Method, .fun = function(x) {
        # index for not NA p-values
        k_pval <- !is.na(x[, "pval"])
        x$FPR_obs001 <- sum(x[, "pval"] < 0.01, na.rm = TRUE) / sum(k_pval)
        x$FPR_obs005 <- sum(x[, "pval"] < 0.05, na.rm = TRUE) / sum(k_pval)
        x$FPR_obs01 <- sum(x[, "pval"] < 0.1, na.rm = TRUE) / sum(k_pval)
        return(x)
    })
    # Compute the mean for each FPR threshold
    df_FPR <- plyr::ddply(.data = df_pval_FPR, .variables = ~ Comparison +
        Method, .fun = function(x) {
        colMeans(x[, c("FPR_obs001", "FPR_obs005", "FPR_obs01")])
    })
    ### QQ and KS ###
    message("3. Computing KS statistics")
    # Create data frame for empirical and theoretical p-values
    # Kolmogorov-Smirnov test
    df_QQ_KS <- plyr::ddply(.data = df_pval, .variables = ~ Comparison + Method,
        .fun = function(x) {
        # Ordered list of p-values
        x <- x[order(x[, "pval"]), ]
        # Index of not NA p-values
        k_pval <- !is.na(x[, "pval"])
        # Theoretical uniform distribution
        x$pval_theoretical[k_pval] <- (seq_len(sum(k_pval))) / (sum(k_pval))
        x$pval_theoretical_rounded[k_pval] <- round((seq_len(sum(k_pval))) /
            (sum(k_pval)), digits = 2)
        x[c("KS","KS_pval")] <- stats::ks.test(
            x = x$pval_theoretical[k_pval],
            y = x$pval[k_pval])[c("statistic","p.value")]
        return(x)
    })
    ### QQ ###
    message("4. Ordering quantiles")
    # For each theoretical p-value, the mean of observed one is computed
    df_QQ <- plyr::ddply(.data = df_QQ_KS, .variables = ~ Method +
        pval_theoretical_rounded, .fun = function(x) {
        pval_mean <- mean(x$pval)
        return(pval_mean)
    })
    names(df_QQ) <- c("Method", "pval_theoretical", "pval_observed")
    ### KS ###
    # Compute the mean for KS statistics and p-values
    df_KS <- plyr::ddply(
        .data = df_QQ_KS, .variables = ~ Comparison + Method,
        .fun = function(x) {
            colMeans(x[, c("KS", "KS_pval")])
        }
    )
    return(list(
        df_pval = df_pval, df_FPR = df_FPR, df_QQ = df_QQ, df_KS =
            df_KS
    ))
} # END - function: createTIEC
