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
    sample_num <- ncol(counts) %/% 2 * 2 # Balanced design sample numerosity
    mock_df <- matrix(NA,nrow = N,ncol = sample_num)
    set.seed(seed)
    for(i in 1:N){ # N random balanced relabellings
        grps <- rep("grp1", sample_num)
        grps[sample(1:sample_num,size = sample_num/2)] <- "grp2"
        mock_df[i,] <- grps
    }
    rownames(mock_df) <- paste0("Comparison",1:1000)
    return(mock_df)
}

############# Differential Abundance Methods #############

### Normalizations ###

#' @title norm_edgeR
#'
#' @export
#' @description
#' Calculate normalization factors from a phyloseq object to scale the raw
#' library sizes. Inherited from edgeR `calcNormaFactors` function.
#' @param physeq A phyloseq object.
#' @inheritParams edgeR::calcNormFactors
#'
#' @return A new column containing the chosen edgeR-based normalization factors
#' is added to the phyloseq `sample_data` slot. The effective library sizes for
#' use in downstream analysis must be multiplied by the normalization factors.
#' @seealso [edgeR::calcNormFactors()] for details.

norm_edgeR <- function(physeq, method = "TMM")
{
    counts <- as(otu_table(physeq), "matrix")
    if (!taxa_are_rows(physeq))
    {
        counts <- t(counts)
    } else {}

    if (method == "upperquartile")
    {
        scaledCounts <- t(otuTab) / colSums(otuTab)
        tmpNF <- apply(scaledCounts, MARGIN = 1L, FUN = function(x)
            quantile(x[x != 0], probs = .75))
        normFacts <- tmpNF/exp(mean(log(tmpNF)))
        method <- "UQ"
    } else {
        normFacts <- edgeR:::calcNormFactors(counts, method = method)
    }# END - ifelse: upperquartile only of non-zero counts

    # VERY IMPORTANT: multiply by library sizes and renormalize.
    # edgeR calculates scaling factors, which still have to be multiplied by
    # library sizes to get to the size factors of effective sequencing depth,
    # i.e. robust estimates of the library sizes.
    # normFacts = normFacts*sample_sums(physeq)
    # normFacts = normFacts/exp(mean(log(normFacts)))

    if (all(is.na(normFacts))) #Resort to proportion normalization in case of
        # failure for all samples
    {
        normFacts = sample_sums(physeq)
    }
    physeq@sam_data@.Data <- c(physeq@sam_data@.Data, list(normFacts))
    aux <- physeq@sam_data@names
    aux[length(aux)] <- paste("NF", method, sep = ".")
    physeq@sam_data@names <- aux
    physeq
}# END - function: normEdgeR
