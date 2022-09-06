#' @title createSplits
#'
#' @export
#' @description
#' Given a phyloseq or TreeSummarizedExperiment object from which the random 
#' splits should be created, this function produces a list of 2 
#' \code{data.frame} objects: \code{Subset1} and \code{Subset2} with as many 
#' rows as the number of splits and as many columns as the half of the number 
#' of samples.
#' @inheritParams get_counts_metadata
#' @param varName name of a factor variable of interest.
#' @param paired name of the unique subject identifier variable. If specified,
#' paired samples will remain in the same split. (default = NULL).
#' @param balanced If \code{TRUE} a balanced design will be created for the
#' splits (default \code{balanced = TRUE}).
#' @param N number of splits to generate.
#'
#' @return A list of 2 \code{data.frame} objects: \code{Subset1} and
#' \code{Subset2} containing \code{N} rows and half of the total number of
#' samples columns. Each cell contains a unique sample identifier.
#' @examples
#' data(ps_plaque_16S)
#' set.seed(123)
#'
#' # Balanced design for repeated measures
# splits_df <- createSplits(
#     object = ps_plaque_16S, varName =
#         "HMP_BODY_SUBSITE", paired = "RSID", balanced = TRUE, N = 100
# )
#'
#' # Balanced design for independent samples
#' splits_df <- createSplits(
#'     object = ps_plaque_16S, varName =
#'         "HMP_BODY_SUBSITE", balanced = TRUE, N = 100
#' )
#'
#' # Unbalanced design
#' splits_df <- createSplits(
#'     object = ps_plaque_16S, varName =
#'         "HMP_BODY_SUBSITE", balanced = FALSE, N = 100
#' )

createSplits <- function(object, assay_name = "counts", varName = NULL, 
    paired = NULL, balanced = TRUE, N = 1000) {
    counts_metadata <- get_counts_metadata(object, assay_name = assay_name)
    metadata <- counts_metadata$metadata
    # Take variable names and levels
    if (is.null(varName)) {
        stop("varName: please supply the name of the variable to perform the", 
            " splitting")
    } else if (!is.element(varName, colnames(metadata))) {
        stop("varName: variable not found")
    } else {
        variable <- metadata[, varName]
        if (!is.factor(variable)) {
            warning("The variable ", varName, " is not a factor.",
                " Coercing to factor.")
            variable <- as.factor(variable)
        }
        var_levels <- levels(variable)
        n_levels <- nlevels(variable)
    }
    # In case of paired samples
    if (!is.null(paired)) {
        # Table for counting paired samples
        paired_IDs_table <- table(metadata[, paired])
        # A sample is considered paired if it is repeated as many times as the
        # others
        paired_IDs <- names(which(table(metadata[, paired]) == 
            max(paired_IDs_table)))
        unpaired <- setdiff(names(paired_IDs_table), paired_IDs)
        # Removing unpaired samples from metadata
        if(length(unpaired) > 0){
            warning(paste(unpaired, collapse = ", "), 
                " samples (", paired, " variable) are not paired.",
                " Removing them...")
            paired_metadata <- which(match(x = metadata[, paired], 
                table = paired_IDs, nomatch = 0) > 0)
            metadata <- metadata[paired_metadata, ]
        }
        num_paired_IDs <- length(paired_IDs)
        # Compute the number of repeated measures for each set of varName
        # and paired variables
        design <- table(metadata[, c(varName, paired)])
        # Count how many reps for each paired ID
        n_reps <- unique(colSums(design))
        # Check whether the reps are inside one group only
        repeated <- ifelse(length(unique(apply(design, 2, unique))) > 1, 
            TRUE, FALSE) 
        if(length(n_reps) > 1){
            stop("The design is too complex to generate paired and", 
                 " balanced subsets. Please set paired = NULL, ", 
                 " balanced = FALSE, or use a simpler design.")
        }
        message("Working on ", num_paired_IDs, " paired samples. ", 
            n_reps, " repeated measures for each ", paired, " ID.")
    } 
    # find indexes for the levels
    names(var_levels) <- var_levels
    # The IDs are the row names
    IDs <- lapply(var_levels, function(level) 
        rownames(metadata[metadata[, varName] == level, ]))
    num_IDs <- lapply(IDs, length)
    # In case of balanced design
    if (balanced) {
        # When the samples are not paired
        if(is.null(paired)){
            # Take as reference the littlest group
            min_length <- min(unlist(num_IDs))
            # Choose the first min_length IDs for each group
            IDs <- lapply(IDs, function(ID) return(ID[seq_len(min_length)]))
            num_IDs <- lapply(IDs, length)
            # We compute the half for each group
            # The first half goes to Subset1, the second half goes to Subset2
            half_s1 <- unlist(num_IDs) %/% 2
            half_s2 <- unlist(num_IDs) - half_s1
            # We prepare the matrices
            Subset1 <- matrix(NA, nrow = N, ncol = sum(half_s1))
            Subset2 <- matrix(NA, nrow = N, ncol = sum(half_s2))
        # When the samples are paired
        } else {
            # Take as reference the littlest group
            min_length <- min(unlist(num_IDs))
            if(repeated){
                min_length <- min_length %/% n_reps
            }
            # Choose the first min_length paired IDs for each group
            IDs <- lapply(IDs, function(ID){
                sub_IDs <- unique(metadata[ID, paired])[seq_len(min_length)]
                sub_index_IDs <- which(match(x = metadata[ID, paired],
                    table = sub_IDs, nomatch = 0) > 0)
                return(ID[sub_index_IDs])
            })
            num_IDs <- lapply(IDs, length)
            # We compute the half for each group
            # The first half goes to Subset1, the second half goes to Subset2
            if(repeated){ # if multiple observation are present for each paired
                # ID
                half_s1 <- (unlist(num_IDs) / n_reps) %/% 2
                half_s2 <- unlist(num_IDs) / n_reps - half_s1
                # We prepare the matrices
                Subset1 <- matrix(NA, nrow = N, ncol = sum(half_s1) * n_reps)
                Subset2 <- matrix(NA, nrow = N, ncol = sum(half_s2) * n_reps)
            } else {
                half_s1 <- unlist(num_IDs) %/% 2
                half_s2 <- unlist(num_IDs) - half_s1
                # We prepare the matrices
                Subset1 <- matrix(NA, nrow = N, ncol = sum(half_s1))
                Subset2 <- matrix(NA, nrow = N, ncol = sum(half_s2))
            }
        }
    # In case of unbalanced design
    } else {
        if(is.null(paired)){
            # We compute the half for each group
            # The first half goes to Subset1, the second half goes to Subset2
            half_s1 <- unlist(num_IDs) %/% 2
            half_s2 <- unlist(num_IDs) - half_s1
            # We prepare the matrices
            Subset1 <- matrix(NA, nrow = N, ncol = sum(half_s1))
            Subset2 <- matrix(NA, nrow = N, ncol = sum(half_s2))
        } else {
            # if multiple observation are present for each paired
            if(repeated){ 
                # ID
                half_s1 <- (unlist(num_IDs) / n_reps) %/% 2
                half_s2 <- unlist(num_IDs) / n_reps - half_s1
                # We prepare the matrices
                Subset1 <- matrix(NA, nrow = N, ncol = sum(half_s1) * n_reps)
                Subset2 <- matrix(NA, nrow = N, ncol = sum(half_s2) * n_reps)
            } else {
                half_s1 <- unlist(num_IDs) %/% 2
                half_s2 <- unlist(num_IDs) - half_s1
                # We prepare the matrices
                Subset1 <- matrix(NA, nrow = N, ncol = sum(half_s1))
                Subset2 <- matrix(NA, nrow = N, ncol = sum(half_s2))
            }
        }
    } 
    # Populating the Subset1 and Subset2 matrices
    for (i in seq_len(N)) {
        # In case of not paired samples, just sample the IDs for each group
        if(is.null(paired)){
            chosenIDs <- mapply(IDs, half_s1 , FUN = function(ID, h){
                list(sample(x = ID, size = h, replace = FALSE))
            })
        # In case of paired samples
        } else {
            # Find the indexes for each paired ID if they are repeated
            if(repeated){
                chosen_paired_IDs <- mapply(IDs, half_s1, FUN = function(ID, h){
                    group_paired_IDs <- unique(metadata[ID, paired])
                    list(sample(x = group_paired_IDs, size = h, 
                                replace = FALSE))
                })
                # Map the chosen paired IDs indexes to metadata rownames
                chosenIDs <- mapply(IDs, chosen_paired_IDs, 
                    FUN = function(ID, i){
                    chosen_index_IDs <- which(match(x = metadata[ID, paired],
                        table = unlist(chosen_paired_IDs), nomatch = 0) > 0)
                    list(ID[chosen_index_IDs])
                })
            # In case of 1 observation for each group
            } else {
                chosen_paired_IDs <- sample(x = unique(metadata[, paired]), 
                    size = unique(half_s1), replace = FALSE)
                # Map the chosen paired IDs indexes to metadata rownames
                chosen_index_IDs <- which(match(x = metadata[, paired],
                    table = unlist(chosen_paired_IDs), nomatch = 0) > 0)
                chosenIDs <- rownames(metadata[chosen_index_IDs, ])
            }
        }
        Subset1[i, ] <- unlist(chosenIDs)
        Subset2[i, ] <- setdiff(unlist(IDs), unlist(chosenIDs))
    }
    rownames(Subset1) <- rownames(Subset2) <- paste0("Comparison", seq_len(N))
    return(list("Subset1" = Subset1, "Subset2" = Subset2))
} # END - function: createSplits

#' @title runSplits
#'
#' @export
#' @importFrom phyloseq prune_samples sample_names filter_taxa
#' @importFrom phyloseq filter_taxa prune_samples sample_data
#' @description
#' Run the differential abundance detection methods on split datasets.
#'
#' @param split_list A list of 2 \code{data.frame} objects: \code{Subset1} and
#' \code{Subset2} produced by the \code{\link{createSplits}} function.
#' @inheritParams runDA
#' @param normalization_list a list object containing the normalization method
#' names and their parameters produced by \code{\link{setNormalizations}}.
#' @inheritParams get_counts_metadata
#'
#' @return A named list containing the results for each method.
#'
#' @examples
#' data(ps_plaque_16S)
#' 
#' # Balanced design
#' my_splits <- createSplits(
#'     object = ps_plaque_16S, varName = "HMP_BODY_SUBSITE", balanced = TRUE,
#'     paired = "RSID", N = 10 # N = 100 suggested
#' )
#' 
#' # Make sure the subject ID variable is a factor
#' phyloseq::sample_data(ps_plaque_16S)[, "RSID"] <- as.factor(
#'     phyloseq::sample_data(ps_plaque_16S)[["RSID"]])
#'
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ RSID + HMP_BODY_SUBSITE, 
#'     coef = "HMP_BODY_SUBSITESupragingival Plaque",
#'     norm = c("TMM", "CSS"))
#'
#' # Set the normalization methods according to the DA methods
#' my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
#'     method = c("TMM", "CSS"))
#'
#' # Run methods on split datasets
#' results <- runSplits(split_list = my_splits, method_list = my_limma,
#'     normalization_list = my_norm, object = ps_plaque_16S)

runSplits <- function(split_list, method_list, normalization_list, object,
    assay_name = "counts", min_counts = 0, min_samples = 0, verbose = TRUE){
    # Create a list with element 1 and 2 corresponding to each subset
    subsets <- seq_along(split_list)
    # Name the elements
    names(subsets) <- names(split_list)
    out <- lapply(X = subsets, FUN = function(subset_number) {
        subset <- split_list[[subset_number]]
        if(verbose)
            message("- Subset", subset_number, "\n")
        # apply -> Comparison1, Comparison2, ..., ComparisonN
        index <- seq_len(nrow(subset))
        apply(X = cbind(index,subset), MARGIN = 1, FUN = function(splits) {
            i <- splits[1]
            splits <- splits[-1]
            if(verbose)
                message("  - Comparison", i, "\n")
            # Splitting
            if(verbose)
                message("    Splitting the samples...")
            # Get counts and metadata
            counts_and_metadata <- get_counts_metadata(object, 
                assay_name = assay_name)
            counts <- counts_and_metadata[[1]]
            is_phyloseq <- counts_and_metadata[[3]]
            if(is_phyloseq){
                # Pruned object
                po <- phyloseq::prune_samples(
                    phyloseq::sample_names(object) %in% splits, object)
                # Keep only present taxa
                if(verbose)
                    message("    Keeping taxa with more than ",
                            min_counts, " counts in more than ", min_samples, 
                            " samples.")
                # Pruned and filtered object
                pfo <- phyloseq::filter_taxa(po, function(x) 
                    sum(x > min_counts) > min_samples, 1)
            } else{
                # Pruned object
                po <- object[, splits]
                if(verbose)
                    message("    Keeping taxa with more than ",
                            min_counts, " counts in more than ", min_samples, 
                            " samples.")
                taxa_to_keep <- apply(counts, 1, function(x) 
                    sum(x > min_counts) > min_samples)
                # Pruned and filtered object
                pfo <- po[taxa_to_keep, ]
            }
            # Adding scaling and normalization factors
            if(!is.null(normalization_list)){
                if(verbose)
                    message("    Computing normalizations...")
                pfo <- runNormalizations(
                    normalization_list = normalization_list,
                    object = pfo, verbose = verbose)
            }
            # Compute weights if necessary
            weights_info <- unlist(lapply(X = method_list, FUN = function(x){
                x[["weights"]]
            }))
            if(sum(weights_info) > 0){
                if(verbose)
                    message("    Computing ZINB weights...")
                weights <- weights_ZINB(pfo, design = ~ 1)
            } else weights <- NULL
            ### DA analysis ###
            if(verbose)
                message("    Differential abundance:")
            runDA(method_list = method_list, object = pfo, weights = weights,
                verbose = verbose)
        })
    })
    return(out)
}

#' @title createConcordance
#'
#' @export
#' @importFrom plyr ddply ldply
#' @importFrom stats runif
#' @description
#' Compute the between and within method concordances comparing the lists of
#' extracted statistics from the outputs of the differential abundance detection
#' methods.
#'
#' @inheritParams extractStatistics
#'
#' @return A long format \code{data.frame} object with several columns:
#' \itemize{
#'     \item{\code{comparison}}{ which indicates the comparison number;}
#'     \item{\code{n_features}}{ which indicates the total number of taxa in
#'     the comparison dataset;}
#'     \item{\code{method1}}{ which contains the first method name;}
#'     \item{\code{method2}}{ which contains the first method name;}
#'     \item{\code{rank}}{;}
#'     \item{\code{concordance}}{ which is defined as the cardinality of the
#'     intersection of the top rank elements of each list, divided by rank, i.e.
#'     , \eqn{(L_{1:rank} \bigcap M_{1:rank})/(rank)}, where L and M represent
#'     the lists of the extracted statistics of method1 and method2
#'     respectively (averaged values between subset1 and subset2).}}
#'
#' @seealso \code{\link{extractStatistics}} and \code{\link{areaCAT}}.
#'
#' @examples
#' data(ps_plaque_16S)
#' 
#' # Balanced design
#' my_splits <- createSplits(
#'     object = ps_plaque_16S, varName = "HMP_BODY_SUBSITE", balanced = TRUE,
#'     paired = "RSID", N = 10 # N = 100 suggested
#' )
#' 
#' # Make sure the subject ID variable is a factor
#' phyloseq::sample_data(ps_plaque_16S)[, "RSID"] <- as.factor(
#'     phyloseq::sample_data(ps_plaque_16S)[["RSID"]])
#'
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ RSID + HMP_BODY_SUBSITE, 
#'     coef = "HMP_BODY_SUBSITESupragingival Plaque",
#'     norm = c("TMM", "CSS"))
#'
#' # Set the normalization methods according to the DA methods
#' my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
#'     method = c("TMM", "CSS"))
#'
#' # Run methods on split datasets
#' results <- runSplits(split_list = my_splits, method_list = my_limma,
#'     normalization_list = my_norm, object = ps_plaque_16S)
#'
#' # Concordance for p-values
#' concordance_pvalues <- createConcordance(
#'     object = results, slot = "pValMat", colName = "rawP", type = "pvalue"
#' )
#'
#' # Concordance for log fold changes
#' concordance_logfc <- createConcordance(
#'     object = results, slot = "statInfo", colName = "logFC", type = "logfc"
#' )
#'
#' # Concordance for log fold changes in the first method and p-values in the
#' # other
#' concordance_logfc_pvalues <- createConcordance(
#'     object = results, slot = c("statInfo", "pValMat"),
#'     colName = c("logFC", "rawP"), type = c("logfc", "pvalue")
#' )

createConcordance <- function(object, slot = "pValMat", colName = "rawP",
    type = "pvalue") {
    data <- lapply(X = object, FUN = function(subset) {
        lapply(X = subset, FUN = function(comparison) {
            c(extractStatistics(
                object = comparison, slot = slot,
                colName = colName, type = type
            ),
            n_features =
                nrow(comparison[[1]][["pValMat"]])
            )
        })
    })
    subset1 <- data[["Subset1"]]
    subset2 <- data[["Subset2"]]
    concordance <- plyr::ldply(seq_len(length(subset1)), function(i){
        comparison1 <- subset1[[i]]
        comparison2 <- subset2[[i]]
        n_features1 <- comparison1[["n_features"]]
        n_features2 <- comparison2[["n_features"]]
        n_features <- min(n_features1, n_features2)
        n_methods <- length(comparison1) - 1
        method_names <- names(comparison1)[seq_len(n_methods)]
        plyr::ldply(seq_len(n_methods), function(j) {
            plyr::ldply(seq_len(n_methods), function(k) {
                if (j != k) { # Compute the Between Methods concordance
                    # For subset1
                    # Make numeric the named vectors and add noise.
                    vec1 <- comparison1[[method_names[j]]] 
                    vec2 <- comparison1[[method_names[k]]]  
                    noise <- stats::runif(length(vec1), 0, 1e-10)
                    conc1 <- CAT(vec1 = vec1 + noise, vec2 = vec2 + noise)
                    conc1 <- data.frame(conc1, "n_features" = n_features1)
                    # For subset2
                    vec1 <- comparison2[[method_names[j]]]
                    vec2 <- comparison2[[method_names[k]]]
                    noise <- stats::runif(length(vec1), 0, 1e-10)
                    conc2 <- CAT(vec1 = vec1 + noise, vec2 = vec2 + noise)
                    conc2 <- data.frame(conc2, "n_features" = n_features2)
                    # Together
                    conc <- rbind(conc1, conc2)
                } else { # Compute the Within Method concordance
                    vec1 <- comparison1[[method_names[j]]]
                    vec2 <- comparison2[[method_names[k]]]
                    noise1 <- stats::runif(length(vec1), 0, 1e-10)
                    noise2 <- stats::runif(length(vec2), 0, 1e-10)
                    conc <- CAT(vec1 = vec1 + noise1, vec2 = vec2 + noise2)
                    conc <- data.frame(conc, "n_features" = n_features)
                }
                conc[, "method1"] <- method_names[j]
                conc[, "method2"] <- method_names[k]
                conc[, "comparison"] <- paste0("Comparison", i)
                return(conc)
            })
        })
    })
    return(concordance)
}

#' @title CAT
#'
#' @export
#' @description
#' For the i top-ranked members of each list, concordance is defined as 
#' \code{length(intersect(vec1[1:i],vec2[1:i]))/i}.
#'
#' @param vec1,vec2 Two numeric vectors, for computing concordance. If these 
#' are numeric vectors with names, the numeric values will be used for sorting 
#' and the names will be used for calculating concordance. Otherwise, they are 
#' assumed to be already-ranked vectors, and the values themselves will be 
#' used for calculating concordance.
#' @param maxrank Optionally specify the maximum size of top-ranked items that 
#' you want to plot.
#'
#' @return a data.frame with two columns: \code{rank} containing the length of 
#' the top lists and \code{concordance} which is the fraction in common that 
#' the two provided lists have in the top \code{rank} items.
#'
#' @seealso \code{\link{createConcordance}}.
#'
#' @examples
#' vec1 <- c("A" = 10, "B" = 5, "C" = 20, "D" = 15)
#' vec2 <- c("A" = 1, "B" = 2, "C" = 3, "D" = 4)
#' 
#' CAT(vec1, vec2)

CAT <- function (vec1, vec2, maxrank = min(length(vec1), length(vec2))) 
{
    if (class(vec1) == "numeric" & class(vec2) == "numeric" & 
        !is.null(names(vec1)) & !is.null(names(vec1))) {
        vec1 <- sort(vec1)
        vec1 <- names(vec1)
        vec2 <- sort(vec2)
        vec2 <- names(vec2)
    }
    if (is.na(maxrank) | is.null(maxrank) | 
        maxrank > min(length(vec1), length(vec2))) {
        maxrank <- min(length(vec1), length(vec2))
    }
    output <- data.frame(rank = 1:maxrank, concordance = NA)
    for (i in 1:nrow(output)) {
        output[i, "concordance"] <- length(
            intersect(vec1[1:i], vec2[1:i]))/i
    }
    return(output)
}

#' @title areaCAT
#'
#' @export
#' @importFrom plyr ddply
#' @importFrom graphics plot abline
#' @description
#' Compute the area between the bisector and the concordance curve.
#'
#' @param concordance A long format \code{data.frame} produced by
#' \link{createConcordance} function.
#' @param plotIt Plot the concordance (default \code{plotIt = FALSE}).
#'
#' @return A long format \code{data.frame} object with several columns:
#' \itemize{
#'     \item{\code{comparison}}{ which indicates the comparison number;}
#'     \item{\code{n_features }}{ which indicates the total number of taxa in
#'     the comparison dataset;}
#'     \item{\code{method1}}{ which contains the first method name;}
#'     \item{\code{method2}}{ which contains the first method name;}
#'     \item{\code{rank}}{;}
#'     \item{\code{concordance}}{ which is defined as the cardinality of the
#'     intersection of the top rank elements of each list, divided by rank, i.e.
#'     , \eqn{(L_{1:rank} \bigcap M_{1:rank})/(rank)}, where L and M represent
#'     the lists of the extracted statistics of method1 and method2
#'     respectively;}
#'     \item{\code{heightOver}}{ which is the distance between the bisector and
#'     the concordance value;}
#'     \item{\code{areaOver}}{ which is the cumulative sum of the
#'     \code{heightOver} value.}}
#'
#' @seealso \code{\link{createConcordance}} and \code{\link{plotConcordance}}
#'
#' @examples
#' data(ps_plaque_16S)
#'
#' # Balanced design for dependent samples
#' my_splits <- createSplits(
#'     object = ps_plaque_16S, varName = "HMP_BODY_SUBSITE", 
#'     balanced = TRUE, paired = "RSID", N = 10 # N = 100 suggested
#' )
#' 
#' # Make sure the subject ID variable is a factor
#' phyloseq::sample_data(ps_plaque_16S)[, "RSID"] <- as.factor(
#'     phyloseq::sample_data(ps_plaque_16S)[["RSID"]])
#'
#' # Initialize some limma based methods
#' my_limma <- set_limma(design = ~ RSID + HMP_BODY_SUBSITE, 
#'     coef = "HMP_BODY_SUBSITESupragingival Plaque",
#'     norm = c("TMM", "CSS"))
#'
#' # Set the normalization methods according to the DA methods
#' my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
#'     method = c("TMM", "CSS"))
#'
#' # Run methods on split datasets
#' results <- runSplits(split_list = my_splits, method_list = my_limma,
#'     normalization_list = my_norm, object = ps_plaque_16S)
#'
#' # Concordance for p-values
#' concordance_pvalues <- createConcordance(
#'     object = results, slot = "pValMat", colName = "rawP", type = "pvalue"
#' )
#'
#' # Add area over the concordance curve
#' concordance_area <- areaCAT(concordance = concordance_pvalues)

areaCAT <- function(concordance, plotIt = FALSE) {
    plyr::ddply(.data = concordance, .variables = ~ comparison + n_features +
        method1 + method2, .fun = function(conc) {
        MaxArea <- unique(conc[, "n_features"])
        estimated <- conc[, "concordance"]
        # y = x values -> bisector
        theoretical <- seq_along(estimated) / unique(conc[, "n_features"])
        if (plotIt) {
            graphics::plot(x = theoretical, y = estimated, type = "l", xlim = c(
                0, 1
            ), ylim = c(0, 1), main = paste0("CAT for ", unique(conc[
                ,
                "method1"
            ]), " & ", unique(conc[, "method2"])))
            graphics::abline(0, 1)
        }
        # Consider the y = x line
        difference <- estimated - theoretical
        HeightOver <- (estimated - theoretical) / MaxArea # rescaling
        HeightOver[difference <= 0] <- 0 # removing negative values
        return(data.frame(
            "rank" = conc[, "rank"],
            "concordance" = conc[, "concordance"], "heightOver" = HeightOver,
            "areaOver" = cumsum(HeightOver)
        ))
    })
}
