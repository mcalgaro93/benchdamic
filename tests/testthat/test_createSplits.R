# Unit test based on createSplits

test_that(desc = "data splits: design for repeated measures", code = {
    data("ps_plaque_16S")
    set.seed(123)
    # Warning for non factor varName variable
    splits_df <- expect_warning(suppressMessages(
        createSplits(object = ps_plaque_16S, varName = "HMP_BODY_SUBSITE", 
            paired = "RSID", balanced = TRUE, N = 10)))
    # Make it a factor
    phyloseq::sample_data(ps_plaque_16S)[,"HMP_BODY_SUBSITE"] <-
        factor(unlist(phyloseq::sample_data(ps_plaque_16S)
        [,"HMP_BODY_SUBSITE"]))
    ### Design for repeated measures
    splits_df <- expect_message(
        createSplits(object = ps_plaque_16S, varName = "HMP_BODY_SUBSITE", 
            paired = "RSID", N = 10))
    expect_equal(nrow(splits_df$Subset1), 10)
    expect_equal(nrow(splits_df$Subset2), 10)
    expect_equal(ncol(splits_df$Subset1),
                 sum(table(phyloseq::sample_data(ps_plaque_16S)[,"RSID"])>1))
    expect_equal(ncol(splits_df$Subset2),
                 sum(table(phyloseq::sample_data(ps_plaque_16S)[,"RSID"])>1))
})

test_that(desc = "data splits: balanced design for independent samples (starting
    from unbalanced groups)", code = {
    data("ps_plaque_16S")
    set.seed(123)
    ### Balanced design for independent samples starting from an unbalanced
    # Remove some sample to create an unbalanced design
    ps_unbalanced <- phyloseq::subset_samples(ps_plaque_16S, c(FALSE,TRUE,TRUE))
    # table(phyloseq::sample_data(ps_unbalanced)[,"HMP_BODY_SUBSITE"])
    splits_df <- expect_warning(createSplits(object = ps_unbalanced, varName =
        "HMP_BODY_SUBSITE", balanced = TRUE, N = 10))
    expect_equal(nrow(splits_df$Subset1), 10)
    expect_equal(nrow(splits_df$Subset2), 10)
    for(i in seq_len(10)){
        t1 <- table(phyloseq::sample_data(ps_unbalanced)[splits_df$Subset1[i,],
                                                         "HMP_BODY_SUBSITE"])
        t2 <- table(phyloseq::sample_data(ps_unbalanced)[splits_df$Subset2[i,],
                                                         "HMP_BODY_SUBSITE"])
        t1_exp <- c(9,9)
        t2_exp <- c(10,10)
        expect_equal(as.vector(t1),t1_exp)
        expect_equal(as.vector(t2),t2_exp)
    }
})

test_that(desc = "data splits: unbalanced design for independent samples
    (starting from unbalanced groups)", code = {
    data("ps_plaque_16S")
    set.seed(123)
    ### Balanced design for independent samples starting from an unbalanced
    # Remove some sample to create an unbalanced design
    ps_unbalanced <- phyloseq::subset_samples(ps_plaque_16S, c(FALSE,TRUE,TRUE))
    # Unbalanced design from unbalanced experiment
    splits_df <- expect_warning(
        createSplits(object = ps_unbalanced, varName = "HMP_BODY_SUBSITE", 
            balanced = FALSE, N = 10))
    expect_equal(nrow(splits_df$Subset1), 10)
    expect_equal(nrow(splits_df$Subset2), 10)
    for(i in seq_len(10)){
        t1 <- table(phyloseq::sample_data(ps_unbalanced)[splits_df$Subset1[i,],
                                                         "HMP_BODY_SUBSITE"])
        t2 <- table(phyloseq::sample_data(ps_unbalanced)[splits_df$Subset2[i,],
                                                         "HMP_BODY_SUBSITE"])
        t1_exp <- c(9,10)
        t2_exp <- c(10,11)
        expect_equal(as.vector(t1),t1_exp)
        expect_equal(as.vector(t2),t2_exp)
    }
})

test_that(desc = "data splits: balanced equals to unbalanced design for 
    repeated measures (starting from unbalanced groups)", code = {
    data("ps_plaque_16S")
    ### Balanced design for independent samples starting from an unbalanced
    # Remove some sample to create an unbalanced design
    ps_unbalanced <- phyloseq::subset_samples(ps_plaque_16S, c(FALSE,TRUE,TRUE))
    # Unbalanced design from unbalanced experiment
    set.seed(123)
    splits_df_unb <- expect_warning(suppressMessages(
        createSplits(object = ps_unbalanced, 
        varName = "HMP_BODY_SUBSITE", paired = "RSID", balanced = FALSE, 
        N = 10)))
    set.seed(123)
    splits_df_b <- expect_warning(suppressMessages(
        createSplits(object = ps_unbalanced, 
        varName = "HMP_BODY_SUBSITE", paired = "RSID", balanced = TRUE, 
        N = 10)))
    expect_equal(splits_df_unb, splits_df_b)
    expect_equal(nrow(splits_df_b$Subset1), 10)
    expect_equal(nrow(splits_df_b$Subset2), 10)
    # Check if each comparison has the right number of samples
    for(i in seq_len(10)){
        # each subset should have 5 paired samples = 10 total samples
        t1 <- table(
            phyloseq::sample_data(ps_unbalanced)[splits_df_unb$Subset1[i,],
            "HMP_BODY_SUBSITE"])
        t2 <- table(
            phyloseq::sample_data(ps_unbalanced)[splits_df_unb$Subset2[i,],
            "HMP_BODY_SUBSITE"])
        t1_exp <- c(5,5)
        t2_exp <- c(5,5)
        expect_equal(as.vector(t1),t1_exp)
        expect_equal(as.vector(t2),t2_exp)
        # each subset should have 2 samples for each RSID
        t1_paired <- table(phyloseq::sample_data(ps_plaque_16S)[
            splits_df_unb$Subset1[i,],"RSID"])>1
        t2_paired <- table(phyloseq::sample_data(ps_plaque_16S)[
            splits_df_unb$Subset2[i,],"RSID"])>1
        expect_equal(as.vector(t1_paired), rep(TRUE, 5))
        expect_equal(as.vector(t2_paired), rep(TRUE, 5))
    }
})

test_that(desc = "data splits: more than 2 groups", code = {
    # 60 subjects
    subjectID <- rep(1:61)
    # 3 observations for each subject
    time <- c("T0", "T1", "T2")
    # subjects belongs to a group
    metadata <- data.frame(
        "subjectID" = rep(subjectID, 3),
        "time" = rep(time, each = 61), 
        "group" = rep(rep(x = c("A", "B", "C"), times = c(11,20,30)), 3), 
        stringsAsFactors = TRUE)
    # Create the phyloseq object
    counts <- matrix(rpois(n = 1000*183, lambda = 5), nrow = 1000, ncol = 183)
    rownames(metadata) <- colnames(counts) <- paste0(
        metadata$time, "_",
        metadata$subjectID, "_",
        metadata$group)
    ps <- phyloseq::phyloseq(
        phyloseq::otu_table(counts, taxa_are_rows = TRUE),
        phyloseq::sample_data(metadata)
    )
    
    # Unbalanced, unpaired
    splits <- createSplits(object = ps, assay_name = "counts", 
        varName = "group", paired = NULL, balanced = FALSE, N = 1)
    t_splits <- table(metadata[splits$Subset1, "group"])
    expect_equal(as.vector(t_splits), c(16,30,45))
    # Balanced, unpaired
    splits_b <- createSplits(object = ps, assay_name = "counts", 
        varName = "group", paired = NULL, balanced = TRUE, N = 1)
    t_splits_b <- table(metadata[splits_b$Subset1, "group"])
    expect_equal(as.vector(t_splits_b), c(16,16,16))
    # Unbalanced, paired
    splits_p <- createSplits(object = ps, assay_name = "counts", 
        varName = "group", paired = "subjectID", balanced = FALSE, N = 1)
    t_splits_p <- table(metadata[splits_p$Subset1, "group"])
    expect_equal(as.vector(t_splits_p), c(15,30,45))
    # Balanced, paired
    splits_b_p <- createSplits(object = ps, assay_name = "counts", 
        varName = "group", paired = "subjectID", balanced = TRUE, N = 1)
    t_splits_p_b <- table(metadata[splits_b_p$Subset1, "group"])
    expect_equal(as.vector(t_splits_p_b), c(15,15,15))
})
