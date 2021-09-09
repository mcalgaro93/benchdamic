# Unit test based on createSplits

test_that(desc = "data splits: design for repeated measures", code = {
    data("ps_plaque_16S")
    set.seed(123)
    # Warning for non factor varName variable
    expect_warning(splits_df <- createSplits(object = ps_plaque_16S, varName =
        "HMP_BODY_SUBSITE", paired = "RSID", balanced = TRUE, N = 10))
    # Make it a factor
    phyloseq::sample_data(ps_plaque_16S)[,"HMP_BODY_SUBSITE"] <-
        factor(unlist(phyloseq::sample_data(ps_plaque_16S)
        [,"HMP_BODY_SUBSITE"]))
    ### Design for repeated measures
    splits_df <- createSplits(object = ps_plaque_16S, varName =
                                  "HMP_BODY_SUBSITE", paired = "RSID", N = 10)
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
    # Warning for non factor varName variable
    expect_warning(splits_df <- createSplits(object = ps_plaque_16S, varName =
        "HMP_BODY_SUBSITE", paired = "RSID", balanced = TRUE, N = 10))
    # Make it a factor
    phyloseq::sample_data(ps_plaque_16S)[,"HMP_BODY_SUBSITE"] <-
        factor(unlist(phyloseq::sample_data(ps_plaque_16S)
                      [,"HMP_BODY_SUBSITE"]))
    ### Balanced design for independent samples starting from an unbalanced
    # Remove some sample to create an unbalanced design
    ps_unbalanced <- phyloseq::subset_samples(ps_plaque_16S, c(FALSE,TRUE,TRUE))
    # table(phyloseq::sample_data(ps_unbalanced)[,"HMP_BODY_SUBSITE"])
    splits_df <- createSplits(object = ps_unbalanced, varName =
                                  "HMP_BODY_SUBSITE", balanced = TRUE, N = 10)
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
    # Warning for non factor varName variable
    expect_warning(splits_df <- createSplits(object = ps_plaque_16S, varName =
        "HMP_BODY_SUBSITE", paired = "RSID", balanced = TRUE, N = 10))
    # Make it a factor
    phyloseq::sample_data(ps_plaque_16S)[,"HMP_BODY_SUBSITE"] <-
        factor(unlist(phyloseq::sample_data(ps_plaque_16S)
                      [,"HMP_BODY_SUBSITE"]))
    ### Balanced design for independent samples starting from an unbalanced
    # Remove some sample to create an unbalanced design
    ps_unbalanced <- phyloseq::subset_samples(ps_plaque_16S, c(FALSE,TRUE,TRUE))
    # Unbalanced design from unbalanced experiment
    splits_df <- createSplits(object = ps_unbalanced, varName =
                                  "HMP_BODY_SUBSITE", balanced = FALSE, N = 10)
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
