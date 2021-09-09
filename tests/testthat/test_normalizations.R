# Unit tests for normalization functions

test_that("Normalization/Scaling factors are directly computable", code = {
    data("ps_stool_16S")
    ncol_before <- ncol(phyloseq::sample_data(ps_stool_16S))
    # Calculate the scaling factors
    ps_stool_16S <- norm_edgeR(object = ps_stool_16S, method = "TMM")
    ps_stool_16S <- norm_DESeq2(object = ps_stool_16S, method = "poscounts")
    ps_stool_16S <- norm_CSS(object = ps_stool_16S, method = "default")
    ps_stool_16S <- norm_TSS(object = ps_stool_16S)
    expect_equal(ncol(phyloseq::sample_data(ps_stool_16S)), ncol_before + 4)
})

test_that("Normalization/Scaling factors are indirectly computable", code = {
    data("ps_stool_16S")
    # Calculate the scaling factors: indirect way
    my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_DESeq2",
        "norm_CSS", "norm_TSS"), method = c("TMM", "poscounts", "default",
        "TSS"))
    ps_new <- runNormalizations(normalization_list = my_norm, ps_stool_16S)
    # Calculate the scaling factors: direct way
    ps_stool_16S <- norm_edgeR(object = ps_stool_16S, method = "TMM")
    ps_stool_16S <- norm_DESeq2(object = ps_stool_16S, method = "poscounts")
    ps_stool_16S <- norm_CSS(object = ps_stool_16S, method = "default")
    ps_stool_16S <- norm_TSS(object = ps_stool_16S)
    expect_equal(phyloseq::sample_data(ps_stool_16S),
        phyloseq::sample_data(ps_new))
})
