# Unit tests for getStatistics and getDA functions

test_that("get statistics: get p-values, no direction", code = {
    # load data
    data("ps_plaque_16S")
    # set normalization
    my_norm <- setNormalizations("norm_edgeR", "TMM")
    ps <- runNormalizations(my_norm, ps_plaque_16S)
    # set method
    my_method <- set_limma(design = ~ 1 + HMP_BODY_SUBSITE, coef = 2,
                           norm = "TMM")
    # Run
    da <- runDA(my_method, object = ps)
    # get p-values, no direction
    t1 <- getStatistics(method = da[[1]], slot = "pValMat", colName = "rawP",
        type = "pvalue", direction = NULL)
    expect_true(is.vector(t1))
})

test_that("get DA/non-DA based on p-value threshold", code = {
    # load data
    data("ps_plaque_16S")
    # set normalization
    my_norm <- setNormalizations("norm_edgeR", "TMM")
    ps <- runNormalizations(my_norm, ps_plaque_16S)
    # set method
    my_method <- set_limma(design = ~ 1 + HMP_BODY_SUBSITE, coef = 2,
        norm = "TMM")
    # Run
    da <- runDA(my_method, object = ps)
    # get DA or non-DA based on p-value threshold
    t1 <- getDA(method = da[[1]], slot = "pValMat", colName = "rawP", type =
        "pvalue", direction = NULL, threshold_pvalue = 0.1, threshold_logfc = 0,
        top = NULL)
    expect_true(is.data.frame(t1))
    expect_true(ncol(t1) == 2)
    expect_true(length(unique(t1[, "DA"])) == 2)
})

test_that("get DA/non-DA based on top 10 values", code = {
    # load data
    data("ps_plaque_16S")
    # set normalization
    my_norm <- setNormalizations("norm_edgeR", "TMM")
    ps <- runNormalizations(my_norm, ps_plaque_16S)
    # set method
    my_method <- set_limma(design = ~ 1 + HMP_BODY_SUBSITE, coef = 2,
                           norm = "TMM")
    # Run
    da <- runDA(my_method, object = ps)
    # get DA or non-DA based on top 10 values
    t1 <- getDA(method = da[[1]], slot = "pValMat", colName = "rawP", type =
        "pvalue", direction = NULL, threshold_pvalue = 1, threshold_logfc = 0,
        top = 10)
    expect_true(is.data.frame(t1))
    expect_true(ncol(t1) == 2)
    expect_true(sum(t1$DA == "DA") == 10)
})

test_that("get p-values and logFC", code = {
    # load data
    data("ps_plaque_16S")
    # set normalization
    my_norm <- setNormalizations("norm_edgeR", "TMM")
    ps <- runNormalizations(my_norm, ps_plaque_16S)
    # set method
    my_method <- set_limma(design = ~ 1 + HMP_BODY_SUBSITE, coef = 2,
                           norm = "TMM")
    # Run
    da <- runDA(my_method, object = ps)
    # get p-values and logFC
    t1 <- getStatistics(method = da[[1]], slot = "pValMat", colName = "rawP",
        type = "pvalue", direction = "logFC")
    expect_true(is.data.frame(t1))
    expect_true(ncol(t1) == 2)
})

test_that("get non-DA/UP/DOWN-Abundant based on p-value threshold", code = {
    # load data
    data("ps_plaque_16S")
    # set normalization
    my_norm <- setNormalizations("norm_edgeR", "TMM")
    ps <- runNormalizations(my_norm, ps_plaque_16S)
    # set method
    my_method <- set_limma(design = ~ 1 + HMP_BODY_SUBSITE, coef = 2,
                           norm = "TMM")
    # Run
    da <- runDA(my_method, object = ps)
    # get non-DA, UP or DOWN-Abundant based on p-value threshold
    t1 <- getDA(method = da[[1]], slot = "pValMat", colName = "rawP", type =
        "pvalue", direction = "logFC", threshold_pvalue = 0.1,
        threshold_logfc = 0, top = NULL)
    expect_true(is.data.frame(t1))
    expect_true(ncol(t1) == 3)
    expect_true(length(unique(t1[, "DA"])) == 3)
})

test_that("get non-DA/UP/DOWN-Abundant based on top 10 values", code = {
    # load data
    data("ps_plaque_16S")
    # set normalization
    my_norm <- setNormalizations("norm_edgeR", "TMM")
    ps <- runNormalizations(my_norm, ps_plaque_16S)
    # set method
    my_method <- set_limma(design = ~ 1 + HMP_BODY_SUBSITE, coef = 2,
                           norm = "TMM")
    # Run
    da <- runDA(my_method, object = ps)
    # get non-DA, UP or DOWN-Abundant based on top 10 values
    t1 <- getDA(method = da[[1]], slot = "pValMat", colName = "rawP", type =
        "pvalue", direction = "logFC", threshold_pvalue = 1,
        threshold_logfc = 0, top = 10)
    expect_true(is.data.frame(t1))
    expect_true(ncol(t1) == 3)
    expect_true(sum(t1$DA != "non-DA") == 10)
})

# Unit tests based on extractDA

test_that("Extract DA features from a list of methods: using top", code = {
    # load data
    data("ps_plaque_16S")
    # set normalization
    my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
        method = c("TMM", "CSS"))
    ps <- runNormalizations(my_norm, ps_plaque_16S)
    # set method
    my_method <- suppressWarnings(set_limma(design = ~ 1 + HMP_BODY_SUBSITE,
        coef = 2, norm = c("TMM", "CSS")))
    # Run
    da <- expect_warning(runDA(my_method, object = ps))
    # Top 10 features (ordered by 'direction') are DA
    t1 = extractDA(object = da, slot = "pValMat", colName = "adjP",
        type = c("pvalue","pvalue"), direction = "logFC", threshold_pvalue = 1,
        threshold_logfc = 0, top = 10)
    expect_equal(class(t1), "list")
    expect_equal(ncol(t1[[2]]), 3)
    expect_equal(length(unique(t1[[2]]$DA)), 3)
})

test_that("Extract DA features from a list of methods: using threshold_pvalue
    and threshold_logfc", code = {
    # load data
    data("ps_plaque_16S")
    # set normalization
    my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
        method = c("TMM", "CSS"))
    ps <- runNormalizations(my_norm, ps_plaque_16S)
    # set method
    my_method <- suppressWarnings(set_limma(design = ~ 1 + HMP_BODY_SUBSITE,
        coef = 2, norm = c("TMM", "CSS")))
    # Run
    da <- expect_warning(runDA(my_method, object = ps))
    # All features with p-value < 0.2 and |logFC| > 0.3 are DA
    t1 = extractDA(object = da, slot = "pValMat", colName = "adjP",
        type = "pvalue", direction = "logFC", threshold_pvalue = 0.2,
        threshold_logfc = 0.3, top = NULL)
    expect_equal(class(t1), "list")
    expect_equal(ncol(t1[[2]]), 3)
    expect_equal(length(unique(t1[[2]]$DA)), 3)
})

test_that("Extract DA features from a list of methods: using threshold_pvalue",
    code = {
    # load data
    data("ps_plaque_16S")
    # set normalization
    my_norm <- setNormalizations(fun = c("norm_edgeR", "norm_CSS"),
                                 method = c("TMM", "CSS"))
    ps <- runNormalizations(my_norm, ps_plaque_16S)
    # set method
    my_method <- suppressWarnings(set_limma(design = ~ 1 + HMP_BODY_SUBSITE,
        coef = 2, norm = c("TMM", "CSS")))
    # Run
    da <- expect_warning(runDA(my_method, object = ps))
    # Feature with p-value < 0.1 are DA (no info about direction)
    t1 = extractDA(object = da, slot = "pValMat", colName = "adjP",
        type = "pvalue", direction = NULL, threshold_pvalue = 0.1,
        threshold_logfc = 0, top = NULL)
    expect_equal(class(t1), "list")
    expect_equal(ncol(t1[[2]]), 2)
    expect_equal(length(unique(t1[[2]]$DA)), 2)
})





