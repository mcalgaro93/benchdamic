# Unit tests based on createMocks

test_that(desc = "mocks are generated as expected", code = {
    data("ps_stool_16S")
    # Generate the patterns for 10 mock comparison for an experiment
    mocks <- createMocks(nsamples = phyloseq::nsamples(ps_stool_16S), N = 10)
    expect_equal(nrow(mocks), 10)
    expect_equal(ncol(mocks), phyloseq::nsamples(ps_stool_16S) %/% 2 * 2)
})
