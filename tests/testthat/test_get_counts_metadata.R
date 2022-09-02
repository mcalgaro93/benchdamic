# Unit tests for input phyloseq and TreeSummarizedExperiment

test_that("Correct extraction of counts and sample data from phyloseq", 
    code = {
    data("ps_stool_16S")
    ps <- ps_plaque_16S
    
    expectations <- function(counts, metadata){
        expect_equal(nrow(counts), phyloseq::ntaxa(ps))
        expect_equal(rownames(counts), phyloseq::taxa_names(ps))
        expect_equal(rownames(metadata), phyloseq::sample_names(ps))
    }
    
    counts_and_metadata <- get_counts_metadata(ps)
    expectations(counts_and_metadata$counts, counts_and_metadata$metadata)
})

test_that("Correct extraction of counts and sample data from t(phyloseq)", 
    code = {
        data("ps_stool_16S")
        ps <- ps_plaque_16S
        # Fit model on the matrix of counts
        otu_tab <- data.frame(phyloseq::otu_table(ps))
        sam_tab <- data.frame(phyloseq::sample_data(ps))
        rownames(sam_tab) <- colnames(otu_tab)
        psT <- phyloseq::phyloseq(phyloseq::otu_table(t(otu_tab),
            taxa_are_rows = FALSE), phyloseq::sample_data(sam_tab))
              
        expectations <- function(counts, metadata){
            expect_equal(nrow(counts), phyloseq::ntaxa(ps))
            expect_equal(rownames(counts), phyloseq::taxa_names(ps))
            expect_equal(rownames(metadata), phyloseq::sample_names(ps))
        }
              
        counts_and_metadata <- get_counts_metadata(ps)
        expectations(counts_and_metadata$counts, counts_and_metadata$metadata)
})

test_that("Correct extraction of counts and sample data from 
    TreeSummarizedExperiment", code = {
    # the count table
    otu_tab <- matrix(rpois(100, 50), nrow = 10)
    rownames(otu_tab) <- paste("F", 1:10, sep = "_")
    colnames(otu_tab) <- paste("S", 1:10, sep = "_")
    # The sample information
    sam_tab <- data.frame(condition = rep(c("control", "trt"), each = 5),
        gender = sample(x = 1:2, size = 10, replace = TRUE))
    rownames(sam_tab) <- colnames(otu_tab)
    
    tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
        assays = list("counts" = otu_tab),
        colData = sam_tab)
    
    expectations <- function(counts, metadata){
        expect_equal(rownames(counts), paste("F", 1:10, sep = "_"))
        expect_equal(rownames(metadata), paste("S", 1:10, sep = "_"))
    }
  
    counts_and_metadata <- get_counts_metadata(tse)
    expectations(counts_and_metadata$counts, counts_and_metadata$metadata)
})
