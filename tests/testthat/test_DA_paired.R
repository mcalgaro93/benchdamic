# Unit tests for DA outputs

test_that("mixMC with multilevel produce pValMat and statInfo", code = {
    data("ps_plaque_16S")
    # Rename
    ps <- ps_plaque_16S
    ps <- phyloseq::filter_taxa(physeq = ps, 
        flist = function(x) sum(x > 0) >= 3, prune = TRUE)
    phyloseq::sample_data(ps)[, "HMP_BODY_SUBSITE"] <- as.factor(
        phyloseq::sample_data(ps)[["HMP_BODY_SUBSITE"]])
    
    # We expect several things from each DA method
    expectations <- function(da, name){
        expect_true(grepl("pValMat",names(da)[1]))
        expect_equal(nrow(da[["pValMat"]]), phyloseq::ntaxa(ps))
        expect_equal(rownames(da[["pValMat"]]), phyloseq::taxa_names(ps))
        expect_true(grepl("statInfo",names(da)[2]))
        expect_true(grepl("name",names(da)[length(names(da))]))
        expect_equal(name, da[["name"]])
    }
    # # DA_mixMC
    # da <- DA_mixMC(object = ps, contrast = c("HMP_BODY_SUBSITE", 
    #     "Supragingival Plaque", "Subgingival Plaque"), ID_variable = "RSID",
    #     verbose = FALSE)
    # expectations(da, name = "mixMC.pc1")
})
