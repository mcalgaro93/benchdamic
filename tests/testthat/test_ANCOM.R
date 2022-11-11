# Unit tests for DA ANCOM variants 

test_that("DA_ANCOM: pValMat and statInfo for single and multi core", code = {
    data("ps_stool_16S")
    ps <- ps_plaque_16S
    set.seed(543)
    group <- sample(x = c("grp1","grp2"), size = phyloseq::nsamples(ps),
                    replace = TRUE)
    phyloseq::sample_data(ps)[, "group"] <- as.factor(group)
    
    expectations <- function(da, name){
        expect_true(grepl("pValMat",names(da)[1]))
        expect_equal(nrow(da[["pValMat"]]), phyloseq::ntaxa(ps))
        expect_equal(rownames(da[["pValMat"]]), phyloseq::taxa_names(ps))
        expect_true(grepl("statInfo",names(da)[2]))
        expect_true(grepl("name",names(da)[length(names(da))]))
        expect_equal(name, da[["name"]])
    }
    # single core
    da_single <- expect_warning(DA_ANCOM(
            object = ps, 
            fix_formula = "group", 
            p_adj_method = "BH", 
            BC = TRUE, 
            n_cl = 1, 
            contrast = c("group", "grp1", "grp2")
        ))
    expectations(da_single, name = "ANCOM.BC")
    
    # Multi-core
    da_multi <- DA_ANCOM(
        object = ps, 
        fix_formula = "group", 
        p_adj_method = "BH", 
        BC = TRUE, 
        n_cl = 2, 
        verbose = FALSE,
        contrast = c("group", "grp1", "grp2")
    )
    expectations(da_multi, name = "ANCOM.BC")
    
    # Check if ANCOM can be parallelized
    # out_bplapply <- BiocParallel::bplapply(
    #     list("ps1" = ps, "ps2" = ps), function(ps_obj){
    #     da_single <- DA_ANCOM(
    #         object = ps_obj, 
    #         fix_formula = "group", 
    #         p_adj_method = "BH", 
    #         BC = TRUE, 
    #         n_cl = 1, verbose = FALSE, 
    #         contrast = c("group", "grp1", "grp2")
    #     )
    #     return(da_single)
    # }, BPPARAM = BiocParallel::MulticoreParam(2))
})
