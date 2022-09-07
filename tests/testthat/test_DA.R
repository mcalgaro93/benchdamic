# Unit tests for DA outputs

test_that("DA methods produce pValMat and statInfo", code = {
    data("ps_stool_16S")
    # Rename
    ps <- ps_plaque_16S
    # Add some normalization/scaling factors
    ps <- norm_edgeR(ps, method = "TMM")
    ps <- norm_DESeq2(ps, method = "poscounts")
    ps <- norm_CSS(ps, method = "CSS")
    ps <- norm_edgeR(ps, method = "none")
    # Setting a random grouping
    set.seed(123)
    group <- sample(x = c("grp1","grp2"), size = phyloseq::nsamples(ps),
        replace = TRUE)
    phyloseq::sample_data(ps)[, "group"] <- as.factor(group)
    # We expect several things from each DA method
    expectations <- function(da, name){
        expect_true(grepl("pValMat",names(da)[1]))
        expect_equal(nrow(da[["pValMat"]]), phyloseq::ntaxa(ps))
        expect_equal(rownames(da[["pValMat"]]), phyloseq::taxa_names(ps))
        expect_true(grepl("statInfo",names(da)[2]))
        expect_true(grepl("name",names(da)[length(names(da))]))
        expect_equal(name, da[["name"]])
    }
    # DA_edgeR
    da <- DA_edgeR(ps, group_name = "group", design = ~ 1 + group, coef = 2,
        norm = "TMM")
    expectations(da, name = "edgeR.TMM")
    # DA_DESeq2
    da <- DA_DESeq2(ps, design = ~ 1 + group, norm = "poscounts", contrast = c(
        "group","grp2","grp1"))
    expectations(da, name = "DESeq2.poscounts")
    # DA_limma
    da <- DA_limma(ps, design = ~ 1 + group, coef = 2, norm = "TMM")
    expectations(da, name = "limma.TMM")
    # DA_metagenomeSeq
    da <- DA_metagenomeSeq(object = ps, design = ~ 1 + group, 
        coef = "groupgrp2", norm = "CSS", model = "fitFeatureModel", 
        verbose = FALSE)
    expectations(da, name = "metagenomeSeq.CSS.fitFeatureModel")
    # da <- DA_metagenomeSeq(object = ps, design = ~ 1 + group, 
    #     coef = "groupgrp2", norm = "CSS", model = "fitZig", 
    #     verbose = FALSE)
    # expectations(da, name = "metagenomeSeq.CSS.fitZig")
    # # DA_corncob
    # da <- DA_corncob(ps, formula = ~ 1 + group, formula_null = ~ 1,
    #                  phi.formula = ~ 1 + group, phi.formula_null = ~ 1 + group,
    #                  test = "Wald", coefficient = "groupgrp2")
    # expectations(da, name = "corncob.Wald")
    # DA_ALDEx2
    da <- DA_ALDEx2(ps, design = "group", mc.samples = 128, test = "t", 
        paired.test = FALSE, denom = "iqlr", 
        contrast = c("group", "grp2", "grp1"), verbose = FALSE)
    expectations(da, name = "ALDEx2.iqlr.t.unpaired")
    # DA_MAST
    da <- DA_MAST(ps, rescale = "median", design = ~ 1 + group, coefficient =
        "groupgrp2", verbose = FALSE)
    expectations(da, name = "MAST.median")
    # DA_Seurat
    da <- DA_Seurat(ps, norm = "LogNormalize", scale.factor = 10000, 
        test = "wilcox", contrast = c("group", "grp2", "grp1"), verbose = FALSE)
    expectations(da, name = "Seurat.LogNormalize.SF10000.wilcox")
    # DA_ANCOM
    da <- DA_ANCOM(object = ps, formula = "group", 
        contrast = c("group", "grp2", "grp1"), BC = TRUE, verbose = FALSE)
    expectations(da, name = "ANCOM.BC")
    # DA_NOISeq
    da <- DA_NOISeq(object = ps, contrast = c("group", "grp2", "grp1"), norm = 
        "tmm", verbose = FALSE)
    expectations(da, name = "NOISeq.tmm")
    # DA_dearseq
    da <- DA_dearseq(object = ps, variables2test = "group", 
        test = "permutation", preprocessed = FALSE, verbose = FALSE)
    expectations(da, name = "dearseq.permutation.1000")
})