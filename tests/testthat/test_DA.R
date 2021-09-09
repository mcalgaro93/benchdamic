# Unit tests for DA outputs

test_that("DA methods produce pValMat and statInfo", code = {
    data("ps_stool_16S")
    ps <- ps_plaque_16S
    ps <- norm_edgeR(ps, method = "TMM")
    ps <- norm_DESeq2(ps, method = "poscounts")
    ps <- norm_CSS(ps, method = "median")
    ps <- norm_edgeR(ps, "none")

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

    da <- DA_edgeR(ps, group_name = "group", design = ~ 1 + group, coef = 2,
        norm = "TMM")
    expectations(da, name = "edgeR.TMM")
    da <- DA_DESeq2(ps, design = ~ 1 + group, norm = "poscounts", contrast = c(
        "group","grp2","grp1"))
    expectations(da, name = "DESeq2.poscounts")
    da <- DA_limma(ps, design = ~ 1 + group, coef = 2, norm = "TMM")
    expectations(da, name = "limma.TMM")
    da <- DA_metagenomeSeq(ps, design = ~ 1 + group, coef = "groupgrp2", norm =
                               "CSSmedian")
    expectations(da, name = "metagenomeSeq.CSSmedian")
    da <- DA_corncob(ps, formula = ~ 1 + group, formula_null = ~ 1,
                     phi.formula = ~ 1 + group, phi.formula_null = ~ 1 + group,
                     test = "Wald", coefficient = "groupgrp2", norm = "none")
    expectations(da, name = "corncob.none.Wald")
    da <- DA_ALDEx2(ps, conditions = group, test = "t", denom = "iqlr", norm =
                        "none")
    expectations(da, name = "ALDEx2.none.iqlr.t")
    da <- DA_MAST(ps, rescale = "median",design = ~ 1 + group, coefficient =
                      "groupgrp2", norm = "none")
    expectations(da, name = "MAST.none.median")
    da <- DA_Seurat(ps, test.use = "wilcox",
        contrast = c("group", "grp2", "grp1"), norm = "none")
    expectations(da, name = "Seurat.none.wilcox")
})
