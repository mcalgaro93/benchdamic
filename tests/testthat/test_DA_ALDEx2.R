# Unit tests for DA ALDEx2 variants 

test_that("DA_ALDEx2: pValMat and statInfo for all tests", code = {
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
    # T 
    ## all, unpaired
    da_t <- expect_warning(DA_ALDEx2(ps, design = "group", mc.samples = 1, 
        test = "t", denom = "all", contrast = c("group", "grp2", "grp1"), 
        paired.test = FALSE))
    expectations(da_t, name = "ALDEx2.all.t.unpaired")
    ## iqlr, paired
    da_t <- expect_warning(DA_ALDEx2(ps, design = "group", mc.samples = 1, 
        test = "t", denom = "iqlr", contrast = c("group", "grp2", "grp1"), 
        paired.test = TRUE))
    expectations(da_t, name = "ALDEx2.iqlr.t.paired")
    ## zero, paired
    da_t <- expect_warning(DA_ALDEx2(ps, design = "group", mc.samples = 2, test = "t", 
        denom = "zero", contrast = c("group", "grp2", "grp1"), 
        paired.test = TRUE))
    expectations(da_t, name = "ALDEx2.zero.t.paired")
    # Wilcox
    ## median, unpaired
    da_w <- expect_warning(DA_ALDEx2(ps, design = "group", mc.samples = 1, 
        test = "wilcox", denom = "median", 
        contrast = c("group", "grp2", "grp1"), paired.test = FALSE))
    expectations(da_w, name = "ALDEx2.median.wilcox.unpaired")
    ## user denom, paired
    da_w <- expect_warning(DA_ALDEx2(ps, design = "group", mc.samples = 1, 
        test = "wilcox", denom = c(1,2,3), 
        contrast = c("group", "grp2", "grp1"), paired.test = TRUE))
    expectations(da_w, name = "ALDEx2.user_denom.wilcox.paired")
    # Kruskal-Wallace
    ## median, kw
    da_kw <- expect_warning(DA_ALDEx2(ps, design = "group", mc.samples = 1, 
        test = "kw", denom = "median", contrast = c("group", "grp2", "grp1"), 
        paired.test = FALSE))
    expectations(da_kw, name = "ALDEx2.median.kw")
    ## all, glm ANOVA, design character
    da_glm <- expect_warning(DA_ALDEx2(ps, design = "~ group", mc.samples = 10, 
        test = "glm", denom = "all", contrast = c("group", "grp2", "grp1")))
    expectations(da_glm, name = "ALDEx2.all.glm")
    
    # New group with 3 levels
    group <- sample(x = c("grp1", "grp2", "grp3"), 
                    size = phyloseq::nsamples(ps),
                    replace = TRUE)
    split <- sample(x = c("split1", "split2"), 
                    size = phyloseq::nsamples(ps),
                    replace = TRUE)
    phyloseq::sample_data(ps)[, "group"] <- as.factor(group)
    phyloseq::sample_data(ps)[, "split"] <- as.factor(split)
    # kw
    ## median, kw
    da_kw <- expect_warning(DA_ALDEx2(ps, design = "group", mc.samples = 1, 
        test = "kw", denom = "median", paired.test = FALSE))
    expectations(da_kw, name = "ALDEx2.median.kw")
    ## median, glm ANOVA like
    da_kw <- expect_warning(DA_ALDEx2(ps, design = "group", mc.samples = 1, 
        test = "kw_glm", denom = "median", paired.test = FALSE))
    expectations(da_kw, name = "ALDEx2.median.kw_glm")
    # glm
    ## all composite design
    da_glm <- expect_warning(DA_ALDEx2(ps, design = "~ split + group", 
        mc.samples = 1, test = "glm", denom = "all", paired.test = FALSE, 
        contrast = c("split", "split2", "split1")))
    expectations(da_glm, name = "ALDEx2.all.glm")
    
    # some errors
    # variable name instead of design if test = "glm"
    expect_error(DA_ALDEx2(ps, design = "split", mc.samples = 1, 
        test = "glm", denom = "all", paired.test = FALSE, 
        contrast = c("split", "split2", "split1")))
    # design instead of variable name if test != "glm"
    expect_error(da_t <- DA_ALDEx2(ps, design = "~ group", mc.samples = 1, 
        test = "t", denom = "all", contrast = c("group", "grp2", "grp1"), 
        paired.test = FALSE))
})
