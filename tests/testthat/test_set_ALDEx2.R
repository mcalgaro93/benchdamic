# Unit tests for set_ALDEx2 variants 

test_that("set_ALDEx2: number of instructions with paired test", code = {
    all_combs <- set_ALDEx2(design = "group", pseudo_count = c(TRUE, FALSE),
        mc.samples = 2, paired.test = c(TRUE, FALSE), 
        denom = c("all", "iqlr", "median", "lvha", "zero"), 
        contrast = c("group", "grp1", "grp2"), 
        test = c("t", "wilcox", "kw", "kw_glm"), expand = TRUE)
    expect_length(all_combs, 60)
})