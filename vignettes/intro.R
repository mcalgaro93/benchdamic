## ----options, include=FALSE, echo=FALSE---------------------------------------
knitr::opts_chunk$set(warning = FALSE, error = FALSE, message = FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  devtools::install_github("mcalgaro93/benchdamic")

## ----load_packs---------------------------------------------------------------
library(benchdamic)
# Data management
library(phyloseq)
library(plyr)
# Graphics
library(ggplot2)
library(cowplot)

## ----dataloading--------------------------------------------------------------
data("ps_stool_16S")
ps_stool_16S

## ----16S_stool_download, eval=FALSE-------------------------------------------
#  ## 16S HMP data download
#  library(HMP16SData)
#  # Extracting V3-V5 16S sequenced regions' count data
#  ps_stool_16S_raw <- V35() %>%
#    subset(select = HMP_BODY_SUBSITE == "Stool" & # Only fecal samples
#      RUN_CENTER == "BI" & # Only sequenced at the BI RUN CENTER
#      SEX == "Male" & # Only male subject
#      VISITNO == 1 & # Only the first visit
#      !duplicated(RSID)) %>% # Duplicated SubjectID removal
#    as_phyloseq()
#  # Remove low depth samples
#  ps_stool_16S_pruned <- prune_samples(sample_sums(ps_stool_16S_raw) >= 10^3, ps_stool_16S_raw)
#  # Remove features with zero counts
#  ps_stool_16S_filtered <- filter_taxa(ps_stool_16S_pruned, function(x) {
#    sum(x > 0) > 0
#  }, 1)
#  # Collapse counts to the genus level
#  ps_stool_16S <- tax_glom(ps_stool_16S_filtered, taxrank = "GENUS")

## -----------------------------------------------------------------------------
example_HURDLE <- fitHURDLE(
  counts = as(otu_table(ps_stool_16S), "matrix"),
  scale = "median"
)
head(example_HURDLE)

## -----------------------------------------------------------------------------
observed_hurdle <- prepareObserved(counts = as(
  otu_table(ps_stool_16S),
  "matrix"
), scale = "median")
head(observed_hurdle)

## -----------------------------------------------------------------------------
head(prepareObserved(counts = as(otu_table(ps_stool_16S), "matrix")))

## -----------------------------------------------------------------------------
head(meanDifferences(
  estimated = example_HURDLE,
  observed = observed_hurdle
))

## ----fitting------------------------------------------------------------------
GOF_stool_16S <- fitModels(
  counts = as(otu_table(ps_stool_16S), "matrix"),
  models = c("NB", "ZINB", "DM", "ZIG", "HURDLE"),
  scale_ZIG = c("median", "default"),
  scale_HURDLE = c("median", "default"),
  verbose = TRUE
)

## ----RMSE_MD------------------------------------------------------------------
RMSE_MD_16S <- plyr::ldply(GOF_stool_16S,
  .fun = function(df) cbind("RMSE" = RMSE(df[, "MD"])),
  .id = "Model"
)
RMSE_MD_16S

## ----RMSE_ZPD-----------------------------------------------------------------
RMSE_ZPD_16S <- plyr::ldply(GOF_stool_16S,
  .fun = function(df) cbind("RMSE" = RMSE(df[, "ZPD"])),
  .id = "Model"
)
RMSE_ZPD_16S

## ----plotGOF_MD, fig.width=15, fig.height=4-----------------------------------
plotMD(
  data = GOF_stool_16S,
  difference = "MD",
  split = TRUE
)

## ----plotGOF_MD_noHurdleDefault, fig.width=15, fig.height=4-------------------
plotMD(
  data = GOF_stool_16S[1:6],
  difference = "MD",
  split = TRUE
)

## ----plotGOF_ZPD, fig.width=15, fig.height=4----------------------------------
plotMD(
  data = GOF_stool_16S[1:6],
  difference = "ZPD",
  split = TRUE
)

## ----plotGOF_MD_collapsed, fig.width=10, fig.height=5-------------------------
plot_grid(plotMD(data = GOF_stool_16S[1:6], difference = "MD", split = FALSE),
  plotMD(data = GOF_stool_16S[1:6], difference = "ZPD", split = FALSE),
  ncol = 2
)

## ----plotGOF_RMSE, fig.width=10,fig.height=7----------------------------------
plot_grid(plotRMSE(GOF_stool_16S, difference = "MD"),
  plotRMSE(GOF_stool_16S, difference = "ZPD"),
  ncol = 2
)

## ----createMocks--------------------------------------------------------------
set.seed(123)
my_mocks <- createMocks(
  nsamples = phyloseq::nsamples(ps_stool_16S),
  N = 10
) # At least N = 1000 is suggested

## ----normalizationManual, eval=FALSE------------------------------------------
#  ps_stool_16S <- norm_edgeR(
#    object = ps_stool_16S,
#    method = "TMM"
#  )
#  ps_stool_16S <- norm_DESeq2(
#    object = ps_stool_16S,
#    method = "poscounts"
#  )
#  ps_stool_16S <- norm_CSS(
#    object = ps_stool_16S,
#    "median"
#  )

## ----setNormalization---------------------------------------------------------
my_normalizations <- setNormalizations(fun = c("norm_edgeR", "norm_DESeq2", "norm_CSS"), method = c("TMM", "poscounts", "median"))
ps_stool_16S <- runNormalizations(normalization_list = my_normalizations, object = ps_stool_16S, verbose = TRUE)

## ----weights------------------------------------------------------------------
zinbweights <- weights_ZINB(
  object = ps_stool_16S,
  K = 0,
  design = "~ 1"
)

## ----set_Methods--------------------------------------------------------------
my_edgeR <- set_edgeR(
    pseudo_count = FALSE,
    group_name = "group", 
    design = ~ group, 
    robust = FALSE, 
    coef = 2, 
    norm = "TMM", 
    weights_logical = c(TRUE, FALSE), 
    expand = TRUE
)
my_DESeq2 <- set_DESeq2(
    pseudo_count = FALSE,
    design = ~ group,
    contrast = c("group", "grp2", "grp1"),
    norm = "poscounts", 
    weights_logical = c(TRUE, FALSE), 
    alpha = 0.05,
    expand = TRUE
)
my_limma <- set_limma(
    pseudo_count = FALSE,
    design = ~ group,
    coef = 2,
    norm = c("TMM", "TMM", "CSSmedian"),
    weights_logical = c(FALSE, TRUE, FALSE), 
    expand = FALSE)

my_methods <- c(my_edgeR, my_DESeq2, my_limma)

## ----runMocks-----------------------------------------------------------------
# Random grouping each time
Stool_16S_mockDA <- runMocks(mocks = my_mocks, method_list = my_methods, object = ps_stool_16S, weights = zinbweights, verbose = FALSE)

## ----customExample, eval=FALSE------------------------------------------------
#  DA_yourMethod <- function(object, parameters) # others
#  {
#      ### your method code ###
#  
#      ### extract important statistics ###
#      vector_of_pval <- NA # contains the p-values
#      vector_of_adjusted_pval <- NA # contains the adjusted p-values
#      name_of_your_features <- NA # contains the OTU, or ASV, or other feature
#                                  # names. Usually extracted from the rownames of
#                                  # the count data
#      vector_of_logFC <- NA # contains the logFCs
#      vector_of_statistics <- NA # contains other statistics
#  
#      ### prepare the output ###
#      pValMat <- data.frame("rawP" = vector_of_pval,
#                            "adjP" = vector_of_adjusted_pval)
#      statInfo <- data.frame("logFC" = vector_of_logFC,
#                             "statistics" = vector_of_statistics)
#      name <- "write.here.the.name"
#      # Be sure that your method hasn't changed the order of the features. If it
#      # happens, you'll need to re-establish the original order.
#      rownames(pValMat) <- rownames(statInfo) <- name_of_your_features
#  
#      # Return the output as a list
#      return(list("pValMat" = pValMat, "statInfo" = statInfo, "name" = name))
#  } # END - function: DA_yourMethod

## ----customExampleInstances, eval=FALSE---------------------------------------
#  my_custom_method <- list(
#      customMethod.1 = list(method = "DA_yourMethod", parameters),
#      customMethod.2 = list(method = "DA_yourMethod", parameters)
#  )

## ----customExampleRun, eval=FALSE---------------------------------------------
#  # Add the custom method instances to the others
#  my_methods <- c(my_edgeR, my_DESeq2, my_limma, my_custom_method)
#  # Run all the methods on a specific data object...
#  runDA(my_methods = my_methods, object = dataObject)
#  # ... Or on the mock datasets to investigate TIEC
#  runMocks(mocks = mock_df, my_methods = my_methods, object = dataObject)

## ----createTIEC---------------------------------------------------------------
TIEC_summary <- createTIEC(Stool_16S_mockDA)

## ----FPRplot------------------------------------------------------------------
cols <- createColors(variable = levels(TIEC_summary$df_pval$Method))
plotFPR(df_FPR = TIEC_summary$df_FPR, cols = cols)

## ----QQplot-------------------------------------------------------------------
plotQQ(df_QQ = TIEC_summary$df_QQ, zoom = c(0, 0.1), cols = cols)

## ----KSplot-------------------------------------------------------------------
plotKS(df_KS = TIEC_summary$df_KS, cols = cols)

## ----dataloading_concordance--------------------------------------------------
data("ps_plaque_16S")

## ----16S_plaque_download, eval=FALSE------------------------------------------
#  ps_plaque_16S_raw <- V35() %>% # Extracting V3-V5 16S sequenced regions' count data
#    subset(select = HMP_BODY_SUBSITE %in% c("Supragingival Plaque", "Subgingival Plaque") & # Only gingival plaque samples
#      RUN_CENTER == "WUGC" & # Only sequenced at the WUCG RUN CENTER
#      SEX == "Male" & # Only male subject
#      VISITNO == 1) %>% # Only the first visit
#    as_phyloseq()
#  
#  # Only paired samples
#  paired <- names(which(table(sample_data(ps_plaque_16S_raw)[, "RSID"]) == 2))
#  ps_plaque_16S_paired <- subset_samples(ps_plaque_16S_raw, RSID %in% paired)
#  
#  # Remove low depth samples
#  ps_plaque_16S_pruned <- prune_samples(sample_sums(ps_plaque_16S_paired) >= 10^3, ps_plaque_16S_paired)
#  
#  # Remove features with zero counts
#  ps_plaque_16S_filtered <- filter_taxa(ps_plaque_16S_pruned, function(x) sum(x > 0) > 0, 1)
#  
#  # Collapse counts to the genus level
#  ps_plaque_16S <- tax_glom(ps_plaque_16S_filtered, taxrank = "GENUS")

## ----createSplits-------------------------------------------------------------
set.seed(123)
sample_data(ps_plaque_16S)$HMP_BODY_SUBSITE <- factor(sample_data(ps_plaque_16S)$HMP_BODY_SUBSITE)

my_splits <- createSplits(
  object = ps_plaque_16S,
  varName = "HMP_BODY_SUBSITE",
  paired = "RSID",
  balanced = TRUE,
  N = 10
) # At least 100 is suggested

## ----set_Methods_noweights----------------------------------------------------
my_edgeR_noWeights <- set_edgeR(group_name = "HMP_BODY_SUBSITE", design = ~ HMP_BODY_SUBSITE, coef = 2, norm = "TMM")
my_DESeq2_noWeights <- set_DESeq2(contrast = c("HMP_BODY_SUBSITE","Supragingival Plaque", "Subgingival Plaque"), design = ~ HMP_BODY_SUBSITE, norm = "poscounts")
my_limma_noWeights <- set_limma(design = ~ HMP_BODY_SUBSITE, coef = 2, norm = c("TMM", "CSSmedian"))

my_methods_noWeights <- c(my_edgeR_noWeights, my_DESeq2_noWeights, my_limma_noWeights)

## ----info_normalizations------------------------------------------------------
str(my_normalizations)

## ----runSplits----------------------------------------------------------------
Plaque_16S_splitsDA <- runSplits(split_list = my_splits, method_list = my_methods_noWeights, normalization_list = my_normalizations, object = ps_plaque_16S, verbose = FALSE)

## ----createConcordance--------------------------------------------------------
concordance <- createConcordance(
    object = Plaque_16S_splitsDA, 
    slot = "pValMat", 
    colName = "rawP", 
    type = "pvalue"
)
head(concordance)

## ----getNames-----------------------------------------------------------------
names(Plaque_16S_splitsDA$Subset1$Comparison1)

## ----logFC_names--------------------------------------------------------------
names(Plaque_16S_splitsDA$Subset1$Comparison1$edgeR.TMM$statInfo)
names(Plaque_16S_splitsDA$Subset1$Comparison1$DESeq2.poscounts$statInfo)
names(Plaque_16S_splitsDA$Subset1$Comparison1$limma.CSSmedian$statInfo)
names(Plaque_16S_splitsDA$Subset1$Comparison1$limma.TMM$statInfo)

## ----alternativeConcordance, eval=FALSE---------------------------------------
#  concordance_alternative <- createConcordance(
#    object = Plaque_16S_splitsDA,
#    slot = "statInfo",
#    colName = c("logFC", "log2FoldChange", "logFC", "logFC"),
#    type = "logfc"
#  )

## ----plotConcordance----------------------------------------------------------
pC <- plotConcordance(concordance = concordance, threshold = 30)
cowplot::plot_grid(plotlist = pC, ncol = 2, align = "h", axis = "tb",
        rel_widths = c(1, 3))

## ----priorKnowledge-----------------------------------------------------------
data("microbial_metabolism")
head(microbial_metabolism)

## ----exampleOfIntegration-----------------------------------------------------
# Extract genera from the phyloseq tax_table slot
genera <- tax_table(ps_plaque_16S)[, "GENUS"]
# Genera as rownames of microbial_metabolism data.frame
rownames(microbial_metabolism) <- microbial_metabolism$Genus
# Match OTUs to their metabolism
priorInfo <- data.frame(genera, "Type" =  microbial_metabolism[genera, "Type"])
unknown_metabolism <- is.na(priorInfo$Type)
priorInfo[unknown_metabolism, "Type"] <- "Unknown"
# Relabel 'F Anaerobic' to 'F_Anaerobic' to remove space
priorInfo$Type <- factor(priorInfo$Type, levels = c("Aerobic","Anaerobic","F Anaerobic","Unknown"), labels = c("Aerobic","Anaerobic","F_Anaerobic","Unknown"))
# Add a more informative names column
priorInfo[, "newNames"] <- paste0(rownames(priorInfo), "|", 
    priorInfo[, "GENUS"])

## ----setNormalizations_enrichment---------------------------------------------
none_normalization <- setNormalizations(fun = "norm_edgeR", method = "none")
my_normalizations_enrichment <- c(my_normalizations, none_normalization)
ps_plaque_16S <- runNormalizations(normalization_list = my_normalizations_enrichment, object = ps_plaque_16S, verbose = FALSE)

## ----set_Methods_enrichment---------------------------------------------------
my_metagenomeSeq <- set_metagenomeSeq(design = ~ HMP_BODY_SUBSITE, coef = 2, norm = "CSSmedian")
my_ALDEx2 <- set_ALDEx2(conditions = "HMP_BODY_SUBSITE", test = "t", norm = "none")
my_corncob <- set_corncob(formula = ~ HMP_BODY_SUBSITE, phi.formula = ~ HMP_BODY_SUBSITE, formula_null = ~ 1, phi.formula_null = ~ HMP_BODY_SUBSITE, test = "Wald", coefficient = "HMP_BODY_SUBSITESupragingival Plaque", norm = "none")
my_MAST <- set_MAST(rescale = "median", design = ~ HMP_BODY_SUBSITE, coefficient = "HMP_BODY_SUBSITESupragingival Plaque", norm = "none")
my_Seurat <- set_Seurat(test.use = "wilcox", contrast = c("HMP_BODY_SUBSITE", "Supragingival Plaque", "Subgingival Plaque"), norm = "none")

my_methods_enrichment <- c(
    my_methods_noWeights, 
    my_metagenomeSeq, 
    my_ALDEx2, 
    my_corncob, 
    my_MAST, 
    my_Seurat
)

## ----runDA_enrichment, message=FALSE------------------------------------------
# Convert to factor
ps_plaque_16S@sam_data$HMP_BODY_SUBSITE <- factor(
    ps_plaque_16S@sam_data$HMP_BODY_SUBSITE
)
# Reference level = "Subgingival Plaque"
ps_plaque_16S@sam_data$HMP_BODY_SUBSITE <- relevel(
    x = ps_plaque_16S@sam_data$HMP_BODY_SUBSITE, 
    ref = "Subgingival Plaque"
)
Plaque_16S_DA <- runDA(method_list = my_methods_enrichment, object = ps_plaque_16S, weights = zinbweights)

## ----info_DA------------------------------------------------------------------
names(Plaque_16S_DA)

## ----createEnrichment---------------------------------------------------------
enrichment <- createEnrichment(
    object = Plaque_16S_DA, priorKnowledge = priorInfo, enrichmentCol = "Type",
    namesCol = "newNames", slot = "pValMat", colName = "adjP", type = "pvalue",
    direction = c("logFC", # edgeR
        "log2FoldChange", # DEseq2
        "logFC", # limma
        "logFC", # limma
        "HMP_BODY_SUBSITESupragingival Plaque", # metagenomeSeq
        "effect", # ALDEx2
        "Estimate", # corncob
        "logFC",# MAST
        "avg_log2FC"), # Seurat
    threshold_pvalue = 0.1,
    threshold_logfc = 0,
    top = NULL,
    alternative = "greater",
    verbose = TRUE
)

## ----plotContingency----------------------------------------------------------
plotContingency(enrichment = enrichment, levels_to_plot = c("Aerobic", "Anaerobic"), method = "metagenomeSeq.CSSmedian")

## ----plotEnrichment, fig.width=7, fig.height=7--------------------------------
plotEnrichment(enrichment = enrichment, enrichmentCol = "Type", levels_to_plot = c("Aerobic", "Anaerobic"))

## ----plotMutualFindings, fig.width=6, fig.height=6----------------------------
plotMutualFindings(enrichment, enrichmentCol = "Type", levels_to_plot = c("Aerobic", "Anaerobic"), n_methods = 1)

## ----createPositives----------------------------------------------------------
positives <- createPositives(
    object = Plaque_16S_DA, priorKnowledge = priorInfo, enrichmentCol = "Type",
    namesCol = "newNames", slot = "pValMat", colName = "rawP", type = "pvalue",
    direction = c("logFC", # edgeR
        "log2FoldChange", # DEseq2
        "logFC", # limma
        "logFC", # limma
        "HMP_BODY_SUBSITESupragingival Plaque", # metagenomeSeq
        "effect", # ALDEx2
        "Estimate", # corncob
        "logFC",# MAST
        "avg_log2FC"), # Seurat
    threshold_pvalue = 1,
    threshold_logfc = 0,
    top = seq.int(from = 0, to = 50, by = 5), 
    alternative = "greater",
    verbose = FALSE,
    TP = list(c("DOWN Abundant", "Anaerobic"), c("UP Abundant", "Aerobic")),
    FP = list(c("DOWN Abundant", "Aerobic"), c("UP Abundant", "Anaerobic"))
)
head(positives)

## ----plotPositives------------------------------------------------------------
plotPositives(positives)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

