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

## ----fitting, warning=FALSE---------------------------------------------------
GOF_stool_16S <- fitModels(
  object = ps_stool_16S,
  models = c("NB", "ZINB", "DM", "ZIG", "HURDLE"),
  scale_ZIG = c("median", "default"),
  scale_HURDLE = c("median", "default")
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
mock_df <- createMocks(
  nsamples = nsamples(ps_stool_16S),
  N = 10
) # At least 1000 is suggested

## ----normalization------------------------------------------------------------
ps_stool_16S <- norm_edgeR(
  object = ps_stool_16S,
  method = "TMM"
)
ps_stool_16S <- norm_DESeq2(
  object = ps_stool_16S,
  method = "poscounts"
)
ps_stool_16S <- norm_CSS(
  object = ps_stool_16S,
  "median"
)
zinbweights <- weights_ZINB(
  object = ps_stool_16S,
  K = 0,
  design = "~ 1"
)

## ----DA_TIEC, message=FALSE---------------------------------------------------
# Random grouping each time
Stool_16S_mockDA <- apply(X = mock_df, MARGIN = 1, FUN = function(x) {
  # Group assignment
  sample_data(ps_stool_16S)$group <- x
  ### DA analysis ###
  returnList <- list()
  returnList <- within(returnList, {
    da.edger <- DA_edgeR(
      object = ps_stool_16S,
      group = unlist(sample_data(ps_stool_16S)[, "group"]),
      design = as.formula("~ group"),
      coef = 2,
      norm = "TMM"
    )
    da.edger.zinb <- DA_edgeR(
      object = ps_stool_16S,
      group = unlist(sample_data(ps_stool_16S)[, "group"]),
      design = as.formula("~ group"),
      coef = 2,
      norm = "TMM",
      weights = zinbweights
    )
    da.deseq <- DA_DESeq2(
      object = ps_stool_16S,
      design = as.formula("~ group"),
      norm = "poscounts",
      contrast = c("group", "grp2", "grp1")
    )
    da.deseq.zinb <- DA_DESeq2(
      object = ps_stool_16S,
      design = as.formula("~ group"),
      norm = "poscounts",
      contrast = c("group", "grp2", "grp1"),
      weights = zinbweights
    )
    da.limma <- DA_limma(
      object = ps_stool_16S,
      design = ~group,
      coef = 2,
      norm = "TMM"
    )
    da.limma.zinb <- DA_limma(
      object = ps_stool_16S,
      design = ~group,
      coef = 2,
      norm = "TMM",
      weights = zinbweights
    )
    da.limma.css <- DA_limma(
      object = ps_stool_16S,
      design = ~group,
      coef = 2,
      norm = "CSSmedian"
    )
  })
  return(returnList)
})

## ----customExample, eval=FALSE------------------------------------------------
#  DA_yourMethod <- function(parameters) # others
#  {
#    ### your method code ###
#  
#    ### extract important variables ###
#    vector_of_pval <- NA # contains the p-values
#    vector_of_adjusted_pval <- NA # contains the adjusted p-values
#    name_of_your_features <- NA # contains the OTU, or ASV, or other feature names. Usually extracted from the rownames of the count data
#  
#    ### prepare the output ###
#    name <- "write.here.the.name"
#    pValMat <- data.frame(
#      "pval" = vector_of_pval,
#      "adjp" = vector_of_adjusted_pval
#    )
#    rownames(pValMat) <- name_of_your_features # Be sure that your method hasn't changed the order of the features. If it happens, you'll need to re-establish the original order.
#    return(list(
#      "pValMat" = pValMat,
#      "name" = name
#    ))
#  } # END - function: DA_yourMethod

## ----createTIEC---------------------------------------------------------------
TIEC_summary <- createTIEC(Stool_16S_mockDA)

## ----FDRplot------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
set.seed(123)
sample_data(ps_plaque_16S)$HMP_BODY_SUBSITE <- factor(sample_data(ps_plaque_16S)$HMP_BODY_SUBSITE)
split_list <- createSplits(
  object = ps_plaque_16S,
  varName = "HMP_BODY_SUBSITE",
  paired = "RSID",
  balanced = TRUE,
  N = 10
) # At least 100 is suggested

## ---- message=FALSE-----------------------------------------------------------
# lapply -> Subset1 and Subset2
Plaque_16S_splitsDA <- lapply(X = split_list, FUN = function(subset) {
  # apply -> Comparison1, Comparison2, ..., ComparisonN
  apply(X = subset, MARGIN = 1, FUN = function(splits) {
    # Splitting
    # cat("Splitting the samples...\n")
    ps <- prune_samples(sample_names(ps_plaque_16S) %in% splits, ps_plaque_16S)
    # Keep only present taxa
    # cat("Removing not present taxa...\n")
    ps <- filter_taxa(ps, function(x) sum(x > 0) > 0, 1)
    # Adding scaling and normalization factors
    # cat("Computing normalizations...\n")
    ps <- norm_edgeR(object = ps, method = "TMM")
    ps <- norm_DESeq2(object = ps, method = "poscounts")
    ps <- norm_CSS(object = ps, method = "median")
    ### DA analysis ###
    # Uncomment the cat code to receive messages
    # cat("Differential abundance:\n")
    returnList <- list()
    returnList <- within(returnList, {
      # cat("edgeR...\n")
      da.edger <- DA_edgeR(
        object = ps,
        group = ps@sam_data$HMP_BODY_SUBSITE,
        design = ~ 1 + HMP_BODY_SUBSITE,
        coef = 2,
        norm = "TMM"
      )
      # cat("DESeq2...\n")
      da.deseq <- DA_DESeq2(
        object = ps,
        design = ~ 1 + HMP_BODY_SUBSITE,
        norm = "poscounts",
        contrast = c(
          "HMP_BODY_SUBSITE",
          "Subgingival Plaque",
          "Supragingival Plaque"
        )
      )
      # cat("limma...\n")
      da.limma <- DA_limma(
        object = ps,
        design = ~ 1 + HMP_BODY_SUBSITE,
        coef = 2,
        norm = "TMM"
      )
      # cat("limma CSS...\n")
      da.limma.css <- DA_limma(
        object = ps,
        design = ~ 1 + HMP_BODY_SUBSITE,
        coef = 2,
        norm = "CSSmedian"
      )
      # cat("\n")
    })
    return(returnList)
  })
})

## -----------------------------------------------------------------------------
concordance <- createConcordance(object = Plaque_16S_splitsDA, slot = "pValMat", colName = "rawP", type = "pvalue")
head(concordance)

## -----------------------------------------------------------------------------
names(Plaque_16S_splitsDA$Subset1$Comparison1)

## -----------------------------------------------------------------------------
names(Plaque_16S_splitsDA$Subset1$Comparison1$da.limma.css$statInfo)
names(Plaque_16S_splitsDA$Subset1$Comparison1$da.limma$statInfo)
names(Plaque_16S_splitsDA$Subset1$Comparison1$da.deseq$statInfo)
names(Plaque_16S_splitsDA$Subset1$Comparison1$da.edger$statInfo)

## -----------------------------------------------------------------------------
concordance_alternative <- createConcordance(
  object = Plaque_16S_splitsDA,
  slot = "statInfo",
  colName = c("logFC", "logFC", "log2FoldChange", "logFC"),
  type = "logfc"
)

## -----------------------------------------------------------------------------
pC <- plotConcordance(concordance = concordance, threshold = 30, cols = cols)

plot_grid(plotlist = list(pC$concordanceDendrogram, pC$concordanceHeatmap), ncol = 2, align = "h", axis = "tb", rel_widths = c(1, 3))

## -----------------------------------------------------------------------------
data("microbial_metabolism")
head(microbial_metabolism)

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
ps_plaque_16S <- norm_edgeR(
  object = ps_plaque_16S,
  method = "TMM"
)
ps_plaque_16S <- norm_edgeR(
  object = ps_plaque_16S,
  method = "none"
)
ps_plaque_16S <- norm_DESeq2(
  object = ps_plaque_16S,
  method = "poscounts"
)
ps_plaque_16S <- norm_CSS(
  object = ps_plaque_16S,
  "median"
)

## ---- message=FALSE-----------------------------------------------------------
# Convert to factor
ps_plaque_16S@sam_data$HMP_BODY_SUBSITE <- factor(
    ps_plaque_16S@sam_data$HMP_BODY_SUBSITE
)
# Reference level = "Subgingival Plaque"
ps_plaque_16S@sam_data$HMP_BODY_SUBSITE <- relevel(
    x = ps_plaque_16S@sam_data$HMP_BODY_SUBSITE, 
    ref = "Subgingival Plaque"
)
Plaque_16S_DA <- list()
Plaque_16S_DA <- within(Plaque_16S_DA, {
    # DA analysis
    cat("edgeR...\n")
    da.edger <- DA_edgeR(
        object = ps_plaque_16S,
        group = ps_plaque_16S@sam_data$HMP_BODY_SUBSITE,
        design = ~ 1 + HMP_BODY_SUBSITE,
        coef = 2,
        norm = "TMM"
    )
    cat("DESeq2...\n")
    da.deseq <- DA_DESeq2(
        object = ps_plaque_16S,
        design = ~ 1 + HMP_BODY_SUBSITE,
        norm = "poscounts",
        contrast = c(
          "HMP_BODY_SUBSITE",
          "Supragingival Plaque",
          "Subgingival Plaque"
        )
    )
    cat("limma...\n")
    da.limma <- DA_limma(
        object = ps_plaque_16S,
        design = ~ 1 + HMP_BODY_SUBSITE,
        coef = 2,
        norm = "TMM"
    )
    cat("limma CSS...\n")
    da.limma.css <- DA_limma(
        object = ps_plaque_16S,
        design = ~ 1 + HMP_BODY_SUBSITE,
        coef = 2,
        norm = "CSSmedian"
    )
    cat("metagenomeSeq CSS...\n")
    da.metagenomeseq.css <- DA_metagenomeSeq(
        object = ps_plaque_16S,
        design = ~ 1 + HMP_BODY_SUBSITE,
        coef = 2,
        norm = "CSSmedian"
    )
    cat("ALDEx2 CSS...\n")
    da.ALDEx2 <- DA_ALDEx2(
        object = ps_plaque_16S,
        conditions = ps_plaque_16S@sam_data$HMP_BODY_SUBSITE,
        test = "t",
        norm = "none"
    )
    cat("corncob Wald...\n")
    da.corconb.Wald <- DA_corncob(
        object = ps_plaque_16S,
        formula = ~ 1 + HMP_BODY_SUBSITE,
        formula_null = ~ 1,
        phi.formula = ~ 1 + HMP_BODY_SUBSITE,
        phi.formula_null = ~ 1 + HMP_BODY_SUBSITE,
        coefficient = "HMP_BODY_SUBSITESupragingival Plaque",
        test = "Wald",
        norm = "none"
    )
    cat("MAST Wald...\n")
    da.MAST <- DA_MAST(
        object = ps_plaque_16S,
        rescale = "median",
        design = ~ 1 + HMP_BODY_SUBSITE,
        coefficient = "HMP_BODY_SUBSITESupragingival Plaque",
        norm = "none"
    )
    cat("Seurat...\n")
    da.Seurat <- DA_Seurat(
        object = ps_plaque_16S,
        test.use = "wilcox",
        contrast = c("HMP_BODY_SUBSITE", "Supragingival Plaque", "Subgingival Plaque"),
        norm = "none"
    )
    cat("\n")
})

## -----------------------------------------------------------------------------
names(Plaque_16S_DA)

## -----------------------------------------------------------------------------
enrichment <- createEnrichment(
    object = Plaque_16S_DA, priorKnowledge = priorInfo, enrichmentCol = "Type",
    namesCol = "newNames", slot = "pValMat", colName = "adjP", type = "pvalue",
    direction = c("avg_logFC", # Seurat
        "logFC", # MAST
        "Estimate", # corncob
        "effect", # ALDEx2
        "HMP_BODY_SUBSITESupragingival Plaque", # metagenomeSeq
        "logFC", # limma
        "logFC", # limma
        "log2FoldChange", # DEseq2
        "logFC"), # edgeR)
    threshold_pvalue = 0.1,
    threshold_logfc = 0,
    top = NULL,
    alternative = "greater",
    verbose = TRUE
)

## -----------------------------------------------------------------------------
plotContingency(enrichment = enrichment, levels_to_plot = c("Aerobic", "Anaerobic", "F_Anaerobic", "Unknown"), method = "metagenomeSeq.CSSmedian")

## ---- fig.width=7, fig.height=7-----------------------------------------------
plotEnrichment(enrichment = enrichment, enrichmentCol = "Type", levels_to_plot = c("Aerobic", "Anaerobic"))

## ---- fig.width=6, fig.height=6-----------------------------------------------
plotMutualFindings(enrichment, enrichmentCol = "Type", levels_to_plot = c("Aerobic", "Anaerobic"), n_methods = 2)

## -----------------------------------------------------------------------------
positives <- createPositives(
    object = Plaque_16S_DA, priorKnowledge = priorInfo, enrichmentCol = "Type",
    namesCol = "newNames", slot = "pValMat", colName = "rawP", type = "pvalue",
    direction = c("avg_logFC", # Seurat
        "logFC", # MAST
        "Estimate", # corncob
        "effect", # ALDEx2
        "HMP_BODY_SUBSITESupragingival Plaque", # metagenomeSeq
        "logFC", # limma
        "logFC", # limma
        "log2FoldChange", # DEseq2
        "logFC"), # edgeR)
    threshold_pvalue = 1,
    threshold_logfc = 0,
    top = seq.int(from = 0, to = 50, by = 5), 
    alternative = "greater",
    verbose = FALSE,
    TP = list(c("DOWN Abundant", "Anaerobic"), c("UP Abundant", "Aerobic")),
    FP = list(c("DOWN Abundant", "Aerobic"), c("UP Abundant", "Anaerobic"))
)
head(positives)

## -----------------------------------------------------------------------------
plotPositives(positives)

