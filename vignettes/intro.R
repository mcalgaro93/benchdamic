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
#  ps_stool_16S_raw <- V35() %>% # Extracting V3-V5 16S sequenced regions' count data
#    subset(select = HMP_BODY_SUBSITE == "Stool" & # Only fecal samples
#             RUN_CENTER == "BI" & # Only sequenced at the BI RUN CENTER
#             SEX == "Male" & # Only male subject
#             VISITNO == 1 & # Only the first visit
#             !duplicated(RSID)) %>% # Duplicated SubjectID removal
#    as_phyloseq()
#  
#  # Remove low depth samples
#  ps_stool_16S_pruned <- prune_samples(sample_sums(ps_stool_16S_raw) >= 10^3, ps_stool_16S_raw)
#  
#  # Remove features with zero counts
#  ps_stool_16S_filtered <- filter_taxa(ps_stool_16S_pruned,function(x) sum(x>0)>0,1)
#  
#  # Collapse counts to the genus level
#  ps_stool_16S <- tax_glom(ps_stool_16S_filtered, taxrank = "GENUS")

## -----------------------------------------------------------------------------
example_HURDLE <- fitHURDLE(counts = as(otu_table(ps_stool_16S),
                                        "matrix"),
                            scale = "median")

head(example_HURDLE)

## -----------------------------------------------------------------------------
observed_hurdle <- prepareObserved(counts = as(otu_table(ps_stool_16S),
                                               "matrix"),
                                   scale = "median")
head(observed_hurdle)

## -----------------------------------------------------------------------------
head(prepareObserved(counts = as(otu_table(ps_stool_16S),"matrix")))

## -----------------------------------------------------------------------------
head(meanDifferences(estimated = example_HURDLE,
                     observed = observed_hurdle))

## ----fitting, warning=FALSE---------------------------------------------------
GOF_stool_16S <- fitModels(object = ps_stool_16S,
                           models = c("NB","ZINB","DM","ZIG","HURDLE"),
                           scale_ZIG = c("median","default"), 
                           scale_HURDLE = c("median","default"))

## ----RMSE_MD------------------------------------------------------------------
RMSE_MD_16S <- plyr::ldply(GOF_stool_16S,
                           .fun = function(df) cbind("RMSE" = RMSE(df[,"MD"])),
                           .id =  "Model")
RMSE_MD_16S

## ----RMSE_ZPD-----------------------------------------------------------------
RMSE_ZPD_16S <- plyr::ldply(GOF_stool_16S, 
                            .fun = function(df) cbind("RMSE" = RMSE(df[,"ZPD"])),
                            .id =  "Model")
RMSE_ZPD_16S

## ----plotGOF_MD, fig.width=15, fig.height=4-----------------------------------
plotMD(data = GOF_stool_16S,
       difference = "MD",
       split = TRUE)

## ----plotGOF_MD_noHurdleDefault, fig.width=15, fig.height=4-------------------
plotMD(data = GOF_stool_16S[1:6],
       difference = "MD",
       split = TRUE)

## ----plotGOF_ZPD, fig.width=15, fig.height=4----------------------------------
plotMD(data = GOF_stool_16S[1:6],
       difference = "ZPD",
       split = TRUE)

## ----plotGOF_MD_collapsed, fig.width=5, fig.height=10-------------------------
plot_grid(plotMD(data = GOF_stool_16S[1:6], difference = "MD", split = FALSE),
          plotMD(data = GOF_stool_16S[1:6], difference = "ZPD", split = FALSE),
          ncol = 1)

## ----plotGOF_RMSE, fig.width=5,fig.height=10----------------------------------
plot_grid(plotRMSE(GOF_stool_16S,difference = "MD"), 
          plotRMSE(GOF_stool_16S,difference = "ZPD"), ncol = 1)

## ----createMocks--------------------------------------------------------------
set.seed(123)
mock_df <- createMocks(nsamples = nsamples(ps_stool_16S),
                       N = 10) # At least 1000 is suggested

## ----normalization------------------------------------------------------------
ps_stool_16S <- norm_edgeR(object = ps_stool_16S,
                           method = "TMM")
ps_stool_16S <- norm_DESeq2(object = ps_stool_16S,
                            method = "poscounts")
ps_stool_16S <- norm_CSS(object = ps_stool_16S,
                         "median")

zinbweights <- weights_ZINB(object = ps_stool_16S, 
                            K = 0,
                            design = "~ 1")

## ----DA_TIEC------------------------------------------------------------------
# Random grouping each time
Stool_16S_mockDA <- apply(X = mock_df, MARGIN = 1, FUN = function(x){
    # Group assignment
    sample_data(ps_stool_16S)$group <- x
    ### DA analysis ###
    returnList = list()
    returnList = within(returnList, {
    da.edger <- DA_edgeR(object = ps_stool_16S,
                         group = unlist(sample_data(ps_stool_16S)[,"group"]),
                         design = as.formula("~ group"),
                         coef = 2,
                         norm = "TMM")
    da.edger.zinb <- DA_edgeR(object = ps_stool_16S,
                         group = unlist(sample_data(ps_stool_16S)[,"group"]),
                         design = as.formula("~ group"),
                         coef = 2,
                         norm = "TMM",
                         weights = zinbweights)
    da.deseq <- DA_DESeq2(object = ps_stool_16S,
                         design = as.formula("~ group"),
                         norm = "poscounts",
                         contrast = c("group","grp2","grp1"))
    da.deseq.zinb <- DA_DESeq2(object = ps_stool_16S,
                         design = as.formula("~ group"),
                         norm = "poscounts",
                         contrast = c("group","grp2","grp1"),
                         weights = zinbweights)
    da.limma <- DA_limma(object = ps_stool_16S,
                         design = ~ group,
                         coef = 2,
                         norm = "TMM")
    da.limma.zinb <- DA_limma(object = ps_stool_16S,
                         design = ~ group,
                         coef = 2,
                         norm = "TMM",
                         weights = zinbweights)
    da.limma.css <- DA_limma(object = ps_stool_16S,
                         design = ~ group,
                         coef = 2,
                         norm = "CSSmedian")
    })
    return(returnList)
})

## ----customExample, eval=FALSE------------------------------------------------
#  DA_yourMethod <- function(parameters) # others
#  {
#      ### your method code ###
#  
#      ### extract important variables ###
#      vector_of_pval = NA # contains the p-values
#      vector_of_adjusted_pval = NA # contains the adjusted p-values
#      name_of_your_features = NA # contains the OTU, or ASV, or other feature names. Usually extracted from the rownames of the count data
#  
#      ### prepare the output ###
#      name = "write.here.the.name"
#      pValMat <- data.frame("pval" = vector_of_pval,
#                            "adjp" = vector_of_adjusted_pval)
#      rownames(pValMat) <- name_of_your_features # Be sure that your method hasn't changed the order of the features. If it happens, you'll need to re-establish the original order.
#  
#      return(list("pValMat" = pValMat,
#                  "name" = name))
#  
#  }# END - function: DA_yourMethod

## ----createTIEC---------------------------------------------------------------
TIEC_summary <- createTIEC(Stool_16S_mockDA)

## ----FDRplot------------------------------------------------------------------
plotFPR(df_FPR = TIEC_summary$df_FPR)

## ----QQplot-------------------------------------------------------------------
plotQQ(df_QQ = TIEC_summary$df_QQ, zoom = c(0,0.1))

## ----KSplot-------------------------------------------------------------------
plotKS(df_KS = TIEC_summary$df_KS)

## ----dataloading_concordance--------------------------------------------------
data("ps_plaque_16S")

## ----16S_plaque_download, eval=FALSE------------------------------------------
#  ps_plaque_16S_raw <- V35() %>% # Extracting V3-V5 16S sequenced regions' count data
#    subset(select = HMP_BODY_SUBSITE %in% c("Supragingival Plaque", "Subgingival Plaque") & # Only gingival plaque samples
#             RUN_CENTER == "WUGC" & # Only sequenced at the WUCG RUN CENTER
#             SEX == "Male" & # Only male subject
#             VISITNO == 1) %>% # Only the first visit
#    as_phyloseq()
#  
#  # Only paired samples
#  paired <- names(which(table(sample_data(ps_plaque_16S_raw)[,"RSID"]) == 2))
#  ps_plaque_16S_paired <- subset_samples(ps_plaque_16S_raw, RSID %in% paired)
#  
#  # Remove low depth samples
#  ps_plaque_16S_pruned <- prune_samples(sample_sums(ps_plaque_16S_paired) >= 10^3, ps_plaque_16S_paired)
#  
#  # Remove features with zero counts
#  ps_plaque_16S_filtered <- filter_taxa(ps_plaque_16S_pruned,function(x) sum(x>0)>0,1)
#  
#  # Collapse counts to the genus level
#  ps_plaque_16S <- tax_glom(ps_plaque_16S_filtered, taxrank = "GENUS")

## -----------------------------------------------------------------------------
set.seed(123)
splits_df <- createSplits(object = ps_plaque_16S, 
                          varName = "HMP_BODY_SUBSITE", 
                          paired = "RSID", 
                          balanced = TRUE, 
                          N = 100)

## ----normalization_plaque-----------------------------------------------------
ps_plaque_16S <- norm_edgeR(object = ps_plaque_16S,
                            method = "TMM")
ps_plaque_16S <- norm_DESeq2(object = ps_plaque_16S,
                             method = "poscounts")
ps_plaque_16S <- norm_CSS(object = ps_plaque_16S,
                          "median")

zinbweights <- weights_ZINB(object = ps_plaque_16S, 
                            K = 0,
                            design = "~ 1 + HMP_BODY_SUBSITE")

## ---- eval=FALSE--------------------------------------------------------------
#  # Random grouping each time
#  Plaque_16S_splitsDA <- apply(X = mock_df, MARGIN = 1, FUN = function(x){
#      # Group assignment
#      sample_data(ps_stool_16S)$group <- x
#      ### DA analysis ###
#      returnList = list()
#      returnList = within(returnList, {
#      da.edger <- DA_edgeR(object = ps_stool_16S,
#                           group = ps_stool_16S@sam_data$group,
#                           design = as.formula("~ group"),
#                           coef = 2,
#                           norm = "TMM")
#  
#      da.deseq <- DA_DESeq2(object = ps_stool_16S,
#                           design = as.formula("~ group"),
#                           norm = "poscounts",
#                           contrast = c("group","grp2","grp1"))
#  
#      da.deseq.zinb <- DA_DESeq2(object = ps_stool_16S,
#                           design = as.formula("~ group"),
#                           norm = "poscounts",
#                           contrast = c("group","grp2","grp1"),
#                           weights = zinbweights)
#  
#      da.limma <- DA_limma(object = ps_stool_16S,
#                           design = ~ group,
#                           coef = 2,
#                           norm = "TMM")
#  
#      da.limma.css <- DA_limma(object = ps_stool_16S,
#                           design = ~ group,
#                           coef = 2,
#                           norm = "CSSmedian")
#      })
#      return(returnList)
#  })

