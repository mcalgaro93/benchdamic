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

## ----16S_data_download, eval=FALSE--------------------------------------------
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

## ----fitting, warning=FALSE---------------------------------------------------
list_16S <- fitModels(counts = ps_stool_16S@otu_table@.Data,
                      models = c("NB","ZINB","DM","ZIG","HURDLE"))

df_16S <- plyr::ldply(list_16S,.id = "Model")

## ----RMSE_MD------------------------------------------------------------------
RMSE_MD_16S <- plyr::ldply(list_16S,.fun = function(df) cbind("RMSE" = RMSE(df[,"MD"])),.id =  "Model")
RMSE_MD_16S

## ----RMSE_ZPD-----------------------------------------------------------------
RMSE_ZPD_16S <- plyr::ldply(list_16S, 
                            .fun = function(df) cbind("RMSE" = RMSE(df[,"ZPD"])),
                            .id =  "Model")
RMSE_ZPD_16S

## ----plotGOF, fig.width=15, fig.height=8--------------------------------------
cowplot::plot_grid(plotlist = list(plotMD(data = df_16S,
                                          difference = "MD",
                                          split = TRUE),
                                   plotMD(data = df_16S,
                                          difference = "ZPD",
                                          split = TRUE)),
                   nrow = 2)

## ----plotGOF_collapsed, fig.width=12, fig.height=5----------------------------
cowplot::plot_grid(plotlist = list(plotMD(data = df_16S, 
                                          difference = "MD", 
                                          split = FALSE),
                                   plotMD(data = df_16S, 
                                          difference = "ZPD", 
                                          split = FALSE)),
                   nrow = 1)

## ----createMocks--------------------------------------------------------------
mock_df <- createMocks(nsamples = nsamples(ps_stool_16S),
                       N = 10, # At least 1000 is suggested
                       seed = 123)

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
                         group = ps_stool_16S@sam_data$group,
                         design = as.formula("~ group"),
                         coef = 2,
                         norm = "TMM")
    da.edger.zinb <- DA_edgeR(object = ps_stool_16S,
                         group = ps_stool_16S@sam_data$group,
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

## ----QQPlot-------------------------------------------------------------------
plotQQ(df_QQ = TIEC_summary$df_QQ)

## -----------------------------------------------------------------------------
plotKS(df_KS = TIEC_summary$df_KS)

