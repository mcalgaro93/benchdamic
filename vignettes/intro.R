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
#  ps_stool_16S_pruned <- prune_samples(sample_sums(ps_stool_16S_raw) >= 10^3, ps_stool_16S_raw)
#  
#  ps_stool_16S <- filter_taxa(ps_stool_16S_pruned,function(x) sum(x>10)>5,1)

## ----fitting_16S, warning=FALSE-----------------------------------------------
list_16S <- fitModels(counts = ps_stool_16S@otu_table@.Data,
                      models = c("NB","ZINB","DM","ZIG","HURDLE"))

df_16S <- plyr::ldply(list_16S,.id = "Model")

## ----RMSE_MD_16S--------------------------------------------------------------
RMSE_MD_16S <- plyr::ldply(list_16S,.fun = function(df) cbind("RMSE" = RMSE(df[,"MD"])),.id =  "Model")
RMSE_MD_16S

## ----RMSE_ZPD_16S-------------------------------------------------------------
RMSE_ZPD_16S <- plyr::ldply(list_16S,.fun = function(df) cbind("RMSE" = RMSE(df[,"ZPD"])),.id =  "Model")
RMSE_ZPD_16S

## ----plot16S, fig.width=15, fig.height=8--------------------------------------
cowplot::plot_grid(plotlist = list(MDPlot(data = df_16S,difference = "MD",split = TRUE),
MDPlot(data = df_16S,difference = "ZPD",split = TRUE)),
nrow = 2)

## ----plot16S_collapsed, fig.width=12, fig.height=5----------------------------
cowplot::plot_grid(plotlist = list(MDPlot(data = df_16S,difference = "MD",split = FALSE),
MDPlot(data = df_16S,difference = "ZPD",split = FALSE)),
nrow = 1)

## -----------------------------------------------------------------------------


