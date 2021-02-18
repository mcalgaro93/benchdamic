#' @title (Data) 33 Stool samples of 16S rRNA (HMP 2012)
#'
#' @description A demonstrative purpose dataset containing microbial abundances
#' for a total of 71 OTUs. The 32 Stool samples belong to the Human Microbiome
#' Project. This particular subset contains the SEX = "Male", RUN_CENTER = "BI",
#' and VISITNO = "1" samples. It is possible to obtain the same dataset after
#' basic filters (remove taxa with zero counts) and collapsing the counts to the
#' genus level; [HMP16Data] Bioconductor package was used to download the data.
#'
#' @docType data
#' @aliases ps_stool_16S
#'
#' @usage data(ps_stool_16S)
#'
#' @format An object of class \code{"phyloseq"}
"ps_stool_16S"
