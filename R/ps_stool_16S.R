#' @title (Data) 33 Stool samples of 16S rRNA (HMP 2012)
#'
#' @description A demonstrative purpose dataset containing microbial abundances
#' for a total of 77 OTUs. The 32 Stool samples belong to the Human Microbiome
#' Project. This particular subset contains the SEX = "Male", RUN_CENTER = "BI",
#' and VISITNO = "1" samples. It is possible to obtain the same dataset after
#' several filters (only taxa present with more than 10 counts in more than 5
#' sample are kept) using the [HMP16Data] Bioconductor package.
#'
#' @docType data
#' @aliases ps_stool_16S
#'
#' @usage data(ps_stool_16S)
#'
#' @format An object of class \code{"phyloseq"}
"ps_stool_16S"
