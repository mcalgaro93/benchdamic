#' @title (Data) 60 Gingival Plaque samples of 16S rRNA (HMP 2012)
#'
#' @description A demonstrative purpose dataset containing microbial abundances
#' for a total of 88 OTUs. The 60 Gingival Plaque paired samples belong to the
#' Human Microbiome Project. This particular subset contains 30 Supragingival
#' and 30 Subgingival Plaque samples from the SEX = "Male", RUN_CENTER = "WUCG",
#' and VISITNO = "1" samples. It is possible to obtain the same dataset after
#' basic filters (remove taxa with zero counts) and collapsing the counts to the
#' genus level; [HMP16Data] Bioconductor package was used to download the data.
#'
#' @docType data
#' @aliases ps_plaque_16S
#'
#' @usage data(ps_plaque_16S)
#'
#' @format An object of class \code{"phyloseq"}
"ps_plaque_16S"
