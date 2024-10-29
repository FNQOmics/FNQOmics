#' Batch - build_OTU_counts analysis of many comparisons
#
#' @description Given a summarized experiment generated using buildSummarized()
#' this function will automatically perform differential abundance (DA)
#' analysis for all possible groups using 3 different methods 1) EdgeR, 2) Voom
#' and 3) ALDEx2 It will also output 3 plots automatically, if the
#' plotting options are selected (see ?diag_plots for more details).
#
#' @param biom A "phyloseq" object with included groups
#'  to be analysed. For format specifications see ?phyloseq E.g.
#'  accessible as "build_OTU_counts$otu". Groups are used to automate colouring of
#'  samples. Default = NULL
#'
#' @param sample_table A sample table to map groups to corresponding samples.  Default = NULL
#'
#' @param taxa taxanomic rank to match ranks. Default = NULL
#'
#' @param verbose Verbosity ON/OFF. Default=FALSE
#'
#' @return A list of all OTU comparisons conducted.
#'
#'
#' @export build_OTU_counts
#'
#' @importFrom phyloseq import_biom import_qiime_sample_data merge_phyloseq tax_table tax_table<-
# @import airway



build_OTU_counts <- function(biom = NULL,
                             sample_table = NULL,
                             taxa = NULL,
                             verbose = FALSE){
  ####///---- check inputs ----\\\###
  if(is.null(biom) | is.null(sample_table)) {
    stop("EITHER a biom or sample_table is not provided. Please provide filenames with full path and rerun.")
  }

  # Import biom and sample table data
  biom <- import_biom(biom)
  samples <- import_qiime_sample_data(sample_table)

  # Merge phyloseq object
  phylo <- merge_phyloseq(biom, samples)

  # Define taxonomy columns

  tax_col <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

  colnames(tax_table(phylo)) <- tax_col

  taxa <-  tax_col

  return (phylo)
}



