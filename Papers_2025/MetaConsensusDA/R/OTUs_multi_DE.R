#' Batch - OTUs_multi_DE analysis of many comparisons
#
#' @description Given a summarized experiment generated using buildSummarized()
#' this function will automatically perform differential expression (DE)
#' analysis for all possible groups using 3 different methods 1) EdgeR, 2) Voom
#' and 3) DEseq2. It will also output 10x diagnostic plots automatically, if the
#' plotting options are selected (see ?diag_plots for more details).
#
#' @param build_OTU_counts_output A "phyloseq" object with included groups
#'  to be analysed. For format specifications see ?phyloseq E.g.
#'  accessible as "build_OTU_counts$otu". Groups are used to automate colouring of
#'  samples. Default = NULL
#' after DE analysis from merged results. Options: TRUE, FALSE. Default = TRUE

#' @param verbose Verbosity ON/OFF. Default=FALSE
#'
#'
#' @return A list of all OTU comparisons conducted.
#'
#' @export OTUs_multi_DE
#'
#' @importFrom ggplot2 ggplot aes geom_bar
#' @importFrom phyloseq phyloseq_to_deseq2 get_variable nsamples otu_table tax_table psmelt
#' @importFrom DESeq2 DESeqDataSet DESeq results
#' @importFrom edgeR DGEList exactTest estimateCommonDisp estimateTagwiseDisp topTags estimateDisp estimateGLMCommonDisp estimateGLMTagwiseDisp calcNormFactors glmFit
#' @importFrom ALDEx2 aldex
# @import airway
#'
#'


OTUs_multi_DE <- function(build_OTU_counts_output = NULL,
                             verbose = FALSE){



      phyloseq_to_edgeR = function(physeq, group, method="RLE", ...){
        # Enforce orientation.
        if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
        x = as(otu_table(physeq), "matrix")
        # Add one to protect against overflow, log(0) issues.
        x = x + 1
        # Check `group` argument
        if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
          # Assume that group was a sample variable name (must be categorical)
          group = get_variable(physeq, group)
        }
        # Define gene annotations (`genes`) as tax_table
        taxonomy = tax_table(physeq, errorIfNULL=FALSE)
        if( !is.null(taxonomy) ){
          taxonomy = data.frame(as(taxonomy, "matrix"))
        }
        # Now turn into a DGEList
        y = DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
        # Calculate the normalization factors
        z = calcNormFactors(y, method=method)
        # Check for division by zero inside `calcNormFactors`
        if( !all(is.finite(z$samples$norm.factors)) ){
          stop("Something wrong with edgeR::calcNormFactors on this data,
               non-finite $norm.factors, consider changing `method` argument")
        }
        # Estimate dispersions
        return(estimateTagwiseDisp(estimateCommonDisp(z)))
      }



    ####///---- check inputs ----\\\###
    if(is.null(build_OTU_counts_output)) {
      stop("A build_OTU_counts_output object is not provided. Please provide filenames with full path and rerun.")
    }


   treat_list <- unique(sample_data(build_OTU_counts_output)$Age_Group)

 # treat_list <- unique(sample_data(build_OTU_counts_output)$Age_Group)
    # Initialize a list to store the results of all comparisons
    all_comparisons_results <- list()

    # Prepare contrast list
    contrast_list <- list()
    DESeq2_OTU_DE_results <- list()
    for (treatment in treat_list) {
      for (other_treatment in treat_list) {
        if (treatment != other_treatment) {

          contrast_list <- c(contrast_list, list(c("Age_Group", treatment, other_treatment)))

          contrast_pair <- paste0(treatment, "_vs_", other_treatment)
          cat("Running comparison:", contrast_pair, "\n")

          # Initialize a list to store the results of each method for this comparison
          comparison_results <- list()

          ######## EdgeR ###########

         print(treatment)
         print(other_treatment)


         # Subset the samples based on the pair of treatments
         filtered_sample_data <- filtered_sample_data <- sample_data(build_OTU_counts_output)[
           sample_data(build_OTU_counts_output)$Age_Group %in% c(treatment, other_treatment), ]

         # Filter OTU table and taxonomy based on filtered sample data
         filtered_OTU_table <- otu_table(build_OTU_counts_output)[, rownames(filtered_sample_data)]
         filtered_tax_table <- tax_table(build_OTU_counts_output)

         # Create a new phyloseq object with the filtered data
         build_OTU_subset <- phyloseq(filtered_OTU_table, filtered_tax_table, filtered_sample_data)

         # Convert the subsetted phyloseq object to an edgeR-compatible object
         test_phylo_reads_edgeR <- phyloseq_to_edgeR(build_OTU_subset, group = "Age_Group")

         # Perform differential abundance analysis
         et <- exactTest(test_phylo_reads_edgeR)

         # Extract top differentially abundant features
         tt <- topTags(et, n = nrow(test_phylo_reads_edgeR$table), adjust.method = "BH", sort.by = "PValue")

         #comparison
         print(tt[1,1])
         print("EdgeR: ")
         print( tt$comparison)

         # Extract results
         EdgeR_OTU_DE_results <- tt$table

         comparison_results$EdgeR <- EdgeR_OTU_DE_results

         write.table(EdgeR_OTU_DE_results, file=paste0(treatment,"_vs_",other_treatment,"_edgeR_DE_results.txt"), quote=F, sep="\t", col.names = TRUE)


          ######## DESeq2 ########
         print(treatment)
         print(other_treatment)


            phylo_reads_collapsed_deseq <- phyloseq_to_deseq2(build_OTU_counts_output, ~Age_Group)

   #         phylo_reads_collapsed_deseq <- DESeq(phylo_reads_collapsed_deseq, test="Wald",fitType="parametric")
            phylo_reads_collapsed_deseq <- DESeq(phylo_reads_collapsed_deseq, sfType = "poscounts")

            DESeq2_OTU_DE_results = results(phylo_reads_collapsed_deseq , contrast=c("Age_Group", treatment, other_treatment), tidy=T, format="DataFrame")

            comparison_results$DESeq2 <- DESeq2_OTU_DE_results

            write.table(DESeq2_OTU_DE_results, file=paste0(treatment,"_vs_",other_treatment,"_DESeq_DE_results.txt"), quote=F, sep="\t", col.names = TRUE)

          ############## ALDEx2 #####

          print(treatment)
          print(other_treatment)

          otu_table <- as.matrix(otu_table(build_OTU_counts_output))

          group_info <- sample_data(build_OTU_counts_output)$Age_Group

          # Filter OTU table and group info for the specified conditions
          condition_indices <- which(group_info %in% c(treatment, other_treatment))

          count_treatment <- sum(sample_data(build_OTU_counts_output)$Age_Group == treatment)

          count_other_treatment <- sum(sample_data(build_OTU_counts_output)$Age_Group == other_treatment)

          filtered_otu <- otu_table[, condition_indices]

          #print(filtered_otu)

          filtered_groups <- group_info[condition_indices]

     #     conds <-  c(rep(treatment, count_treatment), rep(other_treatment, count_other_treatment))

          print(group_info)

          print(filtered_groups)

          ALDEx2_OTU_DE_results <- aldex(reads = filtered_otu,  conditions = filtered_groups, mc.samples = 128, denom = "all", verbose = TRUE, useMC = TRUE, cores = 2, test="t", effect=TRUE, paired.test=FALSE)

          comparison_results$ALDEx2 <- ALDEx2_OTU_DE_results

          #print(ALDEx2_OTU_DE_results$wi.eBH < 0.05)

          write.table(ALDEx2_OTU_DE_results, file=paste0(treatment,"_vs_",other_treatment,"_aldex_DE_results.txt"), quote=FALSE, sep='\t', col.names = NA)

          all_comparisons_results[[contrast_pair]] <- comparison_results

        }
      }
    }

    return(all_comparisons_results)
  }

