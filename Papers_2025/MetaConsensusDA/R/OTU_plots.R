#' Batch - OTU_plots
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
#'
#' @export OTU_plots
#' @importFrom ggplot2 geom_boxplot element_line aes_string labs guide_legend theme_bw guides element_blank theme_minimal scale_color_manual ggplot ggsave aes geom_bar facet_grid scale_fill_manual theme element_text
#' @importFrom utils write.table globalVariables
#' @importFrom phyloseq phyloseq sample_data taxa_are_rows transform_sample_counts tax_glom psmelt estimate_richness ordinate plot_ordination
# @import airway
#'



OTU_plots <- function( build_OTU_counts_output = NULL,
                      verbose = FALSE){

  ####///---- check inputs ----\\\###
  if(is.null(build_OTU_counts_output)) {
    stop("A build_OTU_counts_output object is not provided. Please provide filenames with full path and rerun.")
  }


  ####### scale plots ########
  tax_table(build_OTU_counts_output)

#  colnames(tax_table(build_OTU_counts_output)) <- tax_col

  phylo2 <- transform_sample_counts(build_OTU_counts_output,function(x) x / sum(x))

  phylo2_class <- tax_glom(phylo2,taxrank='Class')

  data_phylo2_class <- psmelt(phylo2_class)

  data_phylo2_class$Class <- as.character(data_phylo2_class$Class)

  data_phylo2_class$Class[data_phylo2_class$Abundance < 0.01] <- "<1% Abundance"

  data_phylo2_class[data_phylo2_class == 0] <- NA

  Count = length(unique(data_phylo2_class$Class))

  unique(data_phylo2_class$Class)

  data_phylo2_class$Class <- factor(data_phylo2_class$Class, levels = c("c__Alphaproteobacteria","c__Bdellovibrionia_A","c__Gammaproteobacteria","c__Clostridia","c__Actinomycetia","c__Bacteroidia","c__Chitinophagia","c__Actinobacteria","c__Fimbriimonadia","c__Verrucomicrobiae","c__Saccharimonadia","c__Bacilli","c__Chthonomonadetes","<1% Abundance"))

  scale_plot <- ggplot(data=data_phylo2_class,aes_string(x="Sample",y="Abundance",fill="Class")) + facet_grid(~Age_Group,scales="free") + geom_bar(aes(), stat="identity", position="stack") +
    scale_fill_manual(values = c("royalblue4", "deepskyblue", "blue","cyan2", "darkorchid","gold1", "forestgreen", "firebrick", "mediumspringgreen", "darkorange1", "deeppink","grey","slategray2", "black"))  +
    theme(  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  write.table(data_phylo2_class, file = "tt_class_table.txt", quote=F, sep="\t", col.names = TRUE)

  ggsave("Adults_scale_plot.pdf", plot = scale_plot, width = 8, height = 6)



  ######## Alpha Diversity #####

  species_physeq <- tax_glom(build_OTU_counts_output, taxrank = "Species")

  species_shannon_diversity <- estimate_richness(species_physeq, measures = "Shannon")

  alpha_diversity <- estimate_richness(build_OTU_counts_output, measures = "Shannon")

  sample_data_df <- data.frame(sample_data(build_OTU_counts_output))

  sample_data_df$SampleID <- rownames(sample_data_df)

  alpha_diversity$SampleID <- rownames(alpha_diversity)

  merged_data <- merge(alpha_diversity, sample_data_df, by = "SampleID")

  Alpha_plot <- ggplot(merged_data, aes_string(x = "Age_Group", y = "Shannon", color = "Age_Group"))

  Alpha_plot2 <-  Alpha_plot  + geom_boxplot() + scale_color_manual(values = c("royalblue4", "deepskyblue", "blue","cyan2", "darkorchid","gold1", "forestgreen", "firebrick", "mediumspringgreen", "darkorange1", "deeppink","grey","slategray2" )) +  theme_minimal() +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          axis.line = element_line(color = "black")) +
    labs(title = "Alpha Diversity: Shannon Index",
         x = "Sample Type",
         y = "Shannon Diversity Index")


  ggsave("Adults_groups_alpha_div_Box_plot.pdf", plot = Alpha_plot2, width = 8, height = 6)


  ####### PCoA / MDS #########


  pcoa <- ordinate(build_OTU_counts_output,"PCoA")

  pcoa_plot <- plot_ordination(build_OTU_counts_output,ordinate(build_OTU_counts_output, "PCoA"),color="Age_Group", label = "SampleID") + theme_bw()

  pcoa_plot2 <- pcoa_plot + theme(panel.grid = element_blank(), panel.border = element_blank(),
                                  axis.line = element_line(color = "black")) +
    guides(colour=guide_legend(override.aes=list(size=5)))  +
    scale_color_manual(values = c("royalblue4", "deepskyblue", "blue","cyan2", "darkorchid","gold1", "forestgreen", "firebrick", "mediumspringgreen", "darkorange1", "deeppink","grey","slategray2", "dodgerblue", "orange2", "maroon","navy"))

  ggsave("Both_groups_mds_plot.pdf", plot = pcoa_plot2, width = 8, height = 6)


  ####### PCoA / MDS #########


  OTU_plot_tables <- list(

    scale_plot_table <- data_phylo2_class,

    AlphaDiv_plot_table <- species_shannon_diversity,

    PCoA_table <- pcoa$vectors
  )

  return(OTU_plot_tables)
}

