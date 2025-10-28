# ConsensusMetaDA

ConsensusMetaDA is an R package for microbiome analysis using multiple algorithms - reaching consensus. No one tool perform best DA anlayes on metaganome data. Hence it needs a consensus approach to obtain more robust differential abundant microbiome. ConsensusMetaDA helps achieve this using five popularly used tools. ConsensusMetaDA uses popular biom format and samples table as input to perform DA along with it generates standard microbiome visualization plots such as Alpha Diversity, Beta Diversity, Rarefaction curve, Bidirectional plot and Scale plot at various taxa levels.

## Installing ConsensusMetaDA

To obtain the original version from github, install devtools in R and use the following:

```R

required packages

# ADAPT R package requires R version 4.4 hence need all R updated to R v4.4.2

# One line per package
# if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# if (!require("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
# if (!require("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
# if (!require("ALDEx2", quietly = TRUE)) BiocManager::install("ALDEx2")
# if (!require("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
# if (!require("ADAPT", quietly = TRUE)) BiocManager::install("ADAPT")
# if (!require("metagenomeSeq", quietly = TRUE)) BiocManager::install("metagenomeSeq")
# 
# if (!require("ConsensusMetaDA", quietly = TRUE)) {
#   if (!require("remotes", quietly = TRUE)) install.packages("remotes")
#   remotes::install_github("kmanoharan01/ConsensusMetaDA")
# }


# Load the packages
library(ggplot2)
library(phyloseq)
library(DESeq2)
library(ALDEx2)
library(edgeR)
library(metagenomeSeq)
library(ADAPT)

```


# ConsensusMetaDA examples

To illustrate functionality of ```ConsensusMetaDA```, we will utilise microbiome biom format data from the ```ecc_saliva``` and annotation libraries as follows. Begin by installing and attaching data from these libraries as follows:

```R
# load ConsensusMetaDA

library(ConsensusMetaDA)

library(ADAPT)

data("ecc_saliva")

otu_mat <- as(otu_table(ecc_saliva), "matrix")
sample_df <- as(sample_data(ecc_saliva), "data.frame")
tax_mat <- as(tax_table(ecc_saliva), "matrix")

# ConsensusMetaDA examples

To illustrate functionality of ```ConsensusMetaDA```, we will utilise microbiome biom format data from the ```ecc_saliva``` and annotation libraries as follows. Begin by installing and attaching data from these libraries as follows:

ecc_saliva_biom <- "../inst//extdata/ecc_saliva_Otu.biom"

ecc_saliva_sample_table <- "../inst//extdata/ecc_saliva_sample_table.txt"

ecc_saliva_tax <- tax_mat

## build
ecc_saliva_build_OTU_counts <- build_OTU_counts(ecc_saliva_biom, 
                                                 ecc_saliva_sample_table, 
                                                 ecc_saliva_tax, 
                                                 abundance_threshold = 10, 
                                                 prevalence_threshold = 0.1, 
                                                 rarity_threshold = 0, 
                                                 variance_threshold = 0 
)

ecc_saliva_build_OTU_counts

## Perform DA

ecc_saliva_DA <- OTUs_multi_DA(ecc_saliva_build_OTU_counts)

## Perform Plots
Creates Five different plots and saves them as pdf files. (1) Alpha diversity - Shannon, (2) Beta diversity – Principal Coordinate Analyses (PcoA), (3) rare faction curve and (4) Scaled plot at taxa level 5) Bi-directional plots. These plots are generated using the ggplot2 package (Wickham, 2016).

ecc_saliva_plots <- OTUs_plots(ecc_saliva_build_OTU_counts)

```

## Contact

For more details, contact Manoharan Kumar:
manoharan.kumar@jcu.edu.au
