# Unload package 'here' if it is already loaded
loaded_packages <- search()[grepl("^package:", search())]
loaded_packages <- sub("^package:", "", loaded_packages)

if ("here" %in% loaded_packages) {
  detach("package:here", unload = TRUE)
}

setwd("/folder/containing/cellranger/output")

# This script requires over 63 Gb of ram for the intDat step

# 1 Load libraries
library(here)
library(scater)
library(Seurat)
library(clustree)
library(dplyr)
library(patchwork)
library(tidyverse)
library(cowplot)
library(sctransform)

# 1. Load Data
# load preprocessed data
# Initialize an empty list to store the 'sce' objects
sce_list <- list()

# This script is batched into two lots of 8 captures because when you o do all
# 16 in one go, the intDat step produced this error:
# 
# Error in .M2C(newTMat(i = c(ij1[, 1], ij2[, 1]), j = c(ij1[, 2], ij2[,  : 
# unable to aggregate TsparseMatrix with 'i' and 'j' slots of length exceeding 2^31-1
#
# This is apparently because the length of slots 'i' and 'j' exceeds a limit of 2^31-1
# which is a maximum limit imposed by the R environment that corresponds to the maximum value of a signed 32-bit integer

# Define a vector of numbers from 1 to 16
numbers <- 1:8

# Define a vector of letters from A to P
letters <- LETTERS[1:8]

# Loop through each number and append letter to cell IDs to denote experiment & avoid duplicate IDs
for (i in numbers) {
  file_path <- here("data", "SCEs", paste("dataset.preprocessed.", i, ".SCE.rds", sep=""))
  sce <- readRDS(file_path)
  
  # Append letter to cell IDs
  colnames(sce) <- paste0(letters[i], "-", colnames(sce))
  
  # Assign the modified sce to the list
  assign(paste0("sce", i), sce)
}

# Create a list to store the sce objects
sce_list <- list(sce1, sce2, sce3, sce4, sce5, sce6, sce7, sce8)

# Identify shared genes
# Extract the ID vectors from each sce object
id_list <- list()

for (i in 1:8) {
  id_list[[i]] <- rowData(get(paste0("sce", i)))$ID
  assign(paste0("id_sce", i), id_list[[i]])
}

# Find the common elements among the ID vectors
shared_genes <- Reduce(intersect, list(id_sce1, id_sce2, id_sce3, id_sce4, 
                                       id_sce5, id_sce6, id_sce7, id_sce8))

# sort & subset each SCE relative to shared genes
# Initialize an empty list to store the 'm' vectors
m_list <- list()

for (i in 1:8) {
  m <- match(shared_genes, rowData(sce_list[[i]])$ID)
  matches <- all(shared_genes == rowData(sce_list[[i]])$ID[m])
  m_list <- c(m_list, list(m))
  assign(paste0("m", i), m)
}

# Initialize an empty list to store the modified 'sce' objects
modified_sce_list <- list()

for (i in 1:8) {
  modified_sce <- sce_list[[i]][get(paste0("m", i)), ]
  modified_sce_list <- c(modified_sce_list, list(modified_sce))
  assign(paste0("sce", i), modified_sce)
}

# Create combined matrix of counts
combo_counts <- cbind(counts(sce1), counts(sce2), counts(sce3), counts(sce4), counts(sce5), counts(sce6), counts(sce7), counts(sce8))  # Combine counts from all four sce objects

# Combine cell metadata
combo_data <- DataFrame(
  bind_rows(
    colData(sce1) %>% 
      data.frame %>% 
      dplyr::select(Barcode, capture, nuclear_fraction, cell_status), 
    colData(sce2) %>% 
      data.frame %>% 
      dplyr::select(Barcode, capture, nuclear_fraction, cell_status),
    colData(sce3) %>%  
      data.frame %>% 
      dplyr::select(Barcode, capture, nuclear_fraction, cell_status),
    colData(sce4) %>%  # Add metadata from sce4
      data.frame %>% 
      dplyr::select(Barcode, capture, nuclear_fraction, cell_status),
    colData(sce5) %>%  
      data.frame %>% 
      dplyr::select(Barcode, capture, nuclear_fraction, cell_status),
    colData(sce6) %>%  
      data.frame %>% 
      dplyr::select(Barcode, capture, nuclear_fraction, cell_status),
    colData(sce7) %>%  
      data.frame %>% 
      dplyr::select(Barcode, capture, nuclear_fraction, cell_status),
    colData(sce8) %>%  
      data.frame %>% 
      dplyr::select(Barcode, capture, nuclear_fraction, cell_status)
  )
)

combo_data$HTO <- as.factor(c(
  as.character(sce1$HTO), 
  as.character(sce2$HTO), 
  as.character(sce3$HTO), 
  as.character(sce4$HTO),
  as.character(sce5$HTO),
  as.character(sce6$HTO),
  as.character(sce7$HTO),
  as.character(sce8$HTO)
))

combo_data$donor <- as.factor(c(
  as.character(sce1$genetic_donor), 
  as.character(sce2$genetic_donor), 
  as.character(sce3$genetic_donor), 
  as.character(sce4$genetic_donor),
  as.character(sce5$genetic_donor),
  as.character(sce6$genetic_donor),
  as.character(sce7$genetic_donor),
  as.character(sce8$genetic_donor)
))

combo_data$experiment <- as.factor(c(
  rep(1, ncol(sce1)), 
  rep(2, ncol(sce2)),
  rep(3, ncol(sce3)),
  rep(4, ncol(sce4)),
  rep(5, ncol(sce5)),
  rep(6, ncol(sce6)),
  rep(7, ncol(sce7)),
  rep(8, ncol(sce8))
))

# Create the combined sce object
sce <- SingleCellExperiment(
  list(counts = combo_counts),
  colData = combo_data,
  rowData = rowData(sce1)  # You can choose to use rowData from any of the original sce objects
)

# 3. QC 
# Identify uninformative genes
# Some useful gene sets
mito_set <- rownames(sce)[which(rowData(sce)$CHR == "MT")]
ribo_set <- grep("^RP(S|L)", rownames(sce), value = TRUE)
# NOTE: A more curated approach for identifying ribosomal protein genes
#       (https://github.com/Bioconductor/OrchestratingSingleCellAnalysis-base/blob/ae201bf26e3e4fa82d9165d8abf4f4dc4b8e5a68/feature-selection.Rmd#L376-L380)
library(msigdbr)
c2_sets <- msigdbr(species = "Homo sapiens", category = "C2")
ribo_set <- union(
  ribo_set,
  c2_sets[c2_sets$gs_name == "KEGG_RIBOSOME", ]$human_gene_symbol)
is_ribo <- rownames(sce) %in% ribo_set
sex_set <- rownames(sce)[any(rowData(sce)$ENSEMBL.SEQNAME %in% c("X", "Y"))]
pseudogene_set <- rownames(sce)[
  any(grepl("pseudogene", rowData(sce)$ENSEMBL.GENEBIOTYPE))]

is_mito <- rownames(sce) %in% mito_set
is_ribo <- rownames(sce) %in% ribo_set
sce <- addPerCellQC(
  sce, 
  subsets = list(Mito = which(is_mito), Ribo = which(is_ribo)))

head(colData(sce)) %>% knitr::kable()

sce$zero_percent <- colSums(counts(sce) == 0)/nrow(sce)
summary(sce$zero_percent)

# Visualise QC metrics
p1 <- plotColData(
  sce,
  "sum",
  x = "capture",
  colour_by = "experiment",
  point_size = 0.5) +
  scale_y_log10() +
  theme(axis.text.x = element_text(size = 6)) +
  annotation_logticks(
    sides = "l",
    short = unit(0.03, "cm"),
    mid = unit(0.06, "cm"),
    long = unit(0.09, "cm"))
p2 <- plotColData(
  sce,
  "detected",
  x = "capture",
  colour_by = "experiment",
  point_size = 0.5) +
  theme(axis.text.x = element_text(size = 6))
p3 <- plotColData(
  sce,
  "subsets_Mito_percent",
  x = "capture",
  colour_by = "experiment",
  point_size = 0.5) +
  theme(axis.text.x = element_text(size = 6))
p4 <- plotColData(
  sce,
  "subsets_Ribo_percent",
  x = "capture",
  colour_by = "experiment",
  point_size = 0.5) +
  theme(axis.text.x = element_text(size = 6))

p5 <- plotColData(
  sce, 
  x = "capture", 
  y = "zero_percent", 
  colour_by = "experiment",
  point_size = 0.5) +
  theme(axis.text.x = element_text(size = 6))

((p1 | p2) / (p3 | p4) / p5) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


p1 <- plotColData(
  sce,
  "sum",
  x = "donor",
  colour_by = "experiment",
  point_size = 0.5) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90)) +
  annotation_logticks(
    sides = "l",
    short = unit(0.03, "cm"),
    mid = unit(0.06, "cm"),
    long = unit(0.09, "cm"))

p2 <- plotColData(
  sce,
  "detected",
  x = "donor",
  colour_by = "experiment",
  point_size = 0.5) +
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 0.5, 
                                   vjust = 1))
p3 <- plotColData(
  sce,
  "subsets_Mito_percent",
  x = "donor",
  colour_by = "experiment",
  point_size = 0.5) +
  theme(axis.text.x = element_text(angle = 90))
p4 <- plotColData(
  sce,
  "subsets_Ribo_percent",
  x = "donor",
  colour_by = "experiment",
  point_size = 0.5) +
  theme(axis.text.x = element_text(angle = 90))

p5 <- plotColData(
  sce, 
  x = "donor", 
  y = "zero_percent", 
  colour_by = "experiment",
  point_size = 0.5) +
  theme(axis.text.x = element_text(angle = 90))

((p1 | p2) / (p3 | p4) / p5) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


plotColData(
  sce, 
  x = "sum", 
  y = "detected",
  other_fields = "capture",
  colour_by = "zero_percent",
  point_size = 0.25, 
  point_alpha = 0.25) +
  facet_wrap(vars(capture), ncol = 2)

plotColData(
  sce, 
  x = "detected", 
  y = "subsets_Mito_percent",
  other_fields = "capture",
  colour_by = "zero_percent",
  point_size = 0.25, 
  point_alpha = 0.25) +
  facet_wrap(vars(capture), ncol = 2)

plotColData(
  sce, 
  x = "sum", 
  y = "subsets_Ribo_percent",
  other_fields = "capture",
  colour_by = "zero_percent",
  point_size = 0.25, 
  point_alpha = 0.25) +
  scale_x_log10() +
  facet_wrap(vars(capture), ncol = 2)

plotColData(
  sce, 
  x = "sum", 
  y = "subsets_Mito_percent",
  other_fields = "capture",
  colour_by = "cell_status",
  point_size = 0.25, 
  point_alpha = 0.25) +
  scale_x_log10() +
  facet_wrap(vars(capture), ncol = 2)

colData(sce) %>%
  data.frame %>%
  ggplot(aes(x = cell_status, y = sum, fill = cell_status)) +
  geom_violin(scale = "count", size = 0.25) +
  scale_y_log10() +
  facet_wrap(~ capture, scales = "free_y")

# Discard uninformative genes & cells
uninformative <- is_mito | is_ribo | rownames(sce) %in% sex_set | rownames(sce) %in% pseudogene_set
sum(uninformative)

junk <- sce$donor %in% c("doublet", "unassigned") | sce$HTO %in% c("Doublet", "Unknown") | (sce$cell_status == "empty_droplet" & sce$experiment == 1)

sceFlt <- sce[!uninformative, !junk]
sceFlt

# Remove low-abundance genes
numCells <- nexprs(sceFlt, byrow = TRUE)
keep <- numCells > 20
sum(keep)

sceFlt <- sceFlt[keep,]
sceFlt

# Convert to Seurat object
counts <- counts(sceFlt)
rownames(counts) <- rowData(sceFlt)$Symbol

seu <- CreateSeuratObject(counts = counts, 
                          meta.data = data.frame(colData(sceFlt)))
seu

# Visualise combined, filtered data
DefaultAssay(seu) <- "RNA"
seu <- FindVariableFeatures(seu) %>% 
  ScaleData() %>% 
  RunPCA(verbose = FALSE, dims = 1:30) %>%
  RunUMAP(verbose = FALSE, dims = 1:30)

DimPlot(seu, group.by = "experiment", combine = FALSE)

DimPlot(seu, split.by = "experiment", combine = FALSE)

# 4. Integrate data
# Normalise the data using SCTransform and integrate across batches/individuals.
# First, import the utility.R script which contains the intDat function. utility.R can
# be found here: https://github.com/Oshlack/paed-cf-cite-seq/blob/submission2/code/utility.R
source(here("code", "utility.R"))
out <- here("data/SCEs/combined.integrated_1.SEU.rds")

if(!file.exists(out)) {
  seuInt_1 <- intDat(seu, split = "donor", type = "RNA")
  saveRDS(seuInt_1, file = out)
} else {
  seuInt_1 <- readRDS(out)
}

###############################################################################
# 1. Load Data
# load preprocessed data
# Initialize an empty list to store the 'sce' objects
sce_list <- list()

# Define a vector of numbers from 1 to 16
numbers <- 9:16

# Define a vector of letters from A to P
letters <- LETTERS[9:16]

# Loop through each number and append letter to cell IDs to denote experiment & avoid duplicate IDs
for (i in 9:16) {
  file_path <- here("data", "SCEs", paste("dataset.preprocessed.", i, ".SCE.rds", sep=""))
  sce <- readRDS(file_path)
  
  # Append letter to cell IDs
  colnames(sce) <- paste0(letters[i - 8], "-", colnames(sce))
  
  # Assign the modified sce to the list
  assign(paste0("sce", i), sce)
}

# Create a list to store the sce objects
sce_list <- list(sce9, sce10, sce11, sce12, sce13, sce14, sce15, sce16)

# Identify shared genes
# Extract the ID vectors from each sce object
id_list <- list()

for (i in 9:16) {
  id_list[[i]] <- rowData(get(paste0("sce", i)))$ID
  assign(paste0("id_sce", i), id_list[[i]])
}

# Find the common elements among the ID vectors
shared_genes <- Reduce(intersect, list(id_sce9, id_sce10, id_sce11, id_sce12, 
                                       id_sce13, id_sce14, id_sce15, id_sce16))

# sort & subset each SCE relative to shared genes
# Initialize an empty list to store the 'm' vectors
m_list <- list()

for (i in 1:8) {
  m <- match(shared_genes, rowData(sce_list[[i]])$ID)
  matches <- all(shared_genes == rowData(sce_list[[i]])$ID[m])
  m_list <- c(m_list, list(m))
  assign(paste0("m", i), m)
}

# Initialize an empty list to store the modified 'sce' objects
modified_sce_list <- list()

for (i in 1:8) {
  modified_sce <- sce_list[[i]][get(paste0("m", i)), ]
  modified_sce_list <- c(modified_sce_list, list(modified_sce))
  assign(paste0("sce", i), modified_sce)
}

# Create combined matrix of counts
combo_counts <- cbind(counts(sce9), counts(sce10), counts(sce11), counts(sce12), counts(sce13), counts(sce14), counts(sce15), counts(sce16))  # Combine counts from all four sce objects

# Combine cell metadata
combo_data <- DataFrame(
  bind_rows(
    colData(sce9) %>%  
      data.frame %>% 
      dplyr::select(Barcode, capture, nuclear_fraction, cell_status),
    colData(sce10) %>%  
      data.frame %>% 
      dplyr::select(Barcode, capture, nuclear_fraction, cell_status),
    colData(sce11) %>%  
      data.frame %>% 
      dplyr::select(Barcode, capture, nuclear_fraction, cell_status),
    colData(sce12) %>%  
      data.frame %>% 
      dplyr::select(Barcode, capture, nuclear_fraction, cell_status),
    colData(sce13) %>%  
      data.frame %>% 
      dplyr::select(Barcode, capture, nuclear_fraction, cell_status),
    colData(sce14) %>%  
      data.frame %>% 
      dplyr::select(Barcode, capture, nuclear_fraction, cell_status),
    colData(sce15) %>%  
      data.frame %>% 
      dplyr::select(Barcode, capture, nuclear_fraction, cell_status),
    colData(sce16) %>%  
      data.frame %>% 
      dplyr::select(Barcode, capture, nuclear_fraction, cell_status)
  )
)

combo_data$HTO <- as.factor(c(
  as.character(sce9$HTO),
  as.character(sce10$HTO),
  as.character(sce11$HTO),
  as.character(sce12$HTO),
  as.character(sce13$HTO),
  as.character(sce14$HTO),
  as.character(sce15$HTO),
  as.character(sce16$HTO)
))

combo_data$donor <- as.factor(c(
  as.character(sce9$genetic_donor),
  as.character(sce10$genetic_donor),
  as.character(sce11$genetic_donor),
  as.character(sce12$genetic_donor),
  as.character(sce13$genetic_donor),
  as.character(sce14$genetic_donor),
  as.character(sce15$genetic_donor),
  as.character(sce16$genetic_donor)
))

combo_data$experiment <- as.factor(c(
  rep(9, ncol(sce9)),
  rep(10, ncol(sce10)),
  rep(11, ncol(sce11)),
  rep(12, ncol(sce12)),
  rep(13, ncol(sce13)),
  rep(14, ncol(sce14)),
  rep(15, ncol(sce15)),
  rep(16, ncol(sce16))
))

# Create the combined sce object
sce <- SingleCellExperiment(
  list(counts = combo_counts),
  colData = combo_data,
  rowData = rowData(sce9)
)

# 3. QC 
# Identify uninformative genes
# Some useful gene sets
mito_set <- rownames(sce)[which(rowData(sce)$CHR == "MT")]
ribo_set <- grep("^RP(S|L)", rownames(sce), value = TRUE)
# NOTE: A more curated approach for identifying ribosomal protein genes
#       (https://github.com/Bioconductor/OrchestratingSingleCellAnalysis-base/blob/ae201bf26e3e4fa82d9165d8abf4f4dc4b8e5a68/feature-selection.Rmd#L376-L380)
library(msigdbr)
c2_sets <- msigdbr(species = "Homo sapiens", category = "C2")
ribo_set <- union(
  ribo_set,
  c2_sets[c2_sets$gs_name == "KEGG_RIBOSOME", ]$human_gene_symbol)
is_ribo <- rownames(sce) %in% ribo_set
sex_set <- rownames(sce)[any(rowData(sce)$ENSEMBL.SEQNAME %in% c("X", "Y"))]
pseudogene_set <- rownames(sce)[
  any(grepl("pseudogene", rowData(sce)$ENSEMBL.GENEBIOTYPE))]

is_mito <- rownames(sce) %in% mito_set
is_ribo <- rownames(sce) %in% ribo_set
sce <- addPerCellQC(
  sce, 
  subsets = list(Mito = which(is_mito), Ribo = which(is_ribo)))

head(colData(sce)) %>% knitr::kable()

sce$zero_percent <- colSums(counts(sce) == 0)/nrow(sce)
summary(sce$zero_percent)

# Visualise QC metrics
p1 <- plotColData(
  sce,
  "sum",
  x = "capture",
  colour_by = "experiment",
  point_size = 0.5) +
  scale_y_log10() +
  theme(axis.text.x = element_text(size = 6)) +
  annotation_logticks(
    sides = "l",
    short = unit(0.03, "cm"),
    mid = unit(0.06, "cm"),
    long = unit(0.09, "cm"))
p2 <- plotColData(
  sce,
  "detected",
  x = "capture",
  colour_by = "experiment",
  point_size = 0.5) +
  theme(axis.text.x = element_text(size = 6))
p3 <- plotColData(
  sce,
  "subsets_Mito_percent",
  x = "capture",
  colour_by = "experiment",
  point_size = 0.5) +
  theme(axis.text.x = element_text(size = 6))
p4 <- plotColData(
  sce,
  "subsets_Ribo_percent",
  x = "capture",
  colour_by = "experiment",
  point_size = 0.5) +
  theme(axis.text.x = element_text(size = 6))

p5 <- plotColData(
  sce, 
  x = "capture", 
  y = "zero_percent", 
  colour_by = "experiment",
  point_size = 0.5) +
  theme(axis.text.x = element_text(size = 6))

((p1 | p2) / (p3 | p4) / p5) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


p1 <- plotColData(
  sce,
  "sum",
  x = "donor",
  colour_by = "experiment",
  point_size = 0.5) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90)) +
  annotation_logticks(
    sides = "l",
    short = unit(0.03, "cm"),
    mid = unit(0.06, "cm"),
    long = unit(0.09, "cm"))

p2 <- plotColData(
  sce,
  "detected",
  x = "donor",
  colour_by = "experiment",
  point_size = 0.5) +
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 0.5, 
                                   vjust = 1))
p3 <- plotColData(
  sce,
  "subsets_Mito_percent",
  x = "donor",
  colour_by = "experiment",
  point_size = 0.5) +
  theme(axis.text.x = element_text(angle = 90))
p4 <- plotColData(
  sce,
  "subsets_Ribo_percent",
  x = "donor",
  colour_by = "experiment",
  point_size = 0.5) +
  theme(axis.text.x = element_text(angle = 90))

p5 <- plotColData(
  sce, 
  x = "donor", 
  y = "zero_percent", 
  colour_by = "experiment",
  point_size = 0.5) +
  theme(axis.text.x = element_text(angle = 90))

((p1 | p2) / (p3 | p4) / p5) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


plotColData(
  sce, 
  x = "sum", 
  y = "detected",
  other_fields = "capture",
  colour_by = "zero_percent",
  point_size = 0.25, 
  point_alpha = 0.25) +
  facet_wrap(vars(capture), ncol = 2)

plotColData(
  sce, 
  x = "detected", 
  y = "subsets_Mito_percent",
  other_fields = "capture",
  colour_by = "zero_percent",
  point_size = 0.25, 
  point_alpha = 0.25) +
  facet_wrap(vars(capture), ncol = 2)

plotColData(
  sce, 
  x = "sum", 
  y = "subsets_Ribo_percent",
  other_fields = "capture",
  colour_by = "zero_percent",
  point_size = 0.25, 
  point_alpha = 0.25) +
  scale_x_log10() +
  facet_wrap(vars(capture), ncol = 2)

plotColData(
  sce, 
  x = "sum", 
  y = "subsets_Mito_percent",
  other_fields = "capture",
  colour_by = "cell_status",
  point_size = 0.25, 
  point_alpha = 0.25) +
  scale_x_log10() +
  facet_wrap(vars(capture), ncol = 2)

colData(sce) %>%
  data.frame %>%
  ggplot(aes(x = cell_status, y = sum, fill = cell_status)) +
  geom_violin(scale = "count", size = 0.25) +
  scale_y_log10() +
  facet_wrap(~ capture, scales = "free_y")

# Discard uninformative genes & cells
uninformative <- is_mito | is_ribo | rownames(sce) %in% sex_set | rownames(sce) %in% pseudogene_set
sum(uninformative)

junk <- sce$donor %in% c("doublet", "unassigned") | sce$HTO %in% c("Doublet", "Unknown") | (sce$cell_status == "empty_droplet" & sce$experiment == 1)

sceFlt <- sce[!uninformative, !junk]
sceFlt

# Remove low-abundance genes
numCells <- nexprs(sceFlt, byrow = TRUE)
keep <- numCells > 20
sum(keep)

sceFlt <- sceFlt[keep,]
sceFlt

# Convert to Seurat object
counts <- counts(sceFlt)
rownames(counts) <- rowData(sceFlt)$Symbol

seu <- CreateSeuratObject(counts = counts, 
                          meta.data = data.frame(colData(sceFlt)))
seu

# Visualise combined, filtered data
DefaultAssay(seu) <- "RNA"
seu <- FindVariableFeatures(seu) %>% 
  ScaleData() %>% 
  RunPCA(verbose = FALSE, dims = 1:30) %>%
  RunUMAP(verbose = FALSE, dims = 1:30)

DimPlot(seu, group.by = "experiment", combine = FALSE)

DimPlot(seu, split.by = "experiment", combine = FALSE)

# 4. Integrate data
# Normalise the data using SCTransform and integrate across batches/individuals.
# First, import the utility.R script which contains the intDat function. utility.R can
# be found here: https://github.com/Oshlack/paed-cf-cite-seq/blob/submission2/code/utility.R
source(here("code", "utility.R"))
out <- here("data/SCEs/combined.integrated_2.SEU.rds")

if(!file.exists(out)) {
  seuInt_2 <- intDat(seu, split = "donor", type = "RNA")
  saveRDS(seuInt_2, file = out)
} else {
  seuInt_2 <- readRDS(out)
}
