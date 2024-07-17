# Unload package 'here' if it is already loaded
loaded_packages <- search()[grepl("^package:", search())]
loaded_packages <- sub("^package:", "", loaded_packages)

if ("here" %in% loaded_packages) {
  detach("package:here", unload = TRUE)
}

setwd("/folder/containing/cellranger/output")

# 1 Load libraries
library(BiocStyle)
library(tidyverse)
library(here)
library(glue)
library(DropletUtils)
library(scran)
library(scater)
library(scuttle)
library(janitor)
library(cowplot)
library(patchwork)
library(scales)
library(DropletQC)
library(EnsDb.Hsapiens.v86)
library(ensembldb)

num_features <- 16

for (i in 1:num_features) {

# 2 Load data
sce <- readRDS(
  here("data", "SCEs", glue::glue("dataset_demultiplexed.{i}.SCE.rds")))

# Some useful colours
hto_colours <- setNames(
  unique(sce$colours$hto_colours),
  unique(names(sce$colours$hto_colours)))
genetic_donor_colours <- setNames(
  unique(sce$colours$genetic_donor_colours),
  unique(names(sce$colours$genetic_donor_colours)))
capture_colours <- setNames(
  unique(sce$colours$capture_colours),
  unique(names(sce$colours$capture_colours)))

p1 <- ggcells(sce) + 
  geom_bar(
    aes(x = genetic_donor, fill = HTO),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  theme_cowplot(font_size = 10) + 
  scale_fill_manual(values = hto_colours)

p2 <- ggcells(sce) + 
  geom_bar(
    aes(x = genetic_donor, fill = capture),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  theme_cowplot(font_size = 10) + 
  scale_fill_manual(values = capture_colours)

p3 <- ggcells(sce) + 
  geom_bar(aes(x = genetic_donor, fill = genetic_donor)) + 
  coord_flip() + 
  ylab("Number of droplets") + 
  theme_cowplot(font_size = 10) + 
  scale_fill_manual(values = genetic_donor_colours) +
  geom_text(stat='count', aes(x = genetic_donor, label=..count..), hjust=1.5, size=2) +
  guides(fill = FALSE)

p1 / p2 / p3 

# Incorporating gene-based annotation
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
# Add chromosome location so we can filter on mitochondrial genes.
location <- mapIds(
  x = EnsDb.Hsapiens.v86,
  # NOTE: Need to remove gene version number prior to lookup.
  keys = rowData(sce)$ID,
  keytype = "GENEID",
  column = "SEQNAME")
rowData(sce)$CHR <- location
# Additional gene metadata from ENSEMBL and NCBI
# NOTE: These columns were customised for this project.
ensdb_columns <- c(
  "GENEBIOTYPE", "GENENAME", "GENESEQSTART", "GENESEQEND", "SEQNAME", "SYMBOL")
names(ensdb_columns) <- paste0("ENSEMBL.", ensdb_columns)
stopifnot(all(ensdb_columns %in% columns(EnsDb.Hsapiens.v86)))
ensdb_df <- DataFrame(
  lapply(ensdb_columns, function(column) {
    mapIds(
      x = EnsDb.Hsapiens.v86,
      keys = rowData(sce)$ID,
      keytype = "GENEID",
      column = column,
      multiVals = "CharacterList")
  }),
  row.names = rowData(sce)$ID)
# NOTE: Can't look up GENEID column with GENEID key, so have to add manually.
ensdb_df$ENSEMBL.GENEID <- rowData(sce)$ID
# NOTE: Homo.sapiens combines org.Hs.eg.db and
#       TxDb.Hsapiens.UCSC.hg19.knownGene (as well as others) and therefore
#       uses entrez gene and RefSeq based data.
library(Homo.sapiens)
# NOTE: These columns were customised for this project.
ncbi_columns <- c(
  # From TxDB: None required
  # From OrgDB
  "ALIAS", "ENTREZID", "GENENAME", "REFSEQ", "SYMBOL")
names(ncbi_columns) <- paste0("NCBI.", ncbi_columns)
stopifnot(all(ncbi_columns %in% columns(Homo.sapiens)))
ncbi_df <- DataFrame(
  lapply(ncbi_columns, function(column) {
    mapIds(
      x = Homo.sapiens,
      keys = rowData(sce)$ID,
      keytype = "ENSEMBL",
      column = column,
      multiVals = "CharacterList")
  }),
  row.names = rowData(sce)$ID)
rowData(sce) <- cbind(rowData(sce), ensdb_df, ncbi_df)
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

head(rowData(sce)) %>%
  knitr::kable()

# 3 Quality control
# Define and visualise the quality control metrics
is_mito <- rownames(sce) %in% mito_set
summary(is_mito)
is_ribo <- rownames(sce) %in% ribo_set
summary(is_ribo)
sce <- addPerCellQC(
  sce, 
  subsets = list(Mito = which(is_mito), Ribo = which(is_ribo)))

p1 <- plotColData(
  sce,
  "sum",
  x = "genetic_donor",
  other_fields = c("capture", "HTO"),
  colour_by = "genetic_donor",
  point_size = 0.5) +
  scale_y_log10() +
  scale_colour_manual(values = genetic_donor_colours, name = "genetic_donor") +
  theme(axis.text.x = element_blank()) +
  annotation_logticks(
    sides = "l",
    short = unit(0.03, "cm"),
    mid = unit(0.06, "cm"),
    long = unit(0.09, "cm"))
p2 <- plotColData(
  sce,
  "detected",
  x = "genetic_donor",
  other_fields = c("capture", "HTO"),
  colour_by = "genetic_donor",
  point_size = 0.5) +
  scale_colour_manual(values = genetic_donor_colours, name = "genetic_donor") +
  theme(axis.text.x = element_blank())
p3 <- plotColData(
  sce,
  "subsets_Mito_percent",
  x = "genetic_donor",
  other_fields = c("capture", "HTO"),
  colour_by = "genetic_donor",
  point_size = 0.5) +
  scale_colour_manual(values = genetic_donor_colours, name = "genetic_donor") +
  theme(axis.text.x = element_blank())
p4 <- plotColData(
  sce,
  "subsets_Ribo_percent",
  x = "genetic_donor",
  other_fields = c("capture", "HTO"),
  colour_by = "genetic_donor",
  point_size = 0.5) +
  scale_colour_manual(values = genetic_donor_colours, name = "genetic_donor") +
  theme(axis.text.x = element_blank())
p1 + p2 + p3 + p4 + plot_layout(guides = "collect", ncol = 2)

p1 + facet_grid(capture ~ HTO) + 
  theme(legend.position="bottom", text = element_text(size=6), axis.text.y = element_text(size=6)) 

p2 + facet_grid(capture ~ HTO) + 
  theme(legend.position="bottom", text = element_text(size=6), axis.text.y = element_text(size=6)) 

p3 + facet_grid(capture ~ HTO) + 
  theme(legend.position="bottom", text = element_text(size=6), axis.text.y = element_text(size=6)) 

p4 + facet_grid(capture ~ HTO) + 
  theme(legend.position="bottom", text = element_text(size=6), axis.text.y = element_text(size=6)) 

# Identify outliers by each metric
sce$batch <- sce$genetic_donor

mito_drop <- isOutlier(
  metric = sce$subsets_Mito_percent, 
  nmads = 3, 
  type = "higher",
  batch = sce$batch,
  subset = !grepl("unassigned", sce$genetic_donor))
mito_drop_df <- data.frame(
  sample = factor(
    colnames(attributes(mito_drop)$thresholds),
    levels(sce$batch)),
  lower = attributes(mito_drop)$thresholds["higher", ])
ribo_drop <- isOutlier(
  metric = sce$subsets_Ribo_percent, 
  nmads = 3, 
  type = "higher",
  batch = sce$batch,
  subset = !grepl("unassigned", sce$genetic_donor))
ribo_drop_df <- data.frame(
  sample = factor(
    colnames(attributes(ribo_drop)$thresholds),
    levels(sce$batch)),
  lower = attributes(ribo_drop)$thresholds["higher", ])

qc_cutoffs_df <- dplyr::inner_join(mito_drop_df, ribo_drop_df, by = "sample")
colnames(qc_cutoffs_df) <- c("batch", "%mito", "%ribo")
inner_join(
  qc_cutoffs_df,
  distinct(as.data.frame(colData(sce)[, c("batch"), drop = FALSE])),
  by = "batch") %>%
  dplyr::select(batch, everything()) %>%
  arrange(batch) %>%
  knitr::kable(caption = "Sample-specific QC metric cutoffs", digits = 1)

sce_pre_QC_outlier_removal <- sce
# TODO: Decide if excluding based on ribo
keep <- !mito_drop
keep <- keep & !is.na(keep)
sce_pre_QC_outlier_removal$keep <- keep
sce <- sce[, keep]
data.frame(
  ByMito = tapply(
    mito_drop, 
    sce_pre_QC_outlier_removal$batch, 
    sum,
    na.rm = TRUE),
  Remaining = as.vector(unname(table(sce$batch))),
  PercRemaining = round(
    100 * as.vector(unname(table(sce$batch))) /
      as.vector(
        unname(
          table(sce_pre_QC_outlier_removal$batch))), 1)) %>%
  tibble::rownames_to_column("batch") %>%
  dplyr::arrange(dplyr::desc(PercRemaining)) %>%
  knitr::kable(
    caption = "Number of samples removed by each QC step and the number of samples remaining.")

# Checking for removal of biologically relevant subpopulations
lost <- calculateAverage(counts(sce_pre_QC_outlier_removal)[, !keep])
kept <- calculateAverage(counts(sce_pre_QC_outlier_removal)[, keep])
library(edgeR)
logged <- cpm(cbind(lost, kept), log = TRUE, prior.count = 2)
logFC <- logged[, 1] - logged[, 2]
abundance <- rowMeans(logged)

is_mito <- rownames(sce) %in% mito_set
is_ribo <- rownames(sce) %in% ribo_set

# Check to see if genes with a log-FC (lost/kept) above 1 might be biologically relevant
high_logfc <- high_logfc[!(rownames(sce)[logFC > 1] %in% mito_set)]
high_logfc <- names(high_logfc)
is_high_logfc <- rownames(sce) %in% high_logfc

par(mfrow = c(1, 1))
plot(
  abundance,
  logFC,
  xlab = "Average count",
  ylab = "Log-FC (lost/kept)",
  pch = 16)
points(
  abundance[is_mito],
  logFC[is_mito],
  col = "dodgerblue",
  pch = 16,
  cex = 1)
points(
  abundance[is_ribo],
  logFC[is_ribo],
  col = "orange",
  pch = 16,
  cex = 1)
# I added this bit too
points(
  abundance[is_high_logfc],
  logFC[is_high_logfc],
  col = "red",
  pch = 16,
  cex = 1)

abline(h = c(-1, 1), col = "red", lty = 2)

# Check whether the cells removed during QC preferentially derive from 
# particular experimental groups
ggcells(sce_pre_QC_outlier_removal) +
  geom_bar(aes(x = genetic_donor, fill = keep)) + 
  ylab("Number of droplets") + 
  theme_cowplot(font_size = 7) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(genetic_donor ~ ., scales = "free_y")

# Compare QC metrics of the discarded and retained droplets
p1 <- plotColData(
  sce_pre_QC_outlier_removal,
  "sum",
  x = "genetic_donor",
  colour_by = "keep",
  point_size = 0.5) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  annotation_logticks(
    sides = "l",
    short = unit(0.03, "cm"),
    mid = unit(0.06, "cm"),
    long = unit(0.09, "cm"))
p2 <- plotColData(
  sce_pre_QC_outlier_removal,
  "detected",
  x = "genetic_donor",
  colour_by = "keep",
  point_size = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3 <- plotColData(
  sce_pre_QC_outlier_removal,
  "subsets_Mito_percent",
  x = "genetic_donor",
  colour_by = "keep",
  point_size = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p4 <- plotColData(
  sce_pre_QC_outlier_removal,
  "subsets_Ribo_percent",
  x = "genetic_donor",
  colour_by = "keep",
  point_size = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1 + p2 + p3 + p4 + plot_layout(guides = "collect")

# Summary
p1 <- plotColData(
  sce,
  "sum",
  x = "genetic_donor",
  other_fields = c("capture", "genetic_donor"),
  colour_by = "genetic_donor",
  point_size = 0.5) +
  scale_y_log10() +
  scale_colour_manual(values = genetic_donor_colours, name = "genetic_donor") +
  theme(axis.text.x = element_blank()) +
  annotation_logticks(
    sides = "l",
    short = unit(0.03, "cm"),
    mid = unit(0.06, "cm"),
    long = unit(0.09, "cm"))
p2 <- plotColData(
  sce,
  "detected",
  x = "genetic_donor",
  other_fields = c("capture", "genetic_donor"),
  colour_by = "genetic_donor",
  point_size = 0.5) +
  scale_colour_manual(values = genetic_donor_colours, name = "genetic_donor") +
  theme(axis.text.x = element_blank())
p3 <- plotColData(
  sce,
  "subsets_Mito_percent",
  x = "genetic_donor",
  other_fields = c("capture", "genetic_donor"),
  colour_by = "genetic_donor",
  point_size = 0.5) +
  scale_colour_manual(values = genetic_donor_colours, name = "genetic_donor") +
  theme(axis.text.x = element_blank())
p4 <- plotColData(
  sce,
  "subsets_Ribo_percent",
  x = "genetic_donor",
  other_fields = c("capture", "genetic_donor"),
  colour_by = "genetic_donor",
  point_size = 0.5) +
  scale_colour_manual(values = genetic_donor_colours, name = "genetic_donor") +
  theme(axis.text.x = element_blank())
p1 + p2 + p3 + p4 + plot_layout(guides = "collect", ncol = 2)

p1 <- ggcells(sce) + 
  geom_bar(
    aes(x = genetic_donor, fill = HTO),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  theme_cowplot(font_size = 10) + 
  scale_fill_manual(values = hto_colours)
p2 <- ggcells(sce) + 
  geom_bar(
    aes(x = genetic_donor, fill = capture),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  theme_cowplot(font_size = 10) + 
  scale_fill_manual(values = capture_colours)
p3 <- ggcells(sce) + 
  geom_bar(aes(x = genetic_donor, fill = genetic_donor)) + 
  coord_flip() + 
  ylab("Number of droplets") + 
  theme_cowplot(font_size = 10) + 
  scale_fill_manual(values = genetic_donor_colours) +
  geom_text(stat='count', aes(x = genetic_donor, label=..count..), hjust=1.5, size=2) +
  guides(fill = FALSE)
p1 / p2 / p3 + plot_layout(guides = "collect")

# 4 Identify further empty droplets using DropletQC
# Calculate nuclear fraction score for every barcode
out <- here("data", "SCEs", glue::glue("dataset_nuclear_fraction_calls.{i}.SCE.rds"))
capture_names <- levels(sce$capture)

if(!file.exists(out)){
  folder <- here()
  
  nf <- lapply(capture_names, function(cn){
    message(cn)
    
    bam <- file.path(folder,
                     glue::glue("dataset_{i}/outs/per_sample_outs/dataset_{i}/count"),
                     "sample_alignments.bam")
    
    bai <- file.path(folder, 
                     glue::glue("dataset_{i}/outs/per_sample_outs/dataset_{i}/count"),
                     "sample_alignments.bam.bai")
    
    nuclear_fraction_tags(bam = bam,
                          bam_index = bai,
                          verbose = TRUE,
                          barcodes = sce[["Barcode"]][sce$Sample == cn])
  })
  names(nf) <- capture_names
  
  dplyr::bind_rows(lapply(nf, 
                          tibble::rownames_to_column, 
                          var = "barcode")) %>%
    dplyr::mutate(cell = paste(rep(seq_along(nf), 
                                   vapply(nf, nrow, 1L)),
                               barcode, sep = "_")) -> nf_dat
  saveRDS(nf_dat, file = out)
  
} else {
  nf_dat <- readRDS(out)
  
}

# Visualise the result
colData(sce) %>% 
  data.frame %>% 
  rownames_to_column(var = "barcode") %>%
  inner_join(nf_dat %>% 
               dplyr::select(barcode, nuclear_fraction)) %>%
  column_to_rownames(var = "barcode") %>%
  DataFrame -> colData(sce)

ed <- lapply(capture_names, function(cn){
  identify_empty_drops(nf_umi = colData(sce) %>% data.frame %>%
                         dplyr::filter(Sample == cn) %>%
                         dplyr::select(nuclear_fraction, sum),
                       include_plot = TRUE,
                       umi_rescue = 10^2.5,
                       nf_rescue = 0.02)
})

# Add DropletQC calls to colData
colData(sce) %>% 
  data.frame %>% 
  rownames_to_column(var = "cell") %>%
  inner_join(ed %>% bind_rows %>%
               rownames_to_column(var = "cell")) %>%
  column_to_rownames(var = "cell") %>%
  DataFrame -> colData(sce)

p <- ggcells(sce, aes(x = nuclear_fraction, y = sum)) +
  aes(colour = cell_status) +
  geom_point(size = 1, alpha = 0.25) +
  scale_y_log10() +
  facet_wrap(~Sample) +
  theme(legend.position = "bottom")
p

plotColData(sce, x = "Sample", 
            y = "sum",
            colour_by = "cell_status") +
  scale_y_log10()

# 5 Save data
out <-  here("data", "SCEs", glue::glue("dataset.preprocessed.{i}.SCE.rds"))

if(!file.exists(out)) saveRDS(sce, out)
}