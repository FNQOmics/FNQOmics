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

num_features <- 16

for (i in 1:num_features) {

# 3 Setting up the data
sce <- readRDS(here("data", "SCEs", glue::glue("dataset_emptyDrops.{i}.SCE.rds")))

# 5 Assign Hashing Tag Oligos (HTOs)
# Preparing HTO data -----------------------------------------------------------
hto_counts <- counts(altExp(sce, "HTO"))

# Demultiplexing HTOs (DropletUtils) -------------------------------------------
confident.min <- 2
hash_stats <- hashedDrops(hto_counts, confident.min = confident.min)

hist(
  hash_stats$LogFC,
  xlab = "Log fold-change from best to second HTO",
  main = "")
abline(v = confident.min, col = "red", lty = 2, lwd = 2)

# Numbers of cells confidently assigned to each hashtag
table(Raw = hash_stats$Best, Confident = hash_stats$Confident)

# Plot HTO logFCs for all cells. Doublets show up in the top-left, singlets in the bottom right.
plot(hash_stats$LogFC, hash_stats$LogFC2, pch=".")
points(hash_stats$LogFC[hash_stats$Confident], 
       hash_stats$LogFC2[hash_stats$Confident], 
       pch=".", 
       col = "red")
abline(v = 2, h = 2, col = "blue", lty = 2)

# Visualise HTO demultiplexing results.
tmp <- SingleCellExperiment(
  assays = list(counts = hto_counts),
  colData = hash_stats)
sf <- librarySizeFactors(counts(tmp))
# NOTE: Kludge to workaround non-positive size factors (corresponding to
#       droplets with zero HTO counts).
sf <- pmax(sf, min(sf[sf > 0]))
sizeFactors(tmp) <- sf
tmp <- logNormCounts(tmp)
tmp$Best <- factor(tmp$Best)
tmp$Confident <- factor(tmp$Confident)
tmp$Doublet <- factor(tmp$Doublet)

plotHeatmap(
  tmp,
  rownames(tmp),
  order_columns_by = c("Best", "Confident", "Doublet"),
  center = TRUE,
  symmetric = TRUE,
  color = hcl.colors(101, "Blue-Red 3"),
  cluster_rows = FALSE)

# Add HTO data to SCE and save
stopifnot(identical(colnames(sce), rownames(hash_stats)))
colData(sce) <- cbind(colData(sce), hash_stats)
sce$HTO <- factor(
  case_when(
    sce$Confident ~ rownames(hto_counts)[sce$Best],
    sce$Doublet ~ "Doublet",
    TRUE ~ "Unknown"),
  levels = c(paste0("Hashtag", c(1:4)), "Doublet", "Unknown"))

dim(sce)

# Add within sample doublet calls to the HTO column
within_sample_doublets <- readRDS(here("data", "SCEs", glue::glue("dataset_experiment2_doublets_{i}.rds")))
sce$HTO <- as.character(sce$HTO)
doublet_indices <- within_sample_doublets$scDblFinder.class == "doublet"
sce$HTO[doublet_indices] <- "Doublet"
sce$HTO <- factor(sce$HTO, levels = c(levels(sce$HTO), "Hashtag1", "Hashtag2", "Hashtag3", "Hashtag4", "Doublet", "Unknown"))

# Plot the number and proportion of droplets assigned to each HTO in 
# each capture and if these were confidently or not confidently assigned.
p1 <- ggcells(sce) + 
  geom_bar(
    aes(x = Best, fill = Confident), 
    position = position_stack(reverse = TRUE)) + 
  coord_flip() +
  ylab("Number of droplets") +
  theme_cowplot(font_size = 7)
p2 <- ggcells(sce) + 
  geom_bar(
    aes(x = Best, fill = Confident), 
    position = position_fill(reverse = TRUE)) + 
  coord_flip() + 
  ylab("Proportion of droplets") + 
  theme_cowplot(font_size = 7)
(p1 + p1 + facet_grid(~capture) + plot_layout(widths = c(1, 2))) / 
  (p2 + p2 + facet_grid(~capture) + plot_layout(widths = c(1, 2))) +
  plot_layout(guides = "collect")

# 6 Demultiplexing cells without genotype reference
capture_names <- levels(sce$capture)
capture_names <- setNames(capture_names, capture_names)

vireo_df <- do.call(
  rbind,
  lapply(capture_names, function(cn) {
    vireo_df <- read.table(
      here("data", "vireo", glue("dataset_{i}"), "donor_ids.tsv"),
      header = TRUE)
    vireo_df$donor_id <- paste0(cn, vireo_df$donor_id)
    captureNumber <- sub(glue::glue("dataset_{i}"), "", cn)
    vireo_df$colname <- paste0(captureNumber, vireo_df$cell)
    # NOTE: Reorder so matches SCE.
    j <- match(colnames(sce)[sce$capture == cn], vireo_df$colname)
    vireo_df <- vireo_df[j, ]
  }))

# Summarise the genotype-based demultiplexing results
knitr::kable(
  tabyl(vireo_df, donor_id) %>%
    adorn_pct_formatting(1),
  caption = "Assignment of droplets to donors using vireo.")

# Assign genetic donor to sce
donor_mapping <- c(
  "Maurice_GEX_Feature_{i}donor0" = "donor_A",
  "Maurice_GEX_Feature_{i}donor1" = "donor_B",
  "Maurice_GEX_Feature_{i}donor2" = "donor_C",
  "Maurice_GEX_Feature_{i}donor3" = "donor_D",
  "Maurice_GEX_Feature_{i}doublet" = "doublet",
  "Maurice_GEX_Feature_{i}unassigned" = "unassigned"
)

# Create the new column using the mapping
genetic_donor_vector <- donor_mapping[as.character(sce$vireo$donor_id)]

# Convert to factor with specified levels
sce$genetic_donor <- factor(genetic_donor_vector,
                            levels = c("donor_A", "donor_B", "donor_C", "donor_D", "doublet", "unassigned"))

# 7 Summary
# Some useful colours
sce$colours <- S4Vectors::make_zero_col_DFrame(ncol(sce))

hto_colours <- setNames(
  palette.colors(nlevels(sce$HTO), "Paired"),
  levels(sce$HTO))
sce$colours$hto_colours <- hto_colours[sce$HTO]

genetic_donor_colours <- setNames(
  c(palette.colors(nlevels(sce$genetic_donor)-1, "Tableau 10"), "#000000"),
  levels(sce$genetic_donor))
sce$colours$genetic_donor_colours <- genetic_donor_colours[sce$genetic_donor]

capture_colours <- setNames(
  palette.colors(nlevels(sce$capture), "Accent"),
  levels(sce$capture))
sce$colours$capture_colours <- capture_colours[sce$capture]

tabyl(
  as.data.frame(colData(sce)[, c("HTO", "genetic_donor")]),
  HTO,
  genetic_donor) %>%
  adorn_title(placement = "combined") %>%
  adorn_totals("both") %>%
  knitr::kable(caption = "Number of droplets assigned to each `HTO`/`donor` combination.")

# Breakdown by capture
p1 <- ggcells(sce) + 
  geom_bar(
    aes(x = capture, fill = HTO),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  theme_cowplot(font_size = 8) + 
  scale_fill_manual(values = hto_colours)
p2 <- ggcells(sce) + 
  geom_bar(
    aes(x = capture, fill = genetic_donor),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  theme_cowplot(font_size = 8) + 
  scale_fill_manual(values = genetic_donor_colours)
p3 <- ggcells(sce) + 
  geom_bar(aes(x = capture, fill = capture)) + 
  coord_flip() + 
  ylab("Number of droplets") + 
  theme_cowplot(font_size = 8) + 
  scale_fill_manual(values = capture_colours) +
  guides(fill = FALSE)
p1 + p2 + p3 + plot_layout(guides = "collect")

# Breakdown by HTO
p1 <- ggcells(sce) + 
  geom_bar(
    aes(x = HTO, fill = genetic_donor),
    position = position_fill(reverse = TRUE)) +
  coord_flip() +
  ylab("Frequency") +
  theme_cowplot(font_size = 8) + 
  scale_fill_manual(values = genetic_donor_colours)
p2 <- ggcells(sce) + 
  geom_bar(aes(x = HTO, fill = HTO)) + 
  coord_flip() + 
  ylab("Number of droplets") + 
  theme_cowplot(font_size = 8) + 
  scale_fill_manual(values = hto_colours) +
  guides(fill = FALSE)
(p1 + p1 + facet_grid(~capture) + plot_layout(widths = c(1, 2))) /
  (p2 + p2 + facet_grid(~capture) + plot_layout(widths = c(1, 2))) +
  plot_layout(guides = "collect")

# Save the data
out <- here("data",
            "SCEs", 
            glue::glue("dataset_demultiplexed.{i}.SCE.rds")
)

if(!file.exists(out)) saveRDS(sce, out)
}
