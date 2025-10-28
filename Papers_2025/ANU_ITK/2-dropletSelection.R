# The following is adapted from https://oshlacklab.com/paed-cf-cite-seq/03_C133_Neeland.emptyDrops.html

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
library(scds)
library(scDblFinder)
library(Matrix)

# This sets the random number generator seed to 42
set.seed(42)

# This sets the scientific notation penalty to 999
options(scipen=999)

# This sets the maximum size of global variables when using the "future" package for parallel and asynchronous programming
options(future.globals.maxSize = 6500 * 1024^2)

num_features <- 16

for (i in 1:num_features) {
  
  # Update the capture names based on the iteration
  capture_names <- paste0("dataset_", i)
  capture_names <- setNames(capture_names, capture_names)

# 2 Load the data
  sce <- readRDS(here("data/SCEs", glue::glue("dataset_{i}.CellRanger.SCE.rds")))

# Preparing HTO data -----------------------------------------------------------
is_hto <- rownames(altExp(sce, "Antibody Capture")) %in%
  paste0("Hashtag", 1:4)
altExp(sce, "HTO") <- altExp(sce, "Antibody Capture")[is_hto, ]
altExp(sce, "ADT") <- altExp(sce, "Antibody Capture")[!is_hto, ]
altExp(sce, "Antibody Capture") <- NULL

expSum <- colSums(counts(sce))
htoSum <- colSums(counts(altExp(sce, "HTO")))
adtSum <- colSums(counts(altExp(sce, "ADT")))

dat <- data.frame(exp = expSum, hto = htoSum, adt = adtSum)

# 3 Identify empty droplets
sce$capture <- factor(sce$Sample)
capture_names <- levels(sce$capture)
capture_names <- setNames(capture_names, capture_names)
empties <- do.call(rbind, lapply(capture_names, function(cn) {
  message(cn)
  empties <- readRDS(
    here("data", "emptyDrops", paste0(cn, ".emptyDrops.rds")))
  empties$capture <- cn
  empties
}))
tapply(
  empties$FDR,
  empties$capture,
  function(x) sum(x <= 0.001, na.rm = TRUE)) %>%
  knitr::kable(
    caption = "Number of non-empty droplets identified using `emptyDrops()` from **DropletUtils**.")

# 3.1 Examine results
par(mfrow = c(2, 2))
lapply(levels(sce$capture), function(s) {
  sce <- sce[, sce$capture == s]
  bcrank <- barcodeRanks(counts(sce))
  
  # Only showing unique points for plotting speed.
  uniq <- !duplicated(bcrank$rank)
  plot(
    x = bcrank$rank[uniq],
    y = bcrank$total[uniq],
    log = "xy",
    xlab = "Rank",
    ylab = "Total UMI count",
    main = s,
    cex.lab = 1.2,
    xlim = c(1, 500000),
    ylim = c(1, 200000))
  abline(h = metadata(bcrank)$inflection, col = "darkgreen", lty = 2)
  abline(h = metadata(bcrank)$knee, col = "dodgerblue", lty = 2)
})

# 4 Remove empty droplets
sce <- sce[, which(empties$FDR <= 0.001)]
sce

# 5 Call within-sample doublets
out <- here("data/SCEs", glue::glue("dataset_experiment2_doublets_{i}.rds"))

if(!file.exists(out)){
  
  sceLst <- sapply(levels(sce$capture), function(cap){
    ## Annotate doublets using scds three step process as run in Demuxafy
    sce1 <- bcds(sce[, sce$capture == cap], 
                 retRes = TRUE, estNdbl = TRUE)
    sce1 <- cxds(sce1, retRes = TRUE, estNdbl = TRUE)
    sce1 <- cxds_bcds_hybrid(sce1, estNdbl = TRUE)
    ## Annotate doublets using scDblFInder with rate estimate from Demuxafy
    sce1 <- scDblFinder(sce1, dbr = ncol(sce1)/1000*0.008)
    sce1
  })  
  
  lapply(sceLst, function(s){
    colData(s) %>% 
      data.frame %>%
      rownames_to_column(var = "cell")
  }) %>% 
    bind_rows() %>% 
    saveRDS(file = out)
} 

# 6 Save data
out <- here("data/SCEs", glue::glue("dataset_emptyDrops.{i}.SCE.rds"))
if (!file.exists(out)) saveRDS(sce, file = out)
}

