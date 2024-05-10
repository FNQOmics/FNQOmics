# Create SingleCellExperiment object from 10X data and identify empty droplets
# using DropletUtils.
# Peter Hickey
# 2021-05-10

# Modified from https://raw.githubusercontent.com/Oshlack/paed-cf-cite-seq/submission2/code/C133_Neeland-dropletutils.R

# Setup ------------------------------------------------------------------------
# Unload package 'here' if it is already loaded
loaded_packages <- search()[grepl("^package:", search())]
loaded_packages <- sub("^package:", "", loaded_packages)

if ("here" %in% loaded_packages) {
  detach("package:here", unload = TRUE)
}

library(DropletUtils)
setwd("/g/data/pq84/single_cell/10X/ANU/scRNA_ITK_IRF5_May23")
library(here)

dir.create(here("data", "SCEs"), recursive = TRUE)
dir.create(here("data", "emptyDrops"))

# Initialize the number of iterations
num_iterations <- 16

# Construct SingleCellExperiment object ----------------------------------------
for (i in 1:num_iterations) {
  capture_names <- paste0("Maurice_GEX_Feature_", i)
  capture_names <- setNames(capture_names, capture_names)
  captures <- setNames(
    here(
      capture_names,
      "outs",
      "multi",
      "count",
      "raw_feature_bc_matrix"),
    capture_names)
  sce <- read10xCounts(samples = captures, col.names = TRUE)
  stopifnot(!anyDuplicated(colnames(sce)))
  sce <- splitAltExps(
    sce,
    rowData(sce)$Type,
    "Gene Expression")
  
  # Identify empty droplets ------------------------------------------------------
  
  set.seed(100)
  list_of_empties <- lapply(capture_names, function(cn) {
    message(cn)
    emptyDrops(counts(sce)[, sce$Sample == cn])
  })
  
  # Check if more permutations are needed; see
  # https://osca.bioconductor.org/quality-control.html#testing-for-empty-droplets
  more_permutations_needed <- sapply(list_of_empties, function(e) {
    table(
      Sig = e$FDR <= 0.001,
      Limited = e$Limited)[1, 2] > 0
  })
  stopifnot(all(!more_permutations_needed))
  
  # Compare `lower` cutoffs ------------------------------------------------------
  
  # Define your lower cutoffs and feature names
  set.seed(100)
  lower <- lower <- c(100, 50, 20)
  names(lower) <- lower
  #names(lower) <- paste0("Maurice_GEX_Feature_", 1:3)
  tmp <- lapply(c(100, 50, 20), function(lower) {
    message(lower)
    lapply(capture_names, function(cn) {
      message(cn)
      emptyDrops(counts(sce)[, sce$Sample == cn], lower = lower)
    })
  })
  tmp2 <- data.frame(
    capture_names = sapply(
      lapply(tmp, "[[", capture_names), function(x) {
        sum(x$FDR < 0.001, na.rm = TRUE)
      }))
  
  rownames(tmp2) <- paste0("emptyDrops(..., lower = ", lower, ")")
  
  # Define the location of the filtered results
  captures_filtered <- setNames(
    here(
      capture_names,
      "outs",
      "per_sample_outs",
      capture_names,
      "count",
      "sample_filtered_feature_bc_matrix"),
    capture_names)
  
  cr <- read10xCounts(captures_filtered)
  
  tmp2 <- rbind(
    tmp2,
    data.frame(
      capture_names = sum(cr$Sample == capture_names),
      row.names = "CellRanger"))
  
  # Save outputs -----------------------------------------------------------------
  output_sce_file <- paste0("Maurice_part2_", i, ".CellRanger.SCE.rds")
  output_emptydrops_file <- paste0("Maurice_part2_", i, "_varying_lower.emptyDrops.rds")
  
  
  saveRDS(
    object = sce,
    file = here(
      "data",
      "SCEs",
      output_sce_file),
    compress = "xz")
  
  for (cn in capture_names) {
    message(cn)
    empties <- list_of_empties[[cn]]
    saveRDS(
      object = empties,
      file = here(
        "data",
        "emptyDrops",
        paste0(cn, ".emptyDrops.rds")),
      compress = "xz")
    writeLines(
      text = sce[["Barcode"]][sce$Sample == cn][which(empties$FDR <= 0.001)],
      con = here(
        "data",
        "emptyDrops",
        paste0(cn, ".barcodes.txt")))
  }
  
  saveRDS(tmp2, file = here("data", "emptyDrops", output_emptydrops_file))
}
