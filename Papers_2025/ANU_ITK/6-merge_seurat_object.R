# Unload package 'here' if it is already loaded
loaded_packages <- search()[grepl("^package:", search())]
loaded_packages <- sub("^package:", "", loaded_packages)

if ("here" %in% loaded_packages) {
  detach("package:here", unload = TRUE)
}

setwd("/folder/containing/cellranger/output")

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
library(harmony)

# Import the two sce objects from the May dataset and the one from the August dataset
seu_1 <- readRDS(here("scRNA_ITK_Aug22", "data", "SCEs", "combined.integrated.SEU.rds"))
seu_2 <- readRDS(here("scRNA_ITK_IRF5_May23", "data", "SCEs", "combined.integrated_1.SEU.rds"))
seu_3 <- readRDS(here("scRNA_ITK_IRF5_May23", "data", "SCEs", "combined.integrated_2.SEU.rds"))

# Add pre/post/healthy $treatment to GEM samples
seu_1$treatment <- as.factor(
  ifelse(
    seu_1$experiment == 1, "healthy_control",
    
    ifelse(
      seu_1$experiment == 2 & seu_1$HTO == "Hashtag2", "pre_treatment",
    
      ifelse(
        seu_1$experiment == 2 & seu_1$HTO == "Hashtag1", "post_treatment",
      
        ifelse(
          seu_1$experiment == 3, "pre_treatment",
        
          ifelse(
            seu_1$experiment == 4, "post_treatment",
            NA
          )
        )
      )
    )
  )
)

# Add cohort data
seu_1$cohort <- as.factor(
  ifelse(
    seu_1$experiment == 1, "HBD",
    ifelse(
      seu_1$experiment != 1, "E42K",
      NA
    )
  )
)

seu_2$cohort <- as.factor(
  ifelse(
    (seu_2$capture %in% c("Maurice_GEX_Feature_2", "Maurice_GEX_Feature_8", "Maurice_GEX_Feature_10", "Maurice_GEX_Feature_16")) &
      seu_2$HTO == "Hashtag2", "HBD",
  
    ifelse(
      (seu_2$capture %in% c("Maurice_GEX_Feature_3", "Maurice_GEX_Feature_5", "Maurice_GEX_Feature_11", "Maurice_GEX_Feature_13")) &
        seu_2$HTO == "Hashtag3", "HBD",
    
      ifelse(
        (seu_2$capture %in% c("Maurice_GEX_Feature_4", "Maurice_GEX_Feature_6", "Maurice_GEX_Feature_12", "Maurice_GEX_Feature_14")) &
          seu_2$HTO == "Hashtag4", "HBD",
      
        ifelse(
          seu_2$capture %in% c("Maurice_GEX_Feature_1", "Maurice_GEX_Feature_7", "Maurice_GEX_Feature_9", "Maurice_GEX_Feature_15") &
            (seu_2$HTO %in% c("Hashtag1", "Hashtag2", "Hashtag3")), "E42K",
        
          ifelse(
            (seu_2$capture %in% c("Maurice_GEX_Feature_2", "Maurice_GEX_Feature_8", "Maurice_GEX_Feature_10", "Maurice_GEX_Feature_16")) &
              (seu_2$HTO %in% c("Hashtag1", "Hashtag3")), "E42K",
          
            ifelse(
              (seu_2$capture %in% c("Maurice_GEX_Feature_2", "Maurice_GEX_Feature_8", "Maurice_GEX_Feature_10", "Maurice_GEX_Feature_16")) &
                seu_2$HTO == "Hashtag4", "T504S",
            
              ifelse(
                (seu_2$capture %in% c("Maurice_GEX_Feature_3", "Maurice_GEX_Feature_5", "Maurice_GEX_Feature_11", "Maurice_GEX_Feature_13")) &
                  (seu_2$HTO %in% c("Hashtag1", "Hashtag2", "Hashtag4")), "T504S",
              
                ifelse(
                  (seu_2$capture %in% c("Maurice_GEX_Feature_4", "Maurice_GEX_Feature_6", "Maurice_GEX_Feature_12", "Maurice_GEX_Feature_14")) &
                    (seu_2$HTO %in% c("Hashtag1", "Hashtag2", "Hashtag3")), "D135Y", "Unknown"
                )
              )
            )
          )
        )
      )
    )
  )
)

seu_3$cohort <- as.factor(
  ifelse(
    (seu_3$capture %in% c("Maurice_GEX_Feature_2", "Maurice_GEX_Feature_8", "Maurice_GEX_Feature_10", "Maurice_GEX_Feature_16")) &
      seu_3$HTO == "Hashtag2", "HBD",
  
    ifelse(
      (seu_3$capture %in% c("Maurice_GEX_Feature_3", "Maurice_GEX_Feature_5", "Maurice_GEX_Feature_11", "Maurice_GEX_Feature_13")) &
        seu_3$HTO == "Hashtag3", "HBD",
    
      ifelse(
        (seu_3$capture %in% c("Maurice_GEX_Feature_4", "Maurice_GEX_Feature_6", "Maurice_GEX_Feature_12", "Maurice_GEX_Feature_14")) &
          seu_3$HTO == "Hashtag4", "HBD",
      
        ifelse(
          (seu_3$capture %in% c("Maurice_GEX_Feature_1", "Maurice_GEX_Feature_7", "Maurice_GEX_Feature_9", "Maurice_GEX_Feature_15") &
            seu_3$HTO %in% c("Hashtag1", "Hashtag2", "Hashtag3")), "E42K",
        
          ifelse(
            (seu_3$capture %in% c("Maurice_GEX_Feature_2", "Maurice_GEX_Feature_8", "Maurice_GEX_Feature_10", "Maurice_GEX_Feature_16")) &
              (seu_3$HTO %in% c("Hashtag1", "Hashtag3")), "E42K",
          
            ifelse(
              (seu_3$capture %in% c("Maurice_GEX_Feature_2", "Maurice_GEX_Feature_8", "Maurice_GEX_Feature_10", "Maurice_GEX_Feature_16")) &
                seu_3$HTO == "Hashtag4", "T504S",
            
              ifelse(
                (seu_3$capture %in% c("Maurice_GEX_Feature_3", "Maurice_GEX_Feature_5", "Maurice_GEX_Feature_11", "Maurice_GEX_Feature_13")) &
                  (seu_3$HTO %in% c("Hashtag1", "Hashtag2", "Hashtag4")), "T504S",
              
                ifelse(
                  (seu_3$capture %in% c("Maurice_GEX_Feature_4", "Maurice_GEX_Feature_6", "Maurice_GEX_Feature_12", "Maurice_GEX_Feature_14")) &
                    (seu_3$HTO %in% c("Hashtag1", "Hashtag2", "Hashtag3")), "D135Y", 
                  
                  ifelse(
                    (seu_3$capture %in% c("Maurice_GEX_Feature_1", "Maurice_GEX_Feature_7", "Maurice_GEX_Feature_9", "Maurice_GEX_Feature_15")) &
                      seu_3$HTO == "Hashtag4", "E42K_unaffected", "Unknown"
                  )
                )
              )
            )
          )
        )
      )
    )
  )
)

# Add patient data
seu_1$patient <- as.factor(
  ifelse(
    seu_1$experiment == 1 &
      seu_1$HTO == "Hashtag1", "HBD046",
    
    ifelse(
      seu_1$experiment == 1 &
        seu_1$HTO == "Hashtag2", "HBD066",
      
      ifelse(
        seu_1$experiment == 1 &
          seu_1$HTO == "Hashtag4", "HBD093",
        ifelse(
          seu_1$experiment != 1, "GEM108",
          NA
        )
      )
    )
  )
)

seu_2$patient <- as.factor(
  ifelse(
    (seu_2$capture %in% c("Maurice_GEX_Feature_1", "Maurice_GEX_Feature_7", "Maurice_GEX_Feature_9", "Maurice_GEX_Feature_15")) &
      seu_2$HTO == "Hashtag1", "GEM108",
  
    ifelse(
      (seu_2$capture %in% c("Maurice_GEX_Feature_1", "Maurice_GEX_Feature_7", "Maurice_GEX_Feature_9", "Maurice_GEX_Feature_15")) &
        seu_2$HTO == "Hashtag2", "PMAI0023",
    
      ifelse(
        (seu_2$capture %in% c("Maurice_GEX_Feature_1", "Maurice_GEX_Feature_7", "Maurice_GEX_Feature_9", "Maurice_GEX_Feature_15")) &
          seu_2$HTO == "Hashtag3", "PMAI0024",
      
        ifelse(
          (seu_2$capture %in% c("Maurice_GEX_Feature_1", "Maurice_GEX_Feature_7", "Maurice_GEX_Feature_9", "Maurice_GEX_Feature_15")) &
            seu_2$HTO == "Hashtag4", "PMAI0025",
        
          ifelse(
            (seu_2$capture %in% c("Maurice_GEX_Feature_2", "Maurice_GEX_Feature_8", "Maurice_GEX_Feature_10", "Maurice_GEX_Feature_16")) &
              seu_2$HTO == "Hashtag1", "PMAI0017",
          
            ifelse(
              (seu_2$capture %in% c("Maurice_GEX_Feature_2", "Maurice_GEX_Feature_8", "Maurice_GEX_Feature_10", "Maurice_GEX_Feature_16")) &
                seu_2$HTO == "Hashtag2", "HBD014",
            
              ifelse(
                (seu_2$capture %in% c("Maurice_GEX_Feature_2", "Maurice_GEX_Feature_8", "Maurice_GEX_Feature_10", "Maurice_GEX_Feature_16")) &
                  seu_2$HTO == "Hashtag3", "PMAI0018",
              
                ifelse(
                  (seu_2$capture %in% c("Maurice_GEX_Feature_2", "Maurice_GEX_Feature_8", "Maurice_GEX_Feature_10", "Maurice_GEX_Feature_16")) &
                    seu_2$HTO == "Hashtag4", "CPI939",
                
                  ifelse(
                    (seu_2$capture %in% c("Maurice_GEX_Feature_3", "Maurice_GEX_Feature_5", "Maurice_GEX_Feature_11", "Maurice_GEX_Feature_13")) &
                      seu_2$HTO == "Hashtag1", "CPI764",
                  
                    ifelse(
                      (seu_2$capture %in% c("Maurice_GEX_Feature_3", "Maurice_GEX_Feature_5", "Maurice_GEX_Feature_11", "Maurice_GEX_Feature_13")) &
                        seu_2$HTO == "Hashtag2", "CPI767",
                    
                      ifelse(
                        (seu_2$capture %in% c("Maurice_GEX_Feature_3", "Maurice_GEX_Feature_5", "Maurice_GEX_Feature_11", "Maurice_GEX_Feature_13")) &
                          seu_2$HTO == "Hashtag3", "HBD002",
                      
                        ifelse(
                          (seu_2$capture %in% c("Maurice_GEX_Feature_3", "Maurice_GEX_Feature_5", "Maurice_GEX_Feature_11", "Maurice_GEX_Feature_13")) &
                            seu_2$HTO == "Hashtag4", "CPI872",
                        
                          ifelse(
                            (seu_2$capture %in% c("Maurice_GEX_Feature_4", "Maurice_GEX_Feature_6", "Maurice_GEX_Feature_12", "Maurice_GEX_Feature_14")) &
                              seu_2$HTO == "Hashtag1", "CPI904",
                          
                            ifelse(
                              (seu_2$capture %in% c("Maurice_GEX_Feature_4", "Maurice_GEX_Feature_6", "Maurice_GEX_Feature_12", "Maurice_GEX_Feature_14")) &
                                seu_2$HTO == "Hashtag2", "CPI874",
                            
                              ifelse(
                                (seu_2$capture %in% c("Maurice_GEX_Feature_4", "Maurice_GEX_Feature_6", "Maurice_GEX_Feature_12", "Maurice_GEX_Feature_14")) &
                                seu_2$HTO == "Hashtag3", "CPI875",
                              
                                ifelse(
                                  (seu_2$capture %in% c("Maurice_GEX_Feature_4", "Maurice_GEX_Feature_6", "Maurice_GEX_Feature_12", "Maurice_GEX_Feature_14")) &
                                    seu_2$HTO == "Hashtag4", "HBD048", "Unknown"
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )
)

seu_3$patient <- as.factor(
  ifelse(
    (seu_3$capture %in% c("Maurice_GEX_Feature_1", "Maurice_GEX_Feature_7", "Maurice_GEX_Feature_9", "Maurice_GEX_Feature_15")) &
      seu_3$HTO == "Hashtag1", "GEM108",
  
    ifelse(
      (seu_3$capture %in% c("Maurice_GEX_Feature_1", "Maurice_GEX_Feature_7", "Maurice_GEX_Feature_9", "Maurice_GEX_Feature_15")) &
        seu_3$HTO == "Hashtag2", "PMAI0023",
    
      ifelse(
        (seu_3$capture %in% c("Maurice_GEX_Feature_1", "Maurice_GEX_Feature_7", "Maurice_GEX_Feature_9", "Maurice_GEX_Feature_15")) &
          seu_3$HTO == "Hashtag3", "PMAI0024",
      
        ifelse(
          (seu_3$capture %in% c("Maurice_GEX_Feature_1", "Maurice_GEX_Feature_7", "Maurice_GEX_Feature_9", "Maurice_GEX_Feature_15")) &
            seu_3$HTO == "Hashtag4", "PMAI0025",
        
          ifelse(
            (seu_3$capture %in% c("Maurice_GEX_Feature_2", "Maurice_GEX_Feature_8", "Maurice_GEX_Feature_10", "Maurice_GEX_Feature_16")) &
              seu_3$HTO == "Hashtag1", "PMAI0017",
          
            ifelse(
              (seu_3$capture %in% c("Maurice_GEX_Feature_2", "Maurice_GEX_Feature_8", "Maurice_GEX_Feature_10", "Maurice_GEX_Feature_16")) &
                seu_3$HTO == "Hashtag2", "HBD014",
            
              ifelse(
                (seu_3$capture %in% c("Maurice_GEX_Feature_2", "Maurice_GEX_Feature_8", "Maurice_GEX_Feature_10", "Maurice_GEX_Feature_16")) &
                  seu_3$HTO == "Hashtag3", "PMAI0018",
              
                ifelse(
                  (seu_3$capture %in% c("Maurice_GEX_Feature_2", "Maurice_GEX_Feature_8", "Maurice_GEX_Feature_10", "Maurice_GEX_Feature_16")) &
                    seu_3$HTO == "Hashtag4", "CPI939",
                
                  ifelse(
                    (seu_3$capture %in% c("Maurice_GEX_Feature_3", "Maurice_GEX_Feature_5", "Maurice_GEX_Feature_11", "Maurice_GEX_Feature_13")) &
                      seu_3$HTO == "Hashtag1", "CPI764",
                  
                    ifelse(
                      (seu_3$capture %in% c("Maurice_GEX_Feature_3", "Maurice_GEX_Feature_5", "Maurice_GEX_Feature_11", "Maurice_GEX_Feature_13")) &
                        seu_3$HTO == "Hashtag2", "CPI767",
                    
                      ifelse(
                        (seu_3$capture %in% c("Maurice_GEX_Feature_3", "Maurice_GEX_Feature_5", "Maurice_GEX_Feature_11", "Maurice_GEX_Feature_13")) &
                          seu_3$HTO == "Hashtag3", "HBD002",
                      
                        ifelse(
                          (seu_3$capture %in% c("Maurice_GEX_Feature_3", "Maurice_GEX_Feature_5", "Maurice_GEX_Feature_11", "Maurice_GEX_Feature_13")) &
                            seu_3$HTO == "Hashtag4", "CPI872",
                        
                          ifelse(
                            (seu_3$capture %in% c("Maurice_GEX_Feature_4", "Maurice_GEX_Feature_6", "Maurice_GEX_Feature_12", "Maurice_GEX_Feature_14")) &
                              seu_3$HTO == "Hashtag1", "CPI904",
                          
                            ifelse(
                              (seu_3$capture %in% c("Maurice_GEX_Feature_4", "Maurice_GEX_Feature_6", "Maurice_GEX_Feature_12", "Maurice_GEX_Feature_14")) &
                                seu_3$HTO == "Hashtag2", "CPI874",
                            
                              ifelse(
                                (seu_3$capture %in% c("Maurice_GEX_Feature_4", "Maurice_GEX_Feature_6", "Maurice_GEX_Feature_12", "Maurice_GEX_Feature_14")) &
                                  seu_3$HTO == "Hashtag3", "CPI875",
                              
                                ifelse(
                                  (seu_3$capture %in% c("Maurice_GEX_Feature_4", "Maurice_GEX_Feature_6", "Maurice_GEX_Feature_12", "Maurice_GEX_Feature_14")) &
                                    seu_3$HTO == "Hashtag4", "HBD048", "Unknown"
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )
)
  
# Add batch info
seu_1$batch <- "Aug"
seu_2$batch <- "May"
seu_3$batch <- "May"

# Identify shared features and filter seurat objects
features_1 <- rownames(seu_1@assays$integrated)
features_2 <- rownames(seu_2@assays$integrated)
features_3 <- rownames(seu_3@assays$integrated)

shared_features <- intersect(intersect(features_1, features_2), features_3)

seu_1_filtered <- subset(seu_1, features = shared_features)
seu_2_filtered <- subset(seu_2, features = shared_features)
seu_3_filtered <- subset(seu_3, features = shared_features)

seu_combo <- merge(seu_1_filtered, y = c(seu_2_filtered, seu_3_filtered), merge.data = TRUE)

# Manually set variable features
VariableFeatures(seu_combo) <- rownames(seu_combo@assays$integrated)

seu_combo <- ScaleData(seu_combo)
seu_combo <- RunPCA(seu_combo, npcs = 30)
DimPlot(seu_combo, group.by = "batch", combine = FALSE)

# Normalise for batch effects with Harmony
options(repr.plot.height = 2.5, repr.plot.width = 6)
seu_combo <- seu_combo %>% RunHarmony("batch", plot_convergene = TRUE)

DimPlot(seu_combo, reduction = "harmony", group.by = "batch", combine = FALSE)

# Perform clustering
seu_combo <- FindNeighbors(seu_combo, reduction = "harmony", dims = 1:30)
seu_combo <- FindClusters(seu_combo)

# 5. Cluster data
# Perform linear dimensional reduction
p1 <- DimPlot(seu_combo, reduction = "harmony", group.by = "cohort")
p2 <- DimPlot(seu_combo, reduction = "harmony", dims = c(1,3), group.by = "cohort")
p3 <- DimPlot(seu_combo, reduction = "harmony", dims = c(2,3), group.by = "cohort")
p4 <- DimPlot(seu_combo, reduction = "harmony", dims = c(3,4), group.by = "cohort")

((p1 | p2) / (p3 | p4)) + plot_layout(guides = "collect") &
  theme(legend.text = element_text(size = 8),
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8))

DimHeatmap(seu_combo, dims = 1:30, cells = 500, balanced = TRUE)

# Determine the dimensionality of the dataset
ElbowPlot(seu_combo, ndims = 30)

# Cluster the cells
dir.create(here("merged"))
out <- here("merged", "merged.clustered.SEU.rds")

if(!file.exists(out)) {
  seu_combo <- FindNeighbors(seu_combo, reduction = "harmony", dims = 1:30)
  seu_combo <- FindClusters(seu_combo, algorithm = 3, 
                            resolution = seq(0.1, 1, by = 0.1))
  seu_combo <- RunUMAP(seu_combo, reduction = "harmony", dims = 1:30)
  saveRDS(seu_combo, file = out)
  
} else {
  seu_combo <- readRDS(out)
  
}

# Visualise clustering at default resolution.
p1 <- DimPlot(seu_combo, reduction = 'umap', label = TRUE, repel = TRUE, group.by = "batch", label.size = 5, pt.size = 1)
p2 <- DimPlot(seu_combo, reduction = 'umap', label = TRUE, repel = TRUE, label.size = 5, pt.size = 1)

p1 + p2

# Clustree resolution visualisation
clustree(seu_combo, prefix = "integrated_snn_res.")

# 7. Annotate data using Azimuth & Human - PBMC reference
# Add Azimuth labels
# Save slimmed down (diet) Seurat object for upload to Azimuth and then append annotations
out <- here("merged", "merged.clustered_diet.SEU.rds")

# Run this a second time after creating combined_azimuth.rds
if(!file.exists(out)){
  DefaultAssay(seu_combo) <- "RNA"
  seuDiet <- DietSeurat(seu_combo, assays = "RNA")
  saveRDS(seuDiet, out)
}

## Once uploading merged.clustered_diet.SEU.rds to Azimuth, running against the PBMC dataset and 
# downloading the results (azimuth_pred.tsv and azimuth_umap.Rds), proceed with the script.
tsv <- as.data.frae(read_tsv(here("merged", "azimuth_pred.tsv")))

azimuth <- readRDS(here("merged", "azimuth_umap.Rds"))
azimuth_combined <- list(azimuth, tsv)
names(azimuth_combined) <- c("umap", "pred.df")

out <- here("merged", "combined_azimuth.rds")
saveRDS(azimuth_combined, out)

seu_combo <- AddAzimuthResults(seu_combo, filename = here("merged", "combined_azimuth.rds"))
seu_combo$predicted.celltype.l2 <- fct_drop(seu_combo$predicted.celltype.l2)

if(any(grepl("predicted.celltype.l2", colnames(seu_combo@meta.data)))) table(seu_combo$predicted.celltype.l2)

# Visualise reference mapping
p1 <- DimPlot(seu_combo, split.by = "cohort", reduction = 'umap.proj', 
        label = TRUE, repel = TRUE, label.size = 3, pt.size = 1,
        group.by = "predicted.celltype.l2") + NoLegend()

p2 <- DimPlot(seu_combo, split.by = "cohort", reduction = 'umap', 
        label = TRUE, repel = TRUE, label.size = 3, pt.size = 1,
        group.by = "predicted.celltype.l2") + NoLegend()

p1 + p2

# Plot QC metrics by label
options(scipen=1)
seu_combo@meta.data %>%
  ggplot(aes(y = predicted.celltype.l2.score,
             x = predicted.celltype.l2,
             fill = predicted.celltype.l2)) +
  geom_violin(scale = "width") +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust = 1)) +
  NoLegend() -> p1

seu_combo@meta.data %>%
  ggplot(aes(y = nCount_RNA,
             x = predicted.celltype.l2,
             fill = predicted.celltype.l2)) +
  geom_violin(scale = "area") +
  scale_y_log10() +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust = 1)) +
  NoLegend() -> p2

seu_combo@meta.data %>%
  ggplot(aes(y = nFeature_RNA,
             x = predicted.celltype.l2,
             fill = predicted.celltype.l2)) +
  geom_violin(scale = "area") +
  scale_y_log10() +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust = 1)) +
  NoLegend() -> p3

(p1 / p2 / p3)

# 10. Save the session
out <- here("merged", "combined.annotated.SEU.rds") 
if(!file.exists(out)) saveRDS(seu_combo, out)
