
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(reshape2)
library(cowplot)
library(umap)
library(Rtsne)
library(pheatmap)
library(ggrepel)
library(ggplot2)
library(SpatialDecon)
library(knitr)
library(stringr)

#BiocManager::install("SpatialDecon")

setwd("/Users/jd116080/myprojects/Matt/GenomX/NanoString_TB/model3/")

datadir2 <- "/Users/jd116080/myprojects/Matt/GenomX/NanoString_TB/"

DCCFiles <- dir(file.path(datadir2, "./dcc"), pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
PKCFiles <- unzip(zipfile = dir(file.path(datadir2, "pkcs"), pattern = ".zip$",full.names = TRUE, recursive = TRUE))
SampleAnnotationFile <- dir(file.path(datadir2, "annotation"), pattern = ".xlsx$", full.names = TRUE, recursive = TRUE)


TBData <- readNanoStringGeoMxSet(dccFiles = DCCFiles, pkcFiles = PKCFiles, 
                                 phenoDataFile = SampleAnnotationFile, 
                                 phenoDataSheet = "TB_Geomx_annoations",                                                  phenoDataDccColName = "Sample_ID", 
                                 protocolDataColNames = c("aoi", "roi"),                                                  experimentDataColNames = c("panel"))


modules <- gsub(".pkc", "", PKCFiles)
modules <- gsub("./", "", modules)

TBData <- shiftCountsOne(TBData,useDALogic = TRUE)
QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 80,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 2,   # Minimum negative control counts (10)
       maxNTCCount = 1000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 100,         # Minimum # of nuclei estimated (100)
       minArea = 5000)         # Minimum segment area (5000)

TBData <- setSegmentQCFlags(TBData,qcCutoffs = QC_params)
QCResults <- protocolData(TBData)[["QCFlags"]] # 19 Fail min negative control counts with default but don't remove as tutorial used 1 instead of 10; I use 2

flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))


col_by <- "TreatmentGroup"
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}
QC_histogram(sData(TBData), "Trimmed (%)", col_by, 80)
QC_histogram(sData(TBData), "Stitched (%)", col_by, 80)
QC_histogram(sData(TBData), "Aligned (%)", col_by, 80)
QC_histogram(sData(TBData), "Saturated (%)", col_by, 50) + labs(title = "Sequencing Saturation (%)",x = "Sequencing Saturation (%)")
QC_histogram(sData(TBData), "area", col_by, 1000, scale_trans = "log10")

negativeGeoMeans <- esBy(negativeControlSubset(TBData), GROUP = "Module", 
                         FUN = function(x) { 
                           assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
                         }) 
protocolData(TBData)[["NegGeoMean"]] <- negativeGeoMeans

negCols <- paste0("NegGeoMean_", modules)
pData(TBData)[, negCols] <- sData(TBData)[["NegGeoMean"]]
for(ann in negCols) {
  plt <- QC_histogram(pData(TBData), ann, col_by, 2, scale_trans = "log10")
  print(plt)
}

pData(TBData) <- pData(TBData)[, !colnames(pData(TBData)) %in% negCols]

# show all NTC values, Freq = # of Segments with a given NTC count:
kable(table(NTC_Count = sData(TBData)$NTC),
      col.names = c("NTC Count", "# of Segments"))
kable(QC_Summary, caption = "QC Summary Table for each Segment")


TBData <- setBioProbeQCFlags(TBData, qcCutoffs = list(minProbeRatio = 0.1, percentFailGrubbs = 20), removeLocalOutliers = TRUE)
ProbeQCResults <- fData(TBData)[["QCFlags"]]
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),Global = sum(ProbeQCResults$GlobalGrubbsOutlier),Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0 & !ProbeQCResults$GlobalGrubbsOutlier))
qc_df
ProbeQCPassed <- subset(TBData, fData(TBData)[["QCFlags"]][,c("LowProbeRatio")] == FALSE & fData(TBData)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)
TBData <- ProbeQCPassed #Removed few target

target_TBData <- aggregateCounts(TBData)
cutoff <- 2
minLOQ <- 2
LOQ <- data.frame(row.names = colnames(target_TBData))
for(module in modules) { vars <- paste0(c("NegGeoMean_", "NegGeoSD_"), module)
if(all(vars[1:2] %in% colnames(pData(target_TBData)))) { LOQ[, module] <- pmax(minLOQ, pData(target_TBData)[, vars[1]] * pData(target_TBData)[, vars[2]] ^ cutoff) } 
}

pData(target_TBData)$LOQ <- LOQ

LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_TBData)$Module == module
  Mat_i <- t(esApply(target_TBData[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_TBData)$TargetName, ]

# Save detection rate information to pheno data
pData(target_TBData)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_TBData)$GeneDetectionRate <-
  pData(target_TBData)$GenesDetected / nrow(target_TBData)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_TBData)$DetectionThreshold <- 
  cut(pData(target_TBData)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
ggplot(pData(target_TBData),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = ROIType)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")


LOQ_Mat <- LOQ_Mat[, colnames(target_TBData)]
fData(target_TBData)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_TBData)$DetectionRate <- fData(target_TBData)$DetectedSegments / nrow(pData(target_TBData))


plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <- unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                                    function(x) {sum(fData(target_TBData)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_TBData))
rownames(plot_detect) <- plot_detect$Freq

ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) + geom_bar(stat = "identity") + geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")), vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue", high = "dodgerblue3", midpoint = 0.65, limits = c(0,1),labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1), expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments", y = "Genes Detected, % of Panel > LOQ")


negativeProbefData <- subset(fData(target_TBData), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_TBData <- target_TBData[fData(target_TBData)$DetectionRate >= 0.1 | fData(target_TBData)$TargetName %in% neg_probes, ] 
dim(target_TBData) #Features  Samples \n 19268       36 -> Remove ~900 genes 

ann_of_interest <- "TreatmentGroup"
Stat_data <- 
  data.frame(row.names = colnames(exprs(target_TBData)),
             Segment = colnames(exprs(target_TBData)),
             Annotation = pData(target_TBData)[, ann_of_interest],
             Q3 = unlist(apply(exprs(target_TBData), 2,
                               quantile, 0.75, na.rm = TRUE)),
             NegProbe = exprs(target_TBData)[neg_probes, ])
Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic", value.name = "Value")

plt1 <- ggplot(Stat_data_m,
               aes(x = Value, fill = Statistic)) +
  geom_histogram(bins = 40) + theme_bw() +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~Annotation, nrow = 1) + 
  scale_fill_brewer(palette = 3, type = "qual") +
  labs(x = "Counts", y = "Segments, #")

plt2 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3, color = Annotation)) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
  geom_point() + guides(color = "none") + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

plt3 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
  geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
  geom_point() + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")

btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                     rel_widths = c(0.43,0.57))
plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))


#Plotting shows good separation for Q3 so proceed with this
target_TBData <- normalize(target_TBData , norm_method = "quant", desiredQuantile = .75, toElt = "q_norm")
target_TBData <- normalize(target_TBData ,norm_method = "neg", fromElt = "exprs",toElt = "neg_norm")

pdf(file="Normalise.pdf")
boxplot(exprs(target_TBData)[,1:36],col = "#9EDAE5", main = "Raw Counts",log = "y", names = 1:36, xlab = "Segment",ylab = "Counts, Raw")
boxplot(assayDataElement(target_TBData[,1:36], elt = "q_norm"),col = "#2CA02C", main = "Q3 Norm Counts",log = "y", names = 1:36, xlab = "Segment",ylab = "Counts, Q3 Normalized") 
dev.off()

custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42
umap_out <- umap(t(log2(assayDataElement(target_TBData , elt = "q_norm"))), config = custom_umap)
pData(target_TBData)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
pdf("Umap.pdf")
ggplot(pData(target_TBData), aes(x = UMAP1, y = UMAP2, color = TreatmentGroup, shape = ROIType)) + geom_point(size = 3) + theme_bw()
dev.off()

#View most variable genes (coefficient of variation)
assayDataElement(object = target_TBData, elt = "log_q") <- assayDataApply(target_TBData, 2, FUN = log, base = 2, elt = "q_norm")
calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- assayDataApply(target_TBData, elt = "log_q", MARGIN = 1, calc_CV)
sort(CV_dat, decreasing = TRUE)[1:5]
write.table(sort(CV_dat, decreasing = TRUE),file="Top100_variable_genes.tsv",sep="\t",row.names = FALSE)

dim(target_TBData)
pData(target_TBData)$testRegion <- factor(pData(target_TBData)$ROIType, c("lesion", "edge.lesion","normal_cd3_edge"))
pData(target_TBData)[["slide"]] <- factor(pData(target_TBData)[["slide name"]])
pData(target_TBData)$testClass <- factor(pData(target_TBData)$TreatmentGroup,c("aCD4","BCG","dos-R1","Mtb","Naive","PE25"))
pData(target_TBData)$combinedClass <- factor(pData(target_TBData)$FullGroup,c("aCD4_edge_lesion","aCD4_lesion","aCD4_normal_cd3_edge","BCG_edge_lesion","BCG_lesion","BCG_normal_cd3_edge","dosR1_edge_lesion","dosR1_lesion","Mtb_edge_lesion","Mtb_lesion","Mtb_normal_cd3_edge","Naïve_normal_cd3_edge","PE25_edge_lesion","PE25_lesion","PE25_normal_cd3_edge"))

assayDataElement(object = target_TBData, elt = "log_q") <- assayDataApply(target_TBData, 2, FUN = log, base = 2, elt = "q_norm")


#target_TBData$testRegion[,"lesion.edge"]
#model3
results <- c()
region  <-  c("edge.lesion") 
ind <- pData(target_TBData)$testRegion == region

mixedOutmc <- mixedModelDE(target_TBData[,ind],
                           elt = "log_q",
                           modelFormula = ~ testClass + (1 | slide   ),
                           groupVar = "testClass",
                           nCores = parallel::detectCores(),
                           multiCore = FALSE)

results <- do.call(rbind, mixedOutmc["lsmeans", ])
tests <- rownames(results)
results <- as.data.frame(results)
results$Contrast <- tests
results$Gene <- unlist(lapply(colnames(mixedOutmc), rep, nrow(mixedOutmc["lsmeans", ][[1]])))
results$FDR <- p.adjust(results$`Pr(>|t|)`, method = "fdr")

results$Subset <- region
results$FDR <- p.adjust(results$`Pr(>|t|)`, method = "fdr")
results <- results[, c("Gene", "Subset", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]




results$Color <- "NS or FC < 0.5"
results$Color[results$`Pr(>|t|)` < 0.05] <- "P < 0.05"
results$Color[results$FDR < 0.05] <- "FDR < 0.05"
results$Color[results$FDR < 0.001] <- "FDR < 0.001"
results$Color[abs(results$Estimate) < 0.5] <- "NS or FC < 0.5"
results$Color <- factor(results$Color,levels = c("NS or FC < 0.5", "P < 0.05","FDR < 0.05", "FDR < 0.001"))


# mixedOutmc_full <- mixedModelDE(target_TBData,
#                            elt = "log_q",
#                            modelFormula = ~ combinedClass + (1 + testRegion | slide),
#                            groupVar = "combinedClass",
#                            nCores = parallel::detectCores(),
#                            multiCore = FALSE)
# 
# full_results <- c()
# full_results <- do.call(rbind, mixedOutmc_full["lsmeans", ])
# tests <- rownames(full_results)
# full_results <- as.data.frame(full_results)
# full_results$Contrast <- tests
# full_results$Gene <- unlist(lapply(colnames(mixedOutmc_full), rep, nrow(mixedOutmc_full["lsmeans", ][[1]])))
# full_results$FDR <- p.adjust(full_results$`Pr(>|t|)`, method = "fdr")
# 
# full_results$Color <- "NS or FC < 0.5"
# full_results$Color[full_results$`Pr(>|t|)` < 0.05] <- "P < 0.05"
# full_results$Color[full_results$FDR < 0.05] <- "FDR < 0.05"
# full_results$Color[full_results$FDR < 0.001] <- "FDR < 0.001"
# full_results$Color[abs(full_results$Estimate) < 0.5] <- "NS or FC < 0.5"
# full_results$Color <- factor(full_results$Color,levels = c("NS or FC < 0.5", "P < 0.05","FDR < 0.05", "FDR < 0.001"))


#Check gene
results[results$Gene =="Stat1",]

head(results$Gene)

ggplot(pData(target_TBData),aes(x = ROIType, fill = ROIType,y = assayDataElement(target_TBData["Stat1", ],elt = "q_norm"))) +
  geom_violin() +
  geom_jitter(width = .2) +
  labs(y = "Stat1 Expression") +
  scale_y_continuous(trans = "log2") +
  facet_wrap(~TreatmentGroup) +
  theme_bw()


#Single contrast from first group (ignores area on slide)
unique(results$Contrast)


pdf(file="Volcano_six_groups.pdf")
for (pair in unique(results$Contrast)) {
  pairs <- str_split_fixed(pair," - ",2)
  x_label = paste("Enriched in ",pairs[,2], "<- log2(FC) -> Enriched in ", pairs[,1])
  results_sub <- subset(results, Contrast %in% pair)
  results_sub$invert_P <- (-log10(results_sub$`Pr(>|t|)`)) * sign(results_sub$Estimate)
  top_A <- results_sub[order(results_sub$invert_P, decreasing = TRUE)[1:5000],] 
  top_B <- results_sub[order(results_sub$invert_P, decreasing = FALSE)[1:5000],]
  top_g <- c(top_A$Gene,top_B$Gene)
  plot <- ggplot(results_sub, aes(x = Estimate, y = -log10(`Pr(>|t|)`),color = Color, label = Gene)) +
    geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
    geom_hline(yintercept = -log10(0.05), lty = "dashed") +
    geom_point() +
    labs(x = x_label,y = "Significance, -log10(P)",color = "Significance") +
    scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",`FDR < 0.05` = "lightblue",`P < 0.05` = "orange2",`NS or FC < 0.5` = "gray"),guide = guide_legend(override.aes = list(size = 4))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
    geom_text_repel(data = subset(results_sub, Gene %in% top_g & FDR < 0.05),size = 4, point.padding = 0.15, color = "black",min.segment.length = .1, box.padding = .2,max.overlaps = 50) +
    theme_bw(base_size = 16) +
    theme(legend.position = "bottom") +
    facet_wrap(~Contrast, scales = "free_y")
  print(plot)
  filename1 = paste(pairs[,1],"_vs_",pairs[,2],"_Enriched_",pairs[,1],".tsv",sep="")
  write.table(top_A,file=filename1,sep="\t",row.names = FALSE)
  filename2 = paste(pairs[,1],"_vs_",pairs[,2],"_Enriched_",pairs[,2],".tsv",sep="")
  write.table(top_B,file=filename2,sep="\t",row.names = FALSE)
}
dev.off()

for (pair in unique(results$Contrast)) {
  pairs <- str_split_fixed(pair," - ",2)
  x_label = paste("Enriched in ",pairs[,2], "<- log2(FC) -> Enriched in ", pairs[,1])
  results_sub <- subset(results, Contrast %in% pair)
  results_sub$invert_P <- (-log10(results_sub$`Pr(>|t|)`)) * sign(results_sub$Estimate)
  top_A <- results_sub[order(results_sub$invert_P, decreasing = TRUE)[1:19268],] 
  filename1 = paste(pairs[,1],"_all_vs_",pairs[,2],"_Enriched_",pairs[,1],".tsv",sep="")
  write.table(top_A,file=filename1,sep="\t",row.names = FALSE)
}



