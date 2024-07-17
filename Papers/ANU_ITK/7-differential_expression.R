# Unload package 'here' if it is already loaded
loaded_packages <- search()[grepl("^package:", search())]
loaded_packages <- sub("^package:", "", loaded_packages)

if ("here" %in% loaded_packages) {
  detach("package:here", unload = TRUE)
}

setwd("/folder/containing/cellranger/output")

library(limma)
library(tidyverse)
library(scater)
library(scran)
library(data.table)
library(edgeR)
library(pheatmap)
library(plyr)
library(dplyr)
library(Seurat)
library(gtable)
library(gridExtra)
library(stringi)
library(ggplot2)
library(ggrepel)
library(foreach)
library(doParallel)
library(MOFA2)
library(MultiAssayExperiment)
library(here)

seu <- readRDS(here("merged", "combined.annotated.SEU.rds"))

################################################
#Differential expression
################################################

#Transform with Sanity
countmatrix <- as.matrix(LayerData(seu, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "countmatrix.txt")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity")
}

expr <- read.table(here("merged", "sanity", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)
cohort <- as.character(seu$cohort)
patient <- as.character(seu$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA
engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/corfit.RData")

cm <- makeContrasts(healthy_D135Y=HBD-D135Y,
                    healthy_E42K=HBD-E42K,
                    healthy_T504S=HBD-T504S,
                    E42K_T504S=E42K-T504S,
                    E42K_E42K_unaffected=E42K-E42K_unaffected,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/limmafit.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# healthy_D135Y
healthy_D135Y <- topTable(fit, coef = "healthy_D135Y", number=Inf, adjust.method = "bonferroni")
sum(healthy_D135Y$adj.P.Val < 0.05)
healthy_D135Y$gene <- sub(".*_", "", rownames(healthy_D135Y))
write.csv(healthy_D135Y, file="/folder/containing/cellranger/output/merged/healthy_D135Y_topTable.csv")

gseatt <- data.frame(gene=healthy_D135Y$gene, t=healthy_D135Y$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/healthy_D135Y_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/healthy_D135Y_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.healthy_D135Y -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/healthy_D135Y_gsea.rnk -scoring_scheme weighted -rpt_label curated.healthy_D135Y -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/healthy_D135Y_gsea.rnk -scoring_scheme weighted -rpt_label GO.healthy_D135Y -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/healthy_D135Y_gsea.rnk -scoring_scheme weighted -rpt_label immune.healthy_D135Y -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/healthy_D135Y_gsea.rnk -scoring_scheme weighted -rpt_label Single_cell.healthy_D135Y -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/ ")

# healthy_E42K
healthy_E42K <- topTable(fit, coef = "healthy_E42K", number=Inf, adjust.method = "bonferroni")
sum(healthy_E42K$adj.P.Val < 0.05)
healthy_E42K$gene <- sub(".*_", "", rownames(healthy_E42K))
write.csv(healthy_E42K, file="/folder/containing/cellranger/output/merged/healthy_E42K_topTable.csv")

gseatt <- data.frame(gene=healthy_E42K$gene, t=healthy_E42K$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/healthy_E42K_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/healthy_E42K_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.healthy_E42K -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/healthy_E42K_gsea.rnk -scoring_scheme weighted -rpt_label curated.healthy_E42K -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/healthy_E42K_gsea.rnk -scoring_scheme weighted -rpt_label GO.healthy_E42K -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/ ")

# healthy_T504S
healthy_T504S <- topTable(fit, coef = "healthy_T504S", number=Inf, adjust.method = "bonferroni")
sum(healthy_T504S$adj.P.Val < 0.05)
healthy_T504S$gene <- sub(".*_", "", rownames(healthy_T504S))
write.csv(healthy_T504S, file="/folder/containing/cellranger/output/merged/healthy_T504S_topTable.csv")

gseatt <- data.frame(gene=healthy_T504S$gene, t=healthy_T504S$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/healthy_T504S_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/healthy_T504S_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.healthy_T504S -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/healthy_T504S_gsea.rnk -scoring_scheme weighted -rpt_label curated.healthy_T504S -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/healthy_T504S_gsea.rnk -scoring_scheme weighted -rpt_label GO.healthy_T504S -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/ ")

# E42K_T504S
E42K_T504S <- topTable(fit, coef = "E42K_T504S", number=Inf, adjust.method = "bonferroni")
sum(E42K_T504S$adj.P.Val < 0.05)
E42K_T504S$gene <- sub(".*_", "", rownames(E42K_T504S))
write.csv(E42K_T504S, file="/folder/containing/cellranger/output/merged/E42K_T504S_topTable.csv")

gseatt <- data.frame(gene=E42K_T504S$gene, t=E42K_T504S$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/E42K_T504S_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/E42K_T504S_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.E42K_T504S -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/E42K_T504S_gsea.rnk -scoring_scheme weighted -rpt_label curated.E42K_T504S -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/E42K_T504S_gsea.rnk -scoring_scheme weighted -rpt_label GO.E42K_T504S -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/ ")

# E42K-E42K_unaffected
E42K-E42K_unaffected <- topTable(fit, coef = "E42K-E42K_unaffected", number=Inf, adjust.method = "bonferroni")
sum(E42K-E42K_unaffected$adj.P.Val < 0.05)
E42K-E42K_unaffected$gene <- sub(".*_", "", rownames(E42K-E42K_unaffected))
write.csv(E42K-E42K_unaffected, file="/folder/containing/cellranger/output/merged/E42K-E42K_unaffected_topTable.csv")

gseatt <- data.frame(gene=E42K-E42K_unaffected$gene, t=E42K-E42K_unaffected$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/E42K-E42K_unaffected_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/E42K-E42K_unaffected_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.E42K-E42K_unaffected -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/E42K-E42K_unaffected_gsea.rnk -scoring_scheme weighted -rpt_label curated.E42K-E42K_unaffected -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/E42K-E42K_unaffected_gsea.rnk -scoring_scheme weighted -rpt_label GO.E42K-E42K_unaffected -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/ ")

##### EXPERIMENT 6 - PROBAND PRE-TREATMENT VS POST-TREATMENT #####
# Create seu_proband by subsetting seu
seu_proband <- subset(seu, subset = treatment %in% c("pre_treatment", "post_treatment", "healthy_control"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_proband, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "proband", "all_cells", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "proband", "all_cells", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "proband", "all_cells", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/proband/all_cells/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/proband/all_cells")
}

expr <- read.table(here("merged", "sanity", "proband", "all_cells", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

#phen <- sub("^(.*?)-.*$", "\\1", seu_proband$treatment)
treatment <- as.character(seu_proband$treatment)
patient <- as.character(seu_proband$patient)

#limma
design <- model.matrix(~0 + treatment)
colnames(design) <- sub("treatment", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA
engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/all_cells/proband_pre_vs_post/corfit_proband_all_cells.RData")

cm <- makeContrasts(pre_post=pre_treatment-post_treatment,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/all_cells/proband_pre_vs_post/limmafit_proband_all_cells.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# pre_post
pre_post <- topTable(fit, coef = "pre_post", number=Inf, adjust.method = "bonferroni")
sum(pre_post$adj.P.Val < 0.05)
pre_post$gene <- sub(".*_", "", rownames(pre_post))
write.csv(pre_post, file="/folder/containing/cellranger/output/merged/all_cells/pre_post_topTable.csv")

gseatt <- data.frame(gene=pre_post$gene, t=pre_post$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/all_cells/proband_pre_vs_post/pre_post_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.pre_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/proband_pre_vs_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label curated.pre_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/proband_pre_vs_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label GO.pre_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/proband_pre_vs_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.pre_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/proband_pre_vs_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.pre_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/proband_pre_vs_post/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "pre_post"],
  logP = -log10(fit$p.value[, "pre_post"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "pre_post"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("GEM108 Pre-treatment vs Post-treatment - All cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "all_cells", "proband_pre_vs_post", "volcanoplot_proband_all_cells.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset CD4 cells and run limma
seu_proband_CD4 <- subset(seu_proband, subset = predicted.celltype.l2 %in% c("CD4 TCM", "CD4 TEM", "CD4 Naive", "CD4 CTL"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_proband_CD4, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "proband", "CD4", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "proband", "CD4", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "proband", "CD4", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/proband/CD4/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/proband/CD4")
}

expr <- read.table(here("merged", "sanity", "pre_post", "CD4", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

#phen <- sub("^(.*?)-.*$", "\\1", seu_proband_CD4$treatment)
treatment <- as.character(seu_proband_CD4$treatment)
patient <- as.character(seu_proband_CD4$patient)

#limma
design <- model.matrix(~0 + treatment)
colnames(design) <- sub("treatment", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA
engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/CD4/proband_pre_vs_post/corfit_proband_CD4.RData")

cm <- makeContrasts(pre_post=pre_treatment-post_treatment,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/CD4/proband_pre_vs_post/limmafit_proband_CD4.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# pre_post
pre_post <- topTable(fit, coef = "pre_post", number=Inf, adjust.method = "bonferroni")
sum(pre_post$adj.P.Val < 0.05)
pre_post$gene <- sub(".*_", "", rownames(pre_post))
write.csv(pre_post, file="/folder/containing/cellranger/output/merged/CD4/proband_pre_vs_post/pre_post_topTable.csv")

gseatt <- data.frame(gene=pre_post$gene, t=pre_post$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/CD4/proband_pre_vs_post/pre_post_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.pre_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/proband_pre_vs_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label curated.pre_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/proband_pre_vs_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label GO.pre_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/proband_pre_vs_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.pre_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/proband_pre_vs_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.pre_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/proband_pre_vs_post/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "pre_post"],
  logP = -log10(fit$p.value[, "pre_post"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "pre_post"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("GEM108 Pre-treatment vs Post-treatment - CD4 cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "CD4", "proband_pre_vs_post", "volcanoplot_proband_CD4.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset CD8 cells and run limma
seu_proband_CD8 <- subset(seu_proband, subset = predicted.celltype.l2 %in% c("CD8 Naive", "CD8 TEM", "CD8 TCM"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_proband_CD8, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "proband", "CD8", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "proband", "CD8", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "proband", "CD8", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/proband/CD8/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/proband/CD8")
}

expr <- read.table(here("merged", "sanity", "proband", "CD8", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

#phen <- sub("^(.*?)-.*$", "\\1", seu_proband_CD8$treatment)
treatment <- as.character(seu_proband_CD8$treatment)
patient <- as.character(seu_proband_CD8$patient)

#limma
design <- model.matrix(~0 + treatment)
colnames(design) <- sub("treatment", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/CD8/proband_pre_vs_post/corfit_proband_CD8.RData")

cm <- makeContrasts(pre_post=pre_treatment-post_treatment,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/CD8/proband_pre_vs_post/limmafit_proband_CD8.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# pre_post
pre_post <- topTable(fit, coef = "pre_post", number=Inf, adjust.method = "bonferroni")
sum(pre_post$adj.P.Val < 0.05)
pre_post$gene <- sub(".*_", "", rownames(pre_post))
write.csv(pre_post, file="/folder/containing/cellranger/output/merged/CD8/proband_pre_vs_post/pre_post_topTable.csv")

gseatt <- data.frame(gene=pre_post$gene, t=pre_post$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/CD8/proband_pre_vs_post/pre_post_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.pre_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/proband_pre_vs_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label curated.pre_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/proband_pre_vs_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label GO.pre_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/proband_pre_vs_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.pre_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/proband_pre_vs_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.pre_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/proband_pre_vs_post/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "pre_post"],
  logP = -log10(fit$p.value[, "pre_post"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "pre_post"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("GEM108 Pre-treatment vs Post-treatment - CD8 cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "CD8", "proband_pre_vs_post", "volcanoplot_proband_CD8.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset B cells cells and run limma
seu_proband_Bcell <- subset(seu_proband, subset = predicted.celltype.l2 %in% c("B naive", "B intermediate", "B memory"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_proband_Bcell, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "proband", "Bcell", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "proband", "Bcell", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "proband", "Bcell", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/proband/Bcell/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/proband/Bcell")
}

expr <- read.table(here("merged", "sanity", "proband", "Bcell", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

#phen <- sub("^(.*?)-.*$", "\\1", seu_proband_Bcell$treatment)
treatment <- as.character(seu_proband_Bcell$treatment)
patient <- as.character(seu_proband_Bcell$patient)

#limma
design <- model.matrix(~0 + treatment)
colnames(design) <- sub("treatment", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/Bcell/proband_pre_vs_post/corfit_proband_Bcell.RData")

cm <- makeContrasts(pre_post=pre_treatment-post_treatment,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/Bcell/proband_pre_vs_post/limmafit_proband_Bcell.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# pre_post
pre_post <- topTable(fit, coef = "pre_post", number=Inf, adjust.method = "bonferroni")
sum(pre_post$adj.P.Val < 0.05)
pre_post$gene <- sub(".*_", "", rownames(pre_post))
write.csv(pre_post, file="/folder/containing/cellranger/output/merged/Bcell/proband_pre_vs_post/pre_post_topTable.csv")

gseatt <- data.frame(gene=pre_post$gene, t=pre_post$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/Bcell/proband_pre_vs_post/pre_post_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.pre_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/proband_pre_vs_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label curated.pre_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/proband_pre_vs_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label GO.pre_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/proband_pre_vs_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.pre_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/proband_pre_vs_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.pre_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/proband_pre_vs_post/ ")

volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "pre_post"],
  logP = -log10(fit$p.value[, "pre_post"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "pre_post"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("GEM108 Pre-treatment vs Post-treatment - B cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "Bcell", "proband_pre_vs_post", "volcanoplot_proband_Bcell.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset Treg cells and run limma
seu_proband_Treg <- subset(seu_proband, subset = predicted.celltype.l2 %in% c("Treg"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_proband_Treg, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "proband", "Treg", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "proband", "Treg", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "proband", "Treg", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/proband/Treg/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/proband/Treg")
}

expr <- read.table(here("merged", "sanity", "proband", "Treg", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

#phen <- sub("^(.*?)-.*$", "\\1", seu_proband_Treg$treatment)
treatment <- as.character(seu_proband_Treg$treatment)
patient <- as.character(seu_proband_Treg$patient)

#limma
design <- model.matrix(~0 + treatment)
colnames(design) <- sub("treatment", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/Treg/proband_pre_vs_post/corfit_proband_Treg.RData")

cm <- makeContrasts(pre_post=pre_treatment-post_treatment,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/Treg/proband_pre_vs_post/limmafit_proband_Treg.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# pre_post
pre_post <- topTable(fit, coef = "pre_post", number=Inf, adjust.method = "bonferroni")
sum(pre_post$adj.P.Val < 0.05)
pre_post$gene <- sub(".*_", "", rownames(pre_post))
write.csv(pre_post, file="/folder/containing/cellranger/output/merged/Treg/proband_pre_vs_post/pre_post_topTable.csv")

gseatt <- data.frame(gene=pre_post$gene, t=pre_post$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/Treg/proband_pre_vs_post/pre_post_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.pre_post_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/proband_pre_vs_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label curated.pre_post_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/proband_pre_vs_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label GO.pre_post_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/proband_pre_vs_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.pre_post_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/proband_pre_vs_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/proband_pre_vs_post/pre_post_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.pre_post_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/proband_pre_vs_post/ ")

volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "pre_post"],
  logP = -log10(fit$p.value[, "pre_post"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "pre_post"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("GEM108 Pre-treatment vs Post-treatment - Treg cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "Treg", "proband_pre_vs_post", "volcanoplot_proband_Treg.pdf"), plot = volcano_plot, width = 16, height = 12)

##### EXPERIMENT 1 - HEALTHY CONTROLS VS D135Y #####
# Create seu_D135Y by subsetting seu
seu_D135Y <- subset(seu, subset = cohort %in% c("D135Y", "HBD"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_D135Y, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "healthy_D135Y", "all_cells", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "healthy_D135Y", "all_cells", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "healthy_D135Y", "all_cells", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/healthy_D135Y/all_cells/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/healthy_D135Y/all_cells")
}

expr <- read.table(here("merged", "sanity", "healthy_D135Y", "all_cells", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

#phen <- sub("^(.*?)-.*$", "\\1", seu_D135Y$cohort)
cohort <- as.character(seu_D135Y$cohort)
patient <- as.character(seu_D135Y$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/all_cells/healthy_D135Y/corfit_healthy_D135Y_all_cells.RData")

cm <- makeContrasts(healthy_D135Y=HBD-D135Y,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/all_cells/healthy_D135Y/limmafit_healthy_D135Y_all_cells.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

healthy_D135Y <- topTable(fit, coef = "healthy_D135Y", number=Inf, adjust.method = "bonferroni")
sum(healthy_D135Y$adj.P.Val < 0.05)
healthy_D135Y$gene <- sub(".*_", "", rownames(healthy_D135Y))
write.csv(healthy_D135Y, file="/folder/containing/cellranger/output/merged/all_cells/healthy_D135Y/topTable_healthy_D135Y_all_cells.csv")

gseatt <- data.frame(gene=healthy_D135Y$gene, t=healthy_D135Y$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/all_cells/healthy_D135Y/healthy_D135Y_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/healthy_D135Y/healthy_D135Y_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.healthy_D135Y -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/healthy_D135Y/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/healthy_D135Y/healthy_D135Y_gsea.rnk -scoring_scheme weighted -rpt_label curated.healthy_D135Y -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/healthy_D135Y/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/healthy_D135Y/healthy_D135Y_gsea.rnk -scoring_scheme weighted -rpt_label GO.healthy_D135Y -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/healthy_D135Y/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/healthy_D135Y/healthy_D135Y_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.healthy_D135Y -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/healthy_D135Y/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/healthy_D135Y/healthy_D135Y_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.healthy_D135Y -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/healthy_D135Y/ ")

volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "healthy_D135Y"],
  logP = -log10(fit$p.value[, "healthy_D135Y"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "healthy_D135Y"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("Healthy controls vs D135Y - All cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "all_cells", "healthy_D135Y", "volcanoplot_healthy_D135Y_all_cells.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset CD4 cells and run limma
seu_D135Y_CD4 <- subset(seu_D135Y, subset = predicted.celltype.l2 %in% c("CD4 TCM", "CD4 TEM", "CD4 Naive", "CD4 CTL"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_D135Y_CD4, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "healthy_D135Y", "CD4", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "healthy_D135Y", "CD4", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "healthy_D135Y", "CD4", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/healthy_D135Y/CD4/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/healthy_D135Y/CD4")
}

expr <- read.table(here("merged", "sanity", "healthy_D135Y", "CD4", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_D135Y_CD4$cohort)
patient <- as.character(seu_D135Y_CD4$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/CD4/healthy_D135Y/corfit_healthy_D135Y_CD4.RData")

cm <- makeContrasts(healthy_D135Y=HBD-D135Y,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/CD4/healthy_D135Y/limmafit_healthy_D135Y_CD4.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# healthy_D135Y
healthy_D135Y <- topTable(fit, coef = "healthy_D135Y", number=Inf, adjust.method = "bonferroni")
sum(healthy_D135Y$adj.P.Val < 0.05)
healthy_D135Y$gene <- sub(".*_", "", rownames(healthy_D135Y))
write.csv(healthy_D135Y, file="/folder/containing/cellranger/output/merged/CD4/healthy_D135Y/healthy_D135Y_CD4_topTable.csv")

gseatt <- data.frame(gene=healthy_D135Y$gene, t=healthy_D135Y$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/CD4/healthy_D135Y/healthy_D135Y_CD4_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/healthy_D135Y/healthy_D135Y_CD4_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.healthy_D135Y_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/healthy_D135Y/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/healthy_D135Y/healthy_D135Y_CD4_gsea.rnk -scoring_scheme weighted -rpt_label curated.healthy_D135Y_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/healthy_D135Y/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/healthy_D135Y/healthy_D135Y_CD4_gsea.rnk -scoring_scheme weighted -rpt_label GO.healthy_D135Y_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/healthy_D135Y/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/healthy_D135Y/healthy_D135Y_CD4_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.healthy_D135Y_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/healthy_D135Y/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/healthy_D135Y/healthy_D135Y_CD4_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.healthy_D135Y_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/healthy_D135Y/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "healthy_D135Y"],
  logP = -log10(fit$p.value[, "healthy_D135Y"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "healthy_D135Y"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("Healthy controls vs D135Y - CD4 cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "CD4", "healthy_D135Y", "volcanoplot_healthy_D135Y_CD4.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset CD8 cells and run limma
seu_D135Y_CD8 <- subset(seu_D135Y, subset = predicted.celltype.l2 %in% c("CD8 Naive", "CD8 TEM", "CD8 TCM"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_D135Y_CD8, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "healthy_D135Y", "CD8", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "healthy_D135Y", "CD8", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "healthy_D135Y", "CD8", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/healthy_D135Y/CD8/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/healthy_D135Y/CD8")
}

expr <- read.table(here("merged", "sanity", "healthy_D135Y", "CD8", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_D135Y_CD8$cohort)
patient <- as.character(seu_D135Y_CD8$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/CD8/healthy_D135Y/corfit_healthy_D135Y_CD8.RData")

cm <- makeContrasts(healthy_D135Y=HBD-D135Y,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/CD8/healthy_D135Y/limmafit_healthy_D135Y_CD8.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# healthy_D135Y
healthy_D135Y <- topTable(fit, coef = "healthy_D135Y", number=Inf, adjust.method = "bonferroni")
sum(healthy_D135Y$adj.P.Val < 0.05)
healthy_D135Y$gene <- sub(".*_", "", rownames(healthy_D135Y))
write.csv(healthy_D135Y, file="/folder/containing/cellranger/output/merged/CD8/healthy_D135Y/healthy_D135Y_CD8_topTable.csv")

gseatt <- data.frame(gene=healthy_D135Y$gene, t=healthy_D135Y$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/CD8/healthy_D135Y/healthy_D135Y_CD8_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/healthy_D135Y/healthy_D135Y_CD8_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.healthy_D135Y_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/healthy_D135Y/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/healthy_D135Y/healthy_D135Y_CD8_gsea.rnk -scoring_scheme weighted -rpt_label curated.healthy_D135Y_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/healthy_D135Y/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/healthy_D135Y/healthy_D135Y_CD8_gsea.rnk -scoring_scheme weighted -rpt_label GO.healthy_D135Y_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/healthy_D135Y/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/healthy_D135Y/healthy_D135Y_CD8_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.healthy_D135Y_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/healthy_D135Y/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/healthy_D135Y/healthy_D135Y_CD8_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.healthy_D135Y_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/healthy_D135Y/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "healthy_D135Y"],
  logP = -log10(fit$p.value[, "healthy_D135Y"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "healthy_D135Y"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("Healthy controls vs D135Y - CD8 cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "CD8", "healthy_D135Y", "volcanoplot_healthy_D135Y_CD8.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset B cells and run limma
seu_D135Y_Bcell <- subset(seu_D135Y, subset = predicted.celltype.l2 %in% c("B naive", "B intermediate", "B memory"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_D135Y_Bcell, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "healthy_D135Y", "Bcell", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "healthy_D135Y", "Bcell", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "healthy_D135Y", "Bcell", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/healthy_D135Y/Bcell/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/healthy_D135Y/Bcell")
}

expr <- read.table(here("merged", "sanity", "healthy_D135Y", "Bcell", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_D135Y_Bcell$cohort)
patient <- as.character(seu_D135Y_Bcell$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/Bcell/healthy_D135Y/corfit_healthy_D135Y_Bcell.RData")

cm <- makeContrasts(healthy_D135Y=HBD-D135Y,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/Bcell/healthy_D135Y/limmafit_healthy_D135Y_Bcell.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# healthy_D135Y
healthy_D135Y <- topTable(fit, coef = "healthy_D135Y", number=Inf, adjust.method = "bonferroni")
sum(healthy_D135Y$adj.P.Val < 0.05)
healthy_D135Y$gene <- sub(".*_", "", rownames(healthy_D135Y))
write.csv(healthy_D135Y, file="/folder/containing/cellranger/output/merged/Bcell/healthy_D135Y/healthy_D135Y_Bcell_topTable.csv")

gseatt <- data.frame(gene=healthy_D135Y$gene, t=healthy_D135Y$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/Bcell/healthy_D135Y/healthy_D135Y_Bcell_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/healthy_D135Y/healthy_D135Y_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.healthy_D135Y_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/healthy_D135Y/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/healthy_D135Y/healthy_D135Y_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label curated.healthy_D135Y_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/healthy_D135Y/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/healthy_D135Y/healthy_D135Y_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label GO.healthy_D135Y_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/healthy_D135Y/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/healthy_D135Y/healthy_D135Y_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.healthy_D135Y_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/healthy_D135Y/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/healthy_D135Y/healthy_D135Y_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.healthy_D135Y_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/healthy_D135Y/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "healthy_D135Y"],
  logP = -log10(fit$p.value[, "healthy_D135Y"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "healthy_D135Y"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("Healthy controls vs D135Y - B cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "Bcell", "healthy_D135Y", "volcanoplot_healthy_D135Y_Bcell.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset Treg cells and run limma
seu_D135Y_Treg <- subset(seu_D135Y, subset = predicted.celltype.l2 %in% c("Treg"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_D135Y_Treg, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "healthy_D135Y", "Treg", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "healthy_D135Y", "Treg", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "healthy_D135Y", "Treg", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/healthy_D135Y/Treg/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/healthy_D135Y/Treg")
}

expr <- read.table(here("merged", "sanity", "healthy_D135Y", "Treg", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_D135Y_Treg$cohort)
patient <- as.character(seu_D135Y_Treg$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/Treg/healthy_D135Y/corfit_healthy_D135Y_Treg.RData")

cm <- makeContrasts(healthy_D135Y=HBD-D135Y,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/Treg/healthy_D135Y/limmafit_healthy_D135Y_Treg.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# healthy_D135Y
healthy_D135Y <- topTable(fit, coef = "healthy_D135Y", number=Inf, adjust.method = "bonferroni")
sum(healthy_D135Y$adj.P.Val < 0.05)
healthy_D135Y$gene <- sub(".*_", "", rownames(healthy_D135Y))
write.csv(healthy_D135Y, file="/folder/containing/cellranger/output/merged/Treg/healthy_D135Y/healthy_D135Y_Treg_topTable.csv")

gseatt <- data.frame(gene=healthy_D135Y$gene, t=healthy_D135Y$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/Treg/healthy_D135Y/healthy_D135Y_Treg_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/healthy_D135Y/healthy_D135Y_Treg_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.healthy_D135Y_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/healthy_D135Y/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/healthy_D135Y/healthy_D135Y_Treg_gsea.rnk -scoring_scheme weighted -rpt_label curated.healthy_D135Y_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/healthy_D135Y/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/healthy_D135Y/healthy_D135Y_Treg_gsea.rnk -scoring_scheme weighted -rpt_label GO.healthy_D135Y_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/healthy_D135Y/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/healthy_D135Y/healthy_D135Y_Treg_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.healthy_D135Y_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/healthy_D135Y/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/healthy_D135Y/healthy_D135Y_Treg_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.healthy_D135Y_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/healthy_D135Y/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "healthy_D135Y"],
  logP = -log10(fit$p.value[, "healthy_D135Y"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "healthy_D135Y"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("Healthy controls vs D135Y - Treg cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "Treg", "healthy_D135Y", "volcanoplot_healthy_D135Y_Treg.pdf"), plot = volcano_plot, width = 16, height = 12)

##### EXPERIMENT 2 - HEALTHY CONTROLS VS E42K #####
# Create seu_D135Y by subsetting seu
seu_E42K <- subset(seu, subset = cohort %in% c("E42K", "HBD"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_E42K, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "healthy_E42K", "all_cells", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "healthy_E42K", "all_cells", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "healthy_E42K", "all_cells", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/healthy_E42K/all_cells/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/healthy_E42K/all_cells")
}

expr <- read.table(here("merged", "sanity", "healthy_E42K", "all_cells", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_E42K$cohort)
patient <- as.character(seu_E42K$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/all_cells/healthy_E42K/corfit_healthy_E42K_all_cells.RData")

cm <- makeContrasts(healthy_E42K=HBD-E42K,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/all_cells/healthy_E42K/limmafit_healthy_E42K_all_cells.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

healthy_E42K <- topTable(fit, coef = "healthy_E42K", number=Inf, adjust.method = "bonferroni")
sum(healthy_E42K$adj.P.Val < 0.05)
healthy_E42K$gene <- sub(".*_", "", rownames(healthy_E42K))
write.csv(healthy_E42K, file="/folder/containing/cellranger/output/merged/all_cells/healthy_E42K/topTable_healthy_E42K_all_cells.csv")

gseatt <- data.frame(gene=healthy_E42K$gene, t=healthy_E42K$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/all_cells/healthy_E42K/healthy_E42K_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/healthy_E42K/healthy_E42K_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.healthy_E42K -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/healthy_E42K/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/healthy_E42K/healthy_E42K_gsea.rnk -scoring_scheme weighted -rpt_label curated.healthy_E42K -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/healthy_E42K/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/healthy_E42K/healthy_E42K_gsea.rnk -scoring_scheme weighted -rpt_label GO.healthy_E42K -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/healthy_E42K/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/healthy_E42K/healthy_E42K_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.healthy_E42K -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/healthy_E42K/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/healthy_E42K/healthy_E42K_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.healthy_E42K -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/healthy_E42K/ ")

volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "healthy_E42K"],
  logP = -log10(fit$p.value[, "healthy_E42K"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "healthy_E42K"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("Healthy controls vs E42K - All cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "all_cells", "healthy_E42K", "volcanoplot_healthy_E42K_all_cells.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset CD4 cells and run limma
seu_E42K_CD4 <- subset(seu_E42K, subset = predicted.celltype.l2 %in% c("CD4 TCM", "CD4 TEM", "CD4 Naive", "CD4 CTL"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_E42K_CD4, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "healthy_E42K", "CD4", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "healthy_E42K", "CD4", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "healthy_E42K", "CD4", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/healthy_E42K/CD4/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/healthy_E42K/CD4")
}

expr <- read.table(here("merged", "sanity", "healthy_E42K", "CD4", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_E42K_CD4$cohort)
patient <- as.character(seu_E42K_CD4$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/CD4/healthy_E42K/corfit_healthy_E42K_CD4.RData")

cm <- makeContrasts(healthy_E42K=HBD-E42K,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/CD4/healthy_E42K/limmafit_healthy_E42K_CD4.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# healthy_E42K
healthy_E42K <- topTable(fit, coef = "healthy_E42K", number=Inf, adjust.method = "bonferroni")
sum(healthy_E42K$adj.P.Val < 0.05)
healthy_E42K$gene <- sub(".*_", "", rownames(healthy_E42K))
write.csv(healthy_E42K, file="/folder/containing/cellranger/output/merged/CD4/healthy_E42K/healthy_E42K_CD4_topTable.csv")

gseatt <- data.frame(gene=healthy_E42K$gene, t=healthy_E42K$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/CD4/healthy_E42K/healthy_E42K_CD4_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/healthy_E42K/healthy_E42K_CD4_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.healthy_E42K_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/healthy_E42K/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/healthy_E42K/healthy_E42K_CD4_gsea.rnk -scoring_scheme weighted -rpt_label curated.healthy_E42K_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/healthy_E42K/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/healthy_E42K/healthy_E42K_CD4_gsea.rnk -scoring_scheme weighted -rpt_label GO.healthy_E42K_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/healthy_E42K/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/healthy_E42K/healthy_E42K_CD4_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.healthy_E42K_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/healthy_E42K/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/healthy_E42K/healthy_E42K_CD4_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.healthy_E42K_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/healthy_E42K/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "healthy_E42K"],
  logP = -log10(fit$p.value[, "healthy_E42K"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "healthy_E42K"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("Healthy controls vs E42K - CD4 cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "CD4", "healthy_E42K", "volcanoplot_healthy_E42K_CD4.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset CD8 cells and run limma
seu_E42K_CD8 <- subset(seu_E42K, subset = predicted.celltype.l2 %in% c("CD8 Naive", "CD8 TEM", "CD8 TCM"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_E42K_CD8, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "healthy_E42K", "CD8", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "healthy_E42K", "CD8", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "healthy_E42K", "CD8", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/healthy_E42K/CD8/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/healthy_E42K/CD8")
}

expr <- read.table(here("merged", "sanity", "healthy_E42K", "CD8", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_E42K_CD8$cohort)
patient <- as.character(seu_E42K_CD8$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/CD8/healthy_E42K/corfit_healthy_E42K_CD8.RData")

cm <- makeContrasts(healthy_E42K=HBD-E42K,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/CD8/healthy_E42K/limmafit_healthy_E42K_CD8.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# healthy_E42K
healthy_E42K <- topTable(fit, coef = "healthy_E42K", number=Inf, adjust.method = "bonferroni")
sum(healthy_E42K$adj.P.Val < 0.05)
healthy_E42K$gene <- sub(".*_", "", rownames(healthy_E42K))
write.csv(healthy_E42K, file="/folder/containing/cellranger/output/merged/CD8/healthy_E42K/healthy_E42K_CD8_topTable.csv")

gseatt <- data.frame(gene=healthy_E42K$gene, t=healthy_E42K$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/CD8/healthy_E42K/healthy_E42K_CD8_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/healthy_E42K/healthy_E42K_CD8_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.healthy_E42K_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/healthy_E42K/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/healthy_E42K/healthy_E42K_CD8_gsea.rnk -scoring_scheme weighted -rpt_label curated.healthy_E42K_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/healthy_E42K/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/healthy_E42K/healthy_E42K_CD8_gsea.rnk -scoring_scheme weighted -rpt_label GO.healthy_E42K_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/healthy_E42K/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/healthy_E42K/healthy_E42K_CD8_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.healthy_E42K_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/healthy_E42K/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/healthy_E42K/healthy_E42K_CD8_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.healthy_E42K_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/healthy_E42K/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "healthy_E42K"],
  logP = -log10(fit$p.value[, "healthy_E42K"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "healthy_E42K"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("Healthy controls vs E42K - CD8 cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "CD8", "healthy_E42K", "volcanoplot_healthy_E42K_CD8.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset B cells cells and run limma
seu_E42K_Bcell <- subset(seu_E42K, subset = predicted.celltype.l2 %in% c("B naive", "B intermediate", "B memory"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_E42K_Bcell, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "healthy_E42K", "Bcell", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "healthy_E42K", "Bcell", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "healthy_E42K", "Bcell", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/healthy_E42K/Bcell/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/healthy_E42K/Bcell")
}

expr <- read.table(here("merged", "sanity", "healthy_E42K", "Bcell", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_E42K_Bcell$cohort)
patient <- as.character(seu_E42K_Bcell$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/Bcell/healthy_E42K/corfit_healthy_E42K_Bcell.RData")

cm <- makeContrasts(healthy_E42K=HBD-E42K,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/Bcell/healthy_E42K/limmafit_healthy_E42K_Bcell.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# healthy_E42K
healthy_E42K <- topTable(fit, coef = "healthy_E42K", number=Inf, adjust.method = "bonferroni")
sum(healthy_E42K$adj.P.Val < 0.05)
healthy_E42K$gene <- sub(".*_", "", rownames(healthy_E42K))
write.csv(healthy_E42K, file="/folder/containing/cellranger/output/merged/Bcell/healthy_E42K/healthy_E42K_Bcell_topTable.csv")

gseatt <- data.frame(gene=healthy_E42K$gene, t=healthy_E42K$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/Bcell/healthy_E42K/healthy_E42K_Bcell_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/healthy_E42K/healthy_E42K_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.healthy_E42K_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/healthy_E42K/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/healthy_E42K/healthy_E42K_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label curated.healthy_E42K_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/healthy_E42K/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/healthy_E42K/healthy_E42K_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label GO.healthy_E42K_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/healthy_E42K/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/healthy_E42K/healthy_E42K_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.healthy_E42K_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/healthy_E42K/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/healthy_E42K/healthy_E42K_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.healthy_E42K_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/healthy_E42K/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "healthy_E42K"],
  logP = -log10(fit$p.value[, "healthy_E42K"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "healthy_E42K"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("Healthy controls vs E42K - B cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "Bcell", "healthy_E42K", "volcanoplot_healthy_E42K_Bcell.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset Treg cells and run limma
seu_E42K_Treg <- subset(seu_E42K, subset = predicted.celltype.l2 %in% c("Treg"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_E42K_Treg, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "healthy_E42K", "Treg", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "healthy_E42K", "Treg", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "healthy_E42K", "Treg", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/healthy_E42K/Treg/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/healthy_E42K/Treg")
}

expr <- read.table(here("merged", "sanity", "healthy_E42K", "Treg", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_E42K_Treg$cohort)
patient <- as.character(seu_E42K_Treg$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/Treg/healthy_E42K/corfit_healthy_E42K_Treg.RData")

cm <- makeContrasts(healthy_E42K=HBD-E42K,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/Treg/healthy_E42K/limmafit_healthy_E42K_Treg.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# healthy_E42K
healthy_E42K <- topTable(fit, coef = "healthy_E42K", number=Inf, adjust.method = "bonferroni")
sum(healthy_E42K$adj.P.Val < 0.05)
healthy_E42K$gene <- sub(".*_", "", rownames(healthy_E42K))
write.csv(healthy_E42K, file="/folder/containing/cellranger/output/merged/Treg/healthy_E42K/healthy_E42K_Treg_topTable.csv")

gseatt <- data.frame(gene=healthy_E42K$gene, t=healthy_E42K$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/Treg/healthy_E42K/healthy_E42K_Treg_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/healthy_E42K/healthy_E42K_Treg_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.healthy_E42K_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/healthy_E42K/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/healthy_E42K/healthy_E42K_Treg_gsea.rnk -scoring_scheme weighted -rpt_label curated.healthy_E42K_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/healthy_E42K/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/healthy_E42K/healthy_E42K_Treg_gsea.rnk -scoring_scheme weighted -rpt_label GO.healthy_E42K_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/healthy_E42K/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/healthy_E42K/healthy_E42K_Treg_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.healthy_E42K_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/healthy_E42K/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/healthy_E42K/healthy_E42K_Treg_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.healthy_E42K_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/healthy_E42K/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "healthy_E42K"],
  logP = -log10(fit$p.value[, "healthy_E42K"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "healthy_E42K"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("Healthy controls vs E42K - Treg cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "Treg", "healthy_E42K", "volcanoplot_healthy_E42K_Treg.pdf"), plot = volcano_plot, width = 16, height = 12)

##### EXPERIMENT 3 - HEALTHY CONTROLS VS T504S #####
# This subsets out the T504S cohort because limma is just taking too long to do the whole lot in one go
# Create seu_T504S by subsetting seu
seu_T504S <- subset(seu, subset = cohort %in% c("T504S", "HBD"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_T504S, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "healthy_T504S", "all_cells", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "healthy_T504S", "all_cells", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "healthy_T504S", "all_cells", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/healthy_T504S/all_cells/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/healthy_T504S/all_cells")
}

expr <- read.table(here("merged", "sanity", "healthy_T504S", "all_cells", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_T504S$cohort)
patient <- as.character(seu_T504S$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/all_cells/healthy_T504S/corfit_healthy_T504S_all_cells.RData")

cm <- makeContrasts(healthy_T504S=HBD-T504S,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/all_cells/healthy_T504S/limmafit_healthy_T504S_all_cells.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

healthy_T504S <- topTable(fit, coef = "healthy_T504S", number=Inf, adjust.method = "bonferroni")
sum(healthy_T504S$adj.P.Val < 0.05)
healthy_T504S$gene <- sub(".*_", "", rownames(healthy_T504S))
write.csv(healthy_T504S, file="/folder/containing/cellranger/output/merged/all_cells/healthy_T504S/topTable_healthy_T504S_all_cells.csv")

gseatt <- data.frame(gene=healthy_T504S$gene, t=healthy_T504S$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/all_cells/healthy_T504S/healthy_T504S_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/healthy_T504S/healthy_T504S_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.healthy_T504S -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/healthy_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/healthy_T504S/healthy_T504S_gsea.rnk -scoring_scheme weighted -rpt_label curated.healthy_T504S -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/healthy_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/healthy_T504S/healthy_T504S_gsea.rnk -scoring_scheme weighted -rpt_label GO.healthy_T504S -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/healthy_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/healthy_T504S/healthy_T504S_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.healthy_T504S -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/healthy_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/healthy_T504S/healthy_T504S_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.healthy_T504S -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/healthy_T504S/ ")

volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "healthy_T504S"],
  logP = -log10(fit$p.value[, "healthy_T504S"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "healthy_T504S"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("Healthy controls vs T504S - All cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "all_cells", "healthy_T504S", "volcanoplot_healthy_T504S_all_cells.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset CD4 cells and run limma
seu_T504S_CD4 <- subset(seu_T504S, subset = predicted.celltype.l2 %in% c("CD4 TCM", "CD4 TEM", "CD4 Naive", "CD4 CTL"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_T504S_CD4, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "healthy_T504S", "CD4", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "healthy_T504S", "CD4", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "healthy_T504S", "CD4", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/healthy_T504S/CD4/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/healthy_T504S/CD4")
}

expr <- read.table(here("merged", "sanity", "healthy_T504S", "CD4", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_T504S_CD4$cohort)
patient <- as.character(seu_T504S_CD4$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/CD4/healthy_T504S/corfit_healthy_T504S_CD4.RData")

cm <- makeContrasts(healthy_T504S=HBD-T504S,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/CD4/healthy_T504S/limmafit_healthy_T504S_CD4.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# healthy_T504S
healthy_T504S <- topTable(fit, coef = "healthy_T504S", number=Inf, adjust.method = "bonferroni")
sum(healthy_T504S$adj.P.Val < 0.05)
healthy_T504S$gene <- sub(".*_", "", rownames(healthy_T504S))
write.csv(healthy_T504S, file="/folder/containing/cellranger/output/merged/CD4/healthy_T504S/healthy_T504S_CD4_topTable.csv")

gseatt <- data.frame(gene=healthy_T504S$gene, t=healthy_T504S$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/CD4/healthy_T504S/healthy_T504S_CD4_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/healthy_T504S/healthy_T504S_CD4_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.healthy_T504S_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/healthy_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/healthy_T504S/healthy_T504S_CD4_gsea.rnk -scoring_scheme weighted -rpt_label curated.healthy_T504S_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/healthy_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/healthy_T504S/healthy_T504S_CD4_gsea.rnk -scoring_scheme weighted -rpt_label GO.healthy_T504S_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/healthy_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/healthy_T504S/healthy_T504S_CD4_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.healthy_T504S_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/healthy_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/healthy_T504S/healthy_T504S_CD4_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.healthy_T504S_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/healthy_T504S/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "healthy_T504S"],
  logP = -log10(fit$p.value[, "healthy_T504S"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "healthy_T504S"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("Healthy controls vs T504S - CD4 cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "CD4", "healthy_T504S", "volcanoplot_healthy_T504S_CD4.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset CD8 cells and run limma
seu_T504S_CD8 <- subset(seu_T504S, subset = predicted.celltype.l2 %in% c("CD8 Naive", "CD8 TEM", "CD8 TCM"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_T504S_CD8, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "healthy_T504S", "CD8", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "healthy_T504S", "CD8", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "healthy_T504S", "CD8", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/healthy_T504S/CD8/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/healthy_T504S/CD8")
}

expr <- read.table(here("merged", "sanity", "healthy_T504S", "CD8", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_T504S_CD8$cohort)
patient <- as.character(seu_T504S_CD8$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/CD8/healthy_T504S/corfit_healthy_T504S_CD8.RData")

cm <- makeContrasts(healthy_T504S=HBD-T504S,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/CD8/healthy_T504S/limmafit_healthy_T504S_CD8.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# healthy_T504S
healthy_T504S <- topTable(fit, coef = "healthy_T504S", number=Inf, adjust.method = "bonferroni")
sum(healthy_T504S$adj.P.Val < 0.05)
healthy_T504S$gene <- sub(".*_", "", rownames(healthy_T504S))
write.csv(healthy_T504S, file="/folder/containing/cellranger/output/merged/CD8/healthy_T504S/healthy_T504S_CD8_topTable.csv")

gseatt <- data.frame(gene=healthy_T504S$gene, t=healthy_T504S$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/CD8/healthy_T504S/healthy_T504S_CD8_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/healthy_T504S/healthy_T504S_CD8_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.healthy_T504S_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/healthy_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/healthy_T504S/healthy_T504S_CD8_gsea.rnk -scoring_scheme weighted -rpt_label curated.healthy_T504S_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/healthy_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/healthy_T504S/healthy_T504S_CD8_gsea.rnk -scoring_scheme weighted -rpt_label GO.healthy_T504S_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/healthy_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/healthy_T504S/healthy_T504S_CD8_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.healthy_T504S_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/healthy_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/healthy_T504S/healthy_T504S_CD8_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.healthy_T504S_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/healthy_T504S/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "healthy_T504S"],
  logP = -log10(fit$p.value[, "healthy_T504S"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "healthy_T504S"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("Healthy controls vs T504S - CD8 cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "CD8", "healthy_T504S", "volcanoplot_healthy_T504S_CD8.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset B cells cells and run limma
seu_T504S_Bcell <- subset(seu_T504S, subset = predicted.celltype.l2 %in% c("B naive", "B intermediate", "B memory"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_T504S_Bcell, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "healthy_T504S", "Bcell", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "healthy_T504S", "Bcell", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "healthy_T504S", "Bcell", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/healthy_T504S/Bcell/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/healthy_T504S/Bcell")
}

expr <- read.table(here("merged", "sanity", "healthy_T504S", "Bcell", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_T504S_Bcell$cohort)
patient <- as.character(seu_T504S_Bcell$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/Bcell/healthy_T504S/corfit_healthy_T504S_Bcell.RData")

cm <- makeContrasts(healthy_T504S=HBD-T504S,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/Bcell/healthy_T504S/limmafit_healthy_T504S_Bcell.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# healthy_T504S
healthy_T504S <- topTable(fit, coef = "healthy_T504S", number=Inf, adjust.method = "bonferroni")
sum(healthy_T504S$adj.P.Val < 0.05)
healthy_T504S$gene <- sub(".*_", "", rownames(healthy_T504S))
write.csv(healthy_T504S, file="/folder/containing/cellranger/output/merged/Bcell/healthy_T504S/healthy_T504S_Bcell_topTable.csv")

gseatt <- data.frame(gene=healthy_T504S$gene, t=healthy_T504S$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/Bcell/healthy_T504S/healthy_T504S_Bcell_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/healthy_T504S/healthy_T504S_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.healthy_T504S_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/healthy_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/healthy_T504S/healthy_T504S_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label curated.healthy_T504S_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/healthy_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/healthy_T504S/healthy_T504S_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label GO.healthy_T504S_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/healthy_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/healthy_T504S/healthy_T504S_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.healthy_T504S_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/healthy_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/healthy_T504S/healthy_T504S_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.healthy_T504S_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/healthy_T504S/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "healthy_T504S"],
  logP = -log10(fit$p.value[, "healthy_T504S"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "healthy_T504S"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("Healthy controls vs T504S - B cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "Bcell", "healthy_T504S", "volcanoplot_healthy_T504S_Bcell.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset Treg cells and run limma
seu_T504S_Treg <- subset(seu_T504S, subset = predicted.celltype.l2 %in% c("Treg"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_T504S_Treg, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "healthy_T504S", "Treg", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "healthy_T504S", "Treg", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "healthy_T504S", "Treg", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/healthy_T504S/Treg/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/healthy_T504S/Treg")
}

expr <- read.table(here("merged", "sanity", "healthy_T504S", "Treg", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_T504S_Treg$cohort)
patient <- as.character(seu_T504S_Treg$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/Treg/healthy_T504S/corfit_healthy_T504S_Treg.RData")

cm <- makeContrasts(healthy_T504S=HBD-T504S,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/Treg/healthy_T504S/limmafit_healthy_T504S_Treg.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# healthy_T504S
healthy_T504S <- topTable(fit, coef = "healthy_T504S", number=Inf, adjust.method = "bonferroni")
sum(healthy_T504S$adj.P.Val < 0.05)
healthy_T504S$gene <- sub(".*_", "", rownames(healthy_T504S))
write.csv(healthy_T504S, file="/folder/containing/cellranger/output/merged/Treg/healthy_T504S/healthy_T504S_Treg_topTable.csv")

gseatt <- data.frame(gene=healthy_T504S$gene, t=healthy_T504S$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/Treg/healthy_T504S/healthy_T504S_Treg_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/healthy_T504S/healthy_T504S_Treg_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.healthy_T504S_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/healthy_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/healthy_T504S/healthy_T504S_Treg_gsea.rnk -scoring_scheme weighted -rpt_label curated.healthy_T504S_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/healthy_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/healthy_T504S/healthy_T504S_Treg_gsea.rnk -scoring_scheme weighted -rpt_label GO.healthy_T504S_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/healthy_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/healthy_T504S/healthy_T504S_Treg_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.healthy_T504S_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/healthy_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/healthy_T504S/healthy_T504S_Treg_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.healthy_T504S_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/healthy_T504S/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "healthy_T504S"],
  logP = -log10(fit$p.value[, "healthy_T504S"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "healthy_T504S"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("Healthy controls vs T504S - Treg cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "Treg", "healthy_T504S", "volcanoplot_healthy_T504S_Treg.pdf"), plot = volcano_plot, width = 16, height = 12)

##### EXPERIMENT 4 - E42K VS T504S #####
# Create seu_E42K_T504S by subsetting seu
seu_E42K_T504S <- subset(seu, subset = cohort %in% c("E42K", "T504S"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_E42K_T504S, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "E42K_T504S", "all_cells", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "E42K_T504S", "all_cells", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "E42K_T504S", "all_cells", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/E42K_T504S/all_cells/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/E42K_T504S/all_cells")
}

expr <- read.table(here("merged", "sanity", "E42K_T504S", "all_cells", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_E42K_T504S$cohort)
patient <- as.character(seu_E42K_T504S$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/all_cells/E42K_T504S/corfit_E42K_T504S_all_cells.RData")

cm <- makeContrasts(E42K_T504S=E42K-T504S,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/all_cells/E42K_T504S/limmafit_E42K_T504S_all_cells.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

E42K_T504S <- topTable(fit, coef = "E42K_T504S", number=Inf, adjust.method = "bonferroni")
sum(E42K_T504S$adj.P.Val < 0.05)
E42K_T504S$gene <- sub(".*_", "", rownames(E42K_T504S))
write.csv(E42K_T504S, file="/folder/containing/cellranger/output/merged/all_cells/E42K_T504S/topTable_E42K_T504S_all_cells.csv")

gseatt <- data.frame(gene=E42K_T504S$gene, t=E42K_T504S$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/all_cells/E42K_T504S/E42K_T504S_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/E42K_T504S/E42K_T504S_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.E42K_T504S -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/E42K_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/E42K_T504S/E42K_T504S_gsea.rnk -scoring_scheme weighted -rpt_label curated.E42K_T504S -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/E42K_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/E42K_T504S/E42K_T504S_gsea.rnk -scoring_scheme weighted -rpt_label GO.E42K_T504S -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/E42K_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/E42K_T504S/E42K_T504S_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.E42K_T504S -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/E42K_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/E42K_T504S/E42K_T504S_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.E42K_T504S -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/E42K_T504S/ ")

volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "E42K_T504S"],
  logP = -log10(fit$p.value[, "E42K_T504S"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "E42K_T504S"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("E42K vs T504S - All cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "all_cells", "E42K_T504S", "volcanoplot_E42K_T504S_all_cells.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset CD4 cells and run limma
seu_E42K_T504S_CD4 <- subset(seu_E42K_T504S, subset = predicted.celltype.l2 %in% c("CD4 TCM", "CD4 TEM", "CD4 Naive", "CD4 CTL"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_E42K_T504S_CD4, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "E42K_T504S", "CD4", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "E42K_T504S", "CD4", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "E42K_T504S", "CD4", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/E42K_T504S/CD4/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/E42K_T504S/CD4")
}

expr <- read.table(here("merged", "sanity", "E42K_T504S", "CD4", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_E42K_T504S_CD4$cohort)
patient <- as.character(seu_E42K_T504S_CD4$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/CD4/E42K_T504S/corfit_E42K_T504S_CD4.RData")

cm <- makeContrasts(E42K_T504S=E42K-T504S,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/CD4/E42K_T504S/limmafit_E42K_T504S_CD4.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# E42K_T504S
E42K_T504S <- topTable(fit, coef = "E42K_T504S", number=Inf, adjust.method = "bonferroni")
sum(E42K_T504S$adj.P.Val < 0.05)
E42K_T504S$gene <- sub(".*_", "", rownames(E42K_T504S))
write.csv(E42K_T504S, file="/folder/containing/cellranger/output/merged/CD4/E42K_T504S/E42K_T504S_CD4_topTable.csv")

gseatt <- data.frame(gene=E42K_T504S$gene, t=E42K_T504S$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/CD4/E42K_T504S/E42K_T504S_CD4_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/E42K_T504S/E42K_T504S_CD4_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.E42K_T504S_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/E42K_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/E42K_T504S/E42K_T504S_CD4_gsea.rnk -scoring_scheme weighted -rpt_label curated.E42K_T504S_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/E42K_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/E42K_T504S/E42K_T504S_CD4_gsea.rnk -scoring_scheme weighted -rpt_label GO.E42K_T504S_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/E42K_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/E42K_T504S/E42K_T504S_CD4_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.E42K_T504S_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/E42K_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/E42K_T504S/E42K_T504S_CD4_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.E42K_T504S_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/E42K_T504S/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "E42K_T504S"],
  logP = -log10(fit$p.value[, "E42K_T504S"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "E42K_T504S"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("E42K vs T504S - CD4 cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "CD4", "E42K_T504S", "volcanoplot_E42K_T504S_CD4.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset CD8 cells and run limma
seu_E42K_T504S_CD8 <- subset(seu_E42K_T504S, subset = predicted.celltype.l2 %in% c("CD8 Naive", "CD8 TEM", "CD8 TCM"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_E42K_T504S_CD8, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "E42K_T504S", "CD8", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "E42K_T504S", "CD8", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "E42K_T504S", "CD8", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/E42K_T504S/CD8/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/E42K_T504S/CD8")
}

expr <- read.table(here("merged", "sanity", "E42K_T504S", "CD8", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_E42K_T504S_CD8$cohort)
patient <- as.character(seu_E42K_T504S_CD8$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/CD8/E42K_T504S/corfit_E42K_T504S_CD8.RData")

cm <- makeContrasts(E42K_T504S=E42K-T504S,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/CD8/E42K_T504S/limmafit_E42K_T504S_CD8.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# E42K_T504S
E42K_T504S <- topTable(fit, coef = "E42K_T504S", number=Inf, adjust.method = "bonferroni")
sum(E42K_T504S$adj.P.Val < 0.05)
E42K_T504S$gene <- sub(".*_", "", rownames(E42K_T504S))
write.csv(E42K_T504S, file="/folder/containing/cellranger/output/merged/CD8/E42K_T504S/E42K_T504S_CD8_topTable.csv")

gseatt <- data.frame(gene=E42K_T504S$gene, t=E42K_T504S$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/CD8/E42K_T504S/E42K_T504S_CD8_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/E42K_T504S/E42K_T504S_CD8_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.E42K_T504S_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/E42K_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/E42K_T504S/E42K_T504S_CD8_gsea.rnk -scoring_scheme weighted -rpt_label curated.E42K_T504S_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/E42K_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/E42K_T504S/E42K_T504S_CD8_gsea.rnk -scoring_scheme weighted -rpt_label GO.E42K_T504S_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/E42K_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/E42K_T504S/E42K_T504S_CD8_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.E42K_T504S_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/E42K_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/E42K_T504S/E42K_T504S_CD8_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.E42K_T504S_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/E42K_T504S/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "E42K_T504S"],
  logP = -log10(fit$p.value[, "E42K_T504S"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "E42K_T504S"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("E42K vs T504S - CD8 cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "CD8", "E42K_T504S", "volcanoplot_E42K_T504S_CD8.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset B cells cells and run limma
seu_E42K_T504S_Bcell <- subset(seu_E42K_T504S, subset = predicted.celltype.l2 %in% c("B naive", "B intermediate", "B memory"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_E42K_T504S_Bcell, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "E42K_T504S", "Bcell", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "E42K_T504S", "Bcell", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "E42K_T504S", "Bcell", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/E42K_T504S/Bcell/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/E42K_T504S/Bcell")
}

expr <- read.table(here("merged", "sanity", "E42K_T504S", "Bcell", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_E42K_T504S_Bcell$cohort)
patient <- as.character(seu_E42K_T504S_Bcell$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/Bcell/E42K_T504S/corfit_E42K_T504S_Bcell.RData")

cm <- makeContrasts(E42K_T504S=E42K-T504S,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/Bcell/E42K_T504S/limmafit_E42K_T504S_Bcell.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# E42K_T504S
E42K_T504S <- topTable(fit, coef = "E42K_T504S", number=Inf, adjust.method = "bonferroni")
sum(E42K_T504S$adj.P.Val < 0.05)
E42K_T504S$gene <- sub(".*_", "", rownames(E42K_T504S))
write.csv(E42K_T504S, file="/folder/containing/cellranger/output/merged/Bcell/E42K_T504S/E42K_T504S_Bcell_topTable.csv")

gseatt <- data.frame(gene=E42K_T504S$gene, t=E42K_T504S$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/Bcell/E42K_T504S/E42K_T504S_Bcell_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/E42K_T504S/E42K_T504S_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.E42K_T504S_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/E42K_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/E42K_T504S/E42K_T504S_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label curated.E42K_T504S_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/E42K_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/E42K_T504S/E42K_T504S_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label GO.E42K_T504S_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/E42K_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/E42K_T504S/E42K_T504S_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.E42K_T504S_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/E42K_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/E42K_T504S/E42K_T504S_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.E42K_T504S_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/E42K_T504S/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "E42K_T504S"],
  logP = -log10(fit$p.value[, "E42K_T504S"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "E42K_T504S"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("E42K vs T504S - B cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "Bcell", "E42K_T504S", "volcanoplot_E42K_T504S_Bcell.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset Treg cells and run limma
seu_E42K_T504S_Treg <- subset(seu_E42K_T504S, subset = predicted.celltype.l2 %in% c("Treg"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_E42K_T504S_Treg, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "E42K_T504S", "Treg", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "E42K_T504S", "Treg", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "E42K_T504S", "Treg", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/E42K_T504S/Treg/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/E42K_T504S/Treg")
}

expr <- read.table(here("merged", "sanity", "E42K_T504S", "Treg", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_E42K_T504S_Treg$cohort)
patient <- as.character(seu_E42K_T504S_Treg$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/Treg/E42K_T504S/corfit_E42K_T504S_Treg.RData")

cm <- makeContrasts(E42K_T504S=E42K-T504S,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/Treg/E42K_T504S/limmafit_E42K_T504S_Treg.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# E42K_T504S
E42K_T504S <- topTable(fit, coef = "E42K_T504S", number=Inf, adjust.method = "bonferroni")
sum(E42K_T504S$adj.P.Val < 0.05)
E42K_T504S$gene <- sub(".*_", "", rownames(E42K_T504S))
write.csv(E42K_T504S, file="/folder/containing/cellranger/output/merged/Treg/E42K_T504S/E42K_T504S_Treg_topTable.csv")

gseatt <- data.frame(gene=E42K_T504S$gene, t=E42K_T504S$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/Treg/E42K_T504S/E42K_T504S_Treg_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/E42K_T504S/E42K_T504S_Treg_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.E42K_T504S_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/E42K_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/E42K_T504S/E42K_T504S_Treg_gsea.rnk -scoring_scheme weighted -rpt_label curated.E42K_T504S_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/E42K_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/E42K_T504S/E42K_T504S_Treg_gsea.rnk -scoring_scheme weighted -rpt_label GO.E42K_T504S_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/E42K_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/E42K_T504S/E42K_T504S_Treg_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.E42K_T504S_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/E42K_T504S/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/E42K_T504S/E42K_T504S_Treg_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.E42K_T504S_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/E42K_T504S/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "E42K_T504S"],
  logP = -log10(fit$p.value[, "E42K_T504S"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "E42K_T504S"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("E42K vs T504S - Treg cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "Treg", "E42K_T504S", "volcanoplot_E42K_T504S_Treg.pdf"), plot = volcano_plot, width = 16, height = 12)

##### EXPERIMENT 5 - E42K AFFECTED VS E42K UNAFFECTED #####
# Create seu_E42K_E42KU by subsetting seu
seu_E42K_E42KU <- subset(seu, subset = cohort %in% c("E42K", "E42K_unaffected"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_E42K_E42KU, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "E42K_E42KU", "all_cells", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "E42K_E42KU", "all_cells", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "E42K_E42KU", "all_cells", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/E42K_E42KU/all_cells/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/E42K_E42KU/all_cells")
}

expr <- read.table(here("merged", "sanity", "E42K_E42KU", "all_cells", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_E42K_E42KU$cohort)
patient <- as.character(seu_E42K_E42KU$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/all_cells/E42K_E42KU/corfit_E42K_E42KU_all_cells.RData")

cm <- makeContrasts(E42K_E42KU=E42K-E42K_unaffected,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/all_cells/E42K_E42KU/limmafit_E42K_E42KU_all_cells.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

E42K_E42KU <- topTable(fit, coef = "E42K_E42KU", number=Inf, adjust.method = "bonferroni")
sum(E42K_E42KU$adj.P.Val < 0.05)
E42K_E42KU$gene <- sub(".*_", "", rownames(E42K_E42KU))
write.csv(E42K_E42KU, file="/folder/containing/cellranger/output/merged/all_cells/E42K_E42KU/topTable_E42K_E42KU_all_cells.csv")

gseatt <- data.frame(gene=E42K_E42KU$gene, t=E42K_E42KU$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/all_cells/E42K_E42KU/E42K_E42KU_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/E42K_E42KU/E42K_E42KU_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.E42K_E42KU -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/E42K_E42KU/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/E42K_E42KU/E42K_E42KU_gsea.rnk -scoring_scheme weighted -rpt_label curated.E42K_E42KU -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/E42K_E42KU/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/E42K_E42KU/E42K_E42KU_gsea.rnk -scoring_scheme weighted -rpt_label GO.E42K_E42KU -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/E42K_E42KU/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/E42K_E42KU/E42K_E42KU_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.E42K_E42KU -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/E42K_E42KU/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/E42K_E42KU/E42K_E42KU_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.E42K_E42KU -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/E42K_E42KU/ ")

volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "E42K_E42KU"],
  logP = -log10(fit$p.value[, "E42K_E42KU"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "E42K_E42KU"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("E42K vs E42K unaffected - All cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "all_cells", "E42K_E42KU", "volcanoplot_E42K_E42KU_all_cells.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset CD4 cells and run limma
seu_E42K_E42KU_CD4 <- subset(seu_E42K_E42KU, subset = predicted.celltype.l2 %in% c("CD4 TCM", "CD4 TEM", "CD4 Naive", "CD4 CTL"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_E42K_E42KU_CD4, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "E42K_E42KU", "CD4", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "E42K_E42KU", "CD4", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "E42K_E42KU", "CD4", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/E42K_E42KU/CD4/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/E42K_E42KU/CD4")
}

expr <- read.table(here("merged", "sanity", "E42K_E42KU", "CD4", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_E42K_E42KU_CD4$cohort)
patient <- as.character(seu_E42K_E42KU_CD4$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/CD4/E42K_E42KU/corfit_E42K_E42KU_CD4.RData")

cm <- makeContrasts(E42K_E42KU=E42K-E42K_unaffected,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/CD4/E42K_E42KU/limmafit_E42K_E42KU_CD4.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# E42K_E42KU
E42K_E42KU <- topTable(fit, coef = "E42K_E42KU", number=Inf, adjust.method = "bonferroni")
sum(E42K_E42KU$adj.P.Val < 0.05)
E42K_E42KU$gene <- sub(".*_", "", rownames(E42K_E42KU))
write.csv(E42K_E42KU, file="/folder/containing/cellranger/output/merged/CD4/E42K_E42KU/E42K_E42KU_CD4_topTable.csv")

gseatt <- data.frame(gene=E42K_E42KU$gene, t=E42K_E42KU$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/CD4/E42K_E42KU/E42K_E42KU_CD4_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/E42K_E42KU/E42K_E42KU_CD4_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.E42K_E42KU_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/E42K_E42KU/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/E42K_E42KU/E42K_E42KU_CD4_gsea.rnk -scoring_scheme weighted -rpt_label curated.E42K_E42KU_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/E42K_E42KU/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/E42K_E42KU/E42K_E42KU_CD4_gsea.rnk -scoring_scheme weighted -rpt_label GO.E42K_E42KU_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/E42K_E42KU/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/E42K_E42KU/E42K_E42KU_CD4_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.E42K_E42KU_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/E42K_E42KU/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/E42K_E42KU/E42K_E42KU_CD4_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.E42K_E42KU_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/E42K_E42KU/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "E42K_E42KU"],
  logP = -log10(fit$p.value[, "E42K_E42KU"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "E42K_E42KU"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("E42K vs E42K unaffected - CD4 cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "CD4", "E42K_E42KU", "volcanoplot_E42K_E42KU_CD4.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset CD8 cells and run limma
seu_E42K_E42KU_CD8 <- subset(seu_E42K_E42KU, subset = predicted.celltype.l2 %in% c("CD8 Naive", "CD8 TEM", "CD8 TCM"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_E42K_E42KU_CD8, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "E42K_E42KU", "CD8", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "E42K_E42KU", "CD8", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "E42K_E42KU", "CD8", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/E42K_E42KU/CD8/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/E42K_E42KU/CD8")
}

expr <- read.table(here("merged", "sanity", "E42K_E42KU", "CD8", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_E42K_E42KU_CD8$cohort)
patient <- as.character(seu_E42K_E42KU_CD8$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/CD8/E42K_E42KU/corfit_E42K_E42KU_CD8.RData")

cm <- makeContrasts(E42K_E42KU=E42K-E42K_unaffected,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/CD8/E42K_E42KU/limmafit_E42K_E42KU_CD8.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# E42K_E42KU
E42K_E42KU <- topTable(fit, coef = "E42K_E42KU", number=Inf, adjust.method = "bonferroni")
sum(E42K_E42KU$adj.P.Val < 0.05)
E42K_E42KU$gene <- sub(".*_", "", rownames(E42K_E42KU))
write.csv(E42K_E42KU, file="/folder/containing/cellranger/output/merged/CD8/E42K_E42KU/E42K_E42KU_CD8_topTable.csv")

gseatt <- data.frame(gene=E42K_E42KU$gene, t=E42K_E42KU$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/CD8/E42K_E42KU/E42K_E42KU_CD8_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/E42K_E42KU/E42K_E42KU_CD8_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.E42K_E42KU_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/E42K_E42KU/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/E42K_E42KU/E42K_E42KU_CD8_gsea.rnk -scoring_scheme weighted -rpt_label curated.E42K_E42KU_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/E42K_E42KU/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/E42K_E42KU/E42K_E42KU_CD8_gsea.rnk -scoring_scheme weighted -rpt_label GO.E42K_E42KU_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/E42K_E42KU/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/E42K_E42KU/E42K_E42KU_CD8_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.E42K_E42KU_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/E42K_E42KU/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/E42K_E42KU/E42K_E42KU_CD8_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.E42K_E42KU_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/E42K_E42KU/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "E42K_E42KU"],
  logP = -log10(fit$p.value[, "E42K_E42KU"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "E42K_E42KU"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("E42K vs E42K unaffected - CD8 cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "CD8", "E42K_E42KU", "volcanoplot_E42K_E42KU_CD8.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset B cells cells and run limma
seu_E42K_E42KU_Bcell <- subset(seu_E42K_E42KU, subset = predicted.celltype.l2 %in% c("B naive", "B intermediate", "B memory"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_E42K_E42KU_Bcell, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "E42K_E42KU", "Bcell", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "E42K_E42KU", "Bcell", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "E42K_E42KU", "Bcell", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/E42K_E42KU/Bcell/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/E42K_E42KU/Bcell")
}

expr <- read.table(here("merged", "sanity", "E42K_E42KU", "Bcell", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_E42K_E42KU_Bcell$cohort)
patient <- as.character(seu_E42K_E42KU_Bcell$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/Bcell/E42K_E42KU/corfit_E42K_E42KU_Bcell.RData")

cm <- makeContrasts(E42K_E42KU=E42K-E42K_unaffected,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/Bcell/E42K_E42KU/limmafit_E42K_E42KU_Bcell.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# E42K_E42KU
E42K_E42KU <- topTable(fit, coef = "E42K_E42KU", number=Inf, adjust.method = "bonferroni")
sum(E42K_E42KU$adj.P.Val < 0.05)
E42K_E42KU$gene <- sub(".*_", "", rownames(E42K_E42KU))
write.csv(E42K_E42KU, file="/folder/containing/cellranger/output/merged/Bcell/E42K_E42KU/E42K_E42KU_Bcell_topTable.csv")

gseatt <- data.frame(gene=E42K_E42KU$gene, t=E42K_E42KU$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/Bcell/E42K_E42KU/E42K_E42KU_Bcell_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/E42K_E42KU/E42K_E42KU_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.E42K_E42KU_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/E42K_E42KU/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/E42K_E42KU/E42K_E42KU_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label curated.E42K_E42KU_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/E42K_E42KU/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/E42K_E42KU/E42K_E42KU_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label GO.E42K_E42KU_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/E42K_E42KU/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/E42K_E42KU/E42K_E42KU_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.E42K_E42KU_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/E42K_E42KU/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/E42K_E42KU/E42K_E42KU_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.E42K_E42KU_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/E42K_E42KU/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "E42K_E42KU"],
  logP = -log10(fit$p.value[, "E42K_E42KU"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "E42K_E42KU"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("E42K vs E42K unaffected - B cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "Bcell", "E42K_E42KU", "volcanoplot_E42K_E42KU_Bcell.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset Treg cells and run limma
seu_E42K_E42KU_Treg <- subset(seu_E42K_E42KU, subset = predicted.celltype.l2 %in% c("Treg"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_E42K_E42KU_Treg, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "E42K_E42KU", "Treg", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "E42K_E42KU", "Treg", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "E42K_E42KU", "Treg", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/E42K_E42KU/Treg/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/E42K_E42KU/Treg")
}

expr <- read.table(here("merged", "sanity", "E42K_E42KU", "Treg", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

cohort <- as.character(seu_E42K_E42KU_Treg$cohort)
patient <- as.character(seu_E42K_E42KU_Treg$patient)

#limma
design <- model.matrix(~0 + cohort)
colnames(design) <- sub("cohort", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/Treg/E42K_E42KU/corfit_E42K_E42KU_Treg.RData")

cm <- makeContrasts(E42K_E42KU=E42K-E42K_unaffected,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/Treg/E42K_E42KU/limmafit_E42K_E42KU_Treg.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# E42K_E42KU
E42K_E42KU <- topTable(fit, coef = "E42K_E42KU", number=Inf, adjust.method = "bonferroni")
sum(E42K_E42KU$adj.P.Val < 0.05)
E42K_E42KU$gene <- sub(".*_", "", rownames(E42K_E42KU))
write.csv(E42K_E42KU, file="/folder/containing/cellranger/output/merged/Treg/E42K_E42KU/E42K_E42KU_Treg_topTable.csv")

gseatt <- data.frame(gene=E42K_E42KU$gene, t=E42K_E42KU$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/Treg/E42K_E42KU/E42K_E42KU_Treg_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/E42K_E42KU/E42K_E42KU_Treg_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.E42K_E42KU_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/E42K_E42KU/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/E42K_E42KU/E42K_E42KU_Treg_gsea.rnk -scoring_scheme weighted -rpt_label curated.E42K_E42KU_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/E42K_E42KU/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/E42K_E42KU/E42K_E42KU_Treg_gsea.rnk -scoring_scheme weighted -rpt_label GO.E42K_E42KU_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/E42K_E42KU/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/E42K_E42KU/E42K_E42KU_Treg_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.E42K_E42KU_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/E42K_E42KU/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/E42K_E42KU/E42K_E42KU_Treg_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.E42K_E42KU_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/E42K_E42KU/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "E42K_E42KU"],
  logP = -log10(fit$p.value[, "E42K_E42KU"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "E42K_E42KU"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("E42K vs E42K unaffected - Treg Cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "Treg", "E42K_E42KU", "volcanoplot_E42K_E42KU_Treg.pdf"), plot = volcano_plot, width = 16, height = 12)

##### EXPERIMENT 7 - HBD VS PRE #####
# Create seu_HBD_pre by subsetting seu
seu_HBD_pre <- subset(seu, subset = treatment %in% c("pre_treatment", "healthy_control"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_HBD_pre, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "HBD_pre", "all_cells", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "HBD_pre", "all_cells", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "HBD_pre", "all_cells", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/HBD_pre/all_cells/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/HBD_pre/all_cells")
}

expr <- read.table(here("merged", "sanity", "HBD_pre", "all_cells", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

treatment <- as.character(seu_HBD_pre$treatment)
patient <- as.character(seu_HBD_pre$patient)

#limma
design <- model.matrix(~0 + treatment)
colnames(design) <- sub("treatment", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/all_cells/HBD_pre/corfit_HBD_pre_all_cells.RData")

cm <- makeContrasts(HBD_pre=healthy_control-pre_treatment,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/all_cells/HBD_pre/limmafit_HBD_pre_all_cells.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

HBD_pre <- topTable(fit, coef = "HBD_pre", number=Inf, adjust.method = "bonferroni")
sum(HBD_pre$adj.P.Val < 0.05)
HBD_pre$gene <- sub(".*_", "", rownames(HBD_pre))
write.csv(HBD_pre, file="/folder/containing/cellranger/output/merged/all_cells/HBD_pre/topTable_HBD_pre_all_cells.csv")

gseatt <- data.frame(gene=HBD_pre$gene, t=HBD_pre$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/all_cells/HBD_pre/HBD_pre_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/HBD_pre/HBD_pre_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.HBD_pre -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/HBD_pre/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/HBD_pre/HBD_pre_gsea.rnk -scoring_scheme weighted -rpt_label curated.HBD_pre -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/HBD_pre/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/HBD_pre/HBD_pre_gsea.rnk -scoring_scheme weighted -rpt_label GO.HBD_pre -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/HBD_pre/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/HBD_pre/HBD_pre_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.HBD_pre -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/HBD_pre/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/HBD_pre/HBD_pre_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.HBD_pre -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/HBD_pre/ ")

volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "HBD_pre"],
  logP = -log10(fit$p.value[, "HBD_pre"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "HBD_pre"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("healthy controls vs pre-treatment - All cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "all_cells", "HBD_pre", "volcanoplot_HBD_pre_all_cells.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset CD4 cells and run limma
seu_HBD_pre_CD4 <- subset(seu_HBD_pre, subset = predicted.celltype.l2 %in% c("CD4 TCM", "CD4 TEM", "CD4 Naive", "CD4 CTL"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_HBD_pre_CD4, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "HBD_pre", "CD4", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "HBD_pre", "CD4", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "HBD_pre", "CD4", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/HBD_pre/CD4/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/HBD_pre/CD4")
}

expr <- read.table(here("merged", "sanity", "HBD_pre", "CD4", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

treatment <- as.character(seu_HBD_pre_CD4$treatment)
patient <- as.character(seu_HBD_pre_CD4$patient)

#limma
design <- model.matrix(~0 + treatment)
colnames(design) <- sub("treatment", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/CD4/HBD_pre/corfit_HBD_pre_CD4.RData")

cm <- makeContrasts(HBD_pre=healthy_control-pre_treatment,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/CD4/HBD_pre/limmafit_HBD_pre_CD4.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# HBD_pre
HBD_pre <- topTable(fit, coef = "HBD_pre", number=Inf, adjust.method = "bonferroni")
sum(HBD_pre$adj.P.Val < 0.05)
HBD_pre$gene <- sub(".*_", "", rownames(HBD_pre))
write.csv(HBD_pre, file="/folder/containing/cellranger/output/merged/CD4/HBD_pre/HBD_pre_CD4_topTable.csv")

gseatt <- data.frame(gene=HBD_pre$gene, t=HBD_pre$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/CD4/HBD_pre/HBD_pre_CD4_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/HBD_pre/HBD_pre_CD4_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.HBD_pre_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/HBD_pre/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/HBD_pre/HBD_pre_CD4_gsea.rnk -scoring_scheme weighted -rpt_label curated.HBD_pre_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/HBD_pre/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/HBD_pre/HBD_pre_CD4_gsea.rnk -scoring_scheme weighted -rpt_label GO.HBD_pre_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/HBD_pre/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/HBD_pre/HBD_pre_CD4_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.HBD_pre_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/HBD_pre/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/HBD_pre/HBD_pre_CD4_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.HBD_pre_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/HBD_pre/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "HBD_pre"],
  logP = -log10(fit$p.value[, "HBD_pre"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "HBD_pre"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("healthy controls vs pre-treatment - CD4 cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "CD4", "HBD_pre", "volcanoplot_HBD_pre_CD4.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset CD8 cells and run limma
seu_HBD_pre_CD8 <- subset(seu_HBD_pre, subset = predicted.celltype.l2 %in% c("CD8 Naive", "CD8 TEM", "CD8 TCM"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_HBD_pre_CD8, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "HBD_pre", "CD8", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "HBD_pre", "CD8", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "HBD_pre", "CD8", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/HBD_pre/CD8/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/HBD_pre/CD8")
}

expr <- read.table(here("merged", "sanity", "HBD_pre", "CD8", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

treatment <- as.character(seu_HBD_pre_CD8$treatment)
patient <- as.character(seu_HBD_pre_CD8$patient)

#limma
design <- model.matrix(~0 + treatment)
colnames(design) <- sub("treatment", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/CD8/HBD_pre/corfit_HBD_pre_CD8.RData")

cm <- makeContrasts(HBD_pre=healthy_control-pre_treatment,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/CD8/HBD_pre/limmafit_HBD_pre_CD8.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# HBD_pre
HBD_pre <- topTable(fit, coef = "HBD_pre", number=Inf, adjust.method = "bonferroni")
sum(HBD_pre$adj.P.Val < 0.05)
HBD_pre$gene <- sub(".*_", "", rownames(HBD_pre))
write.csv(HBD_pre, file="/folder/containing/cellranger/output/merged/CD8/HBD_pre/HBD_pre_CD8_topTable.csv")

gseatt <- data.frame(gene=HBD_pre$gene, t=HBD_pre$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/CD8/HBD_pre/HBD_pre_CD8_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/HBD_pre/HBD_pre_CD8_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.HBD_pre_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/HBD_pre/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/HBD_pre/HBD_pre_CD8_gsea.rnk -scoring_scheme weighted -rpt_label curated.HBD_pre_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/HBD_pre/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/HBD_pre/HBD_pre_CD8_gsea.rnk -scoring_scheme weighted -rpt_label GO.HBD_pre_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/HBD_pre/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/HBD_pre/HBD_pre_CD8_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.HBD_pre_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/HBD_pre/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/HBD_pre/HBD_pre_CD8_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.HBD_pre_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/HBD_pre/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "HBD_pre"],
  logP = -log10(fit$p.value[, "HBD_pre"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "HBD_pre"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("healthy controls vs pre-treatment - CD8 cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "CD8", "HBD_pre", "volcanoplot_HBD_pre_CD8.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset B cells cells and run limma
seu_HBD_pre_Bcell <- subset(seu_HBD_pre, subset = predicted.celltype.l2 %in% c("B naive", "B intermediate", "B memory"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_HBD_pre_Bcell, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "HBD_pre", "Bcell", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "HBD_pre", "Bcell", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "HBD_pre", "Bcell", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/HBD_pre/Bcell/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/HBD_pre/Bcell")
}

expr <- read.table(here("merged", "sanity", "HBD_pre", "Bcell", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

treatment <- as.character(seu_HBD_pre_Bcell$treatment)
patient <- as.character(seu_HBD_pre_Bcell$patient)

#limma
design <- model.matrix(~0 + treatment)
colnames(design) <- sub("treatment", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/Bcell/HBD_pre/corfit_HBD_pre_Bcell.RData")

cm <- makeContrasts(HBD_pre=healthy_control-pre_treatment,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/Bcell/HBD_pre/limmafit_HBD_pre_Bcell.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# HBD_pre
HBD_pre <- topTable(fit, coef = "HBD_pre", number=Inf, adjust.method = "bonferroni")
sum(HBD_pre$adj.P.Val < 0.05)
HBD_pre$gene <- sub(".*_", "", rownames(HBD_pre))
write.csv(HBD_pre, file="/folder/containing/cellranger/output/merged/Bcell/HBD_pre/HBD_pre_Bcell_topTable.csv")

gseatt <- data.frame(gene=HBD_pre$gene, t=HBD_pre$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/Bcell/HBD_pre/HBD_pre_Bcell_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/HBD_pre/HBD_pre_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.HBD_pre_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/HBD_pre/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/HBD_pre/HBD_pre_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label curated.HBD_pre_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/HBD_pre/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/HBD_pre/HBD_pre_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label GO.HBD_pre_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/HBD_pre/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/HBD_pre/HBD_pre_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.HBD_pre_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/HBD_pre/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/HBD_pre/HBD_pre_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.HBD_pre_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/HBD_pre/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "HBD_pre"],
  logP = -log10(fit$p.value[, "HBD_pre"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "HBD_pre"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("healthy controls vs pre-treatment - B cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "Bcell", "HBD_pre", "volcanoplot_HBD_pre_Bcell.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset Treg cells and run limma
seu_HBD_pre_Treg <- subset(seu_HBD_pre, subset = predicted.celltype.l2 %in% c("Treg"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_HBD_pre_Treg, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "HBD_pre", "Treg", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "HBD_pre", "Treg", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "HBD_pre", "Treg", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/HBD_pre/Treg/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/HBD_pre/Treg")
}

expr <- read.table(here("merged", "sanity", "HBD_pre", "Treg", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

treatment <- as.character(seu_HBD_pre_Treg$treatment)
patient <- as.character(seu_HBD_pre_Treg$patient)

#limma
design <- model.matrix(~0 + treatment)
colnames(design) <- sub("treatment", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/Treg/HBD_pre/corfit_HBD_pre_Treg.RData")

cm <- makeContrasts(HBD_pre=healthy_control-pre_treatment,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/Treg/HBD_pre/limmafit_HBD_pre_Treg.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# HBD_pre
HBD_pre <- topTable(fit, coef = "HBD_pre", number=Inf, adjust.method = "bonferroni")
sum(HBD_pre$adj.P.Val < 0.05)
HBD_pre$gene <- sub(".*_", "", rownames(HBD_pre))
write.csv(HBD_pre, file="/folder/containing/cellranger/output/merged/Treg/HBD_pre/HBD_pre_Treg_topTable.csv")

gseatt <- data.frame(gene=HBD_pre$gene, t=HBD_pre$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/Treg/HBD_pre/HBD_pre_Treg_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/HBD_pre/HBD_pre_Treg_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.HBD_pre_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/HBD_pre/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/HBD_pre/HBD_pre_Treg_gsea.rnk -scoring_scheme weighted -rpt_label curated.HBD_pre_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/HBD_pre/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/HBD_pre/HBD_pre_Treg_gsea.rnk -scoring_scheme weighted -rpt_label GO.HBD_pre_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/HBD_pre/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/HBD_pre/HBD_pre_Treg_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.HBD_pre_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/HBD_pre/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/HBD_pre/HBD_pre_Treg_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.HBD_pre_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/HBD_pre/ ")

# This plot doesn't work because I couldn't generate fit
# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "HBD_pre"],
  logP = -log10(fit$p.value[, "HBD_pre"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "HBD_pre"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 2) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("healthy controls vs pre-treatment - Treg Cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "Treg", "HBD_pre", "volcanoplot_HBD_pre_Treg.pdf"), plot = volcano_plot, width = 16, height = 12)

##### EXPERIMENT 8 - HBD VS POST #####
# Create seu_HBD_post by subsetting seu
seu_HBD_post <- subset(seu, subset = treatment %in% c("post_treatment", "healthy_control"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_HBD_post, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "HBD_post", "all_cells", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "HBD_post", "all_cells", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "HBD_post", "all_cells", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/HBD_post/all_cells/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/HBD_post/all_cells")
}

expr <- read.table(here("merged", "sanity", "HBD_post", "all_cells", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

treatment <- as.character(seu_HBD_post$treatment)
patient <- as.character(seu_HBD_post$patient)

#limma
design <- model.matrix(~0 + treatment)
colnames(design) <- sub("treatment", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/all_cells/HBD_post/corfit_HBD_post_all_cells.RData")

cm <- makeContrasts(HBD_post=healthy_control-post_treatment,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/all_cells/HBD_post/limmafit_HBD_post_all_cells.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

HBD_post <- topTable(fit, coef = "HBD_post", number=Inf, adjust.method = "bonferroni")
sum(HBD_post$adj.P.Val < 0.05)
HBD_post$gene <- sub(".*_", "", rownames(HBD_post))
write.csv(HBD_post, file="/folder/containing/cellranger/output/merged/all_cells/HBD_post/topTable_HBD_post_all_cells.csv")

gseatt <- data.frame(gene=HBD_post$gene, t=HBD_post$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/all_cells/HBD_post/HBD_post_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/HBD_post/HBD_post_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.HBD_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/HBD_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/HBD_post/HBD_post_gsea.rnk -scoring_scheme weighted -rpt_label curated.HBD_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/HBD_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/HBD_post/HBD_post_gsea.rnk -scoring_scheme weighted -rpt_label GO.HBD_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/HBD_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/HBD_post/HBD_post_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.HBD_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/HBD_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/all_cells/HBD_post/HBD_post_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.HBD_post -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/all_cells/HBD_post/ ")

volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "HBD_post"],
  logP = -log10(fit$p.value[, "HBD_post"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "HBD_post"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("healthy controls vs post-treatment - All cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "all_cells", "HBD_post", "volcanoplot_HBD_post_all_cells.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset CD4 cells and run limma
seu_HBD_post_CD4 <- subset(seu_HBD_post, subset = predicted.celltype.l2 %in% c("CD4 TCM", "CD4 TEM", "CD4 Naive", "CD4 CTL"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_HBD_post_CD4, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "HBD_post", "CD4", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "HBD_post", "CD4", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "HBD_post", "CD4", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/HBD_post/CD4/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/HBD_post/CD4")
}

expr <- read.table(here("merged", "sanity", "HBD_post", "CD4", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

treatment <- as.character(seu_HBD_post_CD4$treatment)
patient <- as.character(seu_HBD_post_CD4$patient)

#limma
design <- model.matrix(~0 + treatment)
colnames(design) <- sub("treatment", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/CD4/HBD_post/corfit_HBD_post_CD4.RData")

cm <- makeContrasts(HBD_post=healthy_control-post_treatment,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/CD4/HBD_post/limmafit_HBD_post_CD4.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# HBD_post
HBD_post <- topTable(fit, coef = "HBD_post", number=Inf, adjust.method = "bonferroni")
sum(HBD_post$adj.P.Val < 0.05)
HBD_post$gene <- sub(".*_", "", rownames(HBD_post))
write.csv(HBD_post, file="/folder/containing/cellranger/output/merged/CD4/HBD_post/HBD_post_CD4_topTable.csv")

gseatt <- data.frame(gene=HBD_post$gene, t=HBD_post$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/CD4/HBD_post/HBD_post_CD4_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/HBD_post/HBD_post_CD4_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.HBD_post_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/HBD_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/HBD_post/HBD_post_CD4_gsea.rnk -scoring_scheme weighted -rpt_label curated.HBD_post_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/HBD_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/HBD_post/HBD_post_CD4_gsea.rnk -scoring_scheme weighted -rpt_label GO.HBD_post_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/HBD_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/HBD_post/HBD_post_CD4_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.HBD_post_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/HBD_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD4/HBD_post/HBD_post_CD4_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.HBD_post_CD4 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD4/HBD_post/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "HBD_post"],
  logP = -log10(fit$p.value[, "HBD_post"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "HBD_post"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("healthy controls vs post-treatment - CD4 cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "CD4", "HBD_post", "volcanoplot_HBD_post_CD4.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset CD8 cells and run limma
seu_HBD_post_CD8 <- subset(seu_HBD_post, subset = predicted.celltype.l2 %in% c("CD8 Naive", "CD8 TEM", "CD8 TCM"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_HBD_post_CD8, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "HBD_post", "CD8", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "HBD_post", "CD8", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "HBD_post", "CD8", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/HBD_post/CD8/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/HBD_post/CD8")
}

expr <- read.table(here("merged", "sanity", "HBD_post", "CD8", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

treatment <- as.character(seu_HBD_post_CD8$treatment)
patient <- as.character(seu_HBD_post_CD8$patient)

#limma
design <- model.matrix(~0 + treatment)
colnames(design) <- sub("treatment", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/CD8/HBD_post/corfit_HBD_post_CD8.RData")

cm <- makeContrasts(HBD_post=healthy_control-post_treatment,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/CD8/HBD_post/limmafit_HBD_post_CD8.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# HBD_post
HBD_post <- topTable(fit, coef = "HBD_post", number=Inf, adjust.method = "bonferroni")
sum(HBD_post$adj.P.Val < 0.05)
HBD_post$gene <- sub(".*_", "", rownames(HBD_post))
write.csv(HBD_post, file="/folder/containing/cellranger/output/merged/CD8/HBD_post/HBD_post_CD8_topTable.csv")

gseatt <- data.frame(gene=HBD_post$gene, t=HBD_post$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/CD8/HBD_post/HBD_post_CD8_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/HBD_post/HBD_post_CD8_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.HBD_post_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/HBD_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/HBD_post/HBD_post_CD8_gsea.rnk -scoring_scheme weighted -rpt_label curated.HBD_post_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/HBD_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/HBD_post/HBD_post_CD8_gsea.rnk -scoring_scheme weighted -rpt_label GO.HBD_post_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/HBD_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/HBD_post/HBD_post_CD8_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.HBD_post_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/HBD_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/CD8/HBD_post/HBD_post_CD8_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.HBD_post_CD8 -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/CD8/HBD_post/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "HBD_post"],
  logP = -log10(fit$p.value[, "HBD_post"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "HBD_post"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("healthy controls vs post-treatment - CD8 cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "CD8", "HBD_post", "volcanoplot_HBD_post_CD8.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset B cells cells and run limma
seu_HBD_post_Bcell <- subset(seu_HBD_post, subset = predicted.celltype.l2 %in% c("B naive", "B intermediate", "B memory"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_HBD_post_Bcell, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "HBD_post", "Bcell", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "HBD_post", "Bcell", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "HBD_post", "Bcell", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/HBD_post/Bcell/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/HBD_post/Bcell")
}

expr <- read.table(here("merged", "sanity", "HBD_post", "Bcell", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

treatment <- as.character(seu_HBD_post_Bcell$treatment)
patient <- as.character(seu_HBD_post_Bcell$patient)

#limma
design <- model.matrix(~0 + treatment)
colnames(design) <- sub("treatment", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/Bcell/HBD_post/corfit_HBD_post_Bcell.RData")

cm <- makeContrasts(HBD_post=healthy_control-post_treatment,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/Bcell/HBD_post/limmafit_HBD_post_Bcell.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# HBD_post
HBD_post <- topTable(fit, coef = "HBD_post", number=Inf, adjust.method = "bonferroni")
sum(HBD_post$adj.P.Val < 0.05)
HBD_post$gene <- sub(".*_", "", rownames(HBD_post))
write.csv(HBD_post, file="/folder/containing/cellranger/output/merged/Bcell/HBD_post/HBD_post_Bcell_topTable.csv")

gseatt <- data.frame(gene=HBD_post$gene, t=HBD_post$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/Bcell/HBD_post/HBD_post_Bcell_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/HBD_post/HBD_post_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.HBD_post_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/HBD_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/HBD_post/HBD_post_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label curated.HBD_post_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/HBD_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/HBD_post/HBD_post_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label GO.HBD_post_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/HBD_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/HBD_post/HBD_post_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.HBD_post_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/HBD_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Bcell/HBD_post/HBD_post_Bcell_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.HBD_post_Bcell -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Bcell/HBD_post/ ")

# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "HBD_post"],
  logP = -log10(fit$p.value[, "HBD_post"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "HBD_post"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("healthy controls vs post-treatment - B cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "Bcell", "HBD_post", "volcanoplot_HBD_post_Bcell.pdf"), plot = volcano_plot, width = 16, height = 12)

# Subset Treg cells and run limma
seu_HBD_post_Treg <- subset(seu_HBD_post, subset = predicted.celltype.l2 %in% c("Treg"))

#Transform with Sanity
countmatrix <- as.matrix(Seurat::GetAssayData(seu_HBD_post_Treg, assay = "RNA"))
forsanity <- rbind(colnames(countmatrix), countmatrix)
forsanity <- cbind(c("GeneID", rownames(countmatrix)), forsanity)

# Check whether you need to generate the table and if so, do so
out <- here("merged", "sanity", "HBD_post", "Treg", "countmatrix")
if(!file.exists(out)) {
  write.table(forsanity, file=here("merged", "sanity", "HBD_post", "Treg", "countmatrix.txt"), sep="\t", row.names = F, col.names = F, quote=F)
}

# Check if you need to run Sanity and if so, do so
out <- here("merged", "sanity", "HBD_post", "Treg", "likelihood.txt")
if(!file.exists(out)) {
  system("/path/to/Sanity -n 28 -e 1 -f /folder/containing/cellranger/output/merged/sanity/HBD_post/Treg/countmatrix.txt -d /folder/containing/cellranger/output/merged/sanity/HBD_post/Treg")
}

expr <- read.table(here("merged", "sanity", "HBD_post", "Treg", "log_transcription_quotients.txt"), row.names = 1, header = T, stringsAsFactors = F)
expr <- data.matrix(expr)
expr <- expr-min(expr)

treatment <- as.character(seu_HBD_post_Treg$treatment)
patient <- as.character(seu_HBD_post_Treg$patient)

#limma
design <- model.matrix(~0 + treatment)
colnames(design) <- sub("treatment", "", colnames(design))
#Do corfit in parallel
library(statmod)
library(parallel)
ngenes <- nrow(expr)
narrays <- ncol(expr)
nbeta <- ncol(design)
QR <- qr(design)
MaxBlockSize <- max(table(patient))
weights <- getEAWP(expr)$weights
Array <- patient
nafun <- function(e) NA

engine <- function(i) {
  y <- drop(expr[i, ])
  o <- is.finite(y)
  A <- factor(Array[o])
  nobs <- sum(o)
  nblocks <- length(levels(A))
  if (nobs > (nbeta + 2L) && nblocks > 1L && nblocks < 
      (nobs - 1L)) {
    y <- y[o]
    X <- design[o, , drop = FALSE]
    Z <- model.matrix(~0 + A)
    if (!is.null(weights)) {
      w <- drop(weights[i, ])[o]
      s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                             X, Z, w, only.varcomp = TRUE, maxit = 20)$varcomp), 
                    error = nafun)
    }
    else s <- tryCatch(suppressWarnings(statmod::mixedModel2Fit(y, 
                                                                X, Z, only.varcomp = TRUE, maxit = 20)$varcomp), 
                       error = nafun)
    if (!is.na(s[1])) 
      s[2]/sum(s)
  }
}

# Set up parallel backend
num_cores <- 2

# Register the parallel backend
registerDoParallel(num_cores)

# Define the function that will be applied in parallel
calculate_rho <- function(i) {
  engine(i)
}

# Run the calculation in parallel
rho <- foreach(i = 1:ngenes, .combine = "c") %dopar% {
  calculate_rho(i)
}

# Stop the parallel backend
stopImplicitCluster()

rhomax <- 0.99
rhomin <- 1/(1 - MaxBlockSize) + 0.01
m <- min(rho, 0, na.rm = TRUE)
if (m < rhomin) 
  rho[rho < rhomin] <- rhomin
m <- max(rho, 0, na.rm = TRUE)
if (m > rhomax) 
  rho[rho > rhomax] <- rhomax
arho <- atanh(rho)
mrho <- tanh(mean(arho, trim = 0.15, na.rm = TRUE))
corfit <- list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
save(corfit, file="/folder/containing/cellranger/output/merged/Treg/HBD_post/corfit_HBD_post_Treg.RData")

cm <- makeContrasts(HBD_post=healthy_control-post_treatment,
                    levels=design)

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
save(fit, file="/folder/containing/cellranger/output/merged/Treg/HBD_post/limmafit_HBD_post_Treg.RData")
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)

# HBD_post
HBD_post <- topTable(fit, coef = "HBD_post", number=Inf, adjust.method = "bonferroni")
sum(HBD_post$adj.P.Val < 0.05)
HBD_post$gene <- sub(".*_", "", rownames(HBD_post))
write.csv(HBD_post, file="/folder/containing/cellranger/output/merged/Treg/HBD_post/HBD_post_Treg_topTable.csv")

gseatt <- data.frame(gene=HBD_post$gene, t=HBD_post$t)
gseatt <- gseatt[!duplicated(gseatt$gene),]
gseatt <- gseatt[order(gseatt$t, decreasing = T),]
write.table(gseatt, file="/folder/containing/cellranger/output/merged/Treg/HBD_post/HBD_post_Treg_gsea.rnk", sep="\t", row.names = F, col.names = F, quote=F)

system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/h.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/HBD_post/HBD_post_Treg_gsea.rnk -scoring_scheme weighted -rpt_label hallmark.HBD_post_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/HBD_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c2.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/HBD_post/HBD_post_Treg_gsea.rnk -scoring_scheme weighted -rpt_label curated.HBD_post_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/HBD_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c5.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/HBD_post/HBD_post_Treg_gsea.rnk -scoring_scheme weighted -rpt_label GO.HBD_post_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/HBD_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c7.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/HBD_post/HBD_post_Treg_gsea.rnk -scoring_scheme weighted -rpt_label immunologic.HBD_post_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/HBD_post/ ")
system("/path/to/GSEA/gsea-cli.sh GSEAPreranked -gmx /path/to/GSEA/databases/c8.all.v2023.1.Hs.symbols.gmt -norm meandiv -nperm 10000 -rnk /folder/containing/cellranger/output/merged/Treg/HBD_post/HBD_post_Treg_gsea.rnk -scoring_scheme weighted -rpt_label cell_type_signature.HBD_post_Treg -create_svgs false -make_sets true -plot_top_x 40 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /folder/containing/cellranger/output/merged/Treg/HBD_post/ ")

# This plot doesn't work because I couldn't generate fit
# Data frame for volcano plot
volcano_df <- data.frame(
  gene = rownames(fit$coefficients),
  logFC = fit$coefficients[, "HBD_post"],
  logP = -log10(fit$p.value[, "HBD_post"])
)

# Get the top 50 genes with the most significant adjusted p-value
top_genes <- head(volcano_df[order(fit$p.value[, "HBD_post"]), ], 50)

# Create the volcano plot
volcano_plot <- ggplot(data = volcano_df, aes(x = logFC, y = logP)) +
  geom_point(color = ifelse(volcano_df$logP > -log10(0.01), "red", "black"), size = 6) +
  geom_text_repel(data = top_genes, aes(label = gene), box.padding = 0.5, point.padding = 0.5, size = 3) +
  xlab("Log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("healthy controls vs post-treatment - Treg Cells") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(), panel.background = element_blank()) +
  theme(axis.line = element_line(color = "black"))

ggsave(here("merged", "Treg", "HBD_post", "volcanoplot_HBD_post_Treg.pdf"), plot = volcano_plot, width = 16, height = 12)
