## Readme #####

In house piple line for short reads Bulk RNA Seq Analyses.


## This pipleline #####

This pipeline test on Human, Mouse and non-model organisms


## Run through steps ####

#Sanity check fastq file size, etc

Step 1. Fastq QC raw reads

Step 2. Trimming the low quality  and adapter 

Step 3. Fastq QC clean reads

Step 3. STAR alignment

Step 4. bam Indexing

Step 5. Differential expressino analyses using ConsensusDE

Step 6a. Enrichment analyses using GSEA for model organisms (Human and Mouse)

Step 6b. KEGG/GO for Non-model organisms

Step 7. Over representation analyses for Human and Mouse

