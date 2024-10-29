# MetaConsensusDA

MetaConsensusDA is an R package for microbiome analysis using multiple algorithms - reaching consensus.

## Installing consensusDE

To obtain the original version from github, install devtools in R and use the following:

```R

required packages

# Install ggplot2 
install.packages("ggplot2")

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")
BiocManager::install("DESeq2")
BiocManager::install("ALDEx2")
BiocManager::install("edgeR")

# Load the packages
library(ggplot2)
library(phyloseq)
library(DESeq2)
library(ALDEx2)
library(edgeR)

```

## Examples

To run consensusDE, load the library and follow the examples in the vignette.

```R
library(MetaConsensusDA)
```

## Contact

For more details, contact Manoharan Kumar:
manoharan.kumar@jcu.edu.au
# MetaConsensusDA
