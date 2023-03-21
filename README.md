# Multi tissue met@tlas
## Title: Multi-tissue single-cell transcriptomics of immune cells uncovered tissue-specific features 

## Author: Maria Tsagiopoulou

## Graphical summary
![Picture2](https://user-images.githubusercontent.com/19466299/203803295-3c1bf56c-c81b-4bbb-bede-6bfa66f33f50.png)


## Abstract


## Data
The data has been deposited in five levels of organization, from raw to processed data:

- raw data. All the new generated fastq files have been downloaded from HCA. Full table with the information regarding the used project is available on the supplementary table 1 of the published paper
- matrices. All the counts table have been deposited in [Zenodo](https://zenodo.org/record/XX)


## Prerequisites
The packages needed to be installed, in order to run the project are:

### from CRAN
```
install.packages(c("tidyverse", "ggplot2",  "stringi", "boot"))
```
### from Bioconductor
```
BiocManager::install(c("limma",  "ComplexHeatmap", "sva", "DESeq2"))
```
### https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/TonsilAtlas/tree/main/scRNA-seq/gene_regulatory_networks/pyScenic
