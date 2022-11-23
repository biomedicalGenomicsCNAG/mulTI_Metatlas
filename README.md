# mulTI_Metalas
## Title

## Author: Maria Tsagiopoulou

## Graphical summary
![graphImage](https://user-images.githubusercontent.com/19466299/179958319-6a34d3c4-536c-41bf-8fc0-1d1184d33220.png)

## Abstract


## Data
The data has been deposited in five levels of organization, from raw to processed data:

- raw data. All the new generated fastq files have been downloaded from HCA
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
### from anaconda
#### necessary for 
