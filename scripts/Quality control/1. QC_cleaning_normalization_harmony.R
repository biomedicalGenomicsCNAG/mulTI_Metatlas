library(data.table)
library(scater)
library(Seurat)
library(ggpubr)
library(tidyverse)
library(loomR)
library(harmony)

# Data
# download the object from zenodo https://zenodo.org/record/7756209
s_1 <- readRDS('immune_combined.rds')


#calculate the percentage of mitochondrial expression
immune.combined[["percent.mt"]] <- PercentageFeatureSet(immune.combined, pattern = "^MT-")



##############################################
#Low quality cells c: 
#- low library size, 
#- low library complexity
#- high fraction of mitochondrial expression 


# Thresholds
VlnPlot(immune.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

lib_complexity_df <- immune.combined@meta.data %>%
  dplyr::group_by(Tissue, Sample_ID) %>%
  dplyr::summarise(
    mean_n_features = round(mean(nFeature_RNA), 2),
    sd_n_features = round(sd(nFeature_RNA), 2)) %>%
  dplyr::arrange(Tissue)
print(lib_complexity_df)


lib_size_per_tissue <- immune.combined@meta.data %>%
  ggplot(aes(Tissue, nFeature_RNA, fill = Tissue)) +
  geom_violin() +
  scale_y_log10() +
  labs(title = "", x = "", y = "Number of Detected Genes") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    legend.title = element_blank()
  ) 
lib_size_per_tissue



lib_size_hist <- immune.combined@meta.data %>%
  ggplot(aes(x=nCount_RNA)) +
  geom_histogram() + geom_vline(xintercept = min_lib_size, linetype = "dashed", color = "red") +
  geom_vline(xintercept = max_lib_size, linetype = "dashed", color = "red") +
  scale_x_continuous(limits = c(0, 30000))
lib_size_hist


n_genes_hist  <- immune.combined@meta.data %>%
  ggplot(aes(x=nFeature_RNA)) +
  geom_histogram() +
  geom_vline(xintercept = min_n_genes, linetype = "dashed", color = "red") +
  geom_vline(xintercept = max_n_genes, linetype = "dashed", color = "red")
n_genes_hist 


pct_mt_hist <- immune.combined@meta.data %>%
  ggplot(aes(x=percent.mt)) +
  geom_histogram() +
  geom_vline(xintercept = max_pct_mt, linetype = "dashed", color = "red") +
  scale_x_continuous(limits = c(0, 100))
pct_mt_hist


#choose cut offs based on the previous distribution
min_lib_size <- 1000
max_lib_size <- 20000
min_n_genes <- 200
max_n_genes <- 6000
max_pct_mt <- 15
min_cells <- 5

#Subset empty droplets and lysed cells
is_low_quality <- 
  immune.combined$nCount_RNA < min_lib_size |
  immune.combined$nCount_RNA > max_lib_size |
  immune.combined$nFeature_RNA < min_n_genes |
  immune.combined$nFeature_RNA > max_n_genes |
  immune.combined$percent.mt > max_pct_mt

metadata_before_qc <- immune.combined@meta.data
table(is_low_quality)
immune.combined$keep_cells <- !is_low_quality
Idents(immune.combined) <- "keep_cells"
immune.combined
immune.combined <- subset(immune.combined, idents = TRUE)
immune.combined
metadata_after_qc <- immune.combined@meta.data


qc_before <- metadata_before_qc %>%
  group_by(Sample_ID) %>%
  summarise(num_cells_before_qc = n())
qc_after <- metadata_after_qc %>%
  group_by(Sample_ID) %>%
  summarise(
    num_cells_after_qc = n(),
    average_library_size = mean(nCount_RNA),
    average_num_detected_genes = mean(nFeature_RNA),
    average_mitochondrial_fraction = mean(percent.mt)
  )
qc_table <- left_join(qc_before, qc_after, by = "Sample_ID")
DT::datatable(qc_table)

#filter out genes
n_cells <- Matrix::rowSums(immune.combined[["RNA"]]@counts > 0)
gene_qc <- n_cells %>% 
  as.data.frame() %>% 
  ggplot(aes(n_cells)) + 
  geom_histogram(bins = 100, alpha = 0.75) +
  scale_x_log10("Number of cells") +
  theme_bw() 
gene_qc +
  geom_vline(xintercept = min_cells, linetype = "dashed", color = "red")

kept_genes <- rownames(immune.combined)[n_cells > min_cells]
table(n_cells > min_cells)
rm(is_low_quality);rm(n_cells);rm(metadata_before_qc);rm(gene_qc)

immune.combined <- subset(immune.combined, features = kept_genes)

############EXPORT 
saveRDS(immune.combined, file = paste0("immune_combine_filtered.rds")) #change this


# Downstream analysis
immune.combined <- NormalizeData(immune.combined)
immune.combined <- FindVariableFeatures(immune.combined, selection.method = "vst", nfeatures = 2000)
immune.combined <- ScaleData(immune.combined)
immune.combined <- RunPCA(immune.combined)
immune.combined <- RunUMAP(immune.combined, dims = 1:10)

immune.combined <- FindNeighbors(immune.combined, dims = 1:10)
immune.combined <- FindClusters(immune.combined, resolution = 0.3)

immune.combined <- RunHarmony(immune.combined, group.by.vars = c("Sample_ID","Project"))
immune.combined <- RunUMAP(immune.combined, reduction = "harmony", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "harmony", dims = 1:30) %>% FindClusters()
DimPlot(immune.combined, group.by = c("Project", "Tissue", "Sample_ID"), ncol = 2)
DimPlot(immune.combined, label = TRUE)
