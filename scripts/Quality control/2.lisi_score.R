# library("devtools")
# devtools::install_github("immunogenomics/lisi")

library(lisi)
library(data.table)
library(Seurat)
library(ggpubr)
library(tidyverse)


immune.combined_un <- readRDS('object.rds') #upload object before harmony correction
saveRDS(immune.combined_un@reductions$pca@cell.embeddings, "uncorrected_pca_coordinates.rds")
write_csv(as.data.frame(immune.combined_un@reductions$pca@cell.embeddings), "uncorrected_pca_coordinates.csv", col_names = TRUE)
write_csv(rownames_to_column(immune.combined_un@meta.data, var = "cell_barcode"), "metadata.csv", col_names = TRUE)
gc()

immune.combined <- readRDS('object_harmony.rds') #upload object after harmony correction
saveRDS(immune.combined@reductions$harmony@cell.embeddings, "harmony_corrected_pca_coordinates.rds")
write_csv(as.data.frame(immune.combined@reductions$pca@cell.embeddings), "harmony_corrected_pca_coordinates.csv", col_names = TRUE)
write_csv(rownames_to_column(immune.combined@meta.data, var = "cell_barcode"), "harmony_metadata.csv", col_names = TRUE)
gc()
# Load data

metadata <- read_csv("harmony_metadata.csv", col_names = TRUE)
uncorrected_coords <- readRDS("uncorrected_pca_coordinates.rds")
harmony_coords <- readRDS("harmony_corrected_pca_coordinates.rds")

uncorrected_coords <- subset(uncorrected_coords, rownames(uncorrected_coords) %in% rownames(harmony_coords))



meta_data <- data.frame(source = metadata$Sex)
rownames(meta_data) <- metadata$cell_barcode
if (all(rownames(meta_data) == rownames(uncorrected_coords))) {
  print("row names are equal")
} else {
  warning("row names are not equal!")
}

#devtools::install_github("immunogenomics/lisi")
library(lisi)


dim_red_mats <- list(
  uncorrected_coords,
  harmony_coords
)
names(dim_red_mats) <- c("uncorrected", "Harmony")
lisi_scores <- purrr::map(dim_red_mats, function(mat) {
  scores <- compute_lisi(
    X = mat,
    meta_data = meta_data,
    label_colnames = "source"
  )
})
lisi_scores <- bind_rows(lisi_scores, .id = "correction")
head(lisi_scores)


# Load data frame
library(ggridges)

ggplot(df, aes(x = depth, y = color)) +
  geom_density_ridges()
# Plot LISI
sorted_corrections <- c("uncorrected", "Harmony")
palette <- c("#999999",  "#612c63")
lisi_scores_gg <- lisi_scores %>%
  mutate(correction = factor(correction, levels = rev(sorted_corrections))) %>%
  ggplot(aes(source, correction, fill = correction)) +
  #geom_violin() +
  geom_density_ridges() + 
  scale_fill_manual(values = palette) +
  labs(x = "iLISI", y = "") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(color = "black", size = 13),
    axis.text.y = element_text(color = "black", size = 12),
    axis.text.x = element_text(size = 11)
  )  


lisi_scores_gg <-lisi_scores %>%
  mutate(correction = factor(correction, levels = rev(sorted_corrections))) %>%
  ggplot(aes(source, correction, fill = correction)) +
  geom_boxplot()  +
  scale_fill_manual(values = palette) +
  labs(x = "iLISI", y = "") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(color = "black", size = 13),
    axis.text.y = element_text(color = "black", size = 12),
    axis.text.x = element_text(size = 11)
  )  


write_csv(lisi_scores, "lisi_score_tissue.csv", col_names = TRUE)

# Save
ggsave(
  filename = "iLISI_Atlas_tissue.pdf",
  plot = lisi_scores_gg,
  width = 14,
  height = 8,
  units = "cm"
)
