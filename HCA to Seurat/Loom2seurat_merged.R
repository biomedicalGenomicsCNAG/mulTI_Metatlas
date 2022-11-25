# Load packages
library(Seurat)
library(ggpubr)
library(tidyverse)
library(loomR)
library(data.table)
library(rlist)

#paths
working_path = "" #add this path
save_obj_path = "" #add this path
save_merged_path ="" #add this path

#change this file
annotation_rds<- fread(paste0(working_path,"hca-manifest-XXXXX.tsv")) #add the manifest file form HCA


my_list  = list()


try( for (i in 1:nrow(annotation_rds)) {
  
  
  C1 <- connect(filename =   paste0(working_path,annotation_rds$bundle_uuid[i], "/", annotation_rds$file_name[i]),
                mode = "r+",  skip.validate = TRUE)
  
  s1 <- as.Seurat(C1)
  s1 <- AddMetaData(s1, metadata=annotation_rds$specimen_from_organism.organ[i], col.name="Tissue")
  s1 <- AddMetaData(s1, metadata=annotation_rds$project.project_core.project_short_name[i], col.name="Project")
  s1 <- AddMetaData(s1, metadata=annotation_rds$library_preparation_protocol.library_construction_approach[i], col.name="Library_prep")
  s1 <- AddMetaData(s1,  metadata=annotation_rds$bundle_uuid[i], col.name="Sample_ID")
  
  s1 <- subset(s1 , subset=PTPRC>1)
  saveRDS(s1, file = paste0(save_obj_path,annotation_rds$bundle_uuid[i], "_CD45.rds"))
  my_list [[ annotation_rds$bundle_uuid[i] ]] = s1 
  rm(s1, C1)
},
    silent = TRUE)


for (i in names(my_list)) {
  my_list[[i]] <- RenameCells(my_list[[i]],
                                         add.cell.id = i)
}


immune.combined <- merge(my_list [[ annotation_rds$bundle_uuid[1] ]],list.remove(my_list, annotation_rds$bundle_uuid[1]))
rm(my_list)

saveRDS(immune.combined, file = paste0(save_merged_path,"immune_combined.rds"))
