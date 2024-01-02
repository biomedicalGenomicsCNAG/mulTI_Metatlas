library(data.table)
library(Seurat)
library(ggpubr)
library(tidyverse)
library(UCell)




# expression genes including in TNFa via NFkB
VlnPlot(sub1, features = c("TNF","NFKBIA","NFKBID","NFKB1","NFKB2", "REL", "RELA", "RELB"), group.by = "Tissue")
Idents(sub)= sub$Tissue
sub1= subset(sub, subset = Tissue == totest$Var1)
unique(sub1$Tissue)
NFKB_story= c("TNF","NFKBIA","NFKBID","NFKB1","NFKB2", "REL", "RELA", "RELB")

DotPlot(sub1, features = NFKB_story) + RotatedAxis()  + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")


# activity of NFkB target genes
markers$NFKB <- as.data.frame(fread('NFkBtargets_regulonsAnalysis.txt'))$gene

target <- AddModuleScore_UCell(
  obj = sub1,
  features = markers
)


target@meta.data %>% group_by(Tissue) %>% 
  summarise(count = n()) 

targets_score=target@meta.data %>% group_by(Tissue) %>%
summarise_at(vars(NFKB_UCell), funs(mean(., na.rm=TRUE)))

summary(targets_score)

signature.names <- paste0(names(markers), "_UCell")
VlnPlot(target, features = signature.names, group.by = "Tissue", cols = col, sort = T) + 
  geom_hline(yintercept=0.5484, linetype="dashed", color = "black", size=1.2)



#overlap of target genes and diff expression
markers <- FindAllMarkers(sub1, only.pos = TRUE, min.pct = 0.5)
markers = markers %>% filter(avg_logFC>1)
markers1 = markers %>% filter(markers$gene %in% NFKB_reg)
markers_plot = merge(markers1, cluster.averages, by.y= "genes", by.x="gene")
markers_plot=markers_plot[!duplicated(markers_plot$gene),]

zmat <- t(scale(t(markers_plot[,-c(1:7)]), center=TRUE, scale=TRUE))
row.names(zmat) = markers_plot[,1]
library(pals)
palette = as.vector(polychrome(14))
palette = c("#5A5156", "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32",
            "#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C",
            "#8A5F00", "#DEA0FD", "#AA0DFE", "#F8A19F")

library(ComplexHeatmap)
col = setNames(palette, 
               c( "prostate gland", "eye","heart" , 
                  "lung" , "skeletal muscle organ", 
                  "esophagus","blood" , "brain","kidney" ,
                  "liver", "colon", 
                  "bone marrow", "spleen","thymus"))
column_ha = HeatmapAnnotation( 
  Clusters = totest$Var1,
  col = list(Clusters = col))


Heatmap(as.matrix(zmat), clustering_distance_columns = "euclidean",
        top_annotation = column_ha,show_column_names = F,
        row_names_gp = gpar(fontsize = 6),
        clustering_method_columns = "complete", show_row_names = F)

# overlap of cytokines KEGG pathway with cytokines driven by NFkB
chem = fread("geneset_KEGG_2.txt") #cytokines pathway from KEGG
diff = as.data.frame(chem)
names(diff) = "diff"
cytokines = intersect(markers_plot$gene, diff$diff) 

cyto_plot = markers_plot %>% filter(markers_plot$gene %in% cytokines)

zmat <- t(scale(t(cyto_plot[,-c(1:7)]), center=TRUE, scale=TRUE))
row.names(zmat) = cyto_plot[,1]

Heatmap(as.matrix(zmat), clustering_distance_columns = "euclidean",
        top_annotation = column_ha,show_column_names = F,
        row_names_gp = gpar(fontsize = 9),
        clustering_method_columns = "complete", show_row_names = T)

DotPlot(sub1, features = cytokines) + RotatedAxis()  + coord_flip() +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")

