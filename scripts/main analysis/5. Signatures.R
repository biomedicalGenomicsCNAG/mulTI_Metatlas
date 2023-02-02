#automatic signatures across the tissues

sub = #add the rds file of interest

cells<-table(sub$Tissue)
cells = as.data.frame(cells)
cols = cells$Var1
total<-matrix(0, nrow = 1, ncol = 2)
total = as.data.frame(total)
names(total)= c("features","V2")
totest = cells %>%  filter(Freq>= 40)


# average values
Idents(sub)= sub$Tissue
cluster.averages <- AverageExpression(sub, group.by = "Tissue")
head(cluster.averages[["RNA"]][, 1:5])
cluster.averages <- as.data.frame(cluster.averages[["RNA"]])
cluster.averages$genes = row.names(cluster.averages)
head(cluster.averages)
cluster.averages = cluster.averages %>% select(c(totest$Var1, "genes"))


#heatmap
for ( i in 1:nrow(totest))
{
  markers <- FindMarkers(sub, ident.1 = totest$Var1[i], group.by = 'Tissue', min.pct = 0.5, only.pos=T)
  features<-head(markers, n = 50)
  features<- markers %>% filter(avg_log2FC>1) %>% arrange(-p_val_adj)
  features <- na.omit(features)
  features <- features[!grepl("^RP", row.names(features)),]
  
  features<- as.data.frame(row.names(features))
  names(features) = "genes"
  
  mat = merge(cluster.averages, features, by= "genes")
  
  library(ComplexHeatmap)
  zmat <- t(scale(t(mat[,-1]), center=TRUE, scale=TRUE))
  row.names(zmat) = mat[,1]
  
  
  p = Heatmap(as.matrix(zmat), clustering_distance_columns = "euclidean",
              clustering_method_columns = "complete", show_row_names = T,
              column_title = paste(totest$Var1[i]))
  
  
  ggarrange(table.p,
            ncol = 1, nrow = 1,
            heights = c(0.5, 1))  %>%
    ggexport(filename = paste0("B cells_level1_heatmap_allgenes","_",totest$Var1[i], ".pdf"),  width = 15, height = 9)
  
  dev.off()
  
}
