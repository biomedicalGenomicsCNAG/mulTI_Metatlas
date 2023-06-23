###PORPROTIONS

#Data
## sub= subset of the main Seurat object in one cell type (e.g. subset on B cells compartment) 

library(ComplexHeatmap)
library(gplots)
library(tidyverse)
library(data.table)

all_tissues= c( "prostate gland", "eye","heart" , 
                "lung" , "skeletal muscle organ", 
                "esophagus","blood" , "brain","kidney" ,
                "liver", "colon", 
                "bone marrow", "spleen","thymus")
all_tissues = as.data.frame(all_tissues)
all_tissues$test= paste(0)

total_cells<- as.data.frame(table(sub$Tissue))
cells1<-table(sub$Tissue, Idents(sub))
cells1[cells1 < 40] <- 0
total_cells<- as.data.frame(table(sub$Tissue))

cells<-table(sub$Tissue, Idents(sub))
cells[cells < 40] <- 0
cells= as.data.frame(cells)
cells = merge(all_tissues, cells, by.x= "all_tissues", by.y = "Var1", all = T)
cells= cells[,-2]
cells[cells < 40] <- 0.1
cells= spread(cells,Var2,Freq)
row.names(cells)=cells$all_tissues
cells= cells[,-1]
cells[is.na(cells)] <- 0
cells= cells[,-ncol(cells)]
cells[cells < 40] <- 0

sign_table<-matrix(0, nrow = nrow(cells), ncol = ncol(cells))
sign_table= as.data.frame(sign_table)


#statistical analysis - hypergeometric distribution
for(c in 1: ncol(cells))
{
  
  for (r in 1: nrow(cells))
  {
    x=1-phyper(cells[r,c]-1, sum(cells[,c]), sum(cells)- sum(cells[,c]), sum(cells[r,])) 
    sign_table[r,c]= paste(as.numeric(x))
    
  }
  sign_table[,c]= as.numeric(sign_table[,c])
}

row.names(sign_table) = row.names(cells)
colnames(sign_table) = colnames(cells)

sign_table[sign_table==0] <- 0.00001
mat= -log10(sign_table)

annotation = as.data.frame(colnames(cells))
names(annotation)= "Name"
annotation = as.data.frame(annotation)
unique(sub$Tissue)



#create annotation and heatmap
column_ha = HeatmapAnnotation(cells_subpopulation = anno_barplot(colSums(cells)), 
                              Clusters = annotation$Name,
                              col = list(Clusters = cols))
row_ha = rowAnnotation(cells_tissue = anno_barplot(rowSums(cells)))


Heatmap(as.matrix(mat), column_names_gp = gpar(fontsize = 10), name= "-log10(pvalue)",
        top_annotation = column_ha, right_annotation = row_ha, show_column_names = F,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(as.matrix(mat)[i, j] > 1.30103) {
            grid.text("*", x, y)
          }
        })



column_ha = HeatmapAnnotation(cells = anno_barplot(colSums(cells)))
row_ha = rowAnnotation(cells = anno_barplot(rowSums(cells)))


Heatmap(as.matrix(mat), column_names_gp = gpar(fontsize = 10), name= "-log10(pvalue)",
        top_annotation = column_ha,right_annotation = row_ha, show_column_names = T, cluster_rows = F,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(as.matrix(mat)[i, j] > 1.30103) {
            grid.text("*", x, y)
          }
        })


-log10(0.05)
