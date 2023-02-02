library(msigdbr)
library(fgsea)
library(data.table)
library(Seurat)
library(ggpubr)
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(readr)


sub = #add the rds file of interest

cells<-table(sub$Tissue)
cells = as.data.frame(cells)
cols = cells$Var1
total<-matrix(0, nrow = 1, ncol = 2)
total = as.data.frame(total)
names(total)= c("features","V2")
totest = cells %>%  filter(Freq>= 40)

i=1

markers <- FindMarkers(sub, ident.1 = totest$Var1[i], group.by = 'Tissue', min.pct = 0.5)
original_gene_list <- markers$avg_log2FC
names(original_gene_list) <- row.names(markers)

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

#Retrieve human H (hallmark) gene set
msigdbr_df <- msigdbr(species = "human", category = "H")
head(msigdbr_df)

# fixing format to work with fgsea
pathwaysH = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

# run fgsea enrichment
fgseaRes <- fgsea(pathways=pathwaysH, gene_list)

gos = fgseaRes %>% filter(padj<0.05) 
gos$tissue = totest$Var1[i]
gos= as.data.frame(gos)
table= gos

#now run in a loop
for ( i in 2:nrow(totest)){
  markers <- FindMarkers(sub, ident.1 = totest$Var1[i], group.by = 'Tissue', min.pct = 0.5)
  
  original_gene_list <- markers$avg_log2FC
  names(original_gene_list) <- row.names(markers)
  
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  # fixing format to work with fgsea
  pathwaysH = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
  
  # run fgsea enrichment
  fgseaRes <- fgsea(pathways=pathwaysH, gene_list)
  
  gos = fgseaRes %>% filter(padj<0.05) 
  gos= as.data.frame(gos)
  
  if (nrow(gos)==0){
    print("only one sub")
  } 
  
  else{
    gos$tissue =totest$Var1[i]
    table = rbind(table, gos)
  } 
  
}

table= as.data.frame(table)
unique(sub$annotation_level2)

ggplot(table, aes(x=pathway, y=ES , label=ES, fill= tissue)) + 
  geom_bar(stat='identity', width=.4,position="dodge", col= "gray66")  +
  coord_flip() +  xlab("GSEA")+
  ylab("-log10(adj.pvalue)") + scale_fill_manual(values=col) + theme_classic() + facet_grid(~ tissue)


# Writing mtcars data to a txt file
table$leadingEdge= as.character(table$leadingEdge)
write_tsv(table, path = "gsea.txt")
