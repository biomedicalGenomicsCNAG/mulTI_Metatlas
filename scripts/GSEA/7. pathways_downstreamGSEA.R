# frequency of pathways
library(data.table)
library(Seurat)
library(ggpubr)
library(tidyverse)


t1= fread("B cells_gsea.txt")
t1$Subpopulation = paste("B cells")
t2= fread("myeloid cells_gsea.txt")
t2$Subpopulation = paste("myeloid cells")
t_all= rbind(t1,t2)
t4= fread("progenitor cells_gsea.txt")
t4$Subpopulation = paste("progenitor cells")
t5= fread("T-NK cells_gsea.txt")
t5$Subpopulation = paste("T-NK cells")
t_all= rbind(t_all,t3)
t_all= rbind(t_all,t4)
t_all= rbind(t_all,t5)

colnames(t_all)
sum_tiss<-t_all %>% 
  dplyr::group_by(pathway,tissue) %>% 
  dplyr::summarise(count_observedWITHINTissue = n()) 

sum_table1<-sum_tiss %>% 
  dplyr::group_by(pathway) %>% 
  dplyr::summarise(count_alltissues = n())


sum_sub<-t_all %>% 
  dplyr::group_by(pathway,Subpopulation) %>% 
  dplyr::summarise(count_withinTissue = n()) 

sum_table2<-sum_sub %>% 
  dplyr::group_by(pathway) %>% 
  dplyr::summarise(count_allsub = n())

t_all$barcode = paste(t_all$pathway, t_all$Subpopulation)
sum_sub$barcode = paste(sum_sub$pathway, sum_sub$Subpopulation)
t_all_full= merge(t_all, sum_sub, by="barcode",all.x=T)
t_all_full= merge(t_all_full, sum_table2, by.x="pathway.x",by.y="pathway",all.x=T)

t_all_full_uniqueSubp = t_all_full %>% filter(count_allsub==1) %>% filter(count_withinTissue==1)
t_all_full_morethanonetissue = t_all_full %>% filter(count_withinTissue>3)


####col tissues
palette = c("#5A5156", "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32",
            "#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C",
            "#8A5F00", "#DEA0FD", "#AA0DFE", "#F8A19F")

col = setNames(palette, 
               c( "prostate gland", "eye","heart" , 
                  "lung" , "skeletal muscle organ", 
                  "esophagus","blood" , "brain","kidney" ,
                  "liver", "colon", 
                  "bone marrow", "spleen","thymus"))

# plot in more than 1 tissue
t_all_full_morethanonetissue$count = paste(1)
ggplot(t_all_full_morethanonetissue, aes(fill=tissue, y=ES, x=pathway.x)) + 
  geom_bar(stat="identity",position = "dodge",width = 0.7, color = "grey40") + scale_fill_manual(values=col) +
  labs( x ="Regulons", y = "ES", fill = "Tissue/organ") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + coord_flip()+ theme_classic() +
  facet_grid(~Subpopulation.x)

# plot in unique detection in one tissue
t_all_full_uniqueSubp$count = paste(1)
ggplot(t_all_full_uniqueSubp, aes(fill=tissue, y=ES, x=pathway.x)) + 
  geom_bar(stat="identity",width = 0.7, color = "grey40") + scale_fill_manual(values=col) +
  labs( x ="Regulons", y = "Count", fill = "Tissue/organ") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + coord_flip()+ theme_classic() +
  facet_grid(~Subpopulation.x) 

