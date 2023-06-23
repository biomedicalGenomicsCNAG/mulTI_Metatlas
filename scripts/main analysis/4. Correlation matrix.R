library(corrplot)
library(RColorBrewer)
library(tidyverse)

#Data
# sub= subset of the main Seurat object in one cell type (e.g. subset on B cells compartment)

cluster.averages <- AverageExpression(sub, group.by = c("Tissue"))# average values


cluster.averages <- as.data.frame(cluster.averages[["RNA"]])
cluster.averages$genes = row.names(cluster.averages)
list<-as.data.frame(table(sub$Tissue))
list = list %>% filter(Freq>=40)
cluster.averages = cluster.averages %>% select(c(list$Var1, "genes"))

last= ncol(cluster.averages)
M <-cor(cluster.averages[,-as.numeric(last)])
corrplot(M, type="upper", order="hclust",tl.col = 'black',
         col=brewer.pal(n=8, name="RdYlBu"), is.corr = FALSE,col.lim=c(min(M), max(M)))

testRes = cor.mtest(cluster.averages[,-last], conf.level = 0.95)

corrplot(M, p.mat = testRes$p, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', 
         tl.col = 'black',
         col=brewer.pal(n=8, name="RdYlBu"),
         number.cex = 0.8, order = 'AOE', diag=FALSE, is.corr = FALSE,col.lim=c(min(M), max(M)))
