#cluster1
library(pheatmap)
library(openxlsx)
setwd("/home/data/Rseq/P101SC1717020113/DEG/diffGenes/cluster")
mat=read.table("cluster1_heatmap.csv",header=TRUE,row.names=1,sep=",",check.names = F)
dim(mat)
h=pheatmap(mat,cluster_rows = T,scale = "row",angle_col = 45,clustering_method = "average",legend = T,fontsize = 5,fontsize_row = 6.5,fontsize_col = 10,color = colorRampPalette(rev(c("red","white","blue")))(102))
ggsave("cluster1_heatmap.png", plot=h,device = "png",path="/home/data/Rseq/P101SC1717020113/DEG/diffGenes/cluster")
ggsave("cluster1_heatmap.pdf", plot=h,device = "pdf",path="/home/data/Rseq/P101SC1717020113/DEG/diffGenes/cluster")