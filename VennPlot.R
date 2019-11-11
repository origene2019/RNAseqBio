#VENN1
library(VennDiagram)
library(openxlsx)
library(ggplot2)
getwd()
setwd("/home/data/Rseq/P101SC1717020113/DEG/diffGenes")
MDHABR_vs_MAABR=read.xlsx("MDHABR_vs_MAABR_diffGenes.xlsx",colNames = TRUE)
FDHABR_vs_FAABR=read.xlsx("FDHABR_vs_FAABR_diffGenes.xlsx",colNames = TRUE)
MDHABR_vs_FDHABR=read.xlsx("MDHABR_vs_FDHABR_diffGenes.xlsx",colNames = TRUE)
MAABR_vs_FAABR=read.xlsx("MAABR_vs_FAABR_diffGenes.xlsx",colNames = TRUE)

#生成list
input=list(MDHABR_vs_MAABR=MDHABR_vs_MAABR$ensembl_gene_id,FDHABR_vs_FAABR=FDHABR_vs_FAABR$ensembl_gene_id,MDHABR_vs_FDHABR=MDHABR_vs_FDHABR$ensembl_gene_id,MAABR_vs_FAABR=MAABR_vs_FAABR$ensembl_gene_id)
#绘图
p<-venn.diagram(input,filename = NULL,margin =0.1,fill=c("#76E0D5","#F26875","#E5FF8D","#4093A4"),cat.cex = 0.7)
ggsave("Venn1.png", plot=p,device = "png",path="/home/data/Rseq/P101SC1717020113/DEG/diffGenes/venn")
ggsave("Venn1.pdf", plot=p,device = "pdf",path="/home/data/Rseq/P101SC1717020113/DEG/diffGenes/venn")