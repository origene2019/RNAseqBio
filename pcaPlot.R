#PCA1
library(gmodels)
library(ggplot2)
library(ggrepel)
setwd("/home/data/Rseq/P101SC1717020113/Quant")
expr=read.table("FDHABR_fpkm_VS_FAABR_fpkm.csv",header = T,row.names = 1,sep = ",")
data=t(expr)
dim(data)
data.pca=fast.prcomp(data,scale=F,center=T)
#p=autoplot(data.pca)+theme_classic()+ggtitle("PCA plot")
#print(p)
#plot(data.pca)
a=summary(data.pca)
tmp=a[4]$importance
pro1=as.numeric(sprintf("%.3f",tmp[2,1]))*100
pro2=as.numeric(sprintf("%.3f",tmp[2,1]))*100
pc=as.data.frame(a$x)

pc$group=c(rep("FDHABR_fpkm",4),rep("FAABR_fpkm",4))
pc$names=rownames(pc)

xlab=paste("PC1(",pro1,"%)",sep="")
ylab=paste("PC2(",pro2,"%)",sep="")
pca=ggplot(pc,aes(PC1,PC2))+
  geom_point(size=3,aes(shape=group,color=group)) 
pca2=pca+geom_text_repel(aes(label=names),size=3,vjust="inward",hjust="inward")+labs(x=xlab,y=ylab,title="PCA")
pca2=pca2+geom_hline(yintercept=1,linetype=4,color="grey")+geom_vline(xintercept=1,linetype=4,color="grey")+theme_bw()
ggsave("pca1.png",plot=pca2,device = "png",path="/home/data/Rseq/P101SC1717020113/Quant",width=10,height=8)
ggsave("pca1.pdf",plot=pca2,device = "pdf",path="/home/data/Rseq/P101SC1717020113/Quant",width=10,height=8)