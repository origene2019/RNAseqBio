#PCA1
library(gmodels)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(reshape)

pca_dir <- "/home/data/Rseq/results361111/resultsAll/pca"
fpkm_dir <- "/home/data/Rseq/results361111/resultsAll/fpkm"
grpTyp <- read.delim2('/home/data/Rseq/results361111/resultsAll/grpTyp.txt', stringsAsFactors = FALSE)
fpkm_tab <- dir(fpkm_dir, pattern = '_fpkm.xlsx$')
options(stringsAsFactors = FALSE)
fpkm1 <- fpkm_tab[[1]]

pcaPlot <- function(fpkm1){
  expr <- read.xlsx(paste0(fpkm_dir, '/', fpkm1), rowNames = TRUE)
  data <- t(expr)
  rownames <- rownames(data)
  data <- apply(data, 2, as.numeric)
  row.names(data) <- rownames
  dim(data)
  
  data.pca <- fast.prcomp(data,scale=F,center=T)
  a <- summary(data.pca)
  tmp <- a[4]$importance
  pro1 <- as.numeric(sprintf("%.3f",tmp[2,1]))*100
  pro2 <- as.numeric(sprintf("%.3f",tmp[2,1]))*100
  pc <- as.data.frame(a$x)
  x <- str_replace(rownames(data), 'Control_|Treat_', '')
  x <- left_join(as.data.frame(x),grpTyp,by='x')
  pc$group <- x$z
  pc$names <- rownames
  
  nm <- paste0(unique(x$z)[2], '_vs_', unique(x$z)[1])
  xlab <- paste("PC1(",pro1,"%)",sep="")
  ylab <- paste("PC2(",pro2,"%)",sep="")
  pca <- ggplot(pc,aes(PC1,PC2)) + geom_point(size=3,aes(shape=group,color=group)) 
  pca <- pca+geom_text_repel(aes(label=names),size=3,vjust="inward",hjust="inward")+labs(x=xlab,y=ylab,title="PCA",subtitle = nm)
  pca2 <- pca+geom_hline(yintercept=1,linetype=4,color="grey")+geom_vline(xintercept=1,linetype=4,color="grey")+theme_bw()
  ggsave(paste0(nm, "_pca.png"),plot=pca2,device = "png",path=pca_dir,width=10,height=8)
  ggsave(paste0(nm, "_pca.pdf"),plot=pca2,device = "pdf",path=pca_dir,width=10,height=8)
}

library(pbmcapply)
mc <- getOption("mc.cores", detectCores(logical = F)-2)
res <- pbmclapply(fpkm_tab, pcaPlot, mc.preschedule = FALSE, mc.cleanup = FALSE, mc.cores = mc)



