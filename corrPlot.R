
if(!require(reshape2)) install.packages("reshape2")
if(!require(ggplot2)) install.packages("ggplot2")
library(reshape2)
library(ggplot2)

corr_plot <- function(cor, title_nm){
  melted_cor <- melt(round(cor,3))
  p <- ggplot(data = melted_cor, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "white")
  p <- p + scale_fill_gradient2(low = "white", high = "blue",  midpoint = floor(min(melted_cor$value)*10)/10, limit = c(floor(min(melted_cor$value)*10)/10,1), space = "Lab", name="Pearson\nCorrelation")
  p <- p + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))
  p <- p + geom_text(aes(Var1, Var2, label = value), color = "black", size = 4) + labs(title = title_nm)
  corr_p <- p + theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = 'right')
  return(corr_p)
}


if(!require(openxlsx)) install.packages("openxlsx")
if(!require(stringr)) install.packages("stringr")
if(!require(dplyr)) install.packages("dplyr")
library(openxlsx)
library(stringr)
library(dplyr)

options(stringsAsFactors=FALSE)
grpTab <- read.delim('/home/data/Rseq/results361111/resultsAll/grpTab.txt', stringsAsFactors = FALSE)
colnames(grpTab) <- c('grp_no', 'treatment', 'control')
grpTyp <- read.delim('/home/data/Rseq/results361111/resultsAll/grpTyp.txt', stringsAsFactors = FALSE)

#1.载入基因表达量文件，添加列名
counts_dir <- '/home/data/Rseq/results361111/resultsAll/counts'
outcnt_files <- sort(dir(counts_dir, pattern = '_counts.txt'))
outcnt_flnms <- str_replace(outcnt_files,'_counts.txt','')

countdat_lst <- list()
for (nm in outcnt_flnms) {
  count_data <- read.delim(paste0(counts_dir, '/', nm, '_counts.txt'), header=FALSE, stringsAsFactors = FALSE)
  countdat_lst[[nm]] <- count_data
}

# 多个数据集同时merge
multimerge <- function(dat=list(), by="V1"){
  if(length(dat)<2)return(as.data.frame(dat))
  mergedat<-dat[[1]]
  names(mergedat)[2] <- names(dat)[1]
  for(i in 2:length(dat)){
    datlst <- dat[[i]]
    names(datlst)[2] <- names(dat)[i]
    mergedat<-merge(mergedat, datlst, by="V1")
  }
  return(mergedat)
}


grpMerge_lst <- list()
for (i in 1:nrow(grpTab)) {
  nn <- grpTab[i,]
  
  bbc <- subset(grpTyp, grpTyp$z == nn$control)
  grpContr_merge <- multimerge(countdat_lst[bbc$x])
  ali_nmc <- c('gene_id', paste0('Control_', colnames(grpContr_merge)[-1]))
  colnames(grpContr_merge) <- ali_nmc
  
  bbt <- grpTyp[which(grpTyp$z == nn$treatment),]
  grpTree_merge <- multimerge(countdat_lst[bbt$x])
  ali_nmt <- c('gene_id', paste0('Treat_', colnames(grpTree_merge)[-1]))
  colnames(grpTree_merge) <- ali_nmt
  
  grp_merge <- merge(grpContr_merge, grpTree_merge, by="gene_id", sort = FALSE)
  #删除前五行
  grpMerge_lst[[nn$grp_no]] <- grp_merge[-1:-5,]
}

# 相关矩阵画图函数
corr_plot <- function(cor, title_nm){
  melted_cor <- melt(round(cor,3))
  p <- ggplot(data = melted_cor, aes(x=melted_cor[,1], y=melted_cor[,2], fill=value)) + geom_tile(color = "white")
  p <- p + scale_fill_gradient2(low = "white", high = "blue",  midpoint = round(min(melted_cor$value)*10)/10, limit = c(round(min(melted_cor$value)*10)/10,1), space = "Lab", name="Pearson\nCorrelation")
  p <- p + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))
  p <- p + geom_text(aes(melted_cor[,1], melted_cor[,2], label = value), color = "black", size = 4) + labs(title = title_nm)
  corr_p <- p + theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = 'right')
  return(corr_p)
}

outfile_dir <- "/home/data/Rseq/results361111/resultsAll/corrPlot"
sampl_corr <- function(gp){
  gp_nm <- paste0(grpTab[which(grpTab$grp_no==gp),2] , '_vs_' , grpTab[which(grpTab$grp_no==gp),3],'_corr')
  mycounts <- grpMerge_lst[[gp]]
  cor <- cor(mycounts[,-1], method="spearman")
  corr_p <- corr_plot(cor,gp_nm)
  ggsave(paste0(outfile_dir, '/', gp_nm , '.png'), plot = corr_p, limitsize = FALSE, width = 8, height = 8)
  
}

library(pbmcapply)
mc <- getOption("mc.cores", detectCores(logical = F)-2)
res <- pbmclapply(grpTab$grp_no, sampl_corr, mc.preschedule = FALSE, mc.cleanup = FALSE, mc.cores = mc)


