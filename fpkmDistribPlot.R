
# 分组情况表
grpTab <- read.delim2('/home/data/Rseq/results361111/resultsAll/grpTab.txt', stringsAsFactors = FALSE)
#中文字段名：差异比较分析	处理	参考
colnames(grpTab) <- c('grp_no', 'treatment', 'control')
grpTyp <- read.delim2('/home/data/Rseq/results361111/resultsAll/grpTyp.txt', stringsAsFactors = FALSE)

outfile_path <- '/home/data/Rseq/results361111/resultsAll'
fpkm_dir <- '/home/data/Rseq/results361111/resultsAll/fpkm'   
diffgenes_dir <- "/home/data/Rseq/results361111/resultsAll/Differential"
fpkmFile <- dir(fpkm_dir,recursive=T,full.names = TRUE,pattern = '_fpkm$')

if(!require(dplyr)) install.packages("dplyr")
if(!require(stringr)) install.packages("stringr")
if(!require(reshape)) install.packages("reshape")
if(!require(ggplot2)) install.packages("ggplot2")
library(dplyr)
library(stringr)
library(reshape)
library(ggplot2)
options(stringsAsFactors=FALSE)

fpkm_lst <- list()
for (i in 1:length(fpkmFile)) {
  nm <- unlist(strsplit(fpkmFile[i], "[/]"))
  nm <- str_replace(nm[length(nm)], '_fpkm', '')
  fpkm <- read.delim2(fpkmFile[i], stringsAsFactors = FALSE)
  fpkm <- fpkm[!startsWith(fpkm$Gene.ID,'novel'),c("Gene.ID","FPKM")]
  fpkm <- fpkm[!duplicated(fpkm$Gene.ID),]
  fpkm_lst[[nm]] <- fpkm
}

# 多个数据集同时merge
multimerge <- function(dat=list(), by="Gene.ID"){
  if(length(dat)<2)return(as.data.frame(dat))
  mergedat<-dat[[1]]
  names(mergedat)[2] <- names(dat)[1]
  for(i in 2:length(dat)){
    datlst <- dat[[i]]
    names(datlst)[2] <- names(dat)[i]
    mergedat<-merge(mergedat, datlst, by="Gene.ID", all = TRUE)
  }
  return(mergedat)
}

fpkmMerge <- multimerge(fpkm_lst)
library(openxlsx)
write.xlsx(fpkmMerge, paste0(fpkm_dir,'/fpkmMerge.xlsx'))

# 实验组fpkm列表
for (i in 1:nrow(grpTab)) {
  nn <- grpTab[i,]
  cnt_nm <- paste0(nn$treatment, '_vs_', nn$control)
  
  bbc <- subset(grpTyp, grpTyp$z == nn$control)
  grpContr_merge <- multimerge(fpkm_lst[bbc$x])
  ali_nmc <- c('gene_id', paste0('Control_', colnames(grpContr_merge)[-1]))
  colnames(grpContr_merge) <- ali_nmc
  
  bbt <- grpTyp[which(grpTyp$z == nn$treatment),]
  grpTree_merge <- multimerge(fpkm_lst[bbt$x])
  ali_nmt <- c('gene_id', paste0('Treat_', colnames(grpTree_merge)[-1]))
  colnames(grpTree_merge) <- ali_nmt
  
  grp_merge <- merge(grpContr_merge, grpTree_merge, by="gene_id", sort = FALSE)
  #删除前五行
  #grpMerge_lst[[nn$grp_no]] <- grp_merge
  #write.xlsx(grp_merge, paste0(fpkm_dir, '/', cnt_nm, '_fpkm.xlsx'))
  
  diffgenes <- read.xlsx(paste0(diffgenes_dir, '/', cnt_nm, "_diffGenes.xlsx"))
  write.xlsx(grp_merge[which(grp_merge$gene_id %in% diffgenes$ensembl_gene_id),], paste0(diffgenes_dir, '/', cnt_nm, '_diffGenes_fpkm.xlsx'))
}




#library(ggpubr)
library(openxlsx)

fpkmlt <- melt(fpkmMerge,id.vars="Gene.ID", variable_name ="sample_no")
fpkmlt$sample_no <- str_replace(fpkmlt$sample_no, '_fpkm', '')
fpkmlt$log2fpkm1 <- ifelse(as.numeric(as.character(fpkmlt$value))==0, NA, log2(as.numeric(as.character(fpkmlt$value)) + 1))
fpkmlt <- left_join(fpkmlt, grpTyp, by = c('sample_no' = 'x'))

outplot_dir <- '/home/data/Rseq/results361111/resultsAll/fpkm/violinPlot'
#n <- '组合2'
for (n in grpTab$grp_no) {
  grp <- subset(grpTab, grpTab$grp_no==n)
  
  treat <- subset(fpkmlt, fpkmlt$z == grp$treatment)
  treat$Group <- 'Treat'
  contr <- subset(fpkmlt, fpkmlt$z == grp$control)
  contr$Group <- 'Control'
  fpkdf <- rbind(treat,contr)
  
  plt_nm <- paste0(grp[1,2], '_vs_', grp[1,3])
  xsize <- ifelse(max(nchar(fpkdf$sample_no))>=8,60,0)
  p <- ggplot(fpkdf, aes(sample_no, log2fpkm1, fill=Group, colour=Group))+geom_violin(alpha=0.2, width=1)+xlab("") + ylab("log2(FPKM+1)")+labs(title = paste0("FPKM distribution\n", plt_nm))
  #给小提琴图加个箱线图，显示出上、下四分位数与中位数
  p <- p + geom_boxplot(alpha=0.4, width=0.15, outlier.colour=NA) 
  p <- p + theme_light()
  qplot <- p + theme(axis.text.x = element_text(hjust = 1, angle = xsize), panel.grid.major =element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5))
  ggsave(paste0(outplot_dir, '/', plt_nm , '_FPKM.png'), plot = qplot, limitsize = FALSE, width = 8, height = 8)
  ggsave(paste0(outplot_dir, '/', plt_nm , '_FPKM.pdf'), plot = qplot, limitsize = FALSE, width = 8, height = 8)
  print(n)
}



