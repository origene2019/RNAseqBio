
library(dplyr)
library(reshape)

#查找路径下所有排序后的bam文件
outfile_path <- '/home/data/Rseq/P101SC1717020113/DEG/alterSplic'
#生成b1.txt, b2.txt
fpkm_dir <- '/home/data/Rseq/P101SC1717020113/'    #sort排序后的bam文件存储路径
fpkmFile <- dir(fpkm_dir,recursive=T,full.names = TRUE,pattern = '_fpkm$')

fpkm_lst <- list()
for (i in 1:length(fpkmFile)) {
  nm <- unlist(strsplit(fpkmFile[i], "[/]"))
  nm <- nm[length(nm)]
  fpkm <- read.delim2(fpkmFile[i])
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









#install.packages("ggpubr")
library("ggpubr")
library("openxlsx")
library(stringr)
library(dplyr)

#fpkmMerge <- read.xlsx('/home/data/Rseq/P101SC1717020113/fpkmMerge.xlsx')

#grpTab <- read.delim('/home/data/Rseq/P101SC1717020113/DEG/grpTab.txt')
#colnames(grpTab) <- c('grp_no', 'treatment', 'control')
#grpTyp <- read.delim('/home/data/Rseq/P101SC1717020113/DEG/grpTyp.txt')

fpkmlt <- melt(fpkmMerge,id.vars="Gene.ID", variable_name ="sample_no")
fpkmlt$sample_no <- str_replace(fpkmlt$sample_no, '_fpkm', '')
fpkmlt$log2fpkm1 <- ifelse(as.numeric(fpkmlt$value)==0, NA, log2(as.numeric(fpkmlt$value) + 1))
fpkmlt <- left_join(fpkmlt, grpTyp, by = c('sample_no' = 'x'))


outplot_dir <- '/home/data/Rseq/P101SC1717020113/DEG/image/fpkmDistrib'
#n <- '组合2'
for (n in grpTab$grp_no) {
  grp <- subset(grpTab, grpTab$grp_no==n)
  
  treat <- subset(fpkmlt, fpkmlt$z == grp$treatment)
  treat$Group <- 'Treat'
  contr <- subset(fpkmlt, fpkmlt$z == grp$control)
  contr$Group <- 'Control'
  fpkdf <- rbind(treat,contr)
  
  plt_nm <- paste0(grp[1,2], '_VS_', grp[1,3])
  p <- ggplot(fpkdf, aes(sample_no, log2fpkm1, fill=Group, colour=Group))+geom_violin(alpha=0.2, width=1)+xlab("") + ylab("log2(FPKM+1)")+labs(title="FPKM distribution", subtitle = plt_nm)+theme(plot.title = element_text(hjust = 0.5))
  #给小提琴图加个箱线图，显示出上、下四分位数与中位数
  p <- p + geom_boxplot(alpha=0.4, width=0.1, outlier.colour=NA) 
  p <- p + theme_light()
  qplot <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5))
  
  ggsave(paste0(outplot_dir, '/', plt_nm , '_FPKM.png'), plot = qplot, limitsize = FALSE, width = 8, height = 8)
  
}

#p<-p+theme_light()
#p<-p+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5))
#p<-p+theme(axis.ticks.x = element_blank())
#p


