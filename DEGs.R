
if(!require(openxlsx)) install.packages("openxlsx")
if(!require(stringr)) install.packages("stringr")
if(!require(dplyr)) install.packages("dplyr")
library(openxlsx)
library(stringr)
library(dplyr)

options(stringsAsFactors=FALSE)
# 分组情况表
grpTab <- read.delim('/home/data/Rseq/results361111/resultsAll/grpTab.txt', stringsAsFactors = FALSE)
#中文字段名：差异比较分析	处理	参考
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


if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
library(BiocManager)
if(!require(DESeq2)) BiocManager::install('DESeq2')
if(!require(biomaRt)) BiocManager::install('biomaRt')
library(DESeq2)
library(biomaRt)
#library(curl)
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(ggrepel)) install.packages("ggrepel")
library(ggplot2)
library(ggrepel)


outplot_dir <- "/home/data/Rseq/results361111/resultsAll/Differential/volcano"
diffcounts_dir <- "/home/data/Rseq/results361111/resultsAll/Differential"
outcounts_dir <- diffcounts_dir

volcanoPlot <- function(gp){
  gp_nm <- paste0(grpTab[which(grpTab$grp_no==gp),2] , '_vs_' , grpTab[which(grpTab$grp_no==gp),3])
  mycounts <- grpMerge_lst[[gp]]
  #这里有个x，需要去除，先把第一列当作行名来处理
  rownames(mycounts)<-mycounts[,1]
  #把带X的列删除
  mycounts<-mycounts[,-1]
  head(mycounts)
  # 这一步很关键，要明白condition这里是因子，不是样本名称；小鼠数据有对照组和处理组，各两个重复
  tn <- nrow(grpTyp[which(grpTyp$z == grpTab[which(grpTab$grp_no==gp),2]),])
  cn <- nrow(grpTyp[which(grpTyp$z == grpTab[which(grpTab$grp_no==gp),3]),])
  condition <- factor(c(rep("control", cn),rep("treat", tn)), levels = c("control","treat"))
  condition
  #colData也可以自己在excel做好另存为.csv格式，再导入即可
  colData <- data.frame(row.names = colnames(mycounts), condition)
  colData
  
  #2 构建dds对象,开始DESeq流程    注释：dds=DESeqDataSet Object
  dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
  dds <- DESeq(dds)
  dds
  
  #3 总体结果查看
  #res <- results(dds, contrast=c("condition", "control", "treat"))或下面命令
  res <- results(dds)
  res <- res[order(res$pvalue),]
  head(res)
  summary(res)
  #所有结果先进行输出
  table(res$padj<0.05)
  
  #4 提取差异表达genes（DEGs）并进行gene symbol注释
  ## 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于0或者小于-0的差异表达基因。
  #5 用bioMart对差异表达基因进行注释
  mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
  mms_symbols <- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),
                       filters = 'ensembl_gene_id', values = row.names(mycounts), mart = mart)
  head(mms_symbols)

  mycounts <- cbind(ensembl_gene_id=row.names(mycounts), mycounts)
  res <- as.data.frame(res)
  res$ensembl_gene_id <- row.names(res)
  diff_name <- left_join(left_join(mycounts, res, by="ensembl_gene_id"), mms_symbols, by="ensembl_gene_id")
  
  diff_name$sig<-"no"
  diff_name$sig[(diff_name$log2FoldChange > 0)&(diff_name$padj < 0.05)] <- "up"
  diff_name$sig[(diff_name$log2FoldChange < 0)&(diff_name$padj < 0.05)] <- "down"
  
  #画图  到此为止就完成了RNA-seq的数据处理流程，下一步就是用pheatmap绘制热图了
  diff <- as.data.frame(diff_name[,c('ensembl_gene_id','log2FoldChange','padj','sig')])
  diff <- subset(diff, diff$sig %in% c('up', 'down'))
  
  #思路还是一样的，我们先定义一个向量，里面是我们要显示的基因
  if(unique(diff$sig)=="up"){
    color <- "red"
  }else if(unique(diff$sig)=="down"){
      color <- "green"
  } else{
    color <- c("green", "red")
  }
  p <- ggplot(diff)+geom_vline(xintercept = 0,linetype="dashed",color="grey") + geom_hline(yintercept = -log10(0.05),linetype="dashed",color="grey")
  p <- p + geom_point(aes(log2FoldChange,-log10(padj),color=sig))+scale_color_manual(values=color)
  p <- p + xlim(-max(abs(diff$log2FoldChange)), max(abs(diff$log2FoldChange))) + theme_light()
  p <- p + ggtitle(gp_nm) + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5))
  ggsave(paste0(outplot_dir, '/', gp_nm , '_DEseq.png'), plot = p, width = 8, height = 8)
  ggsave(paste0(outplot_dir, '/', gp_nm , '_DEseq.pdf'), plot = p, width = 8, height = 8)
  #差异基因列表
  write.xlsx(subset(diff_name, diff_name$sig %in% c('up', 'down')),paste0(diffcounts_dir, '/', gp_nm, '_diffGenes.xlsx'))
  #write.xlsx(subset(diff_name, diff_name$sig=='up'),paste0(diffcounts_dir, '/', gp_nm, '_diffGenesUp.xlsx'))
  #write.xlsx(subset(diff_name, diff_name$sig=='down'),paste0(diffcounts_dir, '/', gp_nm, '_diffGenesDown.xlsx'))
  #read counts列表
  write.xlsx(diff_name, paste0(outcounts_dir, '/', gp_nm, '_Genes.xlsx'))
}

library(pbmcapply)
mc <- getOption("mc.cores", detectCores(logical = F)-2)
res <- pbmclapply(grpTab$grp_no, volcanoPlot, mc.preschedule = FALSE, mc.cleanup = FALSE, mc.cores = mc)

