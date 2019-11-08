#BiocManager::install('DESeq2')
library(openxlsx)
library(stringr)
library(dplyr)

options(stringsAsFactors=FALSE)
path <- "/home/data/Rseq/P101SC1717020113/DEG"
setwd(path)
# 分组情况表
grpTab <- read.delim('/home/data/Rseq/P101SC1717020113/DEG/grpTab.txt', stringsAsFactors = FALSE)
#中文字段名：差异比较分析	处理	参考
colnames(grpTab) <- c('grp_no', 'treatment', 'control')
grpTyp <- read.delim('/home/data/Rseq/P101SC1717020113/DEG/grpTyp.txt', stringsAsFactors = FALSE)


#1.载入基因表达量文件，添加列名
counts_dir <- '/home/data/Rseq/P101SC1717020113/'
outcnt_files <- sort(dir(counts_dir, pattern = '_count.txt'))
outcnt_flnms <- str_replace(outcnt_files,'_count.txt','')

countdat_lst <- list()
for (nm in outcnt_flnms) {
  count_data <- read.delim(paste0(counts_dir, nm, '_count.txt'), header=FALSE, stringsAsFactors = FALSE)
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




library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(biomaRt)
library(curl)

for(gp in grpTab$grp_no){
  gp_nm <- paste0(grpTab[which(grpTab$grp_no==gp),2] , '_vs_' , grpTab[which(grpTab$grp_no==gp),3])
  mycounts <- grpMerge_lst[[gp]]
  #这里有个x，需要去除，先把第一列当作行名来处理
  rownames(mycounts)<-mycounts[,1]
  #把带X的列删除
  mycounts<-mycounts[,-1]
  head(mycounts)
  # 这一步很关键，要明白condition这里是因子，不是样本名称；小鼠数据有对照组和处理组，各两个重复
  condition <- factor(c(rep("control",4),rep("treat",4)), levels = c("control","treat"))
  condition
  #colData也可以自己在excel做好另存为.csv格式，再导入即可
  colData <- data.frame(row.names=colnames(mycounts), condition)
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
  ## 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
  diff_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 0)
  #或  diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
  dim(diff_gene_deseq2)
  head(diff_gene_deseq2)
  
  
  #5 用bioMart对差异表达基因进行注释
  mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
  mms_symbols <- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),
                       filters = 'ensembl_gene_id', values = row.names(mycounts), mart = mart)
  
  
  #5.2 合并数据:res结果+mms_symbols合并成一个文件
  head(diff_gene_deseq2)
  head(mms_symbols)
  
  ensembl_gene_id <- row.names(mycounts)
  mycounts <- cbind(ensembl_gene_id, mycounts)
  diff_name <- left_join(mycounts, mms_symbols, by="ensembl_gene_id")
  
  diff_name$sig<-"no"
  #然后设定一个取值范围，这里和我们画的线一致，Log2FC>2或<-2，且padj<0.05的我们认为是上调或者下调
  diff_name$sig[(diff_name$log2FoldChange > 0)&(diff_name$padj < 0.05)] <- "up"
  diff_name$sig[(diff_name$log2FoldChange < 0)&(diff_name$padj < 0.05)] <- "down"
  
  
  #画图  到此为止就完成了RNA-seq的数据处理流程，下一步就是用pheatmap绘制热图了
  
  diff <- as.data.frame(diff_name[,c('ensembl_gene_id','log2FoldChange','padj','sig')])
  diff <- subset(diff, diff$ensembl_gene_id %in% row.names(diff_gene_deseq2))
  
  #思路还是一样的，我们先定义一个向量，里面是我们要显示的基因
  p <- ggplot(diff)+geom_vline(xintercept = 0,linetype="dashed",color="grey")+geom_hline(yintercept = -log10(0.05),linetype="dashed",color="grey")+geom_point(aes(log2FoldChange,-log10(padj),color=sig))+scale_color_manual(values=c("green", "red"))
  p <- p + xlim(-max(abs(diff$log2FoldChange)), max(abs(diff$log2FoldChange)))
  #outcnt_plot <- p+geom_text_repel(aes(log2FoldChange,-log10(padj),label=name))
  p <- p + ggtitle(gp_nm)
  if(!'image' %in% list.files()) system('mkdir image')
  ggsave(paste0('image/', gp_nm , '_DEseq.png'), plot = p, width = 8, height = 8)
  

  diffcount_xls <- right_join(as.data.frame(cbind("ensembl_gene_id" = row.names(mycounts), mycounts)), as.data.frame(diff_name[diff_name$ensembl_gene_id %in% row.names(diff_gene_deseq2)]), by = 'ensembl_gene_id')
  #head(diffcount_xls)
  if(!'diffGenes' %in% list.files()) system('mkdir diffGenes')
  write.xlsx(diffcount_xls,paste0(path, '/diffGenes/', gp_nm, '_diffGenes.xlsx'))
  
  counts_xls <- right_join(as.data.frame(cbind("ensembl_gene_id" = row.names(mycounts), mycounts)), as.data.frame(diff_name), by = 'ensembl_gene_id')
  write.xlsx(counts_xls,paste0(path, '/diffGenes/Cnt/', gp_nm, '_counts.xlsx'))

  }