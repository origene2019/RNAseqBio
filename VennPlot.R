#VENN1

outfile_dir <- "/home/data/Rseq/results361111/resultsAll/venn"
diffgenes_dir <- "/home/data/Rseq/results361111/resultsAll/Differential"
grpVenn <- read.delim2('/home/data/Rseq/results361111/resultsAll/grpVenn.txt')
grpTab <- read.delim2('/home/data/Rseq/results361111/resultsAll/grpTab.txt', stringsAsFactors = FALSE)
#中文字段名：差异比较分析	处理	参考
colnames(grpTab) <- c('grp_no', 'treatment', 'control')

if(!require(VennDetail)) BiocManager::install("VennDetail")
if(!require(VennDiagram)) BiocManager::install("VennDiagram")
if(!require(openxlsx)) install.packages("openxlsx")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(purrr)) install.packages("purrr")
if(!require(dplyr)) install.packages("dplyr")
library(VennDetail)
library(VennDiagram)
library(openxlsx)
library(ggplot2)
library(purrr)
library(dplyr)

vennPlot <- function(i){
  venn_nodf <- as.data.frame(na.omit(t(grpVenn[i,-1])))
  venn_nodf <- venn_nodf[which(!venn_nodf[,1]==''),]
  venn_nodf <- grpTab[which(grpTab$grp_no %in% venn_nodf),]
  nm <- paste0(venn_nodf$treatment, '_vs_', venn_nodf$control)
  input <- list()
  dat_list <- list()
  for (gn in nm) {
    diffdat <- read.xlsx(paste0(diffgenes_dir, '/', gn, '_diffGenes.xlsx'))
    input[[gn]] <- diffdat$ensembl_gene_id
    dat_list[[gn]] <- diffdat
  }
  
  res <- venndetail(input, sep = "~")
  summary(res) #show overlap 'details' of all subsets
  result <- VennDetail::result(res, wide = TRUE)
  head(result)
  
  multimerge <- function(dat=list()){
    if(length(dat)<2){
      dat1 <- as.data.frame(dat[[1]])
      names(dat1)[-1] <- paste0(names(dat)[1], '-', str_replace(names(dat1)[-1],"Control_|Treat_",""))
      return(dat1)
    }
    mergedat<-dat[[1]]
    names(mergedat)[-1] <- paste0(names(dat)[1], '-', str_replace(names(mergedat)[-1],"Control_|Treat_",""))
    for(i in 2:length(dat)){
      datlst <- dat[[i]]
      names(datlst)[-1] <- paste0(names(dat)[i], '-', str_replace(names(datlst)[-1],"Control_|Treat_",""))
      mergedat<-merge(mergedat, datlst, by="ensembl_gene_id", all = TRUE)
    }
    return(mergedat)
  }
  
  venn_grp <- unique(res$Subset)
  for (k in 1:length(venn_grp)) {
    grp_res <- getSet(res, venn_grp[k]) # get unique elements in A
    sub_nm <- str_split(venn_grp[k],"~")
    sub_lst <- dat_list[sub_nm[[1]]]
    venn_dat <- multimerge(sub_lst)
    venn_dat <- subset(venn_dat, venn_dat[,1] %in% grp_res$Detail)
    if(!paste0('Venn', i) %in% dir(outfile_dir)) system(paste0('mkdir ', outfile_dir, '/Venn',i))
    write.xlsx(venn_dat,paste0(outfile_dir,'/Venn', i, '/', nrow(grp_res), '_', venn_grp[k], '.xlsx'))
  }
  
  
  #绘图
  p <- venn.diagram(input,filename = NULL,margin =0.1,fill=c("#76E0D5","#F26875","#E5FF8D","#4093A4","#8B8B00","#9ACD32")[1:length(nm)],cat.cex = 0.7)
  ggsave(paste0("Venn",i,".png"), plot=p, device = "png",path=outfile_dir)
  ggsave(paste0("Venn",i,".pdf"), plot=p, device = "pdf",path=outfile_dir)
  
}

library(pbmcapply)
mc <- getOption("mc.cores", detectCores(logical = F)-2)
res <- pbmclapply(1:nrow(grpVenn), vennPlot, mc.preschedule = FALSE, mc.cleanup = FALSE, mc.cores = mc)
