#cluster1

diffFpkm_dir <- "/home/data/Rseq/results361111/resultsAll/Differential"
outfile_dir <- "/home/data/Rseq/results361111/resultsAll/fpkm/cluster"
diffFpkm_file <- dir(diffFpkm_dir, pattern = '_diffGenes_fpkm.xlsx')
diffFpkm_file

if(!require(pheatmap)) install.packages("pheatmap")
if(!require(openxlsx)) install.packages("openxlsx")
if(!require(stringr)) install.packages("stringr")
if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)
library(pheatmap)
library(openxlsx)
library(stringr)

heatmapCluster <- function(diffpkm_nm){
  diffpkm_nm2 <- str_replace(diffpkm_nm, '_diffGenes_fpkm.xlsx', '')
  options(stringsAsFactors = FALSE)
  mat <- read.xlsx(paste0(diffFpkm_dir, '/', diffpkm_nm), rowNames = TRUE)
  for(c in colnames(mat)){
    mat[,c] <- as.numeric(mat[,c])
  }
  dim(mat)
  if(nrow(mat)>1){
  h <- pheatmap(mat,cluster_rows = T,scale = "row",angle_col = 45,clustering_method = "average",
                legend = T,fontsize_col = 6,show_rownames=FALSE,main = diffpkm_nm2,
                color = colorRampPalette(rev(c("red","white","blue")))(dim(mat)[1]*dim(mat)[2]))
  ggsave(paste0(diffpkm_nm2, "_heatmap.png"), plot=h,device = "png",path=outfile_dir)
  ggsave(paste0(diffpkm_nm2, "_heatmap.pdf"), plot=h,device = "pdf",path=outfile_dir)
  }
}

library(pbmcapply)
mc <- getOption("mc.cores", detectCores(logical = F)-2)
res <- pbmclapply(diffFpkm_file, heatmapCluster, mc.preschedule = FALSE, mc.cleanup = FALSE, mc.cores = mc)
