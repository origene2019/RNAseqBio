
path <- '/home/origene/lukunpeng/miniconda3/bin'  
source_dir <- '/home/data/Rseq/wangyangxian36/bam'     #bam数据文件路径
outfile_dir <- '/home/data/Rseq/results361111/resultsAll/counts'   
gtf_file <- '/home/data/Rseq/wangyangxian/genome/Mus_musculus.GRCm38.94.gtf'   #gtf文件完整路径
bamSort_file <- dir(source_dir, pattern = '.bam$', all.files=T, recursive=T)
bamSort_file
gtf_file <- globalVariables(gtf_file)

if(!require(stringr)) install.packages("stringr")
library(stringr)

###### 5. 表达水平分析
exprLevel <- function(fl_nm){
  gtf_file <- gtf_file
  sortbam_nm <- str_replace(fl_nm,'_sort.bam|.bam','')
  input_file <- paste0(source_dir, '/', fl_nm)
  output_file <- paste0(outfile_dir, '/', sortbam_nm, '_counts.txt')
  expr_cmd <- paste0(path,'/htseq-count -f bam ', input_file, ' ', gtf_file, ' > ' , output_file)
  
  expr_res <- system(expr_cmd,intern = TRUE)
}

library(pbmcapply)
mc <- getOption("mc.cores", detectCores(logical = F)-2)
res <- mclapply (bamSort_file, exprLevel, mc.preschedule = FALSE, mc.cleanup = FALSE, mc.cores = mc)

