#计算fpkm值
path <- '/home/origene/zhengguosong/anaconda3/envs/rna/bin'  #虚拟环境，htseq等软件的安装环境
infile_dir <- '/home/data/Rseq/wangyangxian36/bam'    #sort排序后的bam文件存储路径
fpkm_dir <- '/home/data/Rseq/results361111/results_all'    #sort排序后的bam文件存储路径
gtf_file <- '/home/data/Rseq/wangyangxian/genome/Mus_musculus.GRCm38.94.gtf'   #gtf文件完整路径

if(!require(stringr)) install.packages("stringr")
library(stringr)

###### 5.2 计算fpkm值
#查找路径下所有排序后的bam文件
bamSort_file <- dir(infile_dir, pattern = '.bam', all.files=T,recursive=T)
bamSort_file

exprLevel <- function(fl_nm){
  sortbam_nm <- str_replace(fl_nm,'-sort.bam|.bam','')
  input_file <- paste0(infile_dir, '/', fl_nm)
  output_fpkm <- paste0(fpkm_dir, '/', sortbam_nm,'_fpkm')
  output_newGtf <- paste0(fpkm_dir, '/', sortbam_nm,'.gtf')
  output_sortFpkm <- paste0(fpkm_dir, '/', sortbam_nm,'_fpkm_sort')
  fpkm_cmd <- paste0(path,'/stringtie ', input_file, ' -G ', gtf_file, ' -A ' , output_fpkm , ' -l novel -o ', output_newGtf , ' && sort ', output_fpkm , ' > ' , output_sortFpkm)
  # -A 输出基因丰度值（FPKM, TPM）
  # -l 新组装的转录本标签（默认STRG）
  # -o 输出新组装的gtf文件
  fpkm_res <- system(fpkm_cmd,intern = TRUE)
}

library(pbmcapply)
mc <- getOption("mc.cores", detectCores(logical = F)-2)
res <- pbmclapply(bamSort_file, exprLevel, mc.preschedule = FALSE, mc.cleanup = FALSE, mc.cores = mc)




