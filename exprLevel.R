library(parallel)

path <- '/home/origene/zhengguosong/anaconda3/envs/rna/bin/'  #虚拟环境，htseq等软件的安装环境
outfile_dir <- '/home/data/Rseq/samtools排序后文件/'    #sort排序后的bam文件存储路径
gtf_file <- '/home/data/Rseq/wangyangxian/genome/Mus_musculus.GRCm38.94.gtf'   #gtf文件完整路径

###### 5. 表达水平分析
#查找路径下所有排序后的bam文件
srcSort_bamFile <- dir(outfile_dir,all.files=T,recursive=T)
bamSort_file <- sort(srcSort_bamFile[endsWith(srcSort_bamFile,'-sort.bam')])
bamSort_file

exprLevel <- function(fl_nm){
  sortbam_nm <- str_replace(fl_nm,'-sort.bam','')
  input_file <- paste0(outfile_dir,fl_nm)
  output_file <- paste0(outfile_dir,sortbam_nm,'_count.txt')
  expr_cmd <- paste0(path,'htseq-count -f bam ', input_file, ' ', gtf_file, ' > ' , output_file)
  #expr_cmd
  system(expr_cmd,intern = TRUE)
}

mc <- getOption("mc.cores", detectCores(logical = F)-2)
res <- mclapply(bamSort_file, exprLevel, mc.preschedule = FALSE, mc.cleanup = FALSE, mc.cores = mc)








