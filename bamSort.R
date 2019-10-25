library(parallel)

path <- '/home/origene/zhengguosong/anaconda3/envs/rna/bin/'  #虚拟环境，samtools等软件的安装环境
source_dir <- '/home/data/Rseq/wangyangxian/'     #存储bam数据文件路径
outfile_dir <- '/home/data/Rseq/自己新建的文件夹/'    #输出文件存储路径(bam-sort文件)

#在路径中查找所有.bam后缀文件
source_file <- dir(source_dir,all.files=T,recursive=T)
bam_file <- sort(source_file[endsWith(source_file,'.bam')])  
bam_file

######  4. 排序(samtools)  对bam文件进行排序，用染色体位置排序！！！
###### 多个bam一起sort最好给个-T，不然sort的临时文件容易串，最后sort的结果可能是错的。不加好像也会自动生成临时文件
bam_sort <- function(sam_file){
  samfl_nm <- str_replace(sam_file,'.bam','')
  infile_nm <- paste0(source_dir, sam_file)
  outfile_nm <- paste0(outfile_dir, samfl_nm, '-sort.bam')
  sortbam_cmd <- paste0(path,'samtools sort -o ', outfile_nm, ' -T ', infile_nm)
  #sortbam_cmd
  
  sortbam_proc <- system(sortbam_cmd, intern = TRUE)
}

#多核并行
mc <- getOption("mc.cores", detectCores(logical = F)-2)   #多核并行，设置核树
res <- mclapply(bam_file, bam_sort, mc.preschedule = FALSE, mc.cleanup = FALSE, mc.cores = mc)












