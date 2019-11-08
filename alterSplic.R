#install.packages('pbmcapply')
library(parallel)
#可变剪切预测
#rmats.py --b1 b1.txt --b2 b2.txt --gtf /home/data/Rseq/wangyangxian/genome/Mus_musculus.GRCm38.94.gtf --od AS -t paired --nthread 5 --readLength 150
#--s1 s1.txt  txt输入文件，包含了Sample1的fastq序列文件
#--s2 s2.txt  txt输入文件，包含了Sample2的fastq序列文件
#--b1 b1.txt  txt输入文件，包含了Sample1的bam序列文件
#--b2 b2.txt  txt输入文件，包含了Sample2的bam序列文件
#--t  readType RNASeq数据类型，'paired'表示双端数据，'single'表示单端数据
#--readLength <int> 测序read长度
#--gtf gtfFile gtf注释文件
#--od  outDir  输出路径

# 分组情况表
grpTab <- read.delim('/home/data/Rseq/P101SC1717020113/DEG/grpTab.txt')
#中文字段名：差异比较分析	处理	参考
colnames(grpTab) <- c('grp_no', 'treatment', 'control')
grpTyp <- read.delim('/home/data/Rseq/P101SC1717020113/DEG/grpTyp.txt')


library(parallel)
#查找路径下所有排序后的bam文件
outfile_path <- '/home/data/Rseq/P101SC1717020113/DEG/alterSplic'
#生成b1.txt, b2.txt
bamfile_dir <- '/home/data/Rseq/P101SC1717020113/'    #sort排序后的bam文件存储路径
srcSort_bamFile <- dir(bamfile_dir,recursive=T,full.names = TRUE)
bamSort_file <- sort(srcSort_bamFile[endsWith(srcSort_bamFile,'-sort.bam')])
bamSort_file

#nn <- grpTab[1,]
alterSplic <- function(gp_no){
  nn <- subset(grpTab, grp_no == gp_no)
  gp_nm <- paste0(nn[,2] , '_vs_' , nn[,3])
  bb1 <- paste0(grpTyp[which(grpTyp$z == nn$control),1],'-sort.bam')
  bb2 <- paste0(grpTyp[which(grpTyp$z == nn$treatment),1],'-sort.bam')
  
  b1 <- list()
  for (b in bb1) {
    b1[b] <- bamSort_file[grep(b,bamSort_file)]
  }
  
  b2 <- list()
  for (b in bb2) {
    b2[b] <- bamSort_file[grep(b,bamSort_file)]
  }
  
  write.table(b1, paste0(outfile_path, "/", gp_nm, "_b1.txt"), row.names = FALSE,col.names = FALSE,quote = FALSE,sep = ",")
  write.table(b2, paste0(outfile_path, "/", gp_nm, "_b2.txt"), row.names = FALSE,col.names = FALSE,quote = FALSE,sep = ",")
  
  rmats_path <- '/home/origene/zhengguosong/software/rMATS.4.0.2/rMATS-turbo-Linux-UCS4'
  gtf_file <- '/home/data/Rseq/wangyangxian/genome/Mus_musculus.GRCm38.94.gtf'
  
  rmats_cmd <- paste0('python2 ',rmats_path, '/rmats.py --b1 ', outfile_path, '/', gp_nm, '_b1.txt --b2 ', outfile_path, '/', gp_nm, '_b2.txt --gtf ', gtf_file, ' --od ', outfile_path, '/AS_', gp_nm,' -t paired --nthread 5 --readLength 150')
  rmats_cmd
  rmats_proc <- system(rmats_cmd,intern = TRUE)
  write.table(rmats_proc,paste0(outfile_path, '/', gp_nm,"_log.txt"))
}

mc <- floor(getOption("mc.cores", detectCores(logical = F)-2)/5)
res <- mclapply(grpTab$grp_no, alterSplic, mc.preschedule = FALSE, mc.cleanup = FALSE, mc.cores = mc)








