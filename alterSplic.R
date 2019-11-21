
#可变剪切预测
grpTab <- read.delim('/home/data/Rseq/results361111/resultsAll/grpTab.txt')
colnames(grpTab) <- c('grp_no', 'treatment', 'control')
grpTyp <- read.delim('/home/data/Rseq/results361111/resultsAll/grpTyp.txt')

rmats_path <- '/home/origene/zhengguosong/software/rMATS.4.0.2/rMATS-turbo-Linux-UCS4'
gtf_file <- '/home/data/Rseq/wangyangxian/genome/Mus_musculus.GRCm38.94.gtf'

outfile_path <- '/home/data/Rseq/results361111/resultsAll/alterSplic'
#生成b1.txt, b2.txt
bamfile_dir <- '/home/data/Rseq/results361111/resultsAll/sortBam'    #按染色体位置pos，sort排序后的bam文件存储路径
bamSort_file <- dir(bamfile_dir, pattern = '_sort.bam', recursive=T,full.names = TRUE)
bamSort_file

if(!require(dplyr)) install.packages("dplyr")
library(dplyr)

#nn <- grpTab[1,]
alterSplic <- function(gp_no){
  nn <- subset(grpTab, grp_no == gp_no)
  gp_nm <- paste0(nn[,2] , '_vs_' , nn[,3])
  bb1 <- paste0(grpTyp[which(grpTyp$z == nn$control),1],'_sort.bam')
  bb2 <- paste0(grpTyp[which(grpTyp$z == nn$treatment),1],'_sort.bam')
  
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
  
  rmats_cmd <- paste0('python2 ',rmats_path, '/rmats.py --b1 ', outfile_path, '/', gp_nm, '_b1.txt --b2 ', outfile_path, '/', gp_nm, '_b2.txt --gtf ', gtf_file, ' --od ', outfile_path, '/AS_', gp_nm,' -t paired --nthread 5 --readLength 150')
  rmats_cmd
  rmats_proc <- system(rmats_cmd, intern = TRUE, wait = FALSE)
  write.table(rmats_proc,paste0(outfile_path, '/', gp_nm,"_log.txt"))
}

k <- 1 #起始任务值
res_no <- length(grpTab$grp_no) #总任务数
n <- floor(getOption("mc.cores", detectCores(logical = F)-2)/5)  #同时执行进程数

repeat { 
  ps_lst <- system('ps x | grep rmats.py', intern = TRUE) #注意这里检索所有的rmats的任务
  thead_n <- length(ps_lst)-1
  addn <- n - thead_n
  if(addn>0){
    res <- alterSplic(grpTab$grp_no[k])
    print(k)
    k <- k+1
    Sys.sleep(sample(30:120,size=1))
  }else{
    Sys.sleep(300)
  }
  
  if(k>res_no) {
    break
  }
}


