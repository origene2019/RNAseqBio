
file_path <- '/home/data/Rseq/results361111/resultsAll/SNP'
f_nm <- dir(file_path,pattern = '_filter_pass.vcf.gz$')
nm <- str_replace(f_nm,'_filter_pass.vcf.gz','')
path <- '/home/origene/lukunpeng/miniconda3/bin'

if(!require(parallel)) install.packages("parallel")
library(parallel)

SNPann <- function(nma){
  SNPann_cmd <- paste0(path, "/snpEff -s ", file_path, '/', nma, ".html GRCm38.86 ", file_path, "/", nma, "_filter_pass.vcf > ", file_path, '/', nma, "_filter_pass_ann.vcf")
  SNPann <- system(SNPann_cmd,  intern = FALSE, wait = FALSE)
}

k <- 1 #起始任务值
res_no <- length(nm) #总任务数
n <- round(detectCores(logical = F)/24*6,0) #同时执行进程数

repeat { 
  ps_lst <- system('ps x | grep filter_pass.vcf.gz', intern = TRUE) #注意这里检索所有的rmats的任务
  thead_n <- length(ps_lst)-1
  addn <- n - thead_n
  if(addn>0){
    res <- SNPann(nm[k])
    print(k)
    k <- k+1
    Sys.sleep(sample(25:80,size=1))
  }else{
    Sys.sleep(300)
  }
  
  if(k>res_no) {
    break
  }
}













