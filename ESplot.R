# (5.2)可变剪切的可视化
#	rmats2sashimiplot（pyhton2环境）需要调用samtools,使用时特别需要注意把python2环境下的samtools升级到最新版本，否则在调用samtools时会调用失败，
#	原因是自动安装的版本太低，服务器没有它需要的libxxx.so.0.0.0库。
#	升级python2中samtools的方法：
# conda activate python2
# conda update samtools
# conda install -c bioconda samtools openssl=1.0

grpTab <- read.delim2('/home/data/Rseq/results361111/resultsAll/grpTab.txt', stringsAsFactors = FALSE)
colnames(grpTab) <- c('grp_no', 'treatment', 'control')
grpTyp <- read.delim2('/home/data/Rseq/results361111/resultsAll/grpTyp.txt', stringsAsFactors = FALSE)
py2path <- "/home/origene/zhengguosong/anaconda3/envs/py2/bin"
as_dir <- "/home/data/Rseq/results361111/resultsAll/alterSplic"
bam_dir <- "/home/data/Rseq/results361111/resultsAll/sortBam"

if(!require(stringr)) install.packages("stringr")
if(!require(readr)) install.packages("readr")
if(!require(dplyr)) install.packages("dplyr")
library(stringr)
library(readr)
library(dplyr)
SEplot <- function(i){
  bb1 <- subset(grpTyp, grpTyp$z==grpTab[i,2])
  bb2 <- subset(grpTyp, grpTyp$z==grpTab[i,3])
  b2 <- paste0(bam_dir, '/', bb1[,1],'_sort.bam',collapse = ",")
  b1 <- paste0(bam_dir, '/', bb2[,1],'_sort.bam',collapse = ",")
  wd_nm <- paste0('AS_',grpTab[i,2], '_vs_', grpTab[i,3])
  if(!'SEplot' %in% list.files(paste0(as_dir))) system(paste0('mkdir ', as_dir, '/SEplot'))
  grp_nm <- paste0(grpTab[i,2], '_vs_', grpTab[i,3])
  plot_path <- paste0(as_dir, '/SEplot/', grp_nm)
  if(!grp_nm %in% list.files(paste0(as_dir, 'SEplot'))) system(paste0('mkdir ', plot_path))
  for(t in 1:5){
  ESt <- c('A3SS.MATS.JCEC.txt','A5SS.MATS.JCEC.txt','MXE.MATS.JCEC.txt','RI.MATS.JCEC.txt','SE.MATS.JCEC.txt')
  ESt_nm <- c('A3SS','A5SS','MXE','RI','SE')
  nm_path <- paste0(plot_path, '/', ESt_nm[t])
  if(!ESt_nm[t] %in% dir(plot_path)) system(paste0('mkdir ', nm_path))
  jcec <- read.delim2(paste0(as_dir,'/',wd_nm,'/', ESt[t]), stringsAsFactors = FALSE)
  jcec$FDR <- as.numeric(jcec$FDR)
  jcec_sig <- subset(jcec, FDR<0.05) 
  jcec_sig$GeneID <- paste0('"',jcec_sig$GeneID,'"')
  jcec_sig$geneSymbol <- paste0('"',jcec_sig$geneSymbol,'"')
  sig_nm <- paste0(as_dir,'/',wd_nm,'/', ESt[t])
  write_delim(jcec_sig,sig_nm,quote_escape = "none",delim = "\t")
  
  jcec_top20 <- jcec %>%
    subset(FDR<0.05) %>%
    arrange(-FDR) %>%
    mutate(row_number = row_number(-FDR)) %>%
    subset(row_number<=20) %>%
    dplyr::select(-row_number)
  jcec_top20$chr <- as.numeric(str_replace(jcec_top20$chr,'chr',''))
  sign20_nm <- paste0(as_dir,'/',wd_nm,'/', ESt[t])
  jcec_top20$GeneID <- paste0('"',jcec_top20$GeneID,'"')
  jcec_top20$geneSymbol <- paste0('"',jcec_top20$geneSymbol,'"')
  write_delim(jcec_top20,sign20_nm,quote_escape = "none",delim = "\t")
  
  #wd_nm2 <- paste0(grpTab[i,2], '_vs_', grpTab[i,3], '_', ESt_nm[t])
  #if(!wd_nm2 %in% list.files(paste0(as_dir, '/SEplot'))) system(paste0('mkdir ', as_dir, '/SEplot/', nm_path))
  
  SEplot_cmd <- paste0(path, "/rmats2sashimiplot --b1 ", b1, " --b2 ", b2, " -t ", ESt_nm[t]," -e ", sign20_nm, " --l1 ", grpTab[i,2], " --l2 ", grpTab[i,3], " --exon_s 1 --intron_s 5 -o ", nm_path)
  SEplot <- system(paste0("export PATH=$PATH:", py2path," && ", SEplot_cmd),  intern = FALSE, wait = TRUE)
  
  pdf_dir <- paste0(nm_path, "/Sashimi_plot")
  pdf_pics <- dir(pdf_dir, pattern = '.pdf$', full.names = TRUE)
  for(p in pdf_pics){
  png_picnm <- str_replace(p,'.pdf','.png')
  png_cmd <- system(paste0("export PATH=$PATH:", py2path," && magick convert -density 300 ", p, " ", png_picnm), intern = FALSE, wait = FALSE)
  }
  }
}

k <- 1 #起始任务值
res_no <- length(grpTab$grp_no) #总任务数
n <- round(detectCores(logical = F)/24*6,0) #同时执行进程数

repeat { 
  ps_lst <- system('ps x | grep rmats2sashimiplot', intern = TRUE) #注意这里检索所有的rmats的任务
  thead_n <- length(ps_lst)-1
  addn <- n - thead_n
  if(addn>0){
    res <- SEplot(k)
    print(k)
    k <- k+1
    Sys.sleep(sample(10:15,size=1))
  }else{
    Sys.sleep(50)
  }
  
  if(k>res_no) {
    break
  }
}



