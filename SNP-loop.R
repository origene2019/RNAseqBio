
# samtool faidx xx.fa	
# bwa index xx.fa	#已完成 
 
## SNP
grpTyp <- read.delim('/home/data/Rseq/results361111/resultsAll/grpTyp.txt')
path <- '/home/origene/zhengguosong/anaconda3/envs/rna/bin'  #虚拟环境，samtools等软件的安装环境
bam_dir <- '/home/data/Rseq/results361111/resultsAll/sortBam'  #按染色体位置POS排序
outsnp_dir <- '/home/data/Rseq/results361111/resultsAll/SNP'
Ensemble_fa <- '/home/data/Rseq/wangyangxian/genome/Mus_musculus_Ensemble_94.fa'

if(!require(parallel)) install.packages("parallel")
if(!require(stringr)) install.packages("stringr")
library(parallel)
library(stringr)

SNP <- function(bam_fl){
  #bam_fl <- "FAA1BR.bam" 
  input_bam <- paste0(bam_dir, '/', bam_fl, '_sort.bam')
  mrkDupbam <- paste0(outsnp_dir, '/', str_replace(bam_fl,'.bam',''), '_marked_duplicates.bam')
  mrkDupMtr <- paste0(outsnp_dir, '/', str_replace(bam_fl,'.bam',''), '_marked_dup_metrics.txt')
  splitMrkbam <- paste0(outsnp_dir, '/', str_replace(bam_fl,'.bam',''), '_split_marked.bam')
  addSplitMrkbam <- paste0(outsnp_dir, '/', str_replace(bam_fl,'.bam',''), '_add_split_marked.bam')
  vcf <- paste0(outsnp_dir, '/', str_replace(bam_fl,'.bam',''), '.vcf.gz')
  vcf_filter <- paste0(outsnp_dir, '/', str_replace(bam_fl,'.bam',''), '_filter.vcf.gz')
  vcf_filterPass <- paste0(outsnp_dir, '/', str_replace(bam_fl,'.bam',''), '_filter_pass.vcf.gz')
  
  snp_cmd <- paste0(
    # (6.1) 标记比对文件中的pcr重复
    path, "/gatk MarkDuplicates -I ", input_bam, " -O ", mrkDupbam, " -M ", mrkDupMtr, " --ASSUME_SORT_ORDER coordinate --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --VALIDATION_STRINGENCY SILENT --CREATE_INDEX true",
    
    " && ",
    
    path, "/samtools index ", mrkDupbam ,
    
    " && ",
    
    # (6.2) 将落在内含子区间的 reads 片段直接切除，并对 MAPQ 进行调整
    path, "/gatk SplitNCigarReads -I ", mrkDupbam, " -O ", splitMrkbam, " -R ", Ensemble_fa, " --create-output-bam-index true",
    
    " && ",
    
    # (6.3) 将落在内含子区间的 reads 片段直接切除，并对 MAPQ 进行调整
    path, "/gatk AddOrReplaceReadGroups -I ", splitMrkbam, " -O ", addSplitMrkbam, " -RGID 4 -RGLB lib1 -RGPL ILLUMINA -RGPU unitl -RGSM 20",
    
    " && ",
    
    path, "/samtools index ", addSplitMrkbam,
    
    " && ",
    
    # (6.4) 检测snp和indel
    path, "/gatk HaplotypeCaller -I ", addSplitMrkbam, " -O ", vcf, " -R ", Ensemble_fa,
    
    " && ",
    
    # (6.5)snp和indel过滤
    path, "/gatk VariantFiltration -V ", vcf, " --filter-expression \"QUAL < 30.0 || QD < 2.0\" --filter-name \"Filter\" -O ", vcf_filter,
    
    " && ",
    
    path, "/gatk SelectVariants -V ", vcf_filter, " --exclude-filtered true -O ", vcf_filterPass 
  )
  #snp_cmd
  rmats_proc <- system(snp_cmd, intern = FALSE, wait = FALSE)
}


#grpTyp2 <- subset(grpTyp,!grpTyp$x %in% c('FAA1BR','FAA1eye'))
#for (n in 2:ceiling(nrow(grpTyp)/8)) {
#  mc <- floor(getOption("mc.cores", detectCores(logical = F)))
#  res <- mclapply(grpTyp[,1], SNP, mc.preschedule = FALSE, mc.set.seed = TRUE, mc.cleanup = FALSE, mc.cores = mc)
#}


#一次并行8个，8个必须同时执行完才能循环下一组
#system2(intern = FALSE, wait = TRUE, env = character())可以设置执行环境
#Sys.which()可查找可执行程序的路径
k <- 1 #起始任务值
res_no <- length(typ2[,1]) #总任务数
n <- round(detectCores(logical = F)/24*6,0) #同时执行进程数

repeat { 
  ps_lst <- system('ps x | grep filter_pass.vcf.gz', intern = TRUE) #注意这里检索所有的rmats的任务
  thead_n <- length(ps_lst)-1
  addn <- n - thead_n
  if(addn>0){
    res <- SNP(typ2[k,1])
    print(k)
    k <- k+1
    Sys.sleep(sample(30:180,size=1))
  }else{
    Sys.sleep(300)
  }
  
  if(k>res_no) {
    break
  }
}

