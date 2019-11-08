#conda 添加频道
#conda config --add channels bioconda
#如何在R中读取linux进程状态（ps）命令的输出
#library(stringr) # has a convenient function for splitting to a fixed length 
#
#raw          <- system("ps aux", intern = TRUE)
#fields       <- strsplit(raw[1], " +")[[1]]
#ps           <- str_split_fixed(raw[-1], " +", n = length(fields))
#colnames(ps) <- fields
#install.packages('pbmcapply')

## SNP
library(parallel)

grpTyp
#grpTyp2 <- subset(grpTyp,!grpTyp$x %in% c('FAA1BR','FAA1eye'))
i <- round(detectCores(logical = F)/24*6,0)

for (n in 2:ceiling(nrow(grpTyp)/8)) {
  
  #n <- 1
  sta <- 1+(n-1)*i
  end <- min(i+(n-1)*i, nrow(grpTyp2))
  grpTyp_2to <- grpTyp2[sta : end, 1]
  
  # samtool faidx xx.fa	
  # bwa index xx.fa	#已完成 
  path <- '/home/origene/lukunpeng/miniconda3/bin/'  #虚拟环境，samtools等软件的安装环境
  bamfl_dir <- '/home/data/Rseq/wangyangxian'
  
  SNP <- function(bam_fl){
    #bam_fl <- "FAA1BR.bam" 
    input_bam <- paste0(bamfl_dir, '/', bam_fl, '.bam')
    mrkDupbam <- paste0(bamfl_dir, '/SNP/', str_replace(bam_fl,'.bam',''), '_marked_duplicates.bam')
    mrkDupMtr <- paste0(bamfl_dir, '/SNP/', str_replace(bam_fl,'.bam',''), '_marked_dup_metrics.txt')
    splitMrkbam <- paste0(bamfl_dir, '/SNP/', str_replace(bam_fl,'.bam',''), '_split_marked.bam')
    Ensemble_fa <- '/home/data/Rseq/wangyangxian/genome/Mus_musculus_Ensemble_94.fa'
    addSplitMrkbam <- paste0(bamfl_dir, '/SNP/', str_replace(bam_fl,'.bam',''), '_add_split_marked.bam')
    vcf <- paste0(bamfl_dir, '/SNP/', str_replace(bam_fl,'.bam',''), '.vcf.gz')
    vcf_filter <- paste0(bamfl_dir, '/SNP/', str_replace(bam_fl,'.bam',''), '_filter.vcf.gz')
    vcf_filterPass <- paste0(bamfl_dir, '/SNP/', str_replace(bam_fl,'.bam',''), '_filter_pass.vcf.gz')
    
    snp_cmd <- paste0(
      # (6.1) 标记比对文件中的pcr重复
      path, "gatk MarkDuplicates -I ", input_bam, " -O ", mrkDupbam, " -M ", mrkDupMtr, " --ASSUME_SORT_ORDER coordinate --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --VALIDATION_STRINGENCY SILENT --CREATE_INDEX true",
      
      " && ",
      
      path, "samtools index ", mrkDupbam ,
      
      " && ",
      
      # (6.2) 将落在内含子区间的 reads 片段直接切除，并对 MAPQ 进行调整
      path, "gatk SplitNCigarReads -I ", mrkDupbam, " -O ", splitMrkbam, " -R ", Ensemble_fa, " --create-output-bam-index true",
      
      " && ",
      
      # (6.3) 将落在内含子区间的 reads 片段直接切除，并对 MAPQ 进行调整
      path, "gatk AddOrReplaceReadGroups -I ", splitMrkbam, " -O ", addSplitMrkbam, " -RGID 4 -RGLB lib1 -RGPL ILLUMINA -RGPU unitl -RGSM 20",
      
      " && ",
      
      path, "samtools index ", addSplitMrkbam,
      
      " && ",
      
      # (6.4) 检测snp和indel
      path, "gatk HaplotypeCaller -I ", addSplitMrkbam, " -O ", vcf, " -R ", Ensemble_fa,
      
      " && ",
      
      # (6.5)snp和indel过滤
      path, "gatk VariantFiltration -V ", vcf, " --filter-expression \"QUAL < 30.0 || QD < 2.0\" --filter-name \"Filter\" -O ", vcf_filter,
      
      " && ",
      
      path, "gatk SelectVariants -V ", vcf_filter, " --exclude-filtered true -O ", vcf_filterPass 
    )
    #snp_cmd
    rmats_proc <- system(snp_cmd, intern = TRUE)
  }
  
  mc <- floor(getOption("mc.cores", detectCores(logical = F)))
  res <- mclapply(grpTyp_2to, SNP, mc.preschedule = FALSE, mc.set.seed = TRUE, mc.cleanup = FALSE, mc.cores = mc)
  
}



#一次并行8个，8个必须同时执行完才能循环下一组
#system2(intern = FALSE, wait = TRUE, env = character())可以设置执行环境
#Sys.which()可查找可执行程序的路径


