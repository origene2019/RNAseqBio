library(parallel)

file_path <- '/home/data/Rseq/wangyangxian/SNP/vcf'
f_nm <- dir(file_path,pattern = '_filter_pass.vcf')
nm <- str_replace(f_nm,'_filter_pass.vcf','')
path <- '/home/origene/lukunpeng/miniconda3/bin/snpEff'

SNPann <- function(nma){
  SNPann_cmd <- paste0(path, " -s ", file_path, '/', nma, ".html GRCm38.86 ", file_path, "/", nma, "_filter_pass.vcf > ", file_path, '/', nma, "_filter_pass_ann.vcf")
  SNPann <- system(SNPann_cmd, intern = TRUE)
}

for (i in 1:floor(length(nm)/8)) {
  nma <- nm[(8*i-7):(8*i)]
  mc <- getOption("mc.cores", detectCores(logical = F))-2
  res <- mclapply(nma, SNPann, mc.preschedule = FALSE, mc.set.seed = TRUE, mc.cleanup = FALSE, mc.cores = mc)
}
















