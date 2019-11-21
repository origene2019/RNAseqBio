if(!require(reshape2)) install.packages("reshape2")
if(!require(ggplot2)) install.packages("ggplot2")
library(reshape2)
library(ggplot2)

corr_plot <- function(cor, title_nm){
  melted_cor <- melt(round(cor,3))
  p <- ggplot(data = melted_cor, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "white")
  p <- p + scale_fill_gradient2(low = "white", high = "blue",  midpoint = floor(min(melted_cor$value)*10)/10, limit = c(floor(min(melted_cor$value)*10)/10,1), space = "Lab", name="Pearson\nCorrelation")
  p <- p + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))
  p <- p + geom_text(aes(Var1, Var2, label = value), color = "black", size = 4) + labs(title = title_nm)
  corr_p <- p + theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = 'right')
  return(corr_p)
}

sampl_corr <- function(gp){
  gp_nm <- paste0(grpTab[which(grpTab$grp_no==gp),2] , '_vs_' , grpTab[which(grpTab$grp_no==gp),3],'_corr')
  mycounts <- grpMerge_lst[[gp]]
  cor <- cor(mycounts[,-1], method="spearman")
  corr_p <- corr_plot(cor,gp_nm)
  ggsave(paste0('image/sample_corr/', gp_nm , '.png'), plot = corr_p, limitsize = FALSE, width = 8, height = 8)
  
}

#mc <- getOption("mc.cores", detectCores(logical = F)-2)
mc <- getOption("mc.cores", 2)
res <- mclapply(grpTab$grp_no, sampl_corr, mc.preschedule = FALSE, mc.cleanup = FALSE, mc.cores = mc)


#for(n in grpTab$grp_no){
#  Tab <- subset(grpTab, grpTab$grp_no==n)
#  nm <- paste0('/home/data/Rseq/P101SC1717020113/DEG/image/sample_corr/countData/', Tab[1,2], '_vs_', Tab[1,3], '_counts.xlsx')
#  counts_dat <- grpMerge_lst[[n]]
#  write.xlsx(counts_dat, nm)

