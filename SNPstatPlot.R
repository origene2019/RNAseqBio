
if(!require(XML)) install.packages("XML")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(dplyr)) install.packages("dplyr")
if(!require(reshape2)) install.packages("reshape2")
library(XML)
library(ggplot2)
library(reshape2)
library(dplyr)
grpTab <- read.delim2('/home/data/Rseq/results361111/resultsAll/grpTab.txt', stringsAsFactors = FALSE)
colnames(grpTab) <- c('grp_no', 'treatment', 'control')
grpTyp <- read.delim2('/home/data/Rseq/results361111/resultsAll/grpTyp.txt', stringsAsFactors = FALSE)
colnames(grpTyp) <- c('Sample','Group')

snp_dir <- "/home/data/Rseq/results361111/resultsAll/SNP"
html_file <- dir(snp_dir,pattern = '.html')

hm_lst5 <- list()
hm_lst6 <- list()
hm_lst9 <- list()
for(i in 1:length(html_file)){                  #建立循环
  html <- readHTMLTable(paste0(snp_dir, '/', html_file[[i]]))                     #读取html中表格的内容
  htmlxls5 <- html[[5]]
  htmlxls6 <- html[[6]] 
  htmlxls9 <- html[[9]] 
  xls5 <- htmlxls5[-1,-2]  
  xls6 <- htmlxls6[-1,-2]
  xls9 <- htmlxls9[-1,-2]
  xls5$sample <- gsub(".html","",html_file[[i]],fixed = TRUE) 
  xls6$sample <- gsub(".html","",html_file[[i]],fixed = TRUE) 
  xls9$sample <- gsub(".html","",html_file[[i]],fixed = TRUE)
  names(xls5)<- c("Type","Count","Percent","Sample")
  names(xls6)<- c("Type","Count","Percent","Sample")
  names(xls9)<- c("Type","Count","Percent","Sample")
  hm_lst5[[i]] <- xls5
  hm_lst6[[i]] <- xls6
  hm_lst9[[i]] <- xls9
}

all_hm5 <- bind_rows(hm_lst5)
all_hm5$Count<-as.numeric(gsub(",","",all_hm5[,2],fixed = TRUE))   #将count中的字符串中的“，”去掉
all_hm5 <- left_join(all_hm5, grpTyp, by=c('Sample'='Sample'))

all_hm6 <- bind_rows(hm_lst6)
all_hm6$Count<-as.numeric(gsub(",","",all_hm6[,2],fixed = TRUE))   #将count中的字符串中的“，”去掉
all_hm6 <- left_join(all_hm6, grpTyp, by=c('Sample'='Sample'))

all_hm9 <- bind_rows(hm_lst9)
all_hm9$Count<-as.numeric(gsub(",","",all_hm9[,2],fixed = TRUE))   #将count中的字符串中的“，”去掉
all_hm9 <- left_join(all_hm9, grpTyp, by=c('Sample'='Sample'))

#准备画图
if(!'snpPlot' %in% list.files(paste0(snp_dir))) system(paste0('mkdir ', snp_dir, '/snpPlot'))
options(scipen=200)
SNPstatPlot <- function(j){
  nm <- paste0(grpTab$treatment[[j]], "_vs_", grpTab$control[[j]])
  subtab5 <- subset(all_hm5,all_hm5$Group %in% c(grpTab$treatment[[j]],grpTab$control[[j]]))
  p5 <- ggplot(subtab5,aes(x=Sample,y=Count,fill=Type)) + theme_bw()
  p5 <- p5 + geom_bar(stat="identity",width = 0.6)+labs(title="Number of effects by impact", subtitle = nm) 
  p5 <- p5 + theme_classic() + theme(plot.title=element_text(size=rel(2)))
  p5 <- p5 + scale_fill_brewer(palette="Accent")
  ggsave(file=paste0(snp_dir, '/snpPlot/', nm, '_impactStat.png'), plot=p5, width=12, height=8)
  ggsave(file=paste0(snp_dir, '/snpPlot/', nm, '_impactStat.pdf'), plot=p5, width=12, height=8)
  
  subdf5 <- dcast(subtab5[,c(1,2,4)],Type~Sample)
  colnames(subdf5)[1] <- "Type (alphabetical order)"
  
  subtab6 <- subset(all_hm6,all_hm6$Group %in% c(grpTab$treatment[[j]],grpTab$control[[j]]))
  p6 <- ggplot(subtab6,aes(x=Sample,y=Count,fill=Type)) + geom_bar(stat="identity",width = 0.6)+labs(title="Number of effects by functional class", subtitle = nm) 
  p6 <- p6 + theme_classic() + theme(plot.title=element_text(size=rel(2)))
  p6 <- p6 + scale_fill_brewer(palette="Accent")
  ggsave(file=paste0(snp_dir, '/snpPlot/', nm, '_functionalStat.png'), plot=p6, width=12, height=8)
  ggsave(file=paste0(snp_dir, '/snpPlot/', nm, '_functionalStat.pdf'), plot=p6, width=12, height=8)
  
  subdf6 <- dcast(subtab6[,c(1,2,4)],Type~Sample)
  colnames(subdf6)[1] <- "Type (alphabetical order)"
  
  subtab9 <- subset(all_hm9,all_hm9$Group %in% c(grpTab$treatment[[j]],grpTab$control[[j]]))
  p9 <- ggplot(subtab9,aes(x=Sample,y=Count,fill=Type))
  p9 <- p9 + geom_bar(stat="identity",width = 0.6)+labs(title="Number of effects by type and region", subtitle = nm) 
  p9 <- p9 + theme_classic() + theme(plot.title=element_text(size=rel(2)))
  colourCount <- length(unique(subtab9$Type))
  getPalette <- colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))
  p9 <- p9 + scale_fill_manual(values = getPalette(colourCount))
  ggsave(file=paste0(snp_dir, '/snpPlot/', nm, '_regionStat.png'), plot=p9, width=12, height=8)
  ggsave(file=paste0(snp_dir, '/snpPlot/', nm, '_regionStat.pdf'), plot=p9, width=12, height=8)
  
  subdf9 <- dcast(subtab9[,c(1,2,4)],Type~Sample)
  colnames(subdf9)[1] <- "Type (alphabetical order)"
  sheets <- list("by impact"=a,"by functional class"=b,"by region"=c)
  write.xlsx(sheets,file =paste0(snp_dir, '/snpPlot/', nm, '_snp.xlsx'))
}
#画图颜色要调整
library(pbmcapply)
mc <- getOption("mc.cores", detectCores(logical = F)-2)
res <- pbmclapply(1:nrow(grpTab), SNPstatPlot, mc.preschedule = FALSE, mc.cleanup = FALSE, mc.cores = mc)
