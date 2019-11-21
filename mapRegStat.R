
if(!require(dplyr)) install.packages("dplyr")
if(!require(openxlsx)) install.packages("openxlsx")
if(!require(reshape)) install.packages("reshape")
if(!require(scales)) install.packages("scales")
library(dplyr)
library(openxlsx)
library(reshape)
library(scales)
regstat_dir <- "/home/data/Rseq/results361111/resultsAll/sortBam"
output_dir <- "/home/data/Rseq/results361111/resultsAll/map"
regstats <- dir(regstat_dir, pattern = '.distrib$')
regstats

maplst <- list()
for(r in regstats){
stanm <- str_replace(r,'_sort.distrib','')
mapsta <- read.table(paste0(stat_dir,'/',r),header = TRUE, sep = "", fill = TRUE, nrows = 10, skip = 4)
region_stat <- as.data.frame(cbind(region=c("exon","intron","intergenic"),Tag_count=c(sum(mapsta[1:3,"Tag_count"]), sum(mapsta[4,"Tag_count"]), sum(mapsta[5:10,"Tag_count"]))))
region_stat$prop <- as.numeric(region_stat$Tag_count)/sum(as.numeric(region_stat$Tag_count))                                      
region_stat$count_prop <- paste0(region_stat$Tag_count,'(',percent(region_stat$prop),')')

base <- ggplot(region_stat, aes(x = "Cotent", y = prop, fill = region)) + geom_bar(stat = 'identity', position = 'stack',width = 1) + coord_polar(theta = "y") 
newlegend <- paste0(region_stat$region, " (", region_stat$Tag_count, ',', percent(region_stat$prop), ')')
base <- base + scale_fill_manual(breaks = region_stat$region, labels = newlegend, values=c("#7CCD7C","#CD96CD","#FFD39B")) + theme_minimal()
base <- base + theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), legend.title = element_blank(), legend.position=c(0.5, -0.05),plot.margin = unit(c(2,3,3,4),"cm")) 
base <- base + labs(title=paste0("Percent of genome regions (", stanm, ")"))
ggsave(file=paste0(output_dir, '/', stanm, '_mapRegion.png'), plot=base, width=12, height=8)
ggsave(file=paste0(output_dir, '/', stanm, '_mapRegion.pdf'), plot=base, width=12, height=8)

region_stat2 <- region_stat[,c(1,4)]
region_stat2 <- rbind(c('sample',stanm), region_stat2) 
region_stat2 <- as.data.frame(t(region_stat2))
colnames(region_stat2) <- region_stat2[1,]
maplst[[stanm]] <- region_stat2[-1,]
}

mapregion <- bind_rows(maplst)
write.xlsx(mapregion,paste0(output_dir,'/align_region.xlsx'))
