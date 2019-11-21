

stat_dir <- "/home/data/Rseq/results361111/resultsAll/sortBam"
stats <- dir(stat_dir, pattern = '.stat$')

if(!require(dplyr)) install.packages("dplyr")
if(!require(openxlsx)) install.packages("openxlsx")
if(!require(reshape)) install.packages("reshape")
if(!require(scales)) install.packages("scales")
library(dplyr)
library(openxlsx)
library(reshape)
library(scales)
mapstalst <- list()
for (s in stats) {
  stanm <- str_replace(s,'_sort_mapping.stat','')
  mapsta <- read.table(paste0(stat_dir,'/',s),header = FALSE, sep = ":", fill = TRUE, nrows = 15, skip = 5)
  mapstat <- mapsta[c(-2,-3,-4,-15),]
  colnm <- c("total_reads","unmapped_reads","multi_map","unique_map","read1_map","read2_map",
             "positive_map","negative_map","splice_map","unsplice_map","proper_map")
  mapstat[,1] <- colnm
  mapstat <- rbind(rbind(mapstat[1,],c("total_map", sum(as.numeric(mapstat[which(mapstat$V1 %in% c("multi_map","unique_map")),2])))),mapstat[2:11,])
  mapstat[1,2] <- sum(as.numeric(mapstat[which(mapstat$V1 %in% c("unmapped_reads","multi_map","unique_map")),2]))
  mapstat$V2 <- as.numeric(mapstat$V2)
  mapstat$V3 <- mapstat[,2]/mapstat[1,2]
  mapstat[12,3] <- mapstat[12,2]/mapstat[2,2]
  mapstat[,stanm] <- paste0(mapstat[,2],'(',percent(mapstat[,3]),')')
  mapstat[1,4] <- mapstat[1,2]
  mapstalst[[stanm]] <- mapstat[,c(1,4)]
}
mapstall <- merge_recurse(mapstalst,by='V1')
allstats <- t(mapstall)
rownm <- row.names(allstats)
colnm <- allstats[1,]
allstats <- as.data.frame(allstats[-1,])
colnames(allstats) <- colnm
allstats <- cbind(sample=rownm[-1],allstats)
write.xlsx(allstats, paste0(stat_dir,'/mapstats.xlsx'),rowNames=FALSE,colNames=TRUE)
