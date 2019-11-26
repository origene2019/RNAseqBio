path <- '/home/data/Rseq/results361111/9.Fusiongene'
file <- dir(path, pattern = 'tsv$')
#install.packages("RCircos")
library(RCircos)
library(stringi)
library(stringr)
library(ggplot2)
library(dplyr)
options(stringsAsFactors = FALSE)
for (t in file) {
  dat <- read.delim2(paste0(path,'/',t))
  if(nrow(dat)>0){
  sco_nm <- stringr::str_replace(t,'Fusiongene.tsv','')
  dat2 <- dat[,c("X.FusionName","LeftLocalBreakpoint","LeftBreakpoint","RightLocalBreakpoint","RightBreakpoint")]
  genes_nm <- as.data.frame(stringi::stri_split_fixed(dat2$X.FusionName,'--',simplify = TRUE))
  colnames(genes_nm) <- c('Gene','gene_nm2')
  chr_pos1 <- as.data.frame(stringi::stri_split_fixed(dat2$LeftBreakpoint,':',simplify = TRUE))
  colnames(chr_pos1) <- c('Chromosome','chr_pos1','mult1')
  chr_pos2 <- as.data.frame(stringi::stri_split_fixed(dat2$RightBreakpoint,':',simplify = TRUE))
  colnames(chr_pos2) <- c('Chromosome.1','chr_pos2','mult2')
  dat3 <- cbind(dat2$X.FusionName,genes_nm,dat2$LeftLocalBreakpoint,chr_pos1,dat2$RightLocalBreakpoint,chr_pos2)
  dat3$chr_staend1 <- as.numeric(as.character(dat3$chr_pos1))+ifelse(dat3$mult1=='-',-1,1)*dat3$`dat2$LeftLocalBreakpoint`
  dat3$chr_staend2 <- as.numeric(as.character(dat3$chr_pos2))+ifelse(dat3$mult1=='-',-1,1)*dat3$`dat2$RightLocalBreakpoint`
  dat4 <- dat3[,c("Gene","gene_nm2","Chromosome","chr_pos1","chr_staend1","Chromosome","chr_pos2","chr_staend2","Chromosome.1")]
  dat4$chromStart <- ifelse(dat4$chr_staend1>as.numeric(as.character(dat4$chr_pos1)),as.numeric(as.character(dat4$chr_pos1)),dat4$chr_staend1)
  dat4$chromEnd <- ifelse(dat4$chr_staend1>as.numeric(as.character(dat4$chr_pos1)),dat4$chr_staend1,as.numeric(as.character(dat4$chr_pos1)))
  dat4$chromStart.1 <- ifelse(dat4$chr_staend2>as.numeric(as.character(dat4$chr_pos2)),as.numeric(as.character(dat4$chr_pos2)),dat4$chr_staend2)
  dat4$chromEnd.1 <- ifelse(dat4$chr_staend2>as.numeric(as.character(dat4$chr_pos2)),dat4$chr_staend2,as.numeric(as.character(dat4$chr_pos2)))
  dat5 <- dat4[,c("Gene","Chromosome","chromStart","chromEnd","gene_nm2","Chromosome.1","chromStart.1","chromEnd.1")]
  
  Link.Data <- dat5[,c("Chromosome","chromStart","chromEnd","Chromosome.1","chromStart.1","chromEnd.1")]
  aa <- dat5[,c("Chromosome","chromStart","chromEnd","Gene")]
  bb <- dat5[,c("Chromosome.1","chromStart.1","chromEnd.1","gene_nm2")]
  colnames(bb) <- colnames(aa)
  Label.Data <- rbind(aa,bb)
  
  data("UCSC.Mouse.GRCm38.CytoBandIdeogram")
  ucsc <- UCSC.Mouse.GRCm38.CytoBandIdeogram
  if(sum(Label.Data$Chromosome=='chrM')>0) {
    subcyto <- Label.Data[which(Label.Data$Chromosome=='chrM'),c(1:3)]
    cyto <- cbind(subcyto, "Band"=paste0('qA',1:nrow(subcyto)),"Stain"=sample(c("gpos100","gpos33","gneg","gpos66"),nrow(subcyto),replace = TRUE))
    cyto$chromStart[1] <- 0
    colnames(cyto) <- colnames(UCSC.Mouse.GRCm38.CytoBandIdeogram)
    ucsc <- rbind(UCSC.Mouse.GRCm38.CytoBandIdeogram,cyto)
  }else{
    cyto <- as.data.frame(cbind("Chromosome"='chrM',"ChromStart"=0,"ChromEnd"=10000, "Band"='qA1',"Stain"='gpos100'))
    colnames(cyto) <- colnames(UCSC.Mouse.GRCm38.CytoBandIdeogram)
    ucsc <- as.data.frame(rbind(UCSC.Mouse.GRCm38.CytoBandIdeogram,cyto))
    ucsc$ChromStart <- as.numeric(ucsc$ChromStart)
    ucsc$ChromEnd <- as.numeric(ucsc$ChromEnd)
  }
  b <- function(...){
    #Gene Labels and connectors on RCircos Plot
    RCircos.Set.Core.Components(cyto.info=ucsc,chr.exclude=NULL,tracks.inside=5,tracks.outside=0)
    RCircos.Set.Plot.Area()
    RCircos.Chromosome.Ideogram.Plot() ##做出最外层
    #RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col=5,track.num=1, side='in', by.fold=1) ##从外到内的第二层
    #RCircos.Histogram.Plot(RCircos.Histogram.Data,data.col = 4,track.num = 3,side = 'in')##第三层
    RCircos.Link.Plot(Link.Data,track.num=4,TRUE)##最中间的线
    name.col <- 4;
    side <- "in";
    track.num <- 1;
    RCircos.Gene.Connector.Plot(Label.Data, track.num, side);
    track.num <- 2;
    RCircos.Gene.Name.Plot(Label.Data, name.col,track.num, side);
  }
  print(t)
  png(paste0(path,'/Circos/',sco_nm,'_circos.png'), width = 800, height = 800)
  print(b())
  dev.off()
  }
}
