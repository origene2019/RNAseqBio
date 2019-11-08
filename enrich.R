

#BiocManager::install('clusterProfiler')
#BiocManager::install('org.Mm.eg.db')
#BiocManager::install('GSEABase')
#BiocManager::install('GEOquery')
#BiocManager::install('ReactomePA')
#BiocManager::install('pathview')
#BiocManager::install('GenomeInfoDbData')
#BiocManager::available()

 library(clusterProfiler)
 library(org.Mm.eg.db)
 library(GSEABase)
 library(GEOquery)
 library(DOSE)
 library(dplyr)
 library(openxlsx)
 library(enrichplot)


i <- 11
cutn="up"
cutn="down"
cutn=NULL

Enrich <- function(i, cutn=NULL){
  for (d in dev.list()) {
    dev.off(d)
  }

diffgn_dir <- '/home/data/Rseq/P101SC1717020113/DEG/diffGenes'
fl_nm <- dir(diffgn_dir, pattern = 'diffGenes')
nm <- str_replace(fl_nm[i], '_diffGenes.xlsx', '')
diff1 <- read.xlsx(paste0(diffgn_dir, '/', fl_nm[i]))
  #  c(NULL,"_up","_down")
if(!is.null(cutn)) {
  diff1 <- subset(diff1, diff1$sig==cutn)  
  cutn <- paste0('_', cutn)
}
if(!paste0(nm,  cutn) %in% list.files(paste0(diffgn_dir, '/enrich'))) system(paste0('mkdir ', diffgn_dir, '/enrich/', nm, cutn))
setwd(paste0(diffgn_dir, '/enrich/', nm, cutn))
# id转换成ENTREZID
eg <- bitr(diff1$ensembl_gene_id, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
#head(eg)

# 02 GO富集分析
# over-representation test：enrichGO()
genelist <- eg$ENTREZID
genelist[duplicated(genelist)]
go <- enrichGO(genelist, OrgDb = "org.Mm.eg.db", ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.1, keyType = 'ENTREZID')
head(go)
dim(go)
dim(go[go$ONTOLOGY=='BP',])
dim(go[go$ONTOLOGY=='CC',])
dim(go[go$ONTOLOGY=='MF',])
write.xlsx(go, paste0(nm, '_goEnrich.xlsx'))
# 进行简单的可视化
goplt1 <- barplot(go, showCategory=20, drop=T, main = nm)
goplt2 <- dotplot(go, showCategory=20, title = nm)

ggsave(paste0(nm, '_goBar.png'), plot = goplt1, limitsize = FALSE, width = 8, height = 8)
ggsave(paste0(nm, '_goBar.pdf'), plot = goplt1, limitsize = FALSE, width = 8, height = 8)
ggsave(paste0(nm, '_goDot.png'), plot = goplt2, limitsize = FALSE, width = 8, height = 8)
ggsave(paste0(nm, '_goDot.pdf'), plot = goplt2, limitsize = FALSE, width = 8, height = 8)

# 还可以绘制GO的网络关系图，但是值得注意的是这里的数据只能是富集一个GO通路（BP、CC或MF）的数据
go.BP <- enrichGO(genelist, OrgDb = "org.Mm.eg.db", ont='BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.1, keyType = 'ENTREZID')
png(paste0(nm, '_go-BPnet.png'), width = 750, height = 800, title = nm)
print(plotGOgraph(go.BP))
dev.off()
pdf(paste0(nm, '_go-BPnet.pdf'), title = nm)
print(plotGOgraph(go.BP))
dev.off()

go.MF <- enrichGO(genelist, OrgDb = "org.Mm.eg.db", ont='MF', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.1, keyType = 'ENTREZID')
png(paste0(nm, '_go-MFnet.png'), width = 750, height = 800, title = nm)
plotGOgraph(go.MF)
dev.off()
pdf(paste0(nm, '_go-MFnet.pdf'), title = nm)
plotGOgraph(go.MF)
dev.off()

go.CC <- enrichGO(genelist, OrgDb = "org.Mm.eg.db", ont='CC', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.1, keyType = 'ENTREZID')
png(paste0(nm, '_go-CCnet.png'), width = 750, height = 800, title = nm)
plotGOgraph(go.CC)
dev.off()
pdf(paste0(nm, '_go-CCnet.pdf'), title = nm)
plotGOgraph(go.CC)
dev.off()

# 03 KEGG通路富集
kegg <- enrichKEGG(genelist, organism = 'mmu', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', 
                   minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.1,use_internal_data = FALSE)
head(kegg)
if(length(kegg)>0 & dim(kegg)[1]>0){
  write.xlsx(kegg, paste0(nm, '_keggEnrich.xlsx'))
  
  keggplt1 <- barplot(kegg, showCategory=20, drop=T, main = nm)
  keggplt2 <- dotplot(kegg, showCategory=30, title = nm)
  keggplt3 <- cnetplot(kegg, circular = TRUE)
  
  ggsave(paste0(nm, '_keggBar.png'), plot = keggplt1, limitsize = FALSE)
  ggsave(paste0(nm, '_keggBar.pdf'), plot = keggplt1, limitsize = FALSE)
  ggsave(paste0(nm, '_keggDot.png'), plot = keggplt2, limitsize = FALSE)
  ggsave(paste0(nm, '_keggDot.pdf'), plot = keggplt2, limitsize = FALSE)
  ggsave(paste0(nm, '_keggCnet.png'), plot = keggplt3, limitsize = FALSE, width = 65, height = 60, units = "cm")
  ggsave(paste0(nm, '_keggCnet.pdf'), plot = keggplt3, limitsize = FALSE, width = 65, height = 60, units = "cm")
  
  if(nrow(kegg)>1){
    png(paste0(nm, '_keggUpset.png'), width = 1200, height = 800, title = nm)
    print(upsetplot(kegg))
    dev.off()
    pdf(paste0(nm, '_keggUpset.pdf'), title = nm)
    print(upsetplot(kegg))
    dev.off()
  }
}


#数据准备
degenes <- left_join(diff1[c('ensembl_gene_id','log2FoldChange')], eg, by=c('ensembl_gene_id'='ENSEMBL'))
geneList <- degenes[c("ENTREZID", "log2FoldChange")]
id_geneList <- geneList[,2]
names(id_geneList) <- geneList[,1]
id_geneList <- sort(id_geneList, decreasing=TRUE)
head(id_geneList)

# 04 GSEA富集分析
## #GSEA富集分析
## gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
## c5 <- read.gmt(gmtfile)
## head(c5)
## ### 先使用基于超几何分布的富集分析
## enrich <- try(enricher(degenes$ENTREZID, TERM2GENE=c5), silent=T)
## #head(enrich)
## gsea <- try(GSEA(id_geneList, TERM2GENE=c5, verbose=FALSE, pvalueCutoff = 0.05), silent=T)
## #head(gsea)
## if(!class(gsea)=="try-error"){
##   write.xlsx(gsea, paste0(nm, '_gseaEnrich.xlsx'))
##   }
##
##
##  #GSEA的GO富集。
##  gsea.go <- try(gseGO(id_geneList, ont = "BP", OrgDb="org.Mm.eg.db", keyType = "ENTREZID"), silent=T)
##  #head(gsea.go)
##  if(!class(gsea.go)=="try-error"){
##    gseaplt1 <- gseaplot2(gsea.go, geneSetID = rownames(gsea.go[1,]), title = nm)
##    #gseaplt2 <- gseaplot2(gsea, geneSetID = "GO:0006959", title = nm)
##    #gseaplt3 <- gseaplot2(gsea, geneSetID = "GO:0030595", title = nm)
##    tmp <- gsea.go@result
##    table(tmp$pvalue<0.05)
##    write.xlsx(tmp, paste0(nm, '_gsea-GO.xlsx'))
##    
##    ggsave(paste0(nm, '_gseaPlt.png'), plot = gseaplt1, limitsize = FALSE, width = 12, height = 8)
##    ggsave(paste0(nm, '_gseaPlt.pdf'), plot = gseaplt1, limitsize = FALSE, width = 12, height = 8)
##    #ggsave(paste0(nm, '_gseaPlt2.png'), plot = gseaplt2, limitsize = FALSE, width = 8, height = 8)
##    #ggsave(paste0(nm, '_gseaPlt3.png'), plot = gseaplt3, limitsize = FALSE, width = 8, height = 8)
##  }
##  
##  #GSEA的KEGG富集，KEGG富集到的某一条通路的可视化：  # Loading required package: org.Hs.eg.db
##  library(pathview)
##  pathview(id_geneList, pathway.id = 'hsa04658',species="hsa", limit=list(gene=max(abs(id_geneList)), cpd=1)) #直接输出到当前工作路径




# DO富集分析 disease ontology
do <- try(enrichDO(id_geneList, ont = "DO", pvalueCutoff = 0.05, pAdjustMethod = "BH",
        minGSSize = 5, maxGSSize = 500, qvalueCutoff = 0.1, readable = FALSE), silent = TRUE)

if(!is.null(do)){
  write.xlsx(do, paste0(nm, '_doEnrich.xlsx'))
  doplt1 <- barplot(do, showCategory=20, drop=T, main = nm)
  doplt2 <- dotplot(do, showCategory=20, title = nm)
  doplt3 <- cnetplot(do, circular=TRUE)
  
  ggsave(paste0(nm, '_doBar.png'), plot = doplt1, limitsize = FALSE, width = 8, height = 8)
  ggsave(paste0(nm, '_doBar.pdf'), plot = doplt1, limitsize = FALSE, width = 8, height = 8)
  ggsave(paste0(nm, '_doDot.png'), plot = doplt2, limitsize = FALSE, width = 8, height = 8)
  ggsave(paste0(nm, '_doDot.pdf'), plot = doplt2, limitsize = FALSE, width = 8, height = 8)
  ggsave(paste0(nm, '_doCnet.png'), plot = doplt3, limitsize = FALSE, width = 8, height = 8)
  ggsave(paste0(nm, '_doCnet.pdf'), plot = doplt3, limitsize = FALSE, width = 8, height = 8)
  
  png(paste0(nm, '_doUpset.png'), width = 750, height = 800, title = nm)
  upsetplot(do)
  dev.off()
  pdf(paste0(nm, '_doUpset.pdf'), title = nm)
  upsetplot(do)
  dev.off()
  
}


# 04 Reactome pathway enrichment analysis
#Pathway enrichment analysis
require(ReactomePA)
reactome <- try(enrichPathway(id_geneList, pvalueCutoff=0.05, readable=T,organism = "mouse"), silent = TRUE)
#head(as.data.frame(reactome))
#write.table(x,file="ReactomePA.xls",quote=F,sep="\t")
if(!is.null(reactome)){
  write.xlsx(reactome, paste0(nm, '_reactomeEnrich.xlsx'))
  reactomePlt1 <- barplot(reactome, showCategory=24)
  reactomePlt2 <- dotplot(reactome, showCategory=24)
  reactomePlt3 <- cnetplot(reactome, categorySize="pvalue", foldChange=gene, circular=TRUE)  #cnetplot(x, categorySize="pvalue", foldChange=geneList)
  
  ggsave(paste0(nm, '_reactomeBar.png'), plot = reactomePlt1, limitsize = FALSE, width = 8, height = 8)
  ggsave(paste0(nm, '_reactomeDot.png'), plot = reactomePlt2, limitsize = FALSE, width = 8, height = 8)
  ggsave(paste0(nm, '_reactomeCnet.png'), plot = reactomePlt3, limitsize = FALSE, width = 8, height = 8)
  
  ggsave(paste0(nm, '_reactomeBar.pdf'), plot = reactomePlt1, limitsize = FALSE, width = 8, height = 8)
  ggsave(paste0(nm, '_reactomeDot.pdf'), plot = reactomePlt2, limitsize = FALSE, width = 8, height = 8)
  ggsave(paste0(nm, '_reactomeCnet.pdf'), plot = reactomePlt3, limitsize = FALSE, width = 8, height = 8)
  
  png(paste0(nm, '_reactomeUpset.png'), width = 750, height = 800, title = nm)
  upsetplot(reactome)
  dev.off()
  pdf(paste0(nm, '_reactomeUpset.pdf'), title = nm)
  upsetplot(reactome)
  dev.off()
  
  png(paste0(nm, '_reactomeEma.png'), width = 750, height = 800, title = nm)
  emapplot(reactome)  #enrichment map
  dev.off()
  pdf(paste0(nm, '_reactomeEma.pdf'), title = nm)
  emapplot(reactome)  #enrichment map
  dev.off()
}



}


for (i in 13:20) {
  Enrich(i,cutn = 'up')
  Enrich(i,cutn = 'down')
  print(i)
}
