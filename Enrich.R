
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
library(BiocManager)
if(!require(clusterProfiler)) BiocManager::install('clusterProfiler')
library(clusterProfiler)
if(!require(org.Mm.eg.db)) BiocManager::install('org.Mm.eg.db')
library(org.Mm.eg.db)
if(!require(GSEABase)) BiocManager::install('GSEABase')
library(GSEABase)
if(!require(GEOquery)) BiocManager::install('GEOquery')
library(GEOquery)
if(!require(DOSE)) BiocManager::install('DOSE')
library(DOSE)
if(!require(topGO)) BiocManager::install('topGO')
library(topGO)
if(!require(enrichplot)) BiocManager::install('enrichplot')
library(enrichplot)
if(!require(dplyr)) install.packages("dplyr")
if(!require(openxlsx)) install.packages("openxlsx")
if(!require(stringr)) install.packages("stringr")
if(!require(ggplot2)) install.packages("ggplot2")
library(dplyr)
library(openxlsx)
library(stringr)
library(ggplot2)

i <- 8
cutn="up"
cutn="down"
cutn=NULL
Enrich <- function(i, cutn=NULL){
  for (d in dev.list()) {
    dev.off(d)
  }
  
  diffgn_dir <- '/home/data/Rseq/P101SC1717020113/DEG/diffGenes'
  fl_nm <- dir(diffgn_dir, pattern = '_diffGenes.xlsx')
  nm <- str_replace(fl_nm[i], '_diffGenes.xlsx', '')
  diff1 <- read.xlsx(paste0(diffgn_dir, '/', fl_nm[i]))
  if(!is.null(cutn)) {
    diff1 <- subset(diff1, diff1$sig==cutn)  
    cutnm <- paste0('_', cutn)  #  c(NULL,"_up","_down")
  }else{cutnm <- cutn}
  if(!paste0(nm,  cutnm) %in% list.files(paste0(diffgn_dir, '/enrich'))) system(paste0('mkdir ', diffgn_dir, '/enrich/', nm, cutnm))
  setwd(paste0(diffgn_dir, '/enrich/', nm, cutnm))
  if(!is.null(diff1) & nrow(diff1)>0){
  eg <- bitr(diff1$ensembl_gene_id, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")  # id转换成ENTREZID
  # 02 GO富集分析
  # over-representation test：enrichGO()
  genelist <- eg$ENTREZID
  #genelist[duplicated(genelist)]
  go <- enrichGO(genelist, OrgDb = "org.Mm.eg.db", ont='ALL',pAdjustMethod = 'BH', 
                 #pvalueCutoff = 0.05, qvalueCutoff = 0.1, 
                 keyType = 'ENTREZID')
  if(!is.null(go)){
  dim(go[go$ONTOLOGY=='BP',])
  dim(go[go$ONTOLOGY=='CC',])
  dim(go[go$ONTOLOGY=='MF',])
  write.xlsx(go@result, paste0(nm, '_goEnrich.xlsx'))
  
  go@result$Description <- ifelse(nchar(go@result$Description)>60, paste0(substr(go@result$Description,1,60),'...'),go@result$Description)
  go_top10 <- go@result %>%
    group_by(ONTOLOGY) %>%
    arrange(-Count) %>%
    mutate(row_number = row_number(-Count)) %>%
    subset(row_number<=20) %>%
    arrange(ONTOLOGY,-Count)
  
  GO_term_order <- factor(as.integer(rownames(go_top10)), labels = go_top10$Description)
  p <- ggplot(data=go_top10, aes(x=GO_term_order, y=Count, fill=ONTOLOGY))
  p <- p + geom_bar(stat='identity', width = 0.8) + theme_bw() + scale_fill_manual(values = c("#8DA1CB", "#FD8D62", "#66C3A5"))
  p <- p + xlab('GO term') + ylab('Num of Genes') + labs(title = paste0('The Most Enriched GO Terms\n', nm))
  p <- p + theme(axis.text.x = element_text(face = 'bold', color = 'gray50', angle = 70, vjust = 1, hjust = 1))
  goplt1 <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5))
  ggsave(paste0(nm, '_goBar.png'), plot = goplt1, limitsize = FALSE, width = 8, height = 8)
  ggsave(paste0(nm, '_goBar.pdf'), plot = goplt1, limitsize = FALSE, width = 8, height = 8)
  
  gnsets_lst <- list()
  gnsets <- go@geneSets
  for (n in 1:length(gnsets)) {
    gene_id <- gnsets[[n]]
    go_nm <- names(gnsets)[[n]]
    go_id <- rep(go_nm,length(gene_id))
    gnsets_lst[[go_nm]] <- as.data.frame(cbind(gene_id,go_id))
  }
  go_df <- bind_rows(gnsets_lst)
  go_diff <- left_join(go_df, eg, by=c('gene_id'='ENTREZID'))
  go_diff <- left_join(go_diff,diff1[c("ensembl_gene_id","sig")],by=c('ENSEMBL'='ensembl_gene_id'))
  go_diff <- left_join(go_diff,go2@result[c("ID","ONTOLOGY","Description","Count")],by=c('go_id'='ID'))
  go_diff <- na.omit(go_diff)
  go_diff <- distinct(go_diff[,c(-1,-3)])
  go_diff <- arrange(go_diff, -Count)
  
  if(!is.null(go_diff) & nrow(go_diff)>0){
  go_on10 <- go_diff %>%
    group_by(ONTOLOGY) %>%
    arrange(-Count) %>%
    mutate(row_number = row_number(-Count)) %>%
    subset(row_number<=10)
  GO_term_order <- factor(as.integer(rownames(go_on10)), labels = go_on10$Description)
  p <- ggplot(data=go_on10, aes(x=go_on10$Description, y=Count, fill=sig))
  p <- p + geom_bar(stat='identity',width=0.5,position="dodge") + theme_bw() 
  p <- p + xlab('GO term') + ylab('Num of Genes') + labs(title = paste0('The Most Enriched GO Terms\n', nm))
  p <- p + theme(axis.text.x = element_text(face = 'bold', color = 'gray50', angle = 70, vjust = 1, hjust = 1))
  p2 <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5))
  goplton <- p2 + facet_wrap(~ ONTOLOGY, scales = "free_x")
  ggsave(paste0(nm, '_goDiff.png'), plot = goplton, limitsize = FALSE, width = 8, height = 8)
  ggsave(paste0(nm, '_goDiff.pdf'), plot = goplton, limitsize = FALSE, width = 8, height = 8)
  }
  
  barplt1 <- barplot(go, showCategory=20, title = paste0('The Most Enriched GO Terms\n',nm))
  ggsave(paste0(nm, '_goBar20.png'), plot = barplt1, limitsize = FALSE, width = 8, height = 8)
  ggsave(paste0(nm, '_goBar20.pdf'), plot = barplt1, limitsize = FALSE, width = 8, height = 8)
  
  goplt2 <- dotplot(go, showCategory=20, title = paste0('The Most Enriched GO Terms\n',nm))
  ggsave(paste0(nm, '_goDot.png'), plot = goplt2, limitsize = FALSE, width = 8, height = 8)
  ggsave(paste0(nm, '_goDot.pdf'), plot = goplt2, limitsize = FALSE, width = 8, height = 8)
  
  # 还可以绘制GO的网络关系图，但是值得注意的是这里的数据只能是富集一个GO通路（BP、CC或MF）的数据
  go.BP <- enrichGO(genelist, OrgDb = "org.Mm.eg.db", ont='BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.1, keyType = 'ENTREZID')
  if(!is.null(go.BP)){
  png(paste0(nm, '_go-BPnet.png'), width = 750, height = 800, title =  paste0('The Most Enriched GO-BP Terms\n', nm))
  print(plotGOgraph(go.BP))
  dev.off()
  pdf(paste0(nm, '_go-BPnet.pdf'), title =  paste0('The Most Enriched GO-BP Terms\n', nm))
  print(plotGOgraph(go.BP))
  dev.off()
  }
  go.MF <- enrichGO(genelist, OrgDb = "org.Mm.eg.db", ont='MF', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.1, keyType = 'ENTREZID')
  if(!is.null(go.MF)){
  png(paste0(nm, '_go-MFnet.png'), width = 750, height = 800, title =  paste0('The Most Enriched GO-MF Terms\n', nm))
  plotGOgraph(go.MF)
  dev.off()
  pdf(paste0(nm, '_go-MFnet.pdf'), title =  paste0('The Most Enriched GO-MF Terms\n', nm))
  plotGOgraph(go.MF)
  dev.off()
  }
  go.CC <- enrichGO(genelist, OrgDb = "org.Mm.eg.db", ont='CC', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.1, keyType = 'ENTREZID')
  if(!is.null(go.CC)){
  png(paste0(nm, '_go-CCnet.png'), width = 750, height = 800, title =  paste0('The Most Enriched GO-CC Terms\n', nm))
  plotGOgraph(go.CC)
  dev.off()
  pdf(paste0(nm, '_go-CCnet.pdf'), title =  paste0('The Most Enriched GO-CC Terms\n', nm))
  plotGOgraph(go.CC)
  dev.off()
  }
  }
  
  # 03 KEGG通路富集
  kegg <- enrichKEGG(genelist, organism = 'mmu', keyType = 'kegg', pAdjustMethod = 'BH')
  if(!is.null(kegg)){
    if(!class(kegg)=="try-error" & nrow(kegg@result)>0){
      write.xlsx(kegg@result, paste0(nm, '_keggEnrich.xlsx'))
      
      keggplt1 <- barplot(kegg, showCategory=20, title = paste0('The Most Enriched KEGG Terms\n',nm))
      ggsave(paste0(nm, '_keggBar.png'), plot = keggplt1, limitsize = FALSE)
      ggsave(paste0(nm, '_keggBar.pdf'), plot = keggplt1, limitsize = FALSE)
      
      if(nrow(kegg)>0){
        keggplt2 <- dotplot(kegg, showCategory=20, title = paste0('The Most Enriched KEGG Terms\n',nm))
        ggsave(paste0(nm, '_keggDot.png'), plot = keggplt2, limitsize = FALSE)
        ggsave(paste0(nm, '_keggDot.pdf'), plot = keggplt2, limitsize = FALSE)
        
        keggplt3 <- cnetplot(kegg, circular = TRUE)
        ggsave(paste0(nm, '_keggCnet.png'), plot = keggplt3, limitsize = FALSE, width = 65, height = 60, units = "cm")
        ggsave(paste0(nm, '_keggCnet.pdf'), plot = keggplt3, limitsize = FALSE, width = 65, height = 60, units = "cm")
        
        upsetpk <- try(upsetplot(kegg),silent = TRUE)
        if(!is.null(upsetpk) & !'try-error' %in% class(upsetpk)){
        ggsave(paste0(nm, '_keggCnet.png'), plot = keggplt3, limitsize = FALSE, width = 65, height = 60, units = "cm")
        ggsave(paste0(nm, '_keggCnet.pdf'), plot = keggplt3, limitsize = FALSE, width = 65, height = 60, units = "cm")
          
        png(paste0(nm, '_keggUpset.png'), width = 1200, height = 800, title =  paste0('The Most Enriched KEGG-Upset Terms\n', nm))
        print(upsetplot(kegg))
        dev.off()
        pdf(paste0(nm, '_keggUpset.pdf'), title =  paste0('The Most Enriched KEGG-Upset Terms\n', nm))
        print(upsetplot(kegg))
        dev.off()
        }
      }
    }
  }
  
  #数据准备
  degenes <- left_join(diff1[c('ensembl_gene_id','log2FoldChange')], eg, by=c('ensembl_gene_id'='ENSEMBL'))
  geneList <- na.omit(degenes[c("ENTREZID", "log2FoldChange")])
  geneList <- arrange(geneList, desc(log2FoldChange))
  head(geneList)
  id_geneList <- geneList$log2FoldChange
  names(id_geneList) <- geneList$ENTREZID
  id_geneList <- sort(id_geneList, decreasing=TRUE)
  head(id_geneList)
  
  # DO富集分析 disease ontology
  do <- try(enrichDO(geneList$ENTREZID, ont = "DO", pAdjustMethod = "BH",
                     # pvalueCutoff = 0.05, qvalueCutoff = 0.1, 
                     readable = FALSE), silent = TRUE)
  if(!is.null(do)){
    if(!"try-error" %in% class(do) & nrow(do@result)>0){
      write.xlsx(do@result, paste0(nm, '_doEnrich.xlsx'))
      
      do@result$Description <- ifelse(nchar(do@result$Description)>60, paste0(substr(do@result$Description,1,60),'...'),do@result$Description)
      doplt1 <- barplot(do, showCategory=20, title = paste0('The Most Enriched DO Terms\n',nm))
      doplt2 <- dotplot(do, showCategory=20, title = paste0('The Most Enriched DO Terms\n',nm))
      doplt3 <- try(cnetplot(do, circular=TRUE),silent = TRUE)
      
      ggsave(paste0(nm, '_doBar.png'), plot = doplt1, limitsize = FALSE, width = 8, height = 8)
      ggsave(paste0(nm, '_doBar.pdf'), plot = doplt1, limitsize = FALSE, width = 8, height = 8)
      ggsave(paste0(nm, '_doDot.png'), plot = doplt2, limitsize = FALSE, width = 8, height = 8)
      ggsave(paste0(nm, '_doDot.pdf'), plot = doplt2, limitsize = FALSE, width = 8, height = 8)
      if(!"try-error" %in% class(doplt3)){
        ggsave(paste0(nm, '_doCnet.png'), plot = doplt3, limitsize = FALSE, width = 8, height = 8)
        ggsave(paste0(nm, '_doCnet.pdf'), plot = doplt3, limitsize = FALSE, width = 8, height = 8)
      }
      
      upset <- try(upsetplot(do),silent = TRUE)
      if("ggplot" %in% class(upset)){
        png(paste0(nm, '_doUpset.png'), width = 1200, height = 800, title =  paste0('The Most Enriched DO-Upset Terms\n', nm))
        print(upsetplot(do))
        dev.off()
        pdf(paste0(nm, '_doUpset.pdf'), title =  paste0('The Most Enriched DO-Upset Terms\n', nm))
        print(upsetplot(do))
        dev.off()
      }
    }
  }
  
  # 04 Reactome pathway enrichment analysis
  #Pathway enrichment analysis
  library(ReactomePA)
  reactome <- try(ReactomePA::enrichPathway(geneList$ENTREZID, readable=T,organism = "mouse"), silent = TRUE)
  if(!is.null(reactome)){
    if(!"try-error" %in% class(reactome) & nrow(reactome@result)>0){
      write.xlsx(as.data.frame(reactome@result), paste0(nm, '_reactomeEnrich.xlsx'))
      
      reactome@result$Description <- ifelse(nchar(reactome@result$Description)>60, paste0(substr(reactome@result$Description,1,60),'...'),reactome@result$Description)
      reactomePlt1 <- barplot(reactome, showCategory=20, title = paste0('The Most Enriched Reactome Terms\n',nm))
      ggsave(paste0(nm, '_reactomeBar.png'), plot = reactomePlt1, limitsize = FALSE, width = 8, height = 8)
      ggsave(paste0(nm, '_reactomeBar.pdf'), plot = reactomePlt1, limitsize = FALSE, width = 8, height = 8)
      
      if(nrow(reactome)>0) {
        reactomePlt2 <- dotplot(reactome, showCategory=20, title=paste0('The Most Enriched Reactome Terms\n',nm))
        ggsave(paste0(nm, '_reactomeDot.png'), plot = reactomePlt2, limitsize = FALSE, width = 8, height = 8)
        ggsave(paste0(nm, '_reactomeDot.pdf'), plot = reactomePlt2, limitsize = FALSE, width = 8, height = 8)
        
        reactomePlt3 <- cnetplot(reactome, categorySize="pvalue",  circular=TRUE)  
        ggsave(paste0(nm, '_reactomeCnet.png'), plot = reactomePlt3, limitsize = FALSE, width = 8, height = 8)
        ggsave(paste0(nm, '_reactomeCnet.pdf'), plot = reactomePlt3, limitsize = FALSE, width = 8, height = 8)
        
        ggsave(paste0(nm, '_reactomeUpset.png'), plot = upsetplot(reactome), limitsize = FALSE, width = 8, height = 8)
        ggsave(paste0(nm, '_reactomeUpset.pdf'), plot = upsetplot(reactome), limitsize = FALSE, width = 8, height = 8)
        ggsave(paste0(nm, '_reactomeEma.png'), plot = emapplot(reactome), limitsize = FALSE, width = 8, height = 8)
        ggsave(paste0(nm, '_reactomeEma.pdf'), plot = emapplot(reactome), limitsize = FALSE, width = 8, height = 8)
      }
    }
      
    }
  }
  
  
  
  # 04 GSEA富集分析
  #GSEA富集分析
  
  diffgn_dir2 <- '/home/data/Rseq/P101SC1717020113/DEG/diffGenes/Cnt'
  fl_nm2 <- dir(diffgn_dir2, pattern = '_counts.xlsx')
  nm <- str_replace(fl_nm2[i], '_counts.xlsx', '')
  diff2 <- read.xlsx(paste0(diffgn_dir2, '/', fl_nm2[i]))
  eg2 <- bitr(diff2[,1], fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  degenes2 <- left_join(diff2[c('ensembl_gene_id','log2FoldChange')], eg2, by=c('ensembl_gene_id'='ENSEMBL'))
  geneList2 <- na.omit(degenes2[c("ENTREZID", "log2FoldChange")])
  geneList2 <- arrange(geneList2, desc(log2FoldChange))
  head(geneList2)
  id_geneList2 <- geneList2$log2FoldChange
  names(id_geneList2) <- geneList2$ENTREZID
  id_geneList2 <- sort(id_geneList2, decreasing=TRUE)
  head(id_geneList2)
  
  
  #GSEA的GO富集。
  gsea.goBP <- try(gseGO(id_geneList2, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont="BP"), silent=T)
  if(!is.null(gsea.goBP)){
    if(!"try-error" %in% class(gsea.goBP) & nrow(gsea.goBP@result)>0){
      write.xlsx(gsea.goBP@result, paste0(nm, '_gseaGO-BP.xlsx'))
      gsea.goBP@result$Description <- ifelse(nchar(gsea.goBP@result$Description)>60, paste0(substr(gsea.goBP@result$Description,1,60),'...'),gsea.goBP@result$Description)
      gseaplt1BP <- gseaplot2(gsea.goBP, geneSetID = rownames(gsea.goBP[1,]), title = paste0('GO-BP PI3K/AKT Signaling Pathway\n',nm))
      ggsave(paste0(nm, '_gseaGO-BP.png'), plot = gseaplt1BP, limitsize = FALSE, width = 12, height = 8)
      ggsave(paste0(nm, '_gseaGO-BP.pdf'), plot = gseaplt1BP, limitsize = FALSE, width = 12, height = 8)
      #ggsave(paste0(nm, '_gseaPlt2.png'), plot = gseaplt2, limitsize = FALSE, width = 8, height = 8)
      #ggsave(paste0(nm, '_gseaPlt3.png'), plot = gseaplt3, limitsize = FALSE, width = 8, height = 8)
    }}
  
  gsea.goMF <- try(gseGO(id_geneList2, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont="MF"), silent=T)
  if(!is.null(gsea.goMF)){
    if(!"try-error" %in% class(gsea.goMF) & nrow(gsea.goMF@result)>0){
      write.xlsx(gsea.goMF@result, paste0(nm, '_gseaGO-MF.xlsx'))
      gsea.goMF@result$Description <- ifelse(nchar(gsea.goMF@result$Description)>60, paste0(substr(gsea.goMF@result$Description,1,60),'...'),gsea.goMF@result$Description)
      gseaplt1MF <- gseaplot2(gsea.goMF, geneSetID = rownames(gsea.goMF[1,]), title = paste0('GO-MF PI3K/AKT Signaling Pathway\n',nm))
      ggsave(paste0(nm, '_gseaGO-MF.png'), plot = gseaplt1MF, limitsize = FALSE, width = 12, height = 8)
      ggsave(paste0(nm, '_gseaGO-MF.pdf'), plot = gseaplt1MF, limitsize = FALSE, width = 12, height = 8)
    }}
  
  gsea.goCC <- try(gseGO(id_geneList2, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont="CC"), silent=T)
  if(!is.null(gsea.goCC)){
    if(!"try-error" %in% class(gsea.goCC) & nrow(gsea.goCC@result)>0){
      write.xlsx(gsea.goCC@result, paste0(nm, '_gseaGO-CC.xlsx'))
      gsea.goCC@result$Description <- ifelse(nchar(gsea.goCC@result$Description)>60, paste0(substr(gsea.goCC@result$Description,1,60),'...'),gsea.goCC@result$Description)
      gseaplt1CC <- gseaplot2(gsea.goCC, geneSetID = rownames(gsea.goCC[1,]), title = paste0('GO-CC PI3K/AKT Signaling Pathway\n',nm))
      ggsave(paste0(nm, '_gseaGO-CC.png'), plot = gseaplt1CC, limitsize = FALSE, width = 12, height = 8)
      ggsave(paste0(nm, '_gseaGO-CC.pdf'), plot = gseaplt1CC, limitsize = FALSE, width = 12, height = 8)
    }}
  
  
  #GSEA的KEGG富集，KEGG富集到的某一条通路的可视化：  # Loading required package: org.Hs.eg.db
  kk <- try(gseKEGG(id_geneList2, organism = 'mmu'),silent = TRUE)
  if(!is.null(kk)){
    if(!"try-error" %in% class(kk) & nrow(kk@result)>0){
      write.xlsx(kk@result, paste0(nm, '_gseaKEGG.xlsx'))
      kk@result$Description <- ifelse(nchar(kk@result$Description)>60, paste0(substr(kk@result$Description,1,60),'...'),kk@result$Description)
      gseapltkk <- gseaplot2(kk, geneSetID = rownames(kk[1,]), title = paste0('KEGG PI3K/AKT Signaling Pathway\n',nm))
      ggsave(paste0(nm, '_gseaKEGG.png'), plot = gseapltkk, limitsize = FALSE, width = 12, height = 8)
      ggsave(paste0(nm, '_gseaKEGG.pdf'), plot = gseapltkk, limitsize = FALSE, width = 12, height = 8)
    }}
  
  dd <- try(gseDO(id_geneList2),silent = TRUE)
  if(!is.null(dd) & nrow(dd)>0){
    if(!"try-error" %in% class(dd) & nrow(dd@result)>0){
      write.xlsx(dd@result, paste0(nm, '_gseaDO.xlsx'))
      dd@result$Description <- ifelse(nchar(dd@result$Description)>60, paste0(substr(dd@result$Description,1,60),'...'),dd@result$Description)
      gseapltdo <- gseaplot2(dd, geneSetID = rownames(dd[1,]), title = paste0('DO PI3K/AKT Signaling Pathway\n',nm))
      ggsave(paste0(nm, '_gseaDO.png'), plot = gseapltdo, limitsize = FALSE, width = 12, height = 8)
      ggsave(paste0(nm, '_gseaDO.pdf'), plot = gseapltdo, limitsize = FALSE, width = 12, height = 8)
    }}
  
  racPathway <- try(gsePathway(id_geneList2, organism = 'mouse'),silent = TRUE)
  if(!is.null(racPathway)){
    if(!"try-error" %in% class(racPathway) & nrow(racPathway@result)>0){
      write.xlsx(racPathway@result, paste0(nm, '_gseaReactomePathway.xlsx'))
      racPathway@result$Description <- ifelse(nchar(racPathway@result$Description)>60, paste0(substr(racPathway@result$Description,1,60),'...'),racPathway@result$Description)
      gseapltrac <- gseaplot2(racPathway, geneSetID = rownames(racPathway[1,]), title = paste0('Reactome PI3K/AKT Signaling Pathway\n',nm))
      ggsave(paste0(nm, '_gseaReactomePathway.png'), plot = gseapltrac, limitsize = FALSE, width = 12, height = 8)
      ggsave(paste0(nm, '_gseaReactomePathway.pdf'), plot = gseapltrac, limitsize = FALSE, width = 12, height = 8)
    }}
  
  
  
  print(i)
}


for (i in 14:20) {
  #Enrich(i)
  #Enrich(i,cutn = 'up')
  Enrich(i,cutn = 'down')
  print(i)
}

