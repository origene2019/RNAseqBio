if (!requireNamespace("BiocManager", quietly = TRUE))install.packages("BiocManager")
library(BiocManager)
pkgs <- rownames(installed.packages())

deg_pkgs <- c("DESeq2","biomaRt")
BiocManager::install(deg_pkgs[!deg_pkgs %in% pkgs])

pkgs <- rownames(installed.packages())
enr_pkgs <- c("clusterProfiler","org.Mm.eg.db","GSEABase","GEOquery","DOSE","topGO")
BiocManager::install(enr_pkgs[!enr_pkgs %in% pkgs])

library(DESeq2)
library(biomaRt)

library(clusterProfiler)
library(org.Mm.eg.db)
library(GSEABase)
library(GEOquery)
library(DOSE)
library(topGO)

if(!require(devtools)) install.packages("devtools")







