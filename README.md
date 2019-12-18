# 转录组分析 TranscripBio

grpTab：Treatment和Control分组关系表  

差异比较分析 | 处理 | 参考  
-------- | ----- | -----
组合1 | MDHABR | MAABR  
组合2 | FDHABR | FAABR  
组合3 | MDHABR | FDHABR  

grpTyp：样本分组明细表  

x | z  
----- | -----
FAA1BR | FAABR  
FAA1eye | FAAEYE  
FAA1LI | FAALI  
FAA1LN | FAALN  

grpVenn：Venn图分组表

组合n | 组合n | 组合n | 组合n | 组合n    
 -------- | ----- | ----- | ----- | -----  
组合比较Venn 1 | 组合1 | 组合2 | 组合3 | 组合4	  
组合比较Venn 2 | 组合8 | 组合9 | 组合10 | 组合11



R包：pbmcapply或parallel，多核并行运算  
Linux命令：whereis查找生信软件执行路径，for example: whereis samtools 返回samtools得执行路径  

## File Tree

resultsAll  
├── alterSplic  
│   ├── AS_OVA_vs_Control  
│   │   ├── A5SS.MATS.JCEC.txt  
│   │   ├── RI.MATS.JCEC.txt  
│   │   ├── SE.MATS.JCEC.txt  
│   │   ├── sign20.SE.MATS.JCEC.txt  
│   │   └── sign.SE.MATS.JCEC.txt  
│   └── ASplot  
│       └── AS_OVA_vs_Control  
│           ├── Sashimi_index  
│           │   └── SE.event.list.txt  
│           │       ├── 11_Snhg1_19_8725014_8725046_+@19_8725014_8725256_+@19_8725222_8725256_+.pdf  
│           │       └── 11_Snhg1_19_8725014_8725046_+@19_8725014_8725256_+@19_8725222_8725256_+.png  
│           └── SE  
│               └── Sashimi_plot  
│                   ├── 9_Ubox5_2_130629916_130630027_-@2_130614889_130616742_-@2_130604485_130604626_-.pdf  
│                   └── 9_Ubox5_2_130629916_130630027_-@2_130614889_130616742_-@2_130604485_130604626_-.png  
├── corrPlot  
│   ├── O_NS398_vs_O_DIC_corr.png  
│   ├── O_NS398_vs_OVA_corr.png  
│   └── OVA_vs_Control_corr.png  
├── counts  
│   ├── OVA2Reye_counts.txt  
│   ├── OVA3Leye_counts.txt   
│   └── OVA3Reye_counts.txt  
├── Differential  
│   │   └── OVA_vs_Control_up  
│   │       ├── GO  
│   │       ├── GSEA  
│   │       ├── KEGG  
│   │       └── Reactome  
│   └── volcano  
│       ├── O_NS398_vs_OVA_DEseq.png  
│       └── OVA_vs_Control_DEseq.png  
├── fpkm  
│   ├── cluster  
│   │   ├── OVA_vs_Control_heatmap.pdf  
│   │   └── OVA_vs_Control_heatmap.png  
│   └── violinPlot  
│       ├── O_NS398_vs_OVA_FPKM.png  
│       └── OVA_vs_Control_FPKM.png  
├── grpTab.txt  
├── grpTyp.txt  
├── grpVenn.txt  
├── map  
│   ├── align_region.xlsx  
│   ├── OVA3Reye_mapRegion.pdf  
│   └── OVA3Reye_mapRegion.png  
├── pca  
│   ├── OVA_vs_Control_pca.pdf  
│   └── OVA_vs_Control_pca.png  
├── PPI  
│   └── OVA_vs_Control  
│       ├── OVA_vs_Control_all.tsv  
│       ├── OVA_vs_Control_diffGenes.xlsx  
│       └── OVA_vs_Control_up.tsv  
├── SNP  
│   └── snpPlot  
├── sortBam  
├── tree.txt  
└── venn  
    ├── Venn1  
    │   ├── 7_O_10526_vs_OVA~O_NS398_vs_OVA.xlsx  
    │   └── 8_O_DIC_vs_OVA.xlsx  
    ├── Venn1.pdf  
    ├── Venn1.png  
    ├── Venn2  
    │   ├── 1_C_10526_vs_Control.xlsx  
    │   └── 5_C_DIC_vs_Control.xlsx  
    ├── Venn2.pdf  
    └── Venn2.png  


#### 需要使用的生信软件，使用conda进行安装及管理 conda install，conda update, conda uninstall  
#### 建议创建单独的虚拟运行环境 conda create -n rna, conda activate rna, conda deactivate rna  
Anaconda官网<https://www.anaconda.com/>  
conda manager 查找软件安装<https://anaconda.org/continuumcrew/conda-manager>  
Linux命令：  
top 查看当前运行的进程信息   
whereis 查找软件路径   


## bamSort.R
#### 使用samtools对bam文件按染色体位置进行排序，未加-n参数
#### 输出文件以'_sort.bam'后缀
##### 参考基因比对结果sam/bam文件格式说明[https://blog.csdn.net/xcaryyz/article/details/79257604]

# 序列比对（RseqMap.R）
## 比对结果统计汇总 mapStat.R
## 比对区域统计汇总 mapRegStat.R


# 数据质控  
mapQC.R   
mapStat.R  
mapRegStat.R  

# 1、定量分析
## 1.1、exprLevel.R
#### 包含排序功能samtools sort，同上
#### 对排序后的bam进行表达水平分析，htseq-count （htseq 0.9.1），这里是按read name排序，不会出现warning
#### 表达量列表 _count.txt

# 2、差异分析 Read Counts
## sampleCorr.R（含读取count值） 或corrPlot.R（已读取count值）
#### 样本间相关性分析，画相关矩阵图

## DEGs.R
DESeq2 1.26.0
#### 差异基因筛选 padj < 0.05 & abs(log2FoldChange) > 0
#### 合并read counts
#### 输出火山图
#### 实验组基因表达量列表 _diffGenles.xlsx
#### 实验组差异基因表达量列表 _counts.xlsx

VennPlot.R
### 韦恩图，差异基因diff count值



# 3、FPKM定量
## 2.1、FPKM.R
stringtie 1.3.3b
##### -A 输出基因丰度值（FPKM, TPM），这里输入的是按染色体位置进行排序的bam，stringtie只认按Pos排序
##### -l 新组装的转录本标签（默认STRG）
##### -o 输出新组装的gtf文件

## 2.3、fpkmDistribPlot.R
#### 合并fpkm
#### fpkm分布图（箱图和小提琴图）

## 2.4、pcaPlot.R
#### 主成分分析 fpkm值

heatmapCluster.R 
### 差异基因聚类，聚类热图，差异基因fpkm值


# 4、富集分析
clusterProfiler 3.14.0
### 根据差异基因进行富集分析
#### org.Hu.eg.db 人类基因库
#### org.Mm.eg.db 老鼠基因库
#### GO、KEGG、GASE、DO、Reactome


# 5、可变剪接
## 5.1 alterSplic.R
rMATS 4.0.2
##### --s1 s1.txt  txt输入文件，包含了Sample1的fastq序列文件
##### --s2 s2.txt  txt输入文件，包含了Sample2的fastq序列文件
##### --b1 b1.txt  txt输入文件，包含了Sample1的bam序列文件
##### --b2 b2.txt  txt输入文件，包含了Sample2的bam序列文件
##### --t  readType RNASeq数据类型，'paired'表示双端数据，'single'表示单端数据，这里采用<paired双端数据>
##### --readLength <int> 测序read长度，<150>
##### --gtf gtfFile gtf注释文件
##### --nthread 线程数 5
##### --od  outDir  输出路径
 
## 5.3 可变剪切可视化
### ESplot.R
####	rmats2sashimiplot（pyhton2环境）需要调用samtools,使用时特别需要注意把python2环境下的samtools升级到最新版本，否则在调用samtools时会调用失败，  
####	原因是自动安装的版本太低，服务器没有它需要的libxxx.so.0.0.0库。  
####	升级python2中samtools的方法：  
conda activate python2  
conda update samtools  
conda install -c bioconda samtools openssl=1.0  

####	rmats2sashimiplot输出的文件为pdf格式，可以用ImageMagick将pdf格式转换为png格式的图片
####	安装好imageMagick后，运行以下代码
convert -density 400 image.pdf imag.png
####	-density 调整图片分辨率

# 6、变异分析
GATK
snpEff
## 6.1、SNP-loop.R
### 变异位点SNP与INDEL检测
##### (6.1) 标记比对文件中的pcr重复，输入bam文件按染色体位置pos排序
##### (6.2) 将落在内含子区间的 reads 片段直接切除，并对 MAPQ 进行调整
##### (6.3) 将落在内含子区间的 reads 片段直接切除，并对 MAPQ 进行调整
##### (6.4) 检测snp和indel
##### (6.5)snp和indel过滤

### 变异位点统计图（待完善）

## 6.2、SNPann.R
### SNP-loop.R  
### 变异基因注释 SNPann.R  
### SNP结果汇总，画统计图 SNPstatPlot.R  




# 7、融合基因  
### fusion.R  
### 融合基因圈图 circoPlot.R  
