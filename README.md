# 转录组分析 TranscripBio

grpTab：Treatment和Control分组关系表  

差异比较分析	处理	参考  
组合1	MDHABR	MAABR  
组合2	FDHABR	FAABR  
组合3	MDHABR	FDHABR  

grpTyp：样本分组明细表  

x	z  
FAA1BR	FAABR  
FAA1eye	FAAEYE  
FAA1LI	FAALI  
FAA1LN	FAALN  

grpVenn：Venn图分组表

	          组合n	组合n	组合n 组合n	组合n  
组合比较Venn 1	组合1	组合2	组合3	组合4	  
组合比较Venn 2	组合8	组合	组合10   



R包：pbmcapply或parallel，多核并行运算
Linux命令：whereis查找生信软件执行路径，for example: whereis samtools 返回samtools得执行路径

## File Tree

.  
├── 0.SupFile  
├── 1.QC  
│   ├── 1.Error  
│   ├── 2.GC  
│   ├── 3.Filter  
│   └── 4.Stat  
├── 2.Mapping  
│   ├── 1.Region   
│   └── 2.Stat  
├── 3.AS  
│   ├── ASlist  
│   │   ├── AS_FDHABR_vs_FAABR  
│   │   └── AS_FDHALN_vs_FAALN  
│   └── ASplot  
│       ├── FAABR_FDHABR  
│       └── FDHALN_vs_FAALN  
├── 4.SNP  
│   ├── SNP_annotation  
│   ├── SNP_site  
│   └── SNP_stat  
├── 5.Quant  
│   ├── Fig.orrelation.analysis  
│   ├── Fig.PCA.analysis  
│   ├── Fig.violin  
│   ├── FPKM  
│   └── Read_count  
├── 6.Differential  
│   ├── Fig.cluster  
│   ├── Fig.venn  
│   ├── Fig.volcano  
│   └── Venn.DEG  
│       ├── venn1  
│       ├── venn2  
│       └── venn3  
├── 7.Enrichment  
│   ├── FDHABR_vs_FAABR  
│   │   ├── GO  
│   │   └── KEGG   
│   └── FDHAEYE_vs_FAAEYE  
│       ├── GO  
│       └── KEGG  
└── 8.PPI  
    ├── FDHABR_vs_FAABR   
    └── FDHALN_vs_FAALN  








## bamSort.R
#### 使用samtools对bam文件按染色体位置进行排序，未加-n参数
#### 输出文件以'_sort.bam'后缀
##### 参考基因比对结果sam/bam文件格式说明[https://blog.csdn.net/xcaryyz/article/details/79257604]

# 1、定量分析
## 1.1、exprLevel.R
#### 包含排序功能samtools sort，同上
#### 对排序后的bam进行表达水平分析，htseq-count （htseq 0.9.1），这里是按read name排序，不会出现warning
#### 表达量列表 _count.txt

## 1.2、sampleCorr.R
#### 样本间相关性分析，画相关矩阵图

# 2、FPKM定量
## 2.1、FPKM.R
stringtie 1.3.3b
##### -A 输出基因丰度值（FPKM, TPM），这里输入的是按染色体位置进行排序的bam，stringtie只认按Pos排序
##### -l 新组装的转录本标签（默认STRG）
##### -o 输出新组装的gtf文件

## 2.3、fpkmDistribPlot.R
#### 合并fpkm
#### fpkm分布图（箱图和小提琴图）

# 3、差异分析 Read Counts
## DEGs.R
DESeq2 1.26.0
#### 差异基因筛选 padj < 0.05 & abs(log2FoldChange) > 0
#### 合并read counts
#### 输出火山图
#### 实验组基因表达量列表 _diffGenles.xlsx
#### 实验组差异基因表达量列表 _counts.xlsx

## 2.4、pcaPlot.R
#### 主成分分析

VennPlot.R
### 韦恩图
heatmapCluster.R
### 差异基因聚类，聚类热图

# 4、富集分析
clusterProfiler 3.14.0
### 根据差异基因进行富集分析
#### org.Hu.eg.db 老鼠基因库
#### org.Mm.eg.db 人类基因库
#### GO、KEGG、GASE、DO、Reactome


# 5、可变剪接
## alterSplic.R
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


# 6、变异分析
GATK
snpEff
## 6.1、SNP-loop.R
### 变异位点SNP与INDEL检测
##### (6.1) 标记比对文件中的pcr重复
##### (6.2) 将落在内含子区间的 reads 片段直接切除，并对 MAPQ 进行调整
##### (6.3) 将落在内含子区间的 reads 片段直接切除，并对 MAPQ 进行调整
##### (6.4) 检测snp和indel
##### (6.5)snp和indel过滤

### 变异位点统计图（待完善）

## 6.2、SNPann.R
#### 变异基因注释





