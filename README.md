# TranscripBio
转录组分析  
grpTab：样本分组明细表  
grpTyp：Treatment和Control分组关系表  

## bamSort.R
#### 使用samtools对bam文件按染色体位置进行排序
#### 输出文件以'-sort.bam'后缀
##### 参考基因比对结果sam/bam文件格式说明[https://blog.csdn.net/xcaryyz/article/details/79257604]

# 1、定量分析
## 1.1、exprLevel.R
#### 包含排序功能samtools sort，同上
#### 对排序后的bam进行表达水平分析，htseq-count （htseq 0.9.1）
#### 表达量 _count.txt


# 2、FPKM定量
## 2.1、FPKM.R
stringtie 1.3.3b
##### -A 输出基因丰度值（FPKM, TPM）
##### -l 新组装的转录本标签（默认STRG）
##### -o 输出新组装的gtf文件

## 2.2、sampleCorr.R
#### 样本间相关性分析，画相关矩阵图

## 2.3、fpkmDistribPlot.R
#### fpkm分布图（箱图和小提琴图）

#### 主成分分析（待完善）


# 3、差异分析
## DEGs.R
DESeq2 1.26.0
#### 差异基因筛选 padj < 0.05 & abs(log2FoldChange) > 0
#### 输出画火山图
#### 样本基因表达量列表 _diffGenles.xlsx
#### 样本差异基因表达量列表 _counts.xlsx

### 韦恩图（待完善）
### 差异基因聚类，聚类热图（待完善）

# 4、富集分析
clusterProfiler 3.14.0
### 根据差异基因进行富集分析
#### org.Hu.eg.db 老鼠基因u库
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

#### 变异位点统计图（待完善）

## 6.2、SNPann.R
#### 变异基因注释





