
path <- '/home/origene/zhengguosong/anaconda3/envs/rna/bin'
ensemble_index <- "/home/data/Rseq/wangyangxian/genome/Mus_musculus_Ensemble_94_index"
fastq_dir <- "/home/data/Rseq/wangyangxian38"
outfile_dir <- "/home/data/Rseq/wangyangxian38/bam"
dir(fastq_dir, pattern = '.fq.gz')

filename=$(cat hisatfile.txt)
for i in ${filename}; do
hisat2 -x /home/data/Rseq/wangyangxian/genome/Mus_musculus_Ensemble_94_index -p 15 -1 ../${i}_1.clean.fq.gz -2 ../${i}_2.clean.fq.gz -S ./${i}.sam;
samtools view --threads 15 -b ${i}.sam > ${i}.bam;
rm ${i}.sam;
samtools sort --threads 15 -n -o ${i}_sort.bam ${i}.bam;
rm ${i}.bam;
done


## hisat2 -x <ht2-idx>  Index filename prefix（参考基因组的索引）
## hisat2 -p threads 设置线程数
## hisat2 -1, -2 输入cleandata1,输入cleandata2
## hisat2 -S 输出sam文件
## samtools view 将sam文件转换为bam文件
## samtools --threads 设置线程数
## samtools -b 输出文件格式为bam文件格式
## samtools sort bam文件排序
## samtools -n 按name排序
## samtools -o 输出文件的文件名