---
title: RNA-seq call variants的GATK4分析流程
description: 自己亲身实践的数据分析流程，希望能对大家有所帮助
categories:
 - Bioimformations
tags: gatk4, RNA-seq
---
# GATK4：RNA-seq call variants流程
### 软件安装及数据准备
##### 1.准备数据
- 建立分析过程中所需要的文件目录，并将**自己需要分析的数据**建立软链接。这里我以自己做过的绿盘鲍RNA-seq数据作为示例，展示完整的从比对到call snp的流程 
 
```sh
mkdir -p ${HOME}/program/heterosis/{bin,data/{genes,genome,index,samples},ASE_output/{align,align_2pass,index_2pass,pretreat,snp}}

#给数据创建软链接
ln -s /public/home/benthic/Data_benthic/Hfu/Hfu_genome/genome_change.fa ${HOME}/program/heterosis/data/genome/Hfu.fa
cp /public/home/benthic/Data_benthic/Hfu/Hfu_genome/HFU_genome_change.gff ${HOME}/program/heterosis/data/genes/Hfu.gff
gffread Hfu.gff -T -o Hfu.gtf #因为基因组里是gff文件，需要转换，为了不破坏原文件，所以这里我把它复制过来转换了
for i in `find /public/home/benthic/Data_benthic/Hfu/DF_mRNA_immunity/NHT150898_TR/clean_data/ -type f -name "*.fq.gz"`;do ln -s ${i} ${HOME}/program/heterosis/data/samples/${i##/*};done
```

##### 2.安装软件
- 实验室集群已经安装了GATK4和STAR，但是如果需要自己安装的话，可以参考<a href="https://github.com/alexdobin/STAR" target="_blank">STAR官网</a>和<a href="https://software.broadinstitute.org/gatk/" target="_blank">GATK官网</a>

### 数据比对
##### 1.STAR第一次比对
- 比对的第一步就是建立索引

```sh
STAR=$(which STAR)
$STAR --runThreadN 24 --runMode genomeGenerate --genomeDir /public/home/benthic/wangxz/program/immu_mRNA_heterosis/data/index \
--genomeFastaFiles /public/home/benthic/wangxz/program/immu_mRNA_heterosis/data/genome/Hfu.fa \
--sjdbGTFfile /public/home/benthic/wangxz/program/immu_mRNA_heterosis/data/genes/Hfu.gtf \
--sjdbOverhang 100
```

- 接下来进入第一次比对

```sh
#配置路径及变量
STAR=$(which STAR)
BASEDIR="${HOME}/program/immu_mRNA_heterosis"
WORKDIR="${BASEDIR}/ASE_output_Hfu"
FASTQLOC="${BASEDIR}/data/samples"
GENOMEIDX="${BASEDIR}/data/index"
#获取样品列表,双端测序必须以“_1.*”"_2.*"为命名规则
reads1=(${FASTQLOC}/*_1.*)
reads1=("${reads1[@]##*/}")
reads2=("${reads1[@]/_1./_2.}")

#切换到工作目录，开始比对
cd $WORKDIR
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: STAR alignment"
for ((i=0; i<=${#reads1[@]}-1; i++ )); do
    sample="${reads1[$i]%%.*}"
    sample="${sample%_*}"

    $STAR --runThreadN 24 --genomeDir ${BASEDIR}/data/index \
    --readFilesIn ${FASTQLOC}/${reads1[$i]} ${FASTQLOC}/${reads2[$i]} \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesCommand zcat \     #这里要注意，如果不加这个参数是没有办法对压缩文件进行比对的
    --outFileNamePrefix $WORKDIR/align/${sample}
done
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> END: STAR alignment"
```

- 结果文件
这里的“SRR3589959”是前面比对时，`--outFileNamePrefix`参数设置的文件头，这里用的是网上的一个例子
>SRR3589959Aligned.sortedByCoord.out.bam  
>SRR3589959Log.final.out  
>SRR3589959Log.out  
>SRR3589959Log.progress.out  
>SRR3589959SJ.out.tab  

##### 2.STAR的二次比对
- 重新建立索引  
第一遍常规比对之后会生成一个SJ.out.tab文件，如上面所提到的SRR3589959SJ.out.tab。第二次比对之前，就需要使用`--sjdbFileChrStartEnd`参数将所有样品的SJ.out.tab文件作为输入的annotated junction进行第二次建index

```sh
STAR=$(which STAR)
BASEDIR="$HOME/program/immu_mRNA_heterosis"
WORKDIR="$BASEDIR/ASE_output_Hfu"
FASTQLOC="$BASEDIR/data/samples"
INDEX2PASS=`ls $WORKDIR/align/*SJ.out.tab`

$STAR --runThreadN 24 --runMode genomeGenerate --genomeDir $WORKDIR/index_2pass \
--genomeFastaFiles $BASEDIR/data/genome/Hfu.fa \
--sjdbGTFfile $BASEDIR/data/genes/Hfu.gtf \
--sjdbFileChrStartEnd $INDEX2PASS   #这里是利用了前面结果里所有的SJ.out.tab文件
--sjdbOverhang 100
```

- 第二次比对

```sh
STAR=$(which STAR)
BASEDIR="$HOME/program/immu_mRNA_heterosis"
WORKDIR="$BASEDIR/ASE_output_Hfu"
FASTQLOC="$BASEDIR/data/samples"
INDEX2PASS=`ls $WORKDIR/align/*SJ.out.tab`

#获取样品列表,双端测序必须以“_1.*”"_2.*"为命名规则
reads1=(${FASTQLOC}/*_1.*)
reads1=("${reads1[@]##*/}")
reads2=("${reads1[@]/_1./_2.}")
#切换到工作目录，开始比对
cd $WORKDIR
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: STAR alignment"
for ((i=0; i<=${#reads1[@]}-1; i++ )); do
    sample="${reads1[$i]%%.*}"
    sample="${sample%_*}"

    $STAR --runThreadN 24 --genomeDir ${WORKDIR}/index_2pass \
    --readFilesIn ${FASTQLOC}/${reads1[$i]} ${FASTQLOC}/${reads2[$i]} \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesCommand zcat \
    --outFileNamePrefix $WORKDIR/align_2pass/${sample}
done
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> END: STAR alignment"
```

### 找到SNP位点
##### GATK4语法  
GATK4相对于3来说，改变了很多，整合了许多软件，因此在使用时的语法上也不同。下面是一个GATK4的简单语法演示  
`gatk [--java-options "jvm args like -Xmx4G go here"] ToolName [GATK args go here]`  
来一个具体的例子  
`gatk --java-options "-Xmx8G" HaplotypeCaller -R reference.fasta -I input.bam -O output.vcf`  
使用GATK3的同学自己上GATK官网查一下语法，还是有很多人用3的。下面就进入正式的call snp环节  

##### 1.bam文件的预处理
- 加上GATK需要的标签  
因为GATK运行时，对文件的格式要求和STAR比对出来的有些不太一样，所以GATK有自己的程序可以给比对过的文件加上标签，使其符合自己的程序要求。

```sh
gatk AddOrReplaceReadGroups -I your_sample.bam -O your_sample_addsorted.bam \
    -SO coordinate -ID your_sample_name -LB rna -PL illumina -PU hiseq -SM your_sample_name
```

- 去除重复的序列

```sh
gatk MarkDuplicates -I your_sample_addsorted.bam -O your_sample_dep.bam \
    --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M your_sample.metrics
```

- GATK使用专门给RNA-Seq开发的SplitNCigarReads工具，将落在外显子上的reads分离出来同时去除N错误碱基（getting rid of Ns but maintaining grouping information），并去除掉落在内含子区域的reads，以减少假阳性变异；*这里注意，在GATK4的版本中，质量标准的转换已经包含在默认参数当中，因此不需要额外指定参数。*

```sh
$GATK SplitNCigarReads -I your_sample_dep.bam -O your_sample_split.bam \
    -R Hfu.fa
```

为了方便我在运行预处理的时候，只写了一个脚本，把三个命令全部贴进去了

```sh
BASEDIR="$HOME/program/immu_mRNA_heterosis"
WORKDIR="${BASEDIR}/ASE_output_Hfu"
BAMLOC="$WORKDIR/align_2pass"
GATK=$(which gatk)
GENOME="$BASEDIR/data/genome/Hfu.fa"

cd $WORKDIR
for i in `ls $BAMLOC/*.bam`;do
    sample=$(basename $i Aligned.sortedByCoord.out.bam)
    $GATK AddOrReplaceReadGroups -I $i -O $WORKDIR/pretreat/${sample}_addsorted.bam \
    -SO coordinate -ID $sample -LB rna -PL illumina -PU hiseq -SM $sample
    $GATK MarkDuplicates -I $WORKDIR/pretreat/${sample}_addsorted.bam -O $WORKDIR/pretreat/${sample}_dep.bam \
    --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M $WORKDIR/pretreat/${sample}.metrics
    $GATK SplitNCigarReads -I $WORKDIR/pretreat/${sample}_dep.bam -O $WORKDIR/pretreat/${sample}_split.bam \
    -R $GENOME
done
```

##### 2.建立参考基因组索引
在进行SNP位点识别之前，需要对参考基因组建立索引

```sh
BASEDIR="$HOME/program/immu_mRNA_heterosis"
GENOME="$BASEDIR/data/genome/Hfu.fa"

cd $BASEDIR/data/genome
samtools faidx Hfu.fa    # 该命令会在Hfu.fa所在目录下创建一个Hfu.fai索引文件
gatk CreateSequenceDictionary -R Hfu.fa -O Hfu.dict
```

##### 3.变异检测
因为我们缺少已知的snp数据库，所以**碱基质量分数重校准（Base quality score recalibration，BQSR)**这一步我们就跳过了，如果想要做这一步的话可以参考<https://ming-lian.github.io/2019/02/08/call-snp>。貌似可以利用samtools和GATK做一个snp数据库  
对于多样本RNA-seq来说，最好可以先对每个样本生成后续分析需要的中间文件GVCF，然后再把所有样本的文件联合起来做joint calling。因此这里我给出的流程是多样本的情况
- 生成中间文件GVCF
对每个样本，都需要生成一个GVCF文件，在进群上需要提交N(=你的样本数)次脚本，所以我生成了一个python脚本用来批量提交类似的pbs脚本【我真的好懒】  
  
先把python脚本提供给大家，这个脚本是通用的，你只要替换中间的命令部分就可以实现linux上所有命令的批量操作啦，详细的参数解释我也写在下面了

```python
#!/public/bioconda/bin/python
import subprocess    #这个是调用外部命令的包
import time    #这个是用来加载时间，下面的程序在每提交一个命令之后会休眠一会，避免任务提交太快

param_file = './job_params.txt'    #这个文件是用来储存每次命令需要更换的部分
'''
比如，这里的程序是做变异检测的，每次需要换的部分就是我输入的样本名，所以我把样本名全都写在这个文件里
job_params.txt
C_D_1
C_D_2
......
F_24h_1
F_24h_2
'''

with open(param_file,'r') as f:
    for line in f:
        sample = line.strip()    #把每行中的变量名赋值给sample这个变量
        qsub_cmd = 'qsub -N {0} -e ./snp/{0}.err -o ./snp/{0}.log -v SAMPLE={0} GATK_callsnp_high.pbs'.format(sample)    
            #.format命令是用来格式化输出的，我就不讲解了。
            # qsub里的命令，-N 指定任务名，-e -o 错误和日志的输出文件，
            #-v 把外部的变量名导入pbs脚本里[这个变量是最重要的，所有你需要在外部定义的变量都需要在这里导入，比如从外部导入每个样本名]
        exit_status = subprocess.call(qsub_cmd, shell=True)    #执行上面的命令
        if exit_status is 1:
            print('Job "{}" failed to submit'.format(qsub_cwd))    #如果提交失败就标准输出错误到屏幕
        time.sleep(1)    #暂停一段时间再提交

print("Done submitting jobs!")
```

下面就是生成GVCF文件的pbs脚本

```sh
#!/bin/bash -x
#PBS -q low
#PBS -j oe
#PBS -l walltime=999:00:00
#PBS -l nodes=1:ppn=12

GATK=$(which gatk)
BASEDIR="$HOME/program/immu_mRNA_heterosis"
WORKDIR="$BASEDIR/ASE_output_Hfu"
GENOME="$BASEDIR/data/genome/Hfu.fa"

cd $WORKDIR
$GATK HaplotypeCaller -R $GENOME -I ${WORKDIR}/pretreat/${SAMPLE}_split.bam \
  --dont-use-soft-clipped-bases -stand-call-conf 20.0 \
  -O ${WORKDIR}/snp/${SAMPLE}.g.vcf -ERC GVCF
```

- 合并多个GVCF文件得到GenomicsDB  
因为这一步的合并每次只能对单个染色体做，所以这里就又利用了批量提交的脚本，对每个染色体分别提交程序。这里只给出pbs脚本，批量提交还是上面的方法

```sh
#!/bin/bash -x
#PBS -q low
#PBS -j oe
#PBS -l walltime=999:00:00
#PBS -l nodes=1:ppn=28

BASEDIR="$HOME/program/immu_mRNA_heterosis"
WORKDIR="${BASEDIR}/ASE_output_Hfu/snp"
GATK=$(which gatk)
GENOME="$BASEDIR/data/genome/Hfu.fa"

cd $WORKDIR
$GATK --java-options "-Xmx4g -Xms4g" GenomicsDBImport --genomicsdb-workspace-path ${CHR}db \
  --reader-threads 28 --sample-name-map sample_map -L ${CHR}    #这里的CHR是用外部变量来替换的，代表的是染色体的编号，比如chr1
```

- 联合检测 joint calling  
对每个染色体分别生成最终的vcf文件

```sh
#!/bin/bash -x
#PBS -q low
#PBS -j oe
#PBS -l walltime=999:00:00
#PBS -l nodes=1:ppn=28

BASEDIR="$HOME/program/immu_mRNA_heterosis"
WORKDIR="${BASEDIR}/ASE_output_Hfu/snp"
GATK=$(which gatk)
GENOME="$BASEDIR/data/genome/Hfu.fa"

cd $WORKDIR
for bed in chr{1..18}
do
$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" GenotypeGVCFs -R $GENOME \
  -V gendb://${bed}db -O final_${bed}.vcf
done
```

- 变异位点过滤 Variant filtering 
 
```sh
#!/bin/bash -x
#PBS -q low
#PBS -N GATK_filter
#PBS -j oe
#PBS -e GATK_filter.err
#PBS -o GATK_filter.log
#PBS -l walltime=999:00:00
#PBS -l nodes=1:ppn=28

BASEDIR="$HOME/program/immu_mRNA_heterosis"
WORKDIR="${BASEDIR}/ASE_output_Hfu/snp"
GATK=$(which gatk)
GENOME="$BASEDIR/data/genome/Hfu.fa"

cd $WORKDIR
$GATK SelectVariants -R $GENOME -V merge.vcf \
  --select-type-to-include SNP -O filter/merge_snp.vcf    #提取所有的snp位点

$GATK SelectVariants -R $GENOME -V merge.vcf \
  --select-type-to-include INDEL -O filter/merge_indel.vcf    #提取所有的indel位点

$GATK VariantFiltration -R $GENOME -V filter/merge_snp.vcf \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"\
  --filter-name "SNP_FILTER" -O filter/merge_filter_snp.vcf    #根据GATK官网推荐的筛选标准对SNP位点进行过滤

$GATK VariantFiltration -R $GENOME -V filter/merge_indel.vcf \
  --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -20.0"\
  --filter-name "INDEL_FILTER" -O filter/merge_filter_indel.vcf    #根据GATK推荐的标准对indel位点进行过滤

$GATK MergeVcfs -I filter/merge_filter_snp.vcf -I filter/merge_filter_indel.vcf \
  -O filter/merge_filter.vcf    #将过滤好的snp位点和indel位点结合在一起

$GATK SelectVariants -R $GENOME -V filter/merge_filter.vcf -O filter/merge_pass.vcf -select "vc.isNotFiltered()"    #挑出所有通筛选标准的位点
```

### 总结
到此为止，snp calling的步骤就结束了，剩下的分析就是根据自己的需要来对vcf文件再加工，或者提取出想要的信息。