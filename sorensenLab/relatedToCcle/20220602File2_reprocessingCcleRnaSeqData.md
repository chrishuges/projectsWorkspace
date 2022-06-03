## Reprocessing some RNAseq data

This document describes the reprocessing of some RNAseq data from a the CCLE database, or DepMap. The data are all located in SRA: PRJNA523380.

### Description

I am going to use [sraDownloader](https://github.com/s-andrews/sradownloader) to get at these data.

I will parse these downloaded raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/).

I will use [STAR](https://github.com/alexdobin/STAR) for the alignment. I am following the instructions in the STAR [manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf).

I will use [Salmon](https://github.com/COMBINE-lab/salmon) to do read quantification. There is a great info page for this tool [here](https://salmon.readthedocs.io/en/latest/).

I will then use [samtools](http://www.htslib.org/) to prepare the files before using [FeatureCounts](http://subread.sourceforge.net/featureCounts.html) (part of [Subread](http://subread.sourceforge.net/)) for quantification.

This is kind of a useful website for a general pipeline, [here](https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html), also [here](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html).

### Data pipeline

First we will move into our working directory and create our shell and snakemake scripts.

```shell
cd /projects/ptx_results/Sequencing/publishedStudies/sequencing20220602_ccleDataReprocessing
touch sraDataProcessingScript.sh
chmod +x sraDataProcessingScript.sh
touch snakefile
```

Edit the contents of the snakefile to include the text below. I use vim for this.

```python
"""
Author: Christopher Hughes
Affiliation: BCCRC
Aim: Workflow for RNA-Seq data
Date: 20211105
"""

###############################
#locations of tools we will use
#BBDUK = "/home/chughes/softwareTools/bbmap-38.90/bbduk.sh"
BBDUK = "/projects/ptx_analysis/chughes/softwareTools/bbmap-38.90/bbduk.sh"
#SALMON = "/home/chughes/softwareTools/salmon-1.5.2/bin/salmon"
SALMON = "/projects/ptx_analysis/chughes/softwareTools/salmon-1.5.2/bin/salmon"
#HISAT2 = "/home/chughes/softwareTools/hisat2-2.2.1/hisat2"
#STAR = "/home/chughes/softwareTools/STAR-2.7.9a/bin/Linux_x86_64/STAR"
STAR = "/projects/ptx_analysis/chughes/softwareTools/STAR-2.7.9a/bin/Linux_x86_64/STAR"
#SAMTOOLS="/home/chughes/softwareTools/samtools-1.12/samtools"
SAMTOOLS="/gsc/software/linux-x86_64-centos7/samtools-1.14/bin/samtools"
#SAMBAMBA="/home/chughes/softwareTools/sambamba-0.8.1/sambamba"
#FEATURECOUNTS="/home/chughes/softwareTools/subread-2.0.3/bin/featureCounts"
FEATURECOUNTS="/projects/ptx_analysis/chughes/softwareTools/subread-2.0.3/bin/featureCounts"
#BAMCOVERAGE="/home/chughes/virtualPython368/bin/bamCoverage"
BAMCOVERAGE="/home/chughes/Virtual_Python383/bin/bamCoverage"


###############################
#locations of our index files
#DATABASE_DIR = "/home/chughes/databases/projectEwsDlg2"
DATABASE_DIR = "/projects/ptx_analysis/chughes/databases/projectEwsDlg2"
STARINDEX = DATABASE_DIR + "/starIndex"
SALMONINDEX = DATABASE_DIR + "/salmonIndex/salmon_index"
GTF = DATABASE_DIR + "/baseGenomeFiles/genome.gtf"
FASTA = DATABASE_DIR + "/baseGenomeFiles/genome.fa"


###############################
#get our list of sample files
SAMPLES, = glob_wildcards("raw/{smp}_1.fastq.gz")
print("There are " + str(len(SAMPLES)) + " total samples to be processed.")
for smp in SAMPLES:
  print("Sample " + smp + " will be processed")


###############################
#processing workflow
rule all:
    input:
      expand("results/{smp}.counts.txt", smp = SAMPLES),
      expand("results/{smp}.sorted.bam.bai", smp = SAMPLES),
      expand("results/{smp}.sorted.bw", smp = SAMPLES),
      expand("quants/{smp}/quant.sf", smp = SAMPLES)

rule bbduk:
  input:
      r1 = "raw/{smp}_1.fastq.gz",
      r2 = "raw/{smp}_2.fastq.gz"
  output:
      ro1 = "results/{smp}_1.clean.fastq.gz",
      ro2 = "results/{smp}_2.clean.fastq.gz"
  message:
      "Processing with BBDuk."
  shell:
      "{BBDUK} in1={input.r1} in2={input.r2} ref=adapters out1={output.ro1} out2={output.ro2} ktrim=r k=23 mink=11 hdist=1 tpe tbo"

rule star:
  input:
      r1 = "results/{smp}_1.clean.fastq.gz",
      r2 = "results/{smp}_2.clean.fastq.gz"
  output:
      "results/{smp}.sorted.bam"
  message:
      "Aligning with STAR."
  shell:
      "{STAR} --runThreadN 16 --genomeDir {STARINDEX} --readFilesIn {input.r1} {input.r2} --readFilesCommand zcat --sjdbGTFfile {GTF} --outStd SAM | {SAMTOOLS} sort -o {output}"

rule bam_indexing:
  input:
      "results/{smp}.sorted.bam"
  output:
      "results/{smp}.sorted.bam.bai"
  message:
      "Indexing BAM with samtools."
  shell:
      "{SAMTOOLS} index {input}"

rule bam_coverage:
  input:
      r1 = "results/{smp}.sorted.bam",
      w2 = "results/{smp}.sorted.bam.bai"
  output:
      "results/{smp}.sorted.bw"
  message:
      "Calculating coverage with deeptools."
  shell:
      "{BAMCOVERAGE} -b {input.r1} -o {output} -p 16 --normalizeUsing BPM"

rule featurecounts:
  input:
      r1 = "results/{smp}.sorted.bam",
      w2 = "results/{smp}.sorted.bam.bai",
      w3 = "results/{smp}.sorted.bw"
  output:
      "results/{smp}.counts.txt"
  message:
      "Counting reads with featureCounts."
  shell:
      "{FEATURECOUNTS} -p --countReadPairs -t exon -g gene_id -a {GTF} -o {output} {input.r1}"

rule salmon:
  input:
      r1 = "results/{smp}_1.clean.fastq.gz",
      r2 = "results/{smp}_2.clean.fastq.gz",
      w3 = "results/{smp}.counts.txt"
  output:
      "quants/{smp}/quant.sf"
  params:
      dir = "quants/{smp}"
  message:
      "Quantifying with salmon."
  shell:
      "{SALMON} quant -i {SALMONINDEX} -l A -p 16 --gcBias --validateMappings -o {params.dir} -1 {input.r1} -2 {input.r2}"
```

Below is the shell script I will use to process these data with snakemake.

```shell
#!/bin/bash

##set the location of software tools and the working directory where files will be stored
##I have two sets of locations here because sometimes I use different servers
#sraDownloader="/home/chughes/softwareTools/sradownloader-3.8/sradownloader"
sraDownloader="/projects/ptx_analysis/chughes/softwareTools/sradownloader-3.8/sradownloader"
#sraCacheLocation="/mnt/Data/chughes/sratoolsRepository"
sraCacheLocation="/projects/ptx_results/Sequencing/sraCache"
#workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211112_ccleSequencingDataReprocessing"
workingDirectory="/projects/ptx_results/Sequencing/publishedStudies/sequencing20220602_ccleDataReprocessing"
eval cd ${workingDirectory}
eval mkdir raw
eval mkdir results
eval mkdir quants

##loop over the accessions
for i in SRR8615265 SRR8615420 SRR8615854 SRR8618314 SRR8615264 SRR8616193 SRR8615630 SRR8616192 SRR8615981 SRR8616059 SRR8615578 SRR8616090 SRR8615774 SRR8615902 SRR8615264 SRR8615437 SRR8615993 SRR8615846 SRR8618317 SRR8616130 SRR8616129 SRR8615682 SRR8618314 SRR8615355 SRR8615414 SRR8615507 SRR8615270 SRR8615559 SRR8615359 SRR8615558 SRR8615949 SRR8616076 SRR8616020 SRR8615637 SRR8615459 SRR8615580 SRR8615664 SRR8615282 SRR8616135 SRR8616032 SRR8615407 SRR8615532 SRR8615803 SRR8616044 SRR8615920 SRR8616053 SRR8615640 SRR8615634 SRR8615238 SRR8615889 SRR8615914 SRR8615430 SRR8616015 SRR8615778 SRR8615805 SRR8616045 SRR8615661 SRR8616131 SRR8616004 SRR8615898 SRR8616054 SRR8615271 SRR8616033 SRR8615341 SRR8616041 SRR8616058 SRR8616197 SRR8615451 SRR8615528 SRR8615756 SRR8615504 SRR8615843 SRR8615758 SRR8616116 SRR8618316 SRR8615494 SRR8615990 SRR8616174 SRR8615895 SRR8616002 SRR8615464 SRR8615825 SRR8615402 SRR8615232 SRR8615812 SRR8615393 SRR8615641 SRR8615511 SRR8615432 SRR8615392 SRR8615284 SRR8615606 SRR8615519 SRR8615799 SRR8616059 SRR8615228 SRR8618301 SRR8615639 SRR8615642 SRR8615458 SRR8615971 SRR8615283 SRR8615856 SRR8615840 SRR8615870 SRR8615788 SRR8615987 SRR8615340 SRR8616189 SRR8616017 SRR8615375 SRR8615834 SRR8616108 SRR8616180 SRR8615584 SRR8615636 SRR8616024 SRR8616102 SRR8616046 SRR8615250 SRR8615806 SRR8615260 SRR8615792 SRR8615939 SRR8615963 SRR8615452 SRR8616177 SRR8616146 SRR8616192 SRR8615980 SRR8615253 SRR8616184 SRR8616019 SRR8615518 SRR8616123 SRR8615822 SRR8615563 SRR8615589 SRR8616190 SRR8615603 SRR8615804 SRR8616081 SRR8615362 SRR8615336 SRR8615404 SRR8615819 SRR8615989 SRR8616106 SRR8615462 SRR8615976 SRR8615736 SRR8615797 SRR8616007 SRR8615480 SRR8615515 SRR8616068 SRR8615611 SRR8615814 SRR8616111 SRR8616158 SRR8615476 SRR8615467 SRR8616179 SRR8615281 SRR8615482 SRR8615905 SRR8615944 SRR8615813 SRR8615894 SRR8618319 SRR8615255 SRR8615844 SRR8616115 SRR8616186 SRR8615577 SRR8615984 SRR8615261 SRR8615258 SRR8616113 SRR8615475 SRR8615429 SRR8615817 SRR8615662 SRR8616188 SRR8615933 SRR8615538 SRR8615878 SRR8615856 SRR8615426 SRR8616001 SRR8615456 SRR8616042 SRR8616110 SRR8615826 SRR8615828 SRR8615339 SRR8618305 SRR8615289 SRR8615665 SRR8615500 SRR8615395 SRR8615810 SRR8615798 SRR8615299 SRR8616055 SRR8615568 SRR8615638 SRR8616199 SRR8615533 SRR8615622 SRR8615630 SRR8616037 SRR8615695 SRR8615745 SRR8616056 SRR8616169 SRR8615648 SRR8616124 SRR8615531 SRR8616215 SRR8616075 SRR8615696 SRR8615871 SRR8616091 SRR8615643 SRR8616201 SRR8615625 SRR8615717 SRR8615450 SRR8615595 SRR8615750 SRR8615674 SRR8615832 SRR8615740 SRR8615688 SRR8616034 SRR8615600 SRR8615495 SRR8615877 SRR8615680 SRR8615378 SRR8616012 SRR8615767 SRR8616132 SRR8615800 SRR8615351 SRR8615836 SRR8615398 SRR8616069 SRR8615715 SRR8615270 SRR8615300 SRR8615739 SRR8615841 SRR8615553 SRR8615704 SRR8615361 SRR8616176 SRR8615723 SRR8615412 SRR8616031 SRR8615425 SRR8616216 SRR8615588 SRR8615708 SRR8615510 SRR8615827 SRR8618318 SRR8616164 SRR8615645 SRR8615281 SRR8615845 SRR8616011 SRR8615838 SRR8615547 SRR8616156 SRR8615597 SRR8615275 SRR8615669 SRR8615348 SRR8615401 SRR8615267 SRR8615438 SRR8615566 SRR8615265 SRR8615754 SRR8616013 SRR8618303 SRR8615354 SRR8615346 SRR8615254 SRR8615874 SRR8615679 SRR8615237 SRR8616218 SRR8615561 SRR8615628 SRR8615716 SRR8615922 SRR8615992 SRR8615445 SRR8616101 SRR8615470 SRR8615483 SRR8615245 SRR8616114 SRR8615594 SRR8616043 SRR8616039 SRR8615692 SRR8615732 SRR8615848 SRR8616143 SRR8615617 SRR8616163 SRR8615876 SRR8615722 SRR8616028 SRR8615678 SRR8615626 SRR8615854 SRR8615272 SRR8615951 SRR8616074 SRR8615236 SRR8615529 SRR8616153 SRR8615440 SRR8615263 SRR8615497 SRR8616079 SRR8615292 SRR8615623 SRR8616064 SRR8615685 SRR8616097 SRR8615431 SRR8615858 SRR8616093 SRR8615646 SRR8615962 SRR8615560 SRR8615360 SRR8615391 SRR8616038 SRR8616141 SRR8615811 SRR8615502 SRR8615959 SRR8615703 SRR8615277 SRR8615489 SRR8615712 SRR8615700 SRR8615225 SRR8615251 SRR8615655 SRR8615644 SRR8616178 SRR8615765 SRR8615490 SRR8615919 SRR8616155 SRR8615747 SRR8615466 SRR8615369 SRR8615235 SRR8615279 SRR8615278 SRR8615921 SRR8615486 SRR8615839 SRR8615668 SRR8615567
do
  printf "Downloading files associated with ${i}."
  eval ${sraDownloader} --outdir ${workingDirectory}/raw ${i}
  ##the file gets renamed upon download, but I just want it to have the SRR id and I can annotate it later
  eval mv ${workingDirectory}/raw/${i}*_1.fastq.gz ${workingDirectory}/raw/${i}_1.fastq.gz
  eval mv ${workingDirectory}/raw/${i}*_2.fastq.gz ${workingDirectory}/raw/${i}_2.fastq.gz
  #eval conda activate snakemake
  eval snakemake --cores 16 --latency-wait 300
  #eval conda deactivate
  eval rm ${workingDirectory}/raw/${i}*.fastq.gz
  eval rm ${sraCacheLocation}/sra/${i}*
  eval rm ${workingDirectory}/results/${i}*.clean.fastq.gz
  eval rm ${workingDirectory}/*.out
  eval rm ${workingDirectory}/*.tab
  eval rm -r ${workingDirectory}/*STAR*
  eval rm ${workingDirectory}/results/${i}*.bam
  eval rm ${workingDirectory}/results/${i}*.bai
done
```
