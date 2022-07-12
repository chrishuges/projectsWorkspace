## Processing some RNAseq data

This document details analysis of the A673 timecourse RNAseq data.

The first thing we will do is parse out the UMI's from the individual reads using the [umi-tools](https://github.com/weng-lab/umitools) package.

I will parse these raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/).

I will use [STAR](https://github.com/alexdobin/STAR) for the alignment. I am following the instructions in the STAR [manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf).

I will use [Salmon](https://github.com/COMBINE-lab/salmon) to do read quantification. There is a great info page for this tool [here](https://salmon.readthedocs.io/en/latest/).

I will then use [samtools](http://www.htslib.org/) to prepare the files before using [FeatureCounts](http://subread.sourceforge.net/featureCounts.html) (part of [Subread](http://subread.sourceforge.net/)) for quantification.

This is kind of a useful website for a general pipeline, [here](https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html), also [here](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html).

The raw data are stored here: `/projects/analysis/analysis34/`. I transferred the data to our own server from that location.

### Data pipeline

First we will move into our working directory and create our shell and snakemake scripts.

```shell
cd /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput
touch rnaDataProcessingScript.sh
chmod +x rnaDataProcessingScript.sh
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
BBDUK = "/home/chughes/softwareTools/bbmap-38.90/bbduk.sh"
#BBDUK = "/projects/ptx_analysis/chughes/softwareTools/bbmap-38.90/bbduk.sh"
SALMON = "/home/chughes/softwareTools/salmon-1.5.2/bin/salmon"
#SALMON = "/projects/ptx_analysis/chughes/softwareTools/salmon-1.5.2/bin/salmon"
#HISAT2 = "/home/chughes/softwareTools/hisat2-2.2.1/hisat2"
STAR = "/home/chughes/softwareTools/STAR-2.7.9a/bin/Linux_x86_64/STAR"
#STAR = "/projects/ptx_analysis/chughes/softwareTools/STAR-2.7.9a/bin/Linux_x86_64/STAR"
SAMTOOLS="/home/chughes/softwareTools/samtools-1.12/samtools"
#SAMTOOLS="/gsc/software/linux-x86_64-centos7/samtools-1.14/bin/samtools"
#SAMBAMBA="/home/chughes/softwareTools/sambamba-0.8.1/sambamba"
FEATURECOUNTS="/home/chughes/softwareTools/subread-2.0.3/bin/featureCounts"
#FEATURECOUNTS="/projects/ptx_analysis/chughes/softwareTools/subread-2.0.3/bin/featureCounts"
BAMCOVERAGE="/home/chughes/virtualPython368/bin/bamCoverage"
#BAMCOVERAGE="/home/chughes/Virtual_Python383/bin/bamCoverage"
UMITOOLS="/home/chughes/virtualPython368/bin/umitools"


###############################
#locations of our index files
DATABASE_DIR = "/home/chughes/databases/projectEwsDlg2"
#DATABASE_DIR = "/projects/ptx_analysis/chughes/databases/projectEwsDlg2"
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

rule umi:
  input:
      r1 = "raw/{smp}_1.fastq.gz",
      r2 = "raw/{smp}_2.fastq.gz"
  output:
      ro1 = "results/{smp}_1.umi.fastq.gz",
      ro2 = "results/{smp}_2.umi.fastq.gz"
  message:
      "Processing with UMI tools."
  shell:
      "{UMITOOLS} reformat_fastq -l {input.r1} -r {input.r2} -L {output.ro1} -R {output.ro2}"

rule bbduk:
  input:
      r1 = "results/{smp}_1.umi.fastq.gz",
      r2 = "results/{smp}_2.umi.fastq.gz"
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
#sraDownloader="/home/chughes/softwareTools/sradownloader-3.8/sradownloader"
#sraDownloader="/projects/ptx_analysis/chughes/softwareTools/sradownloader-3.8/sradownloader"
#sraCacheLocation="/mnt/Data/chughes/sratoolsRepository"
#sraCacheLocation="/projects/ptx_results/Sequencing/sraCache"
workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20210421_a673TimecourseRnaSeqOutput"
eval cd ${workingDirectory}
eval mkdir raw
eval mkdir results
eval mkdir quants

##loop over the accessions
for i in ACAGTG_setA ACTTGA_setA ATCACG_setA CAGATC_setA CGATGT_setA GCCAAT_setA TGACCA_setA TTAGGC_setA ACAGTG_setB ACTTGA_setB ATCACG_setB CAGATC_setB CGATGT_setB GCCAAT_setB TGACCA_setB TTAGGC_setB ACAGTG_setC ACTTGA_setC ATCACG_setC CAGATC_setC CGATGT_setC GCCAAT_setC TGACCA_setC TTAGGC_setC
do
  printf "Processing files associated with ${i}.\n"
  #move the files to the working directory
  eval mv ${workingDirectory}/fastq/${i}_1.fastq.gz ${workingDirectory}/raw/${i}_1.fastq.gz
  eval mv ${workingDirectory}/fastq/${i}_2.fastq.gz ${workingDirectory}/raw/${i}_2.fastq.gz
  #eval conda activate snakemake
  eval snakemake --cores 8 --latency-wait 300
  #eval conda deactivate
  #eval rm ${workingDirectory}/raw/${i}*.fastq.gz
  #eval rm ${workingDirectory}/results/${i}*.clean.fastq.gz
  #eval rm ${workingDirectory}/*.out
  #eval rm ${workingDirectory}/*.tab
  #eval rm -r ${workingDirectory}/*STAR*
  #eval rm ${workingDirectory}/results/${i}*.bam
  #eval rm ${workingDirectory}/results/${i}*.bai
  #move the files back to the fastq directory
  eval mv ${workingDirectory}/raw/${i}_1.fastq.gz ${workingDirectory}/fastq/${i}_1.fastq.gz
  eval mv ${workingDirectory}/raw/${i}_2.fastq.gz ${workingDirectory}/fastq/${i}_2.fastq.gz
done
```

For downstream analysis, file file layout is:

sampleName | cell | treatment | barcode | experiment
| --- | --- | --- | --- | --- |
day0_rep1 | a673 | none | ATCACG | PX1955
day0_rep2 | a673 | none | ATCACG | PX1956
day0_rep3 | a673 | none | ATCACG | PX1957
day7_rep1 | a673 | none | CGATGT | PX1955
day7_rep2 | a673 | none | CGATGT | PX1956
day7_rep3 | a673 | none | CGATGT | PX1957
day9_rep1 | a673 | none | TTAGGC | PX1955
day9_rep2 | a673 | none | TTAGGC | PX1956
day9_rep3 | a673 | none | TTAGGC | PX1957
day10_rep1 | a673 | none | TGACCA | PX1955
day10_rep2 | a673 | none | TGACCA | PX1956
day10_rep3 | a673 | none | TGACCA | PX1957
day11_rep1 | a673 | none | ACAGTG | PX1955
day11_rep2 | a673 | none | ACAGTG | PX1956
day11_rep3 | a673 | none | ACAGTG | PX1957
day14_rep1 | a673 | none | GCCAAT | PX1955
day14_rep2 | a673 | none | GCCAAT | PX1956
day14_rep3 | a673 | none | GCCAAT | PX1957
day17_rep1 | a673 | none | CAGATC | PX1955
day17_rep2 | a673 | none | CAGATC | PX1956
day17_rep3 | a673 | none | CAGATC | PX1957
day22_rep1 | a673 | none | ACTTGA | PX1955
day22_rep2 | a673 | none | ACTTGA | PX1956
day22_rep3 | a673 | none | ACTTGA | PX1957

Now that we have these data, process them downstream with R.
