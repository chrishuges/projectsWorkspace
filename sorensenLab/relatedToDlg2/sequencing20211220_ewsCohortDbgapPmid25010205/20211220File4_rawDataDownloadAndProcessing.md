## Reprocessing some RNAseq data

This document describes the reprocessing of some RNAseq data from a previous publication. Specifically:

"The genomic landscape of the Ewing Sarcoma family of tumors reveals recurrent STAG2 mutation"
PLoS Genetics, 2021, Pubmed ID: 25010205, dbGaP: phs000768.v2.p1

### Description

I am interested in the RNAseq data across the entire cohort of EwS patients. These data are controlled access, so I had to request to be allowed to work with the data. These are the information provided for the RNAseq provided in the manuscript:

PolyA selected RNA libraries were prepared for sequencing on the Illumina HiSeq2000 using TruSeq v3 chemistry according to the manufacturer's protocol (Illumina). RNA sequencing was performed with an average yield of 18.6 Gb per sample. Raw reads were mapped using to ENSEMBL reference (hg19) using TopHat2.0 [60]. Fusion analysis was done using TopHat 2.0 and DeFuse 0.6 [61]. The 3 alternated fusions described were confirmed using RT-PCR using flanking primers and Sanger sequencing of the resultant product.

Expression FPKM results were obtained at both gene and transcript level using CuffLinks 2.1 [62]. The log2 FPKM expression results from TopHat mapping were median-normalized using in-house data from 63 normal tissue samples. Exon level expression was calculated using the formula RPKM = (r * 109)/(f * R), with r being the number of reads mapped to an exon, f being the exon length, and R being the total read count of the sample. Hierarchical clustering was performed on normalized log2 FPKM expression values at the gene level using Euclidean distance and Ward agglomeration method.

For variant detection, samtools (http://samtools.sourceforge.net/) is used to count the number of reads uniquely mapped to a position found as variant in DNA sequencing of the same sample or a position of interest based on a mutation being present in the TCGA (http://cancergenome.nih.gov/) or compared to the reference genome hg19 in genes of interest. If there are reads supporting a variant base then the total reads supporting it are counted and variant allele frequency is calculated.

I am going to use [sraDownloader](https://github.com/s-andrews/sradownloader) to get at these data.

I will parse these downloaded raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/).

I will use [STAR](https://github.com/alexdobin/STAR) for the alignment. I am following the instructions in the STAR [manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf).

I will use [Salmon](https://github.com/COMBINE-lab/salmon) to do read quantification. There is a great info page for this tool [here](https://salmon.readthedocs.io/en/latest/).

I will then use [samtools](http://www.htslib.org/) to prepare the files before using [FeatureCounts](http://subread.sourceforge.net/featureCounts.html) (part of [Subread](http://subread.sourceforge.net/)) for quantification.

I will then use the bamCoverage function from [deeptools](https://deeptools.readthedocs.io/en/develop/) to get coverage estimates across the transcriptome.

### Data pipeline

First we will move into our working directory and create our shell and snakemake scripts.

```shell
cd /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211220_ewsCohortDbgapPmid25010205
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
#working directory
BASE_DIR = "/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211220_ewsCohortDbgapPmid25010205"


###############################
#locations of tools we will use
BBDUK = "/home/chughes/softwareTools/bbmap-38.90/bbduk.sh"
SALMON = "/home/chughes/softwareTools/salmon-1.5.2/bin/salmon"
STAR = "/home/chughes/softwareTools/STAR-2.7.9a/bin/Linux_x86_64/STAR"
SAMTOOLS="/home/chughes/softwareTools/samtools-1.12/samtools"
#SAMBAMBA="/home/chughes/softwareTools/sambamba-0.8.1/sambamba"
FEATURECOUNTS="/home/chughes/softwareTools/subread-2.0.3/bin/featureCounts"
BAMCOVERAGE="/home/chughes/virtualPython368/bin/bamCoverage"


###############################
#locations of our index files
DATABASE_DIR = "/home/chughes/databases/projectEwsDlg2"
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
      expand("results/{smp}.sorted.bam.bai", smp = SAMPLES),
      expand("results/{smp}.sorted.bw", smp = SAMPLES),
      expand("results/{smp}.counts.txt", smp = SAMPLES),      
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
      "{STAR} --runThreadN 8 --genomeDir {STARINDEX} --readFilesIn {input.r1} {input.r2} --readFilesCommand zcat --sjdbGTFfile {GTF} --outStd SAM | {SAMTOOLS} sort -o {output}"

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
      "{BAMCOVERAGE} -b {input.r1} -o {output} -p 8"

rule featurecounts:
  input:
      r1 = "results/{smp}.sorted.bam",
      w2 = "results/{smp}.sorted.bam.bai",
      w3 = "results/{smp}.sorted.bw
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
      "{SALMON} quant -i {SALMONINDEX} -l A -p 8 --gcBias --validateMappings -o {params.dir} -1 {input.r1} -2 {input.r2}"
```

Below is the shell script I will use to process these data with snakemake. 

```shell
#!/bin/bash

##set the location of software tools and the working directory where files will be stored
sraDownloader="/home/chughes/softwareTools/sradownloader-3.8/sradownloader"
sraCacheLocation="/mnt/Data/chughes/sratoolsRepository"
workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211220_ewsCohortDbgapPmid25010205"
eval cd ${workingDirectory}
eval mkdir raw
eval mkdir results
eval mkdir quants

##loop over the accessions
for i in SRR5163{665..757} #was {671..757}
do
  printf "Downloading files associated with ${i}."

  ##prefetch the files from SRA
  eval prefetch --ngc prj_29838.ngc -p ${i}
  eval mv ${i}/${i}*.sra ./${i}.sra
  eval rm -r ${i}
  eval fasterq-dump --ngc prj_29838.ngc -p ${i}.sra
  eval rm ${i}.sra
  eval mv ${i}*.fastq ${workingDirectory}/raw
  eval gzip -v ${workingDirectory}/raw/${i}_1.fastq
  eval gzip -v ${workingDirectory}/raw/${i}_2.fastq
  
  #process data files
  eval snakemake --cores 8 --latency-wait 300

  #clean-up
  eval rm ${workingDirectory}/raw/${i}*
  eval rm ${workingDirectory}/results/${i}*.clean.fastq.gz
  eval rm ${workingDirectory}/*.out
  eval rm ${workingDirectory}/*.tab
  eval rm -r ${workingDirectory}/*STAR*
  eval rm ${workingDirectory}/results/${i}*.bam
  eval rm ${workingDirectory}/results/${i}*.bai
done
```


This was the script I used previously:

```shell
#! /bin/bash
salmonLocation="/home/chughes/softwareTools/salmon-1.5.1/bin/salmon"
indexLocation="/home/chughes/databases/refgenieManualGenomes/hg38/salmonPartialSaIndex_072021/default/"
starLocation="/home/chughes/softwareTools/STAR-2.7.1a/bin/Linux_x86_64/STAR"
referenceLocation="/home/chughes/databases/refgenieManualGenomes/hg38/starIndex_072021/default/"
annotationLocation="/home/chughes/databases/refgenieManualGenomes/hg38/gencodeGtf_072021/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.gtf"
samtoolsLocation="/home/chughes/softwareTools/samtools-1.12/samtools"
rawDataOutputDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencingEwsCohortDbgapPmid25010205/"


#create the directory to store the files
if [ ! -d $rawDataOutputDirectory ]; then
  printf "Data directory $rawDataOutputDirectory does not exist, creating it.\n"
  eval "mkdir $rawDataOutputDirectory"
else
  printf "Data directory $rawDataOutputDirectory exists, moving on.\n"
fi

eval "cd ${rawDataOutputDirectory}"
eval "mkdir ${rawDataOutputDirectory}starResults"

#download the process the files
for i in SRR5163{738..748}
do
  echo $i

  ##prefetch the files from SRA
  eval prefetch --ngc prj_29838.ngc -p ${i}
  eval mv ${i}/${i}*.sra ./${i}.sra
  eval rm -r ${i}
  eval fasterq-dump --ngc prj_29838.ngc -p ${i}.sra
  eval rm ${i}.sra
#  eval rm -r fasterq.tmp*
  eval gzip -v ${i}_1.fastq
  eval gzip -v ${i}_2.fastq

  ##salmon analylsis
  salmonCall="$salmonLocation quant -i $indexLocation -l A -1 ${rawDataOutputDirectory}${i}_1.fastq.gz -2 ${rawDataOutputDirectory}${i}_2.fastq.gz --gcBias --validateMappings -o ${rawDataOutputDirectory}${i}_quant"
  eval $salmonCall

  ##star alignment
  starCall="$starLocation --runThreadN 12 --genomeDir ${referenceLocation} --readFilesIn ${rawDataOutputDirectory}${i}_1.fastq.gz ${rawDataOutputDirectory}${i}_2.fastq.gz --readFilesCommand zcat --sjdbGTFfile ${annotationLocation} --outFileNamePrefix ${rawDataOutputDirectory}starResults/${i}_ --outSAMtype BAM Unsorted"
  eval $starCall

  ##bamsort
  bamSortCall="$samtoolsLocation sort ${rawDataOutputDirectory}/starResults/${i}_Aligned.out.bam -o ${rawDataOutputDirectory}/starResults/${i}.sorted.bam"
  eval $bamSortCall

  ##index bam files
  bamIndexCall="$samtoolsLocation index ${rawDataOutputDirectory}/starResults/${i}.sorted.bam"
  eval $bamIndexCall

  ##calculate coverage on chromosome 11
  coverageCall="bamCoverage -b ${rawDataOutputDirectory}/starResults/${i}.sorted.bam -o ${rawDataOutputDirectory}/starResults/${i}.chr11.bw --binSize 10 --region chr11 --normalizeUsing BPM --smoothLength 30 --centerReads -p 6"
  eval $coverageCall

  ##remove original files
  eval rm ${i}*.fastq.gz
  eval rm ${rawDataOutputDirectory}/starResults/*Aligned*.bam

done
```

This was a script I used previously to get coverage maps.

```shell
#!/bin/bash

##set the location of software tools and the working directory where files will be stored
bamCoverage="/home/chughes/virtualPython368/bin/bamCoverage"
workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211220_ewsCohortDbgapPmid25010205/archive/starResults"
eval cd ${workingDirectory}


##loop over the accessions
for i in *.bam
do
  printf "Processing files associated with ${i}.\n\n"

  ##get coverage
  eval $bamCoverage -b ${i} -o ${i::-3}.bw -p 6
  eval rm ${i}
done
```