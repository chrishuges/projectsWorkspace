## Reprocessing some RNAseq data

This document describes the re-processing of some RNAseq data from neuroblastoma cells published in the paper:

Zeid R, Lawlor MA, Poon E, Reyes JM et al. Enhancer invasion shapes MYCN-dependent transcriptional amplification in neuroblastoma. Nat Genet 2018 Apr;50(4):515-523. PMID: 29379199, GEO: GSE80153

### Description

I will parse these raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/).

I will use [Salmon](https://github.com/COMBINE-lab/salmon) to do read quantification. There is a great info page for this tool [here](https://salmon.readthedocs.io/en/latest/).

### Data pipeline

First I need to create my working directory on the sorensen lab server.

```shell
cd /projects/ptx_results/Sequencing/publishedStudies/sequencing20220711_nbCellLineMycNPmid29379199
touch seqDataProcessingScript.sh
chmod +x seqDataProcessingScript.sh
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
#STAR = "/projects/ptx_analysis/chughes/softwareTools/STAR-2.7.9a/bin/Linux_x86_64/STAR"
#SAMTOOLS="/home/chughes/softwareTools/samtools-1.12/samtools"
#SAMTOOLS="/gsc/software/linux-x86_64-centos7/samtools-1.14/bin/samtools"
#SAMBAMBA="/home/chughes/softwareTools/sambamba-0.8.1/sambamba"
#FEATURECOUNTS="/home/chughes/softwareTools/subread-2.0.3/bin/featureCounts"
#FEATURECOUNTS="/projects/ptx_analysis/chughes/softwareTools/subread-2.0.3/bin/featureCounts"
#BAMCOVERAGE="/home/chughes/virtualPython368/bin/bamCoverage"
#BAMCOVERAGE="/home/chughes/Virtual_Python383/bin/bamCoverage"


###############################
#locations of our index files
#DATABASE_DIR = "/home/chughes/databases/projectEwsDlg2"
DATABASE_DIR = "/projects/ptx_analysis/chughes/databases/projectEwsDlg2"
#STARINDEX = DATABASE_DIR + "/starIndex"
SALMONINDEX = DATABASE_DIR + "/salmonIndex/salmon_index"
GTF = DATABASE_DIR + "/baseGenomeFiles/genome.gtf"
FASTA = DATABASE_DIR + "/baseGenomeFiles/genome.fa"


###############################
#get our list of sample files
SAMPLES, = glob_wildcards("raw/{smp}.fastq.gz")
print("There are " + str(len(SAMPLES)) + " total samples to be processed.")
for smp in SAMPLES:
  print("Sample " + smp + " will be processed")


###############################
#processing workflow
rule all:
    input:
      expand("quants/{smp}/quant.sf", smp = SAMPLES)

rule bbduk:
  input:
      r1 = "raw/{smp}.fastq.gz"
  output:
      ro1 = "results/{smp}.clean.fastq.gz"
  message:
      "Processing with BBDuk."
  shell:
      "{BBDUK} in={input.r1} ref=adapters out={output.ro1} ktrim=r k=23 mink=11 hdist=1"

rule salmon:
  input:
      r1 = "results/{smp}.clean.fastq.gz"
  output:
      "quants/{smp}/quant.sf"
  params:
      dir = "quants/{smp}"
  message:
      "Quantifying with salmon."
  shell:
      "{SALMON} quant -i {SALMONINDEX} -l A -p 16 --gcBias --validateMappings -o {params.dir} -r {input.r1}"
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
#workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToHaifeng/sequencing20220628_nbCellLineRnaSeq"
workingDirectory="/projects/ptx_results/Sequencing/publishedStudies/sequencing20220711_nbCellLineMycNPmid29379199"
eval cd ${workingDirectory}
eval mkdir raw
eval mkdir results
eval mkdir quants

##loop over the accessions
for i in SRR3356{427..446} 
do
  printf "Downloading files associated with ${i}."
  eval ${sraDownloader} --outdir ${workingDirectory}/raw ${i}
  eval mv ${workingDirectory}/raw/${i}_GSM*.fastq.gz ${workingDirectory}/raw/${i}.fastq.gz
  #eval conda activate snakemake
  eval snakemake --cores 16 --latency-wait 300
  #eval conda deactivate
  eval rm ${workingDirectory}/raw/${i}*.fastq.gz
  eval rm ${sraCacheLocation}/sra/${i}*
  eval rm ${workingDirectory}/results/${i}*.clean.fastq.gz
  #eval rm ${workingDirectory}/*.out
  #eval rm ${workingDirectory}/*.tab
  #eval rm -r ${workingDirectory}/*STAR*
  #eval rm ${workingDirectory}/results/${i}*.bam
  #eval rm ${workingDirectory}/results/${i}*.bai
done
```