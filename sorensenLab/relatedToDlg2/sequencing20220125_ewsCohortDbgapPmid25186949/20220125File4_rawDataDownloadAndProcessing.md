## Reprocessing some RNAseq data

This document describes the reprocessing of some RNAseq data from a previous publication. Specifically:

"The genomic landscape of pediatric Ewing sarcoma"
Cancer Discovery, 2014, Pubmed ID: 25186949, dbGaP: phs000804.v1.p1

### Description

I am interested in the RNAseq data across the entire cohort of EwS patients. These data are controlled access, so I had to request to be allowed to work with the data. These are the information provided for the RNAseq provided in the manuscript:

"RNA for each sample was converted into a library of template molecules for sequencing according to the protocol for the Illumina TruSeq unstranded RNA Sample Protocol (61). All libraries were sequenced on the Illumina HiSeq 2000 instrument set to generate 76-bp paired-end reads for WES and RNASeq and 101-bp paired-end reads for WGS."

I am going to use the instructions details [here](https://www.ncbi.nlm.nih.gov/sra/docs/sra-dbgap-download/) to get at these data.

I will parse these downloaded raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/).

I will use [STAR](https://github.com/alexdobin/STAR) for the alignment. I am following the instructions in the STAR [manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf).

I will use [Salmon](https://github.com/COMBINE-lab/salmon) to do read quantification. There is a great info page for this tool [here](https://salmon.readthedocs.io/en/latest/).

I will then use [samtools](http://www.htslib.org/) to prepare the files before using [FeatureCounts](http://subread.sourceforge.net/featureCounts.html) (part of [Subread](http://subread.sourceforge.net/)) for quantification.

I will then use the bamCoverage function from [deeptools](https://deeptools.readthedocs.io/en/develop/) to get coverage estimates across the transcriptome.

### Data pipeline

First we will move into our working directory and create our shell and snakemake scripts.

```shell
cd /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20220125_ewsCohortDbgapPmid25186949
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
BASE_DIR = "/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20220125_ewsCohortDbgapPmid25186949"


###############################
#locations of tools we will use
BBDUK = "/home/chughes/softwareTools/bbmap-38.90/bbduk.sh"
SALMON = "/home/chughes/softwareTools/salmon-1.5.2/bin/salmon"
#HISAT2 = "/home/chughes/softwareTools/hisat2-2.2.1/hisat2"
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
      "{SALMON} quant -i {SALMONINDEX} -l A -p 8 --gcBias --validateMappings -o {params.dir} -1 {input.r1} -2 {input.r2}"
```

Below is the shell script I will use to process these data with snakemake. 

```shell
#!/bin/bash

##set the location of software tools and the working directory where files will be stored
sraDownloader="/home/chughes/softwareTools/sradownloader-3.8/sradownloader"
sraCacheLocation="/mnt/Data/chughes/sratoolsRepository"
workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20220125_ewsCohortDbgapPmid25186949"
eval cd ${workingDirectory}
eval mkdir raw
eval mkdir results
eval mkdir quants

##loop over the accessions
for i in SRR3591303 SRR3772620 SRR3625954 SRR3626996 SRR3627394 SRR3586330 SRR3775145 SRR3537890 SRR3476973 SRR3479294 SRR3691911 SRR3668689 SRR3626203 SRR3631976 SRR3584860 SRR3609632 SRR3741339 SRR3767711 SRR3769004 SRR3770748 SRR3774859 SRR3608717 SRR3679321 SRR3644440 SRR3641619 SRR3642521 SRR3538014 SRR3480650 SRR3485513 SRR3478236 SRR3478277 SRR3489518 SRR3492874
do
  printf "Downloading files associated with ${i}."

  ##prefetch the files from SRA
  eval prefetch --ngc prj_29838.ngc -p ${i}
  eval mv ${sraCacheLocation}/sra/${i}*.sra ${workingDirectory}/raw/${i}.sra
  eval fasterq-dump --ngc prj_29838.ngc -p ${workingDirectory}/raw/${i}.sra
  eval mv ${workingDirectory}/${i}*.fastq ${workingDirectory}/raw
  eval gzip -v ${workingDirectory}/raw/${i}_1.fastq
  eval gzip -v ${workingDirectory}/raw/${i}_2.fastq

  #process data files
  eval snakemake --cores 8 --latency-wait 300

  #clean-up
  eval rm ${workingDirectory}/raw/${i}*.sra
  eval rm ${workingDirectory}/raw/${i}*.fastq.gz
  eval rm ${workingDirectory}/results/${i}*.clean.fastq.gz
  eval rm ${workingDirectory}/*.out
  eval rm ${workingDirectory}/*.tab
  eval rm -r ${workingDirectory}/*STAR*
  eval rm ${workingDirectory}/results/${i}*.bam
  eval rm ${workingDirectory}/results/${i}*.bai
done
```

This is a script I used to get the sra files alongside the processing, just to save time.

```shell
#!/bin/bash

##set the location of software tools and the working directory where files will be stored
#sraDownloader="/home/chughes/softwareTools/sradownloader-3.8/sradownloader"
sraCacheLocation="/mnt/Data/chughes/sratoolsRepository"
workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20220125_ewsCohortDbgapPmid25186949"
eval cd ${workingDirectory}
#eval mkdir raw
#eval mkdir results
#eval mkdir quants

##loop over the accessions
for i in SRR3591303 SRR3772620 SRR3625954 SRR3626996 SRR3627394 SRR3586330 SRR3775145 SRR3537890 SRR3476973 SRR3479294 SRR3691911 SRR3668689 SRR3626203 SRR3631976 SRR3584860 SRR3609632 SRR3741339 SRR3767711 SRR3769004 SRR3770748 SRR3774859 SRR3608717 SRR3679321 SRR3644440 SRR3641619 SRR3642521 SRR3538014 SRR3480650 SRR3485513 SRR3478236 SRR3478277 SRR3489518 SRR3492874
do
  printf "Downloading files associated with ${i}."

  ##prefetch the files from SRA
  eval prefetch --ngc prj_29838.ngc -p ${i}
  #eval mv ${sraCacheLocation}/sra/${i}*.sra ${workingDirectory}/raw/${i}.sra
  #eval fasterq-dump --ngc prj_29838.ngc -p ${workingDirectory}/raw/${i}.sra
  #eval mv ${workingDirectory}/${i}*.fastq ${workingDirectory}/raw
  #eval gzip -v ${workingDirectory}/raw/${i}_1.fastq
  #eval gzip -v ${workingDirectory}/raw/${i}_2.fastq

  #process data files
  #eval snakemake --cores 8 --latency-wait 300

  #clean-up
  #eval rm ${workingDirectory}/raw/${i}*.sra
  #eval rm ${workingDirectory}/raw/${i}*.fastq.gz
  #eval rm ${workingDirectory}/results/${i}*.clean.fastq.gz
  #eval rm ${workingDirectory}/*.out
  #eval rm ${workingDirectory}/*.tab
  #eval rm -r ${workingDirectory}/*STAR*
  #eval rm ${workingDirectory}/results/${i}*.bam
  #eval rm ${workingDirectory}/results/${i}*.bai
done
```
