## Reprocessing some RNAseq data

This document describes the reprocessing of some RNAseq data from a previous publication. Specifically:

"Cancer-Specific Retargeting of BAF Complexes by a Prion-like Domain"
Cell, 2017, Pubmed ID: 28844694, GEO: GSE94278

### Description

I am interested in the RNAseq data where they have knocked down EWS-FLI1 expression. I want to reprocess these to get expression data, .

I am going to use [sraDownloader](https://github.com/s-andrews/sradownloader) to get at these data. I am going to save them in the directory `/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20201204_ewsFli1PrionBoulayPmid28844694/`.

I will parse these downloaded raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/).

I will use [STAR](https://github.com/alexdobin/STAR) for the alignment. I am following the instructions in the STAR [manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf).

I will use [Salmon](https://github.com/COMBINE-lab/salmon) to do read quantification. There is a great info page for this tool [here](https://salmon.readthedocs.io/en/latest/).

For library construction, they don't mention anything about whether the libraries are stranded or not, so we will assume they aren't.

I will then use [samtools](http://www.htslib.org/) to prepare the files before using [FeatureCounts](http://subread.sourceforge.net/featureCounts.html) (part of [Subread](http://subread.sourceforge.net/)) for quantification.

This is kind of a useful website for a general pipeline, [here](https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html), also [here](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html).

### Data pipeline

First we will move into our working directory and create our shell and snakemake scripts.

```shell
cd /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20201204_ewsFli1PrionBoulayPmid28844694
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
BASE_DIR = "/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20201204_ewsFli1PrionBoulayPmid28844694"


###############################
#locations of tools we will use
BBDUK = "/home/chughes/softwareTools/bbmap-38.90/bbduk.sh"
SALMON = "/home/chughes/softwareTools/salmon-1.5.2/bin/salmon"
#HISAT2 = "/home/chughes/softwareTools/hisat2-2.2.1/hisat2"
STAR = "/home/chughes/softwareTools/STAR-2.7.9a/bin/Linux_x86_64/STAR"
SAMTOOLS="/home/chughes/softwareTools/samtools-1.12/samtools"
SAMBAMBA="/home/chughes/softwareTools/sambamba-0.8.1/sambamba"
FEATURECOUNTS="/home/chughes/softwareTools/subread-2.0.3/bin/featureCounts"


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

rule featurecounts:
  input:
      "results/{smp}.sorted.bam"
  output:
      "results/{smp}.counts.txt"
  message:
      "Counting reads with featureCounts."
  shell:
      "{FEATURECOUNTS} -p --countReadPairs -t exon -g gene_id -a {GTF} -o {output} {input}"

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
workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20201204_ewsFli1PrionBoulayPmid28844694"
eval cd ${workingDirectory}
eval mkdir raw
eval mkdir results
eval mkdir quants

##loop over the accessions
for i in SRR52176{67..70}
do
  printf "Downloading files associated with ${i}."
  eval ${sraDownloader} --noena --outdir ${workingDirectory}/raw ${i}
  ##the file gets renamed upon download, but I just want it to have the SRR id and I can annotate it later
  eval mv ${workingDirectory}/raw/${i}*_1.fastq.gz ${workingDirectory}/raw/${i}_1.fastq.gz
  eval mv ${workingDirectory}/raw/${i}*_2.fastq.gz ${workingDirectory}/raw/${i}_2.fastq.gz
  #eval conda activate snakemake
  eval snakemake --cores 8
  #eval conda deactivate
  eval rm ${workingDirectory}/raw/${i}*.fastq.gz
  eval rm ${sraCacheLocation}/sra/${i}*.sra.cache
  eval rm ${workingDirectory}/results/${i}*.clean.fastq.gz
  eval rm ${workingDirectory}/results/${i}_Aligned.out.bam
done
```

This is an alternative script that uses ENA.

```shell
#!/bin/bash

##set the location of software tools and the working directory where files will be stored
ftpLocation="ftp://ftp.sra.ebi.ac.uk/vol1/fastq"
baseAccession="SRR52176"
workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20201204_ewsFli1PrionBoulayPmid28844694"
eval cd ${workingDirectory}
eval mkdir raw
eval mkdir results
eval mkdir quants

##loop over the accessions
for i in {68..70}
do
  printf "Downloading files associated with ${baseAccession}${i}."

  eval wget ${ftpLocation}/${baseAccession:0:6}/00${i: -1}/${baseAccession}${i}/${baseAccession}${i}_1.fastq.gz
  eval wget ${ftpLocation}/${baseAccession:0:6}/00${i: -1}/${baseAccession}${i}/${baseAccession}${i}_2.fastq.gz
  eval mv ${workingDirectory}/${baseAccession}${i}_1.fastq.gz ${workingDirectory}/raw
  eval mv ${workingDirectory}/${baseAccession}${i}_2.fastq.gz ${workingDirectory}/raw
  #eval conda activate snakemake
  eval snakemake --cores 8
  #eval conda deactivate
  eval rm ${workingDirectory}/raw/${baseAccession}${i}*.fastq.gz.fastq.gz
  eval rm ${workingDirectory}/results/${baseAccession}${i}*.clean.fastq.gz
  eval rm ${workingDirectory}/results/${baseAccession}${i}_Aligned.out.bam
done
```