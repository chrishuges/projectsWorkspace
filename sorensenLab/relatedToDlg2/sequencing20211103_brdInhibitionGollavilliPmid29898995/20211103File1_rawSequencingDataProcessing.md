## Reprocessing some RNAseq data

This document describes the reprocessing of some RNAseq data from a previous publication. Specifically:

"BET bromodomain dependency in EWS/ETS driven Ewing Sarcoma"
GEO: 113604, PMID: 29898995

### Description

I am interested in the CHLA10 RNAseq data where they have knocked down EWS-FLI1 expression. I want to reprocess these to get expression data, .

I am going to use [sraDownloader](https://github.com/s-andrews/sradownloader) to get at these data. I am going to save them in the directory `/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995/`. 

I will parse these downloaded raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/). 

I am going to use [HiSAT2](http://daehwankimlab.github.io/hisat2/) to do the alignment. For the index, I created one as described in `..\relatedToDlg2\sequencing20211029_databaseSetupAndSnakemake\20211029File1_downloadingDatabasesAndIndexes.md`. I am going to use the transcriptome aware index so I can capture splicing information. There is some good info on the alignment details [here](http://daehwankimlab.github.io/hisat2/manual/). 

For library construction, it appears the authors used the [TruSeq RNA Library Prep Kit v2](https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/truseq-rna-v2.html). This appears to be non-stranded, so we will leave HiSAT2 with default settings.


I will then use [samtools](http://www.htslib.org/) to prepare the files before using [FeatureCounts](http://subread.sourceforge.net/featureCounts.html) (part of [Subread](http://subread.sourceforge.net/)) for quantification.

This is kind of a useful website for a general pipeline, [here](https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html).

### Data pipeline

First we will move into our working directory and create our shell and snakemake scripts.

```shell
cd /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995
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
Date: 20211103
"""

###############################
#working directory
BASE_DIR = "/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995"


###############################
#locations of tools we will use
BBDUK = "/home/chughes/softwareTools/bbmap-38.90/bbduk.sh"
SALMON = "/home/chughes/softwareTools/salmon-1.5.2/bin/salmon"
HISAT2 = "/home/chughes/softwareTools/hisat2-2.2.1/hisat2"
SAMTOOLS="/home/chughes/softwareTools/samtools-1.12/samtools"
FEATURECOUNTS="/home/chughes/softwareTools/subread-2.0.3/bin/featureCounts"

/home/chughes/databases/projectEwsDlg2/hisat2Index
###############################
#locations of our index files
DATABASE_DIR = "/home/chughes/databases/projectEwsDlg2"
INDEX = DATABASE_DIR + "/hisat2Index/genome_tran"
SALMONINDEX = DATABASE_DIR + "/salmonIndex/salmon_index"
GTF   = DATABASE_DIR + "/baseGenomeFiles/genome.gtf"
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
    input: expand("results/{smp}.counts.txt", smp = SAMPLES)

rule bbduk:
  input:
      r1 = "raw/{smp}_1.fastq.gz"
      r2 = "raw/{smp}_2.fastq.gz"
  output:
      ro1 = "results/{smp}_1.clean.fastq.gz"
      ro2 = "results/{smp}_2.clean.fastq.gz"
  message:
      "Processing with BBDuk."
  shell:
      "{BBDUK} in1={input.r1} in2={input.r2} ref=adapters out1={output.ro1} out2={output.ro2} ktrim=r k=23 mink=11 hdist=1"

rule salmon:
  input:
      r1 = "results/{smp}_1.clean.fastq.gz"
      r2 = "results/{smp}_2.clean.fastq.gz"
  output:
      "quants/{smp}/quant.sf"
  params:
      dir = "quants/{smp}"
  message:
      "Quantifying with salmon."
  shell:
      "{SALMON} quant -i {SALMONINDEX} -l A -p 8 --gcBias --validateMappings -o {params.dir} -1 {input.r1} -2 {input.r2}"

rule hisat2:
  input:
      r1 = "results/{smp}_1.clean.fastq.gz"
      r2 = "results/{smp}_2.clean.fastq.gz"
  output:
      "results/{smp}.unsorted.sam"
  message:
      "Aligning with hisat2."
  shell:
      "{HISAT2} -p 8 -x {INDEX} -1 {input.r1} -2 {input.r2} -S {output}"

rule bam_conversion:
  input:
      "results/{smp}.unsorted.sam"
  output:
      "results/{smp}.unsorted.bam"
  message:
      "SAM to BAM conversion with samtools."
  shell:
      "{SAMTOOLS} view -h -b -o {output} {input}"

rule bam_sorting:
  input:
      "results/{smp}.unsorted.bam"
  output:
      "results/{smp}.sorted.bam"
  message:
      "BAM sorting with sambamba."
  shell:
      "{SAMBAMBA} sort -t 6 -o {output} {input}"

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
```

Below is the shell script I will use to process these data with snakemake.

```shell
#!/bin/bash

##set the location of software tools and the working directory where files will be stored
sraDownloader="/home/chughes/softwareTools/sradownloader-3.8/sradownloader"
sraCacheLocation="/mnt/Data/chughes/sratoolsRepository"
workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211103_brdInhibitionGollavilliPmid29898995"
eval cd ${workingDirectory}
eval mkdir raw
eval mkdir results

##loop over the accessions
for i in SRR7059{715..722}
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
  eval rm ${workingDirectory}/results/${i}.unsorted.sam
  eval rm ${workingDirectory}/results/${i}.unsorted.bam
  eval rm ${workingDirectory}/results/${i}.sorted.bam
  eval rm ${workingDirectory}/results/${i}.sorted.bam.bai
done
```




