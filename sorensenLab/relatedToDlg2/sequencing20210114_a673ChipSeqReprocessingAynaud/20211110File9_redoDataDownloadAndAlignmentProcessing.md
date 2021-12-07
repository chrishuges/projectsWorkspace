## Reprocessing some ChIPseq data

This document describes the reprocessing of some ChIPseq data from a previous publication. Specifically:

"Transcriptional Programs Define Intratumoral Heterogeneity of Ewing Sarcoma at Single-Cell Resolution"
Cell Reports, 2020, Pubmed ID: 32049009, GEO: GSE130025

### Description

I am interested in the ChIPseq data across the cell lines they profiled. There is a metadata file in the server directory detailed below that describes these data.

I am going to use [SRAExplorer](https://sra-explorer.info/#) to get at these data. I am going to save them in the directory `/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211027_grunewaldEwsEncyclopedia/`. I will parse these downloaded raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/). Because the authors originally used Bowtie for the alignment, I will as well. For the index, I created one as described in `..\relatedToDlg2\sequencing20211029_databaseSetupAndSnakemake\20211029File1_downloadingDatabasesAndIndexes.md`. There is a pretty good walkthrough of the alignment process [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/03_align_and_filtering.html). I will then use [samtools](http://www.htslib.org/) and [sambamba](https://lomereiter.github.io/sambamba/) to sort and filter the sam and bam files that are generated, prior to peak calling. I ended up writing this into a script that uses snakemake. Proceed with the analysis as below.

### Data pipeline

First we will move into our working directory and create our shell and snakemake scripts.

```shell
cd /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211110_aynaudScRnaseqPmid32049009
touch sraDataProcessingScript.sh
chmod +x sraDataProcessingScript.sh
touch snakefile
```

Edit the contents of the snakefile to include the text below. I use vim for this.

```shell
"""
Author: Christopher Hughes
Affiliation: BCCRC
Aim: Workflow for ChIP-Seq data
Date: 20211030
"""

###############################
#working directory
BASE_DIR = "/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211110_aynaudScRnaseqPmid32049009"


###############################
#locations of tools we will use
BBDUK = "/home/chughes/softwareTools/bbmap-38.90/bbduk.sh"
BOWTIE2 = "/home/chughes/softwareTools/bowtie2-2.4.4/bowtie2"
SAMTOOLS="/home/chughes/softwareTools/samtools-1.12/samtools"
SAMBAMBA="/home/chughes/softwareTools/sambamba-0.8.1/sambamba"

###############################
#locations of our index files
DATABASE_DIR = "/home/chughes/databases/projectEwsDlg2"
INDEX = DATABASE_DIR + "/bowtie2Index/human_genome"
GTF   = DATABASE_DIR + "/baseGenomeFiles/genome.gtf"
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
    input: expand("results/{smp}.filtered.bam.bai", smp = SAMPLES)

rule bbduk:
  input:
      "raw/{smp}.fastq.gz"
  output:
      "results/{smp}.clean.fastq.gz"
  message:
      "Processing with BBDuk."
  shell:
      "{BBDUK} in={input} ref=adapters out={output} ktrim=r k=23 mink=11 hdist=1"

rule bowtie:
  input:
      "results/{smp}.clean.fastq.gz"
  output:
      "results/{smp}.unsorted.sam"
  message:
      "Aligning with bowtie2."
  shell:
      "{BOWTIE2} -p 8 -q --local -x {INDEX} -U {input} -S {output}"

rule bam_conversion:
  input:
      "results/{smp}.unsorted.sam"
  output:
      "results/{smp}.unsorted.bam"
  message:
      "SAM to BAM conversion with samtools."
  shell:
      "{SAMTOOLS} view -h -S -b -o {output} {input}"

rule bam_sorting:
  input:
      "results/{smp}.unsorted.bam"
  output:
      "results/{smp}.sorted.bam"
  message:
      "BAM sorting with sambamba."
  shell:
      "{SAMBAMBA} sort -t 6 -o {output} {input}"

rule bam_filtering:
  input:
      "results/{smp}.sorted.bam"
  output:
      "results/{smp}.filtered.bam"
  message:
      "Filtering BAM with sambamba."
  shell:
      """{SAMBAMBA} view -h -t 6 -f bam -F "[XS] == null and not unmapped and not duplicate" {input} > {output}"""

rule bam_indexing:
  input:
      "results/{smp}.filtered.bam"
  output:
      "results/{smp}.filtered.bam.bai"
  message:
      "Indexing BAM with samtools."
  shell:
      "{SAMTOOLS} index {input}"

```

ENA stopped working for a bit so I needed to the downloading of the raw files directly from SRA. I did this using [sraDownloader](https://github.com/s-andrews/sradownloader). This is not ideal because I find it is a bit slow from SRA, but at least it works. As you will see in the loop sequence below, there are some accession numbers I skip (EWS24 cell line). This is because these were sequenced using a paired-end method on a NovaSeq, so would need a different analysis pipeline (all of the other data are single-end HiSeq).

```shell
#!/bin/bash

##set the location of software tools and the working directory where files will be stored
sraDownloader="/home/chughes/softwareTools/sradownloader-3.8/sradownloader"
sraCacheLocation="/mnt/Data/chughes/sratoolsRepository"
workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211110_aynaudScRnaseqPmid32049009"
eval cd ${workingDirectory}
eval mkdir raw
eval mkdir results

##loop over the accessions
for i in SRR8832{666..674}
do
  printf "Downloading files associated with ${i}."
  eval ${sraDownloader} --noena --outdir ${workingDirectory}/raw ${i}
  ##the file gets renamed upon download, but I just want it to have the SRR id and I can annotate it later
  eval mv ${workingDirectory}/raw/${i}*.fastq.gz ${workingDirectory}/raw/${i}.fastq.gz
  #eval conda activate snakemake
  eval snakemake --cores 8
  #eval conda deactivate
  eval rm ${workingDirectory}/raw/${i}.fastq.gz
  eval rm ${sraCacheLocation}/sra/${i}.sra.cache
  eval rm ${workingDirectory}/results/${i}.clean.fastq.gz
  eval rm ${workingDirectory}/results/${i}.unsorted.sam
  eval rm ${workingDirectory}/results/${i}.unsorted.bam
  eval rm ${workingDirectory}/results/${i}.sorted.bam
  eval rm ${workingDirectory}/results/${i}.sorted.bam.bai
done
```

Now we can run the analysis, writing a log so we can track it. I am inside the working directory here, so I can just execute it. If you are not, you need to set the absolute path to the script. I also activate the conda environment. I am a complete snakemake n00b, so I haven't figured out the nuts and bolts of making it work with bash script integration, but this does the job for now.

```shell
conda activate snakemake
./sraDataProcessingScript.sh |& tee sraDataProcessingScriptLog.txt
```

Now we move on to peak calling using [MACS](https://github.com/macs3-project/MACS). There is a great walkthrough of this process [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html). Also [here](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) for discussion of [deeptools](https://deeptools.readthedocs.io/en/develop/), specifically bamCoverage.

The first file here is for narrow peak calling on the transcription factors and broad for histone marks. We also do the coverage call here.

```shell
#!/bin/bash
annotationLocation="/home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gtf"
rawDataOutputDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211110_aynaudScRnaseqPmid32049009"
for i in SRR8832669 SRR8832670 SRR8832671 SRR8832672 SRR8832673 SRR8832674
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR8832666.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01
  eval bamCoverage -b ${rawDataOutputDirectory}/results/${i}.filtered.bam -o ${rawDataOutputDirectory}/results/${i}.chr11.bw --binSize 10 --region chr11 --normalizeUsing BPM --smoothLength 30 --extendReads 150 --centerReads -p 6
done

for i in SRR8832667 SRR8832668
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR8832666.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 --broad
  eval bamCoverage -b ${rawDataOutputDirectory}/results/${i}.filtered.bam -o ${rawDataOutputDirectory}/results/${i}.chr11.bw --binSize 10 --region chr11 --normalizeUsing BPM --smoothLength 30 --extendReads 150 --centerReads -p 6
done
```

The rest I can do in R.