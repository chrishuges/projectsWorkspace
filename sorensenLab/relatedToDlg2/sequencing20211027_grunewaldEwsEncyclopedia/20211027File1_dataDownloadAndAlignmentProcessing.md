## Reprocessing some ChIPseq data

This document describes the reprocessing of some ChIPseq data from a previous publication. Specifically:

"The Ewing Sarcoma Cell Line Atlas (ESCLA)"
GEO: GSE176339, https://www.biorxiv.org/content/10.1101/2021.06.08.447518v1

### Description

I am interested in the ChIPseq data across the cell lines they profiled. There is a metadata file in the server directory detailed below that describes these data.

I am going to use [SRAExplorer](https://sra-explorer.info/#) to get at these data. I am going to save them in the directory `/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211027_grunewaldEwsEncyclopedia/`. I will parse these downloaded raw files for quality and adapter sequences using bbduk from the [bbTools package](https://sourceforge.net/projects/bbmap/). There are some good walkthroughs on how to use this package in the packages documentation itself, as well as [here](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/). Because the authors originally used Bowtie for the alignment, I will as well. For the index, I created one as described in `..\relatedToDlg2\sequencing20211029_databaseSetupAndSnakemake\20211029File1_downloadingDatabasesAndIndexes.md`. There is a pretty good walkthrough of the alignment process [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/03_align_and_filtering.html). I will then use [samtools](http://www.htslib.org/) and [sambamba](https://lomereiter.github.io/sambamba/) to sort and filter the sam and bam files that are generated, prior to peak calling. I ended up writing this into a script that uses snakemake. Proceed with the analysis as below.

### Data pipeline

First we will move into our working directory and create our shell and snakemake scripts.

```shell
cd /mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211027_grunewaldEwsEncyclopedia
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
BASE_DIR = "/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211027_grunewaldEwsEncyclopedia"


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
workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211027_grunewaldEwsEncyclopedia"
eval cd ${workingDirectory}
eval mkdir raw
eval mkdir results

##loop over the accessions
for i in SRR1476{{0997..1021},{1027..1087}}
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

Now we move on to peak calling using [MACS](https://github.com/macs3-project/MACS). There is a great walkthrough of this process [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html). Also [here](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) for discussion of deepTools, specifically bamCoverage.

The first file here is for narrow peak calling on the transcription factors.

```shell
#!/bin/bash
annotationLocation="/home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gtf"
rawDataOutputDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211027_grunewaldEwsEncyclopedia"
for i in SRR14760997
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761001.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done

for i in SRR14761002
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761006.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done

for i in SRR14761007
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761011.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done

for i in SRR14761012
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761016.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done

for i in SRR14761017
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761021.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done

for i in SRR14761027
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761031.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done

for i in SRR14761032
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761036.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done

for i in SRR14761037
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761041.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done

for i in SRR14761042
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761046.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done

for i in SRR14761047
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761051.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done

for i in SRR14761052
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761057.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done

for i in SRR14761058
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761062.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done

for i in SRR14761063
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761067.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done

for i in SRR14761068
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761072.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done

for i in SRR14761073
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761077.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done

for i in SRR14761078
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761082.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done

for i in SRR14761083
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761087.filtered.bam -f BAM -g hs -n ${rawDataOutputDirectory}/results/${i} -B -q 0.01 
done
```

This second file is for broad peak calling with the histone marks.

```shell
#!/bin/bash
annotationLocation="/home/chughes/databases/projectEwsDlg2/baseGenomeFiles/genome.gtf"
rawDataOutputDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211027_grunewaldEwsEncyclopedia"
for i in SRR14760998 SRR14760999 SRR14761000
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761001.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done

for i in SRR14761003 SRR14761004 SRR14761005
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761006.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done

for i in SRR14761008 SRR14761009 SRR14761010
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761011.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done

for i in SRR14761013 SRR14761014 SRR14761015
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761016.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done

for i in SRR14761018 SRR14761019 SRR14761020
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761021.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done

for i in SRR14761028 SRR14761029 SRR14761030
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761031.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done

for i in SRR14761033 SRR14761034 SRR14761035
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761036.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done

for i in SRR14761038 SRR14761039 SRR14761040
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761041.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done

for i in SRR14761043 SRR14761044 SRR14761045
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761046.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done

for i in SRR14761048 SRR14761049 SRR14761050
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761051.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done

for i in SRR14761053 SRR14761055 SRR14761056
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761057.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done

for i in SRR14761059 SRR14761060 SRR14761061
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761062.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done

for i in SRR14761064 SRR14761065 SRR14761066
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761067.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done

for i in SRR14761069 SRR14761070 SRR14761071
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761072.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done

for i in SRR14761074 SRR14761075 SRR14761076
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761077.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done

for i in SRR14761079 SRR14761080 SRR14761081
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761082.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done

for i in SRR14761084 SRR14761085 SRR14761086
do
  eval macs3 callpeak -t ${rawDataOutputDirectory}/results/${i}.filtered.bam -c ${rawDataOutputDirectory}/results/SRR14761087.filtered.bam -f BAM -g hs -n ${i} -B -q 0.01 --broad
done
```

Now we are ready for visualization. There is a good walkthrough of visualization [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/10_data_visualization.html). First, I ran [deeptools](https://deeptools.readthedocs.io/en/develop/) on the bam files to create coverage maps:

```shell
#!/bin/bash

##set the location of software tools and the working directory where files will be stored
workingDirectory="/mnt/Data/chughes/projectsRepository/sorensenLab/relatedToDlg2/sequencing20211027_grunewaldEwsEncyclopedia"

##loop over the accessions
for i in SRR1476{{0997..1021},{1027..1053},{1055..1087}}
do
  eval bamCoverage -b ${workingDirectory}/results/${i}.filtered.bam -o ${workingDirectory}/results/${i}.chr11.bw --binSize 10 --region chr11 --normalizeUsing BPM --smoothLength 30 --extendReads 150 --centerReads -p 6
done
```

The rest I can do in R.